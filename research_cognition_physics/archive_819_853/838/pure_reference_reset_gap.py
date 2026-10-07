"""838: exact covariance cost of resetting original pure CAR references.

Numerical matrices calibrate a projection identity applying analytically to
the original 730 pure Cauchy covariance. They are not a spacetime PDE solve.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'pure_reference_reset_gap_results.json'

def analyze(p,inside):
    n=len(p);q=np.diag([1.]*inside+[0.]*(n-inside));rest=np.eye(n)-q
    assert np.linalg.norm(p@p-p)<1e-12
    assert np.linalg.norm(p+p.conj()-np.eye(n))<1e-12
    a=p[:inside,:inside];b=p[:inside,inside:]
    new=q/2+rest@p@rest;delta=new-p
    variance=float(np.trace(a-a@a).real)
    square=float(np.linalg.norm(delta)**2)
    column_error=float(np.linalg.norm(delta[:,:inside].conj().T@delta[:,:inside]-np.eye(inside)/4))
    gram_error=float(np.linalg.norm(a-a@a-b@b.conj().T))
    assert column_error<1e-12 and gram_error<1e-12
    assert abs(square-inside/4-variance)<1e-12
    assert inside/4-1e-12<=square<=inside/2+1e-12
    assert np.linalg.eigvalsh(new).min()>-1e-12 and np.linalg.eigvalsh(new).max()<1+1e-12
    e=np.eye(n)[:,0];v=2j*delta@e
    assert np.linalg.norm(v.imag)<1e-12 and abs(np.linalg.norm(v)-1)<1e-12 and abs(e@v)<1e-12
    before=complex(2j*e@p@v);after=complex(2j*e@new@v)
    assert abs(before-1)<1e-12 and abs(after)<1e-12
    return dict(inside_self_dual_dimension=inside,
                cross_covariance_squared_norm=float(np.linalg.norm(b)**2),
                total_covariance_change_squared_HS_norm=square,
                projection_block_identity_residual=gram_error,
                every_internal_column_half_norm_identity_residual=column_error,
                even_unit_witness_expectation_original=float(before.real),
                even_unit_witness_expectation_reset=float(after.real),
                witness_may_require_region_outside_code=True)

def run():
    n=40;k=20
    j0=np.kron(np.eye(n//2),np.array([[0.,1.],[-1.,0.]]))
    rng=np.random.default_rng(838)
    rows=[]
    for _ in range(5):
        o,_=np.linalg.qr(rng.normal(size=(n,n)))
        p=(np.eye(n)+1j*o@j0@o.T)/2
        rows.append(analyze(p,k))
    reducing=analyze((np.eye(n)+1j*j0)/2,k)
    j_cross=np.block([[np.zeros((k,k)),np.eye(k)],[-np.eye(k),np.zeros((k,k))]])
    correlated=analyze((np.eye(n)+1j*j_cross)/2,k)
    assert abs(reducing['total_covariance_change_squared_HS_norm']-5)<1e-12
    assert abs(correlated['total_covariance_change_squared_HS_norm']-10)<1e-12
    # Independent full Fock check: a pure Gaussian Bell state versus its
    # product of unchanged one-mode marginals. No physical source is renamed.
    x=np.array([[0,1],[1,0]],complex);y=np.array([[0,-1j],[1j,0]]);z=np.diag([1.,-1.]);i=np.eye(2)
    gam=[np.kron(x,i),np.kron(y,i),np.kron(z,x),np.kron(z,y)]
    bell=np.array([1,0,0,1],complex)/np.sqrt(2)
    rho=np.outer(bell,bell.conj());mixed=np.eye(4)/4
    p=np.array([[np.trace(rho@a@b)/2 for b in gam] for a in gam])
    checked=analyze(p,2)
    delta=np.eye(4)/2-p
    v=2j*delta[:,0];gv=sum(v[a]*gam[a] for a in range(4));w=1j*gam[0]@gv
    assert np.linalg.norm(w-w.conj().T)<1e-12 and np.linalg.norm(w@w-np.eye(4))<1e-12
    initial=float(np.trace(rho@w).real);final=float(np.trace(mixed@w).real)
    assert abs(initial-1)<1e-12 and abs(final)<1e-12
    # Its local one-mode density is unchanged, despite the detectable joint
    # quadratic difference. Trace over either mode gives I/2 for both states.
    tensor=rho.reshape(2,2,2,2)
    marginal_error=max(np.linalg.norm(np.einsum('abcb->ac',tensor)-i/2),
                       np.linalg.norm(np.einsum('abad->bd',tensor)-i/2))
    assert marginal_error<1e-12
    return dict(round=838,all_checks_passed=True,fresh_test_groups=1,
                random_pure_self_dual_projection_rows=rows,
                reducing_subspace_endpoint=reducing,
                identical_internal_covariance_cross_only_endpoint=correlated,
                independent_two_mode_Fock_check=dict(
                    local_marginal_difference=float(marginal_error),
                    original_even_witness_expectation=initial,
                    reset_even_witness_expectation=final,projection_check=checked),
                original_730_pure_projection_used_in_analytic_proof=True,
                ten_mode_reset_covariance_squared_HS_interval=[5,10],
                exact_even_bounded_witness_difference=1,
                bound_includes_self_dual_doubling_convention=True,
                original_background_PDE_solved=False,
                specific_existing_local_source_nonzero_proven=False,
                original_fixed_region_contains_full_witness_proven=False,
                all_internal_preparation_routes_excluded=False,
                native_reset_or_compensation_implemented=False,
                scope='For the declared product-rest organized preparation, pure original covariance P is replaced by Q/2+(1-Q)P(1-Q). Projection purity fixes an unavoidable finite-rank covariance change and a bounded even bilinear witness. Logical input variations remain source-blind, but the absolute reference shift cannot be set to zero. This is not a no-go for other correlated preparations or internal compensation.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
