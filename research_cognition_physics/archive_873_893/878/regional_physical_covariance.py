"""878: local gauge representative identity and covariance/product calibration.
Finite Stueckelberg phase-space fixture, not the original gravity PDE.
The continuum trace-class and quasi-equivalence conclusions are analytic.
"""
from pathlib import Path
import argparse,json
import numpy as np
TARGET=Path(__file__).with_name("regional_physical_covariance_results.json")
def sym(a):return (a+a.conj().T)/2
def power(a,t):
    val,vec=np.linalg.eigh(sym(a))
    assert val.min()>-1e-10
    if t<0:assert val.min()>1e-12
    return (vec*np.maximum(val,0)**t)@vec.conj().T
def norm(a):return float(np.linalg.norm(a,2))
def block(a,b):
    z=np.zeros((len(a),len(b)))
    return np.block([[a,z],[z.T,b]])
def run():
    n=32;m=0.7
    eye=np.eye(n)
    d=(np.roll(eye,1,axis=0)-np.roll(eye,-1,axis=0))/2
    kq=np.vstack((d,m*eye))
    k=np.vstack((kq,np.zeros_like(kq)))
    left=np.zeros((n,4*n));left[:,n:2*n]=eye/m
    p= np.eye(2*n)-kq@np.linalg.solve(kq.T@kq,kq.T)
    r=block(p,p)
    j=np.block([[np.zeros((2*n,2*n)),np.eye(2*n)],[-np.eye(2*n),np.zeros((2*n,2*n))]])
    omega=power(eye-d@d,0.5)
    cov=0.5*r@block(block(omega,omega),block(np.linalg.inv(omega),np.linalg.inv(omega)))@r
    groups=[];local_errors=[];wrong_errors=[];tails=[]
    for seeds,interval in [([3,4],range(0,9)),([19,20],range(16,25))]:
        u=[]
        for seed in seeds:
            f=eye[:,seed]
            u.append(np.concatenate((-f/m,np.zeros(3*n))))
            u.append(np.concatenate((np.zeros(2*n),-m*f,d.T@f)))
        u=np.array(u).T;a=r@u
        chi=np.zeros(n);chi[list(interval)]=1
        chi0=np.diag(chi);chi4=np.diag(np.tile(chi,4))
        loc=chi4+(k@chi0-chi4@k)@left
        alpha=np.linalg.solve(kq.T@kq,kq.T@u[:2*n,:])
        v=loc@a
        expected=u-k@chi0@alpha
        # A lattice derivative has one-hop range. Reserve that buffer explicitly;
        # continuum differential operators preserve support exactly.
        collar=np.maximum.reduce((chi,np.roll(chi,1),np.roll(chi,-1)))
        mask=np.tile(collar,4)==0
        local_errors.append(dict(
            quotient_identity=norm(r@v-a),
            commutator_formula=norm(v-expected),
            gauge_constraint=norm(kq.T@v[2*n:,:]),
            outside_collar=norm(v[mask,:]),
            covariance_identity=norm(v.T@cov@v-a.T@cov@a)))
        wrong_errors.append(norm(r@chi4@a-a))
        tails.append(norm(a[mask,:]))
        groups.append((u,a,v))
    a,b=groups[0][1],groups[1][1]
    mu_a=a.T@cov@a;mu_b=b.T@cov@b
    normalize=block(power(mu_a,-0.5),power(mu_b,-0.5))
    both=np.concatenate((a,b),axis=1)
    gram=sym(normalize.T@both.T@cov@both@normalize)
    tau=normalize.T@both.T@j@both@normalize
    s=sym(gram+0.5j*tau)
    s0=sym(np.eye(len(gram))+0.5j*tau)
    delta=s-s0
    cross=gram[:len(mu_a),len(mu_a):]
    angle=norm(cross)
    ps_lhs=float(np.linalg.norm(power(s,0.5)-power(s0,0.5),'fro')**2)
    ps_rhs=float(np.linalg.svd(delta,compute_uv=False).sum())
    covariance=dict(cross_symplectic_defect=norm(a.T@j@b),
        regional_gram_lower_bound=float(np.linalg.eigvalsh(gram).min()),
        predicted_two_region_lower_bound=1-angle,
        cross_correlation_operator_norm=angle,
        positive_two_point_min=float(np.linalg.eigvalsh(s).min()),
        product_two_point_min=float(np.linalg.eigvalsh(s0).min()),
        sqrt_difference_hilbert_schmidt_squared=ps_lhs,
        covariance_difference_trace_norm=ps_rhs)
    assert norm(left@k-eye)<1e-12
    assert norm(r@r-r)<1e-12 and norm(r@k)<1e-12
    assert max(v for row in local_errors for v in row.values())<2e-12
    assert min(wrong_errors)>1e-3 and min(tails)>1e-3
    assert covariance['cross_symplectic_defect']<1e-12
    assert 0<angle<1 and abs(covariance['regional_gram_lower_bound']-(1-angle))<1e-12
    assert min(covariance['positive_two_point_min'],covariance['product_two_point_min'])>0
    assert ps_lhs<=ps_rhs+1e-12
    return dict(round=878,date='2026-10-06',fresh_numbered_groups=1,
        cumulative_numbered_groups=3663,all_checks_passed=True,
        fixture=dict(periodic_sites=n,gauge_mass=m,raw_phase_dimension=4*n,
            physical_phase_dimension=int(round(np.trace(r))),two_regional_modes_each=4),
        left_inverse_defect=norm(left@k-eye),local_representative_checks=local_errors,
        nonlocal_physical_representative_tails=tails,
        omitted_commutator_quotient_defects=wrong_errors,covariance_checks=covariance,
        analytic_scope='Original mixed bosonic free physical CCR in finitely many separated Cauchy regions inside one regular reference patch: original covariance is quasi-equivalent to product of its regional marginals. No global type-I collar, original energy/source estimates or interacting finite-coupling Q/E identification claimed.',
        finite_fixture_proves_continuum_theorem=False,full_goal_completed=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    result=run()
    if args.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
