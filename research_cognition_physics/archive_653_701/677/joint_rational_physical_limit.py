"""677: bounded full-spectrum regulator and whole physical-source convergence.

Finite original two-time candidate; no physical-time/RP/HF identification.
The note proves dominated convergence. Numerical groups audit the actual
dictionary, covariance and nonzero source coefficients, not a Haar integral.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_gauss_boundary_functional as base

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_rational_physical_limit_results.json'
GOLD=(1+np.sqrt(5))/2
ALPHA=np.sqrt(5)/4
LAM=.37


def norm(a):
    return float(np.linalg.norm(a,2))


def soft(eps,mat,lam=LAM):
    jm,jp,m,mb,pair=mat
    n=len(eps);r=n//2;eye=np.eye(n)
    qp=jp@jp.conj().T;qm=jm@jm.conj().T;g5=qp-qm
    pu=(eye-eps)/2;pv=(eye+eps)/2
    d=(eye+g5@eps)/2
    phi=np.block([[jm.conj().T@pv,np.zeros((r,n))],
                  [-jp.T@m@(pu+.5*pv),jp.T]])
    n0=np.block([[m@qp,-d.T],[d,mb@qm]])
    insertion=phi.T@pair@phi
    return dict(D=d,Phi=phi,N=n0+lam*insertion,insertion=insertion,Pu=pu)


def regulate(h,g5,a,layers):
    x=g5@h;n=len(h)
    ha=g5@(a*x)@np.linalg.inv(2*np.eye(n)+a*x)
    herm=norm(ha-ha.conj().T);assert herm<2e-13
    ev,v=np.linalg.eigh((ha+ha.conj().T)/2)
    assert np.max(abs(ev))<1
    values=np.tanh(layers*np.arctanh(ev))
    eps=(v*values)@v.conj().T
    return eps,ha,ev


def source_coeff(q,z):
    k=z.shape[1]
    c=q['Phi'].T@z
    return base.pf(np.block([[q['N'],c],[-c.T,np.zeros((k,k))]]))


def fixture():
    links=np.tile(np.eye(16,dtype=complex),(2,4,1,1))
    for i in range(4):
        links[0,i]=base.prior.rep(*base.prior.group(67322+i,.16))
    for x in range(2):
        links[1,x]=base.prior.rep(*base.prior.group(67326+x,.32))
        links[1,2+x]=base.prior.rep(*base.prior.group(67328+x,.38))
    rng=np.random.default_rng(67321)
    e=rng.normal(size=(4,10));e/=np.linalg.norm(e,axis=1)[:,None]
    return links,e,base.mass.car.PHI.copy()


def original_source_limit():
    links,e,phis=fixture()
    u,v,d,h,gap=base.kernel(links)
    mat=base.fixed_matrices(e,phis)
    g5=mat[1]@mat[1].conj().T-mat[0]@mat[0].conj().T
    target=v@v.conj().T-u@u.conj().T
    exact=soft(target,mat)
    old=base.regular(u,v,d,mat)
    assert norm(exact['N']-old['N'])<1e-12
    rng=np.random.default_rng(67731)
    z=rng.normal(size=(256,4))+1j*rng.normal(size=(256,4))
    z/=np.linalg.norm(z,axis=0)
    w=base.pf(exact['N'])
    sing=float(np.linalg.svd(exact['N'],compute_uv=False)[-1])
    assert sing>1e-6 and abs(w)>1e-30
    coeff=[source_coeff(exact,z[:,:k]) for k in (2,4)]
    deriv=.5*w*np.trace(np.linalg.solve(exact['N'],exact['insertion']))
    rows=[];pnorm=norm(mat[-1]);hnorm=norm(h)
    assert hnorm<=3+1e-12
    for a,layers in ((.125,1024),(.03125,16384),(.0078125,262144),(.001953125,4194304)):
        eps,ha,ev=regulate(h,g5,a,layers)
        q=soft(eps,mat);delta=norm(eps-target)
        perturb=a*hnorm**2/(2*(1-a*hnorm/2))
        gap_b=gap-perturb
        spectral=2*np.exp(-2*layers*min(abs(np.arctanh(ev))))
        analytic_bound=float(spectral+perturb/gap_b) if gap_b>0 else 2.
        assert delta<=analytic_bound+1e-11
        dp=norm(q['Phi']-exact['Phi']);dn=norm(q['N']-exact['N'])
        assert dp<=ALPHA*delta+1e-12
        assert dn<=(.5+2*LAM*GOLD*ALPHA*pnorm)*delta+1e-12
        assert norm(q['Phi'])<=GOLD+1e-12
        assert norm(q['N'])<=2+LAM*GOLD**2*pnorm+1e-12
        weight=base.pf(q['N']);cs=[source_coeff(q,z[:,:k]) for k in (2,4)]
        ds=.5*weight*np.trace(np.linalg.solve(q['N'],q['insertion']))
        rows.append(dict(a=a,L=layers,sign_error=delta,sign_bound=analytic_bound,
            physical_map_error=dp,full_matrix_error=dn,
            scalar_weight=base.old.cpair(weight),scalar_absolute_error=float(abs(weight-w)),
            physical_source_coefficients=[base.old.cpair(c) for c in cs],
            source_absolute_errors=[float(abs(c-b)) for c,b in zip(cs,coeff)],
            mass_derivative=base.old.cpair(ds),mass_derivative_absolute_error=float(abs(ds-deriv)),
            scalar_relative_error_resolved_fixture_only=float(abs(weight/w-1))))
    assert rows[-1]['full_matrix_error']<rows[0]['full_matrix_error']/30
    assert rows[-1]['scalar_absolute_error']<rows[0]['scalar_absolute_error']/10
    assert max(rows[-1]['source_absolute_errors'])<max(rows[0]['source_absolute_errors'])/10
    return dict(original_nonflat_full_group=True,original_16_channels=True,
        original_variable_phi_and_S9=True,Wilson_gap=gap,Wilson_norm=hnorm,
        target_smallest_N_singular=sing,target_weight=base.old.cpair(w),
        target_source_coefficients=[base.old.cpair(c) for c in coeff],
        target_mass_derivative=base.old.cpair(deriv),rows=rows,
        inverse_used_only_for_nonsingular_derivative_numerical_check=True,
        complete_average_convergence_is_analytic_not_sampled=True)


def finite_covariance_and_regulator():
    links,e,phis=fixture();mat=base.fixed_matrices(e,phis)
    _,_,_,h,_=base.kernel(links)
    g5=mat[1]@mat[1].conj().T-mat[0]@mat[0].conj().T
    a=.23;layers=7
    eps,ha,ev=regulate(h,g5,a,layers);q=soft(eps,mat)
    plus=np.linalg.matrix_power(np.eye(256)+ha,layers)
    minus=np.linalg.matrix_power(np.eye(256)-ha,layers)
    direct=(plus-minus)@np.linalg.inv(plus+minus)
    rational_error=norm(eps-direct);assert rational_error<2e-12
    groups=[base.prior.group(67740+i,.24) for i in range(4)]
    lt,et,pt,rs=base.transform(links,e,phis,groups)
    _,_,_,ht,_=base.kernel(lt);mt=base.fixed_matrices(et,pt)
    epst,_,_=regulate(ht,g5,a,layers);qt=soft(epst,mt)
    def four_blocks(blocks):
        return base.old.diag(base.old.diag(*blocks[:2]),base.old.diag(*blocks[2:]))
    r4=four_blocks([np.kron(np.eye(4),r) for r in rs])
    r2=four_blocks([np.kron(np.eye(2),r) for r in rs])
    gx=base.old.diag(r4,r4.conj());gy=base.old.diag(r2,r2.conj())
    errs=dict(epsilon=norm(epst-r4@eps@r4.conj().T),
        physical_source_map=norm(qt['Phi']@gx-gy@q['Phi']),
        full_mass_congruence=norm(gx.T@qt['N']@gx-q['N']),
        skew=norm(q['N']+q['N'].T))
    assert max(errs.values())<3e-12
    rng=np.random.default_rng(67744)
    z=rng.normal(size=(256,2))+1j*rng.normal(size=(256,2))
    z/=np.linalg.norm(z,axis=0)
    c=source_coeff(q,z);ct=source_coeff(qt,gy.conj()@z)
    w=base.pf(q['N']);wt=base.pf(qt['N'])
    assert np.linalg.svd(q['N'],compute_uv=False)[-1]>1e-6
    werr=float(abs(wt/w-1));cerr=float(abs(ct-c)/max(abs(ct),abs(c)))
    assert max(werr,cerr)<3e-10
    defect=norm(q['Pu']@q['Pu']-q['Pu'])
    assert defect>1e-3
    return dict(a=a,L=layers,epsilon_spectral_radius=float(max(abs(np.tanh(layers*np.arctanh(ev))))),
        rational_identity_error=rational_error,covariance_errors=errs,
        scalar_phase_preserving_relative_error=werr,source_phase_preserving_relative_error=cerr,
        nonprojector_defect=defect,finite_projectors_not_exact=True,
        no_exact_chiral_factorization_applied=True,no_finite_RP_claim=True)


def run():
    deps=('research_note_659.md','research_note_669.md','research_note_673.md','research_note_676.md',
          'joint_gauss_boundary_functional.py','round677_drafts/rational_limit_probe_results.json')
    return dict(date='2026-10-02',round=677,tests_run=2,failures=0,errors=0,
        physical_source_limit=original_source_limit(),
        covariance_and_soft_regulator=finite_covariance_and_regulator(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope='At fixed finite box and positive original bosonic heat time: full bounded rational regulator converges coefficientwise and in weighted L1 to673 original unnormalized source functional, without uniform Wilson gap or rank restriction. Full gauge/S9/field averages and finite mass derivatives included analytically. Fifth coordinate, finite regulator positivity, local bulk measure and original HF physical time are not identified. No continuum volume limit, normalized state, quantum GR or spatial-dimension derivation.',
        all_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=677,tests_run=2,all_checks_passed=True,
        convergence=result['physical_source_limit']['rows'],
        covariance=result['covariance_and_soft_regulator'])))
