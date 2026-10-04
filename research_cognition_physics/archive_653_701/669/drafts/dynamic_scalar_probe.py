"""669 entry: actual original scalar growth and nonconstant mass interface.

Finite reflection-paired background diagnostic, not integration of scalar paths.
Original full gauge process and real instruments remain separate requirements.
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;sys.path.insert(0,str(BASE))
import joint_mass_auxiliary_reflection as mass
import joint_operator_domain_completion as domain
old=mass.old


def run():
    rng=np.random.default_rng(66901)
    basis=[]
    for x in np.eye(5):
        p,_=mass.mass_pairing(1,domain.ball(x));basis.append(p)
    growth=[]
    for radius in (.2,1.,4.,12.,30.):
        x=rng.normal(size=5);x*=radius/np.linalg.norm(x)
        phi=domain.ball(x);p,_=mass.mass_pairing(1,phi)
        linear=sum(xi*pi for xi,pi in zip(x,basis))
        error=old.norm(p-linear)/(1+old.norm(p));assert error<4e-13
        u=float(domain.value_x(x));constants=domain.constants()
        tail=radius**2>=constants['t_threshold']
        if tail:assert u>=constants['kappa']*radius**4
        lr2=np.array([.5+.25*np.sin(phi[4]),.5-.25*np.sin(phi[4])])
        assert abs(sum(lr2)-1)<1e-14 and min(lr2)>=.25 and max(lr2)<=.75
        growth.append(dict(radius=radius,mass_linear_identity_relative_error=error,
            original_potential=u,quartic_tail_lower_bound=constants['kappa']*radius**4 if tail else None,
            original_s_record_squared=lr2.tolist()))
    # A genuinely varying phi on nx=2,nt=4, paired by the same time reflection.
    nx,nt=2,4;count=nx*nt;r=32*count
    positive=rng.normal(size=(4,5))*.5
    x=np.concatenate((positive.reshape(2,2,5)[::-1].reshape(4,5),positive))
    phis=domain.ball(x)
    n0,_,select,_=mass.physical(nx,nt,phis[0])
    pair=np.zeros((2*r,2*r),complex)
    for i,phi in enumerate(phis):
        p,_=mass.mass_pairing(1,phi);sl=slice(32*i,32*i+32);sr=slice(r+32*i,r+32*i+32)
        pair[sl,sl]=p[:32,:32];pair[sr,sr]=p[32:,32:]
    sites=[(t,z) for t in range(nt) for z in range(nx)]
    rr=np.eye(count)[[sites.index((nt-1-t,z)) for t,z in sites]]
    a=np.kron(rr,np.eye(32));zero=np.zeros_like(a)
    theta=np.block([[zero,a],[a,zero]])
    refl=old.err(pair+theta.T@pair.conj()@theta);assert refl<1e-13
    pm=select.T@pair@select;new=n0+.37*pm
    pf=old.old.old.pfaffian;ratio=pf(new)/pf(n0)
    assert ratio.real>0 and abs(ratio.imag)/ratio.real<2e-12
    g=select@(-np.linalg.inv(new))@select.T
    idx=np.flatnonzero(np.tile(np.repeat([t>=nt//2 for t,z in sites],32),2))
    gram=(theta@g)[np.ix_(idx,idx)]
    herm=old.err(gram-gram.conj().T);mineig=float(min(np.linalg.eigvalsh((gram+gram.conj().T)/2)))
    assert herm<3e-12 and mineig>-3e-12
    deps=('joint_mass_auxiliary_reflection.py','joint_operator_domain_completion.py',
          'research_note_603.md','research_note_623.md','research_note_624.md','research_note_643.md',
          'research_note_668.md')
    return dict(date='2026-10-02',entry_for_round=669,formal_round_complete=False,
        original_target_growth=growth,original_target_not_flat_Gaussian=True,
        nonconstant_reflection_paired_background=dict(nx=nx,nt=nt,phi=phis.tolist(),
            weight_ratio=old.cpair(ratio),mass_reflection_error=refl,
            reflection_Gram_minimum=mineig,reflection_hermiticity_error=herm),
        scalar_path_measure_not_integrated=True,full_Gauss_or_instrument_equivalence_not_claimed=True,
        dependency_hashes={p:hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in deps},
        all_checks_passed=True)


if __name__=='__main__':
    target=HERE/'dynamic_scalar_probe_results.json';result=run()
    if '--write-results' in sys.argv:
        with target.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(target.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('entry_for_round','formal_round_complete','all_checks_passed')}))
