"""679: equation-of-motion shift and contact cancellation for physical sources.

Exact reflection identities, not reflection positivity. Fixed original matrices;
arbitrary-epsilon fixtures test algebra, not realization by Wilson backgrounds.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_rational_physical_limit as old

base=old.base
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_physical_source_reflection_results.json'


def error(x):
    return float(np.max(abs(x),initial=0))


def reflection_matrices(mat,permutation):
    jm,jp=mat[:2];n=len(jm);r=jm.shape[1]
    ell=np.kron(permutation,np.kron(base.internal.spin.GAMMA[3],np.eye(16)))
    c=jp.T@ell@jm
    ry=np.block([[np.zeros((r,r)),c.T],[c,np.zeros((r,r))]])
    rx=np.block([[np.zeros((n,n)),ell.T],[ell,np.zeros((n,n))]])
    return ell,ry,rx


def algebra(eps,epst,mat,matt,perm):
    jm,jp,m,mb,pair=mat;n=len(m);r=n//2
    ell,ry,rx=reflection_matrices(mat,perm)
    q=old.soft(eps,mat,0);qt=old.soft(epst,matt,0)
    d=q['D'];g5=jp@jp.T-jm@jm.T
    s=np.block([[-.5*jm.T@m,jm.T],[jp.T,.5*jp.T@mb]])
    a=q['N'].conj();f=q['Phi'].conj()
    tilde=ry@qt['Phi']@rx
    direct=np.block([[jm.T,-jm.T@m@(np.eye(n)-.5*d.conj().T)],
                     [np.zeros((r,n)),jp.T@(np.eye(n)-d.conj().T)]])
    contact=s@f.T-f@s.T-s@a@s.T
    errors=dict(M_unitary=error(mb@m-np.eye(n)),M_skew=error(m+m.T),
        M_chiral=error(m@g5-g5@m),M_reflection=error(ell.T@matt[2]@ell-m),
        gamma5_Hermiticity=error(d.conj().T-g5@d@g5),
        D_reflection=error(qt['D']-ell@d.conj().T@ell),
        bare_N0_reflection=error(rx.T@qt['N']@rx+a),
        reflected_Phi_formula=error(tilde-direct),
        source_shift=error(tilde-f-s@a),source_contact=error(contact),
        mass_reflection=error(matt[-1]+ry.T@pair.conj()@ry))
    assert max(errors.values())<3e-12,errors
    singular=float(np.linalg.svd(a,compute_uv=False)[-1])
    covariance=None
    if singular>1e-8:
        covariance=error(tilde@np.linalg.solve(a,tilde.T)-f@np.linalg.solve(a,f.T))
        assert covariance<2e-10
    return dict(errors=errors,N0_smallest_singular=singular,
                covariance_shift_error_on_regular_set=covariance),ry,rx


def coefficients(eps,epst,mat,matt,ry,rx,z,lam):
    q=old.soft(eps,mat,lam);qt=old.soft(epst,matt,lam)
    singular=float(min(np.linalg.svd(q['N'],compute_uv=False)[-1],
                       np.linalg.svd(qt['N'],compute_uv=False)[-1]))
    rows=[]
    for k in (0,2,4,6):
        a=base.pf(q['N']) if not k else old.source_coeff(q,z[:,:k])
        zz=ry.T@z[:,:k].conj()[:,::-1]
        b=base.pf(qt['N']) if not k else old.source_coeff(qt,zz)
        absolute=float(abs(b-a.conjugate()))
        assert absolute<2e-10*max(1.,abs(a),abs(b))
        relative=None
        if singular>1e-8 and max(abs(a),abs(b))>1e-280:
            relative=float(absolute/max(abs(a),abs(b)))
            assert relative<2e-8,(k,relative)
        rows.append(dict(sources=k,original=base.old.cpair(a),
            reflected=base.old.cpair(b),absolute_error=absolute,
            relative_error_resolved_only=relative))
    return dict(lambda_mass=lam,N_smallest_singular=singular,coefficients=rows,
        bare_full_N_error=error(rx.T@qt['N']@rx+q['N'].conj()),
        bare_Phi_error=error(ry@qt['Phi']@rx-q['Phi'].conj()))


def general_algebra():
    rng=np.random.default_rng(67931)
    e=rng.normal(size=(1,10));e/=np.linalg.norm(e,axis=1)[:,None]
    full=base.fixed_matrices(e,base.mass.car.PHI[:1])
    singular_mass=base.fixed_matrices(e,np.zeros_like(base.mass.car.PHI[:1]))
    jm,jp=full[:2];g5=jp@jp.T-jm@jm.T
    ell,_,_=reflection_matrices(full,np.eye(1));t=g5@ell
    raw=rng.normal(size=(64,64))+1j*rng.normal(size=(64,64))
    v,_=np.linalg.qr(raw)
    z=rng.normal(size=(64,6))+1j*rng.normal(size=(64,6))
    z/=np.linalg.norm(z,axis=0)
    cases=[('soft',np.linspace(-.8,.9,64),full),
           ('unequal_projector_30_34',np.r_[-np.ones(30),np.ones(34)],full),
           ('odd_projector_31_33',np.r_[-np.ones(31),np.ones(33)],full),
           ('soft_singular_mass',np.linspace(-.8,.9,64),singular_mass),
           ('singular_massless_kernel',None,full)]
    rows=[]
    for label,ev,mat in cases:
        eps=np.eye(64) if ev is None else (v*ev)@v.conj().T
        epst=-t@eps@t.conj().T
        a,ry,rx=algebra(eps,epst,mat,mat,np.eye(1))
        b=coefficients(eps,epst,mat,mat,ry,rx,z,.37)
        rows.append(dict(case=label,algebra=a,reflection=b))
    assert rows[0]['reflection']['bare_full_N_error']>.01
    assert rows[0]['reflection']['bare_Phi_error']>.1
    return dict(rows=rows,arbitrary_epsilon_not_claimed_as_Wilson_realizations=True,
                singular_cases_not_divided_by_weight=True)


def dynamic_original():
    links,e,phis=old.fixture()
    perm=np.eye(4)[[2,3,0,1]]
    lt=links[:,[2,3,0,1]].copy()
    lt[1]=np.array([x.conj().T for x in links[1]])
    mat=base.fixed_matrices(e,phis)
    matt=base.fixed_matrices(e[[2,3,0,1]],phis[[2,3,0,1]])
    u,v,_,h,_=base.kernel(links);ut,vt,_,ht,_=base.kernel(lt)
    jm,jp=mat[:2];g5=jp@jp.T-jm@jm.T
    rng=np.random.default_rng(67932)
    z=rng.normal(size=(256,6))+1j*rng.normal(size=(256,6))
    z/=np.linalg.norm(z,axis=0)
    rows=[]
    for label,a,layers in [('exact',None,None),('soft',.23,1)]:
        eps=v@v.conj().T-u@u.conj().T if a is None else old.regulate(h,g5,a,layers)[0]
        epst=vt@vt.conj().T-ut@ut.conj().T if a is None else old.regulate(ht,g5,a,layers)[0]
        ident,ry,rx=algebra(eps,epst,mat,matt,perm)
        coeff=coefficients(eps,epst,mat,matt,ry,rx,z,.37)
        rows.append(dict(mode=label,a=a,L=layers,algebra=ident,reflection=coeff))
    return dict(original_nonflat_full_group=True,original_16_channels=True,
        original_variable_phi_and_E=True,rows=rows,
        full_Haar_S9_average_not_numerically_sampled=True,
        Hermiticity_after_full_average_proved_in_note=True)


def run():
    return dict(date='2026-10-02',round=679,tests_run=2,failures=0,errors=0,
        general_algebra=general_algebra(),dynamic_original=dynamic_original(),
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
            for name in ('research_note_661.md','research_note_668.md','research_note_673.md',
                         'research_note_674.md','research_note_677.md','research_note_678.md',
                         'joint_gauss_boundary_functional.py','joint_rational_physical_limit.py',
                         'round679_drafts/physical_reflection_probe_results.json')},
        scope=dict(all_physical_Grassmann_source_reflection=True,
            finite_soft_and_exact=True,zero_weights_and_singular_masses_retained=True,
            complete_physical_form_Hermitian=True,reflection_positivity=False,
            nonzero_normalization=False,original_HF_identification=False,
            continuum_or_quantum_GR=False,old_space_contracts_inherited=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=679,tests_run=2,all_checks_passed=True)))
