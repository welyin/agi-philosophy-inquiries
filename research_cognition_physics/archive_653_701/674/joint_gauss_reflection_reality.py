"""674: exact scalar/mass-source reflection through complementary Pfaffians.

The proof applies without inverse Dirac propagators or divisions by weights.
It establishes Hermiticity/reality of673's full scalar Gauss functional, not RP.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_gauss_boundary_functional as prior
old=prior.old;internal=prior.internal;mass=prior.mass
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_gauss_reflection_reality_results.json'


def pf(a):
    return 1.+0j if len(a)==0 else prior.pf(a)


def ingredients(u,v,mat,lam=.37):
    jm,jp,m,mb,pair=mat;n=len(m);r=n//2
    p=lam*pair[:r,:r];q=lam*pair[r:,r:]
    b=jm.conj()@p@jm.conj().T+jp.conj()@np.linalg.inv(q)@jp.conj().T
    phase=np.linalg.det(np.column_stack((u,v)))
    au=pf(u.T@m@u);av=pf(v.T@m@v)
    bv=pf(v.T@b@v);du=pf(u.T@(-np.linalg.inv(b.conj()))@u)
    pf_b=pf(b);a=pf(q)
    return dict(B=b,M=m,au=au,av=av,bv=bv,du=du,phase=phase,pf_B=pf_b,a=a,
        aux_complement_error=float(abs(au-phase*av.conjugate())),
        mass_complement_error=float(abs(du-bv.conjugate()/(phase.conjugate()*pf_b.conjugate()))
            /max(abs(du),abs(bv),1.)),
        mass_prefactor_error=float(abs(a/pf_b.conjugate()-a.conjugate())/max(abs(a),1e-280)),
        factorized_weight=au*a*bv/phase)


def algebra_check():
    rng=np.random.default_rng(67411)
    z=rng.normal(size=(64,64))+1j*rng.normal(size=(64,64));v0,_=np.linalg.qr(z)
    g5=np.kron(internal.spin.G5,np.eye(16));g0=np.kron(internal.spin.GAMMA[3],np.eye(16));t=g5@g0
    mat=prior.fixed_matrices(np.eye(10)[0][None,:],mass.car.PHI[0][None,:])
    rows=[]
    for ru in (0,30,31,32,34,64):
        u=v0[:,:ru];v=v0[:,ru:];ut=t@v;vt=t@u
        d=(np.eye(64)+g5@(v@v.conj().T-u@u.conj().T))/2
        dt=(np.eye(64)+g5@(vt@vt.conj().T-ut@ut.conj().T))/2
        q=prior.regular(u,v,d,mat);qt=prior.regular(ut,vt,dt,mat)
        f=ingredients(u,v,mat)
        inverse_error=old.err(t.T@f['B']@t+np.linalg.inv(f['B'].conj()))
        eps=max(f['aux_complement_error'],f['mass_complement_error'],f['mass_prefactor_error'],inverse_error)
        assert eps<2e-11
        if ru%2==0:
            relative=float(abs(qt['weight']/q['weight'].conjugate()-1))
            factor=float(abs(q['weight']/f['factorized_weight']-1))
            assert relative<2e-11 and factor<2e-11
        else:
            relative=None;factor=None
            assert max(abs(q['weight']),abs(qt['weight']))<1e-24
        rows.append(dict(ru=ru,rv=64-ru,reflected_ru=64-ru,
            weight=old.cpair(q['weight']),reflected_weight=old.cpair(qt['weight']),
            reflection_relative_when_resolved=relative,factor_relative_when_resolved=factor,
            auxiliary_complement_absolute=f['aux_complement_error'],
            mass_complement_scaled=f['mass_complement_error'],
            mass_prefactor_relative=f['mass_prefactor_error'],mass_inverse_reflection_error=inverse_error))
    # Singular physical mass is a polynomial limit, not a forbidden configuration.
    phi=np.array([0.,0.,0.,0.,.4])
    singular=prior.fixed_matrices(np.eye(10)[0][None,:],phi[None,:])
    u=v0[:,:32];v=v0[:,32:];ut=t@v;vt=t@u
    d=(np.eye(64)+g5@(v@v.conj().T-u@u.conj().T))/2
    dt=g0@d.conj().T@g0
    a=prior.regular(u,v,d,singular);b=prior.regular(ut,vt,dt,singular)
    rel=float(abs(b['weight']/a['weight'].conjugate()-1))
    assert rel<3e-11 and np.linalg.matrix_rank(singular[-1])<64
    return dict(projector_fixtures_not_Wilson_topology_claim=True,
        determinant_fixed_chiral_frame=old.cpair(np.linalg.det(np.column_stack((mat[1],mat[0])))),
        auxiliary_full_Pfaffian=old.cpair(prior.pf(mat[2])),rows=rows,
        singular_original_Higgs_mass=dict(physical_pair_rank=int(np.linalg.matrix_rank(singular[-1])),
            reflection_relative=rel,proof_uses_polynomial_extension_not_inverse_mass=True))


def actual_check():
    e=np.random.default_rng(67321).normal(size=(4,10));e/=np.linalg.norm(e,axis=1)[:,None]
    phis=mass.car.PHI.copy()
    reflect=np.eye(4)[[2,3,0,1]]
    r4=np.kron(np.kron(reflect,np.eye(4)),np.eye(16))
    g5=np.kron(np.kron(np.eye(4),internal.spin.G5),np.eye(16))
    g0=np.kron(np.kron(np.eye(4),internal.spin.GAMMA[3]),np.eye(16))
    ell=r4@g0;t=g5@ell
    mat=prior.fixed_matrices(e,phis);matt=prior.fixed_matrices(e[[2,3,0,1]],phis[[2,3,0,1]])
    rows=[]
    for seed,scale in ((67326,.32),(67350,.9),(67351,.9),(67352,.9)):
        links=np.array([[prior.prior.rep(*prior.prior.group(seed+8*mu+i,scale)) for i in range(4)] for mu in range(2)])
        u,v,d,h,gap=prior.kernel(links)
        reversed_links=links[:,[2,3,0,1]].copy()
        reversed_links[1]=np.array([z.conj().T for z in links[1]])
        ur,vr,dr,hr,gapr=prior.kernel(reversed_links)
        f=ingredients(u,v,mat);ft=ingredients(ur,vr,matt)
        q=prior.regular(u,v,d,mat);qt=prior.regular(ur,vr,dr,matt)
        errors=dict(Wilson_reflection=old.err(hr+t@h@t.conj().T),
            overlap_reflection=old.err(dr-ell@d.conj().T@ell),
            auxiliary_reflection=old.err(t.T@matt[2]@t-mat[2]),
            mass_duality=old.err(t.T@ft['B']@t+np.linalg.inv(f['B'].conj())),
            auxiliary_complement_absolute=f['aux_complement_error'],
            mass_prefactor_relative=f['mass_prefactor_error'])
        assert max(errors.values())<4e-10
        smallest=float(np.linalg.svd(q['N'],compute_uv=False)[-1])
        resolved=smallest>1e-10
        rel=float(abs(qt['weight']/q['weight'].conjugate()-1)) if resolved else None
        if resolved:assert rel<3e-10 and f['mass_complement_error']<4e-10
        # Balance the absolute polynomial values; do not divide by unresolved
        # near-zero Pfaffians. The spectral/Clifford identities remain checked.
        assert u.shape==v.shape
        balance_scale=max(float(np.linalg.norm(f['B'],2)),float(np.linalg.norm(np.linalg.inv(f['B']),2)))
        du_bal=pf(u.T@(-np.linalg.inv(f['B'].conj())/balance_scale)@u)
        bv_bal=pf(v.T@(f['B']/balance_scale)@v)
        balanced=float(abs(du_bal-bv_bal.conjugate()/(f['phase'].conjugate()*f['pf_B'].conjugate())))
        assert balanced<2e-12
        rows.append(dict(seed=seed,gauge_scale=scale,Pfaffian_balance_scale=balance_scale,Wilson_gap=min(gap,gapr),errors=errors,
            weight=old.cpair(q['weight']),reflected_weight=old.cpair(qt['weight']),
            relative_only_when_resolved=rel,near_zero_phase_not_interpreted=not resolved,
            mass_complement_raw_relative_diagnostic=f['mass_complement_error'],
            mass_complement_relative_is_resolved=resolved,
            mass_complement_balanced_absolute_error=balanced,
            finite_matrix_reflection_proof_does_not_rely_on_Pfaffian_ratio=True))
    return dict(rows=rows,full_group_original_fields_and_masses=True,
        no_Haar_sampling_or_artificial_Hermitian_projection=True)


def run():
    deps=('joint_gauss_boundary_functional.py','research_note_643.md','research_note_670.md',
          'research_note_673.md','round674_drafts/unnormalized_source_probe_results.json')
    return dict(date='2026-10-02',round=674,tests_run=2,failures=0,errors=0,
        complementary_Pfaffian_reflection=algebra_check(),actual_full_gauge_reflection=actual_check(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope='Scalar and local mass-source reflection identity for original finite signed weight via complementary Pfaffians, including unequal chiral ranks and singular physical masses by polynomial extension. With673 absolute integrability and full two-boundary Haar measure, the scalar orbit kernel is Hermitian and total weight real. No positivity, nonzero normalization, arbitrary fermionic-source reflection, original Hamiltonian identity, continuum or quantum GR conclusion.',
        all_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=674,tests_run=2,all_checks_passed=True)))
