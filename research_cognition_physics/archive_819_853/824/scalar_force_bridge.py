"""824: original target metric, all 64 mass matrices and legal CAR rotation.

The original curved source identity is analytic. The numerical spinor and
spatial Gram below calibrate it; they do not evolve the original curved PDE.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'802'))
import original_bff_vertex as vertex
sys.path.insert(0,str(HERE.parent/'805'))
import original_past_covariance_action as original
sys.path.insert(0,str(HERE.parent/'819'))
import logical_input_source as packets
TARGET=HERE/'scalar_force_bridge_results.json'

def run():
    nums=vertex.NUMERATORS;Ns=nums[4];C=original.CHARGE
    rng=np.random.default_rng(824);error=0.;missing_metric=0.
    for _ in range(12):
        phi=original.PHI+.06*rng.normal(size=5);F=2-np.dot(phi,phi)/6
        kinv=F*(np.eye(5)-np.outer(phi,phi)/12)
        derivative=[vertex.dmass(phi,e) for e in np.eye(5)]
        for A in range(5):
            actual=sum(kinv[A,B]*derivative[B] for B in range(5))
            error=max(error,float(np.max(abs(actual-np.sqrt(F)*nums[A]))))
        wrong=F*derivative[4]
        missing_metric=max(missing_metric,float(np.max(abs(wrong-np.sqrt(F)*Ns))))
    assert error<1e-12 and missing_metric>1e-4
    assert np.max(abs(C@Ns.conj()@C+Ns))<1e-13
    s=np.zeros(64,complex)
    vals,vecs=np.linalg.eigh(original.GAMMA[2][np.ix_([30,31],[30,31])])
    s[[30,31]]=vecs[:,-1]
    nu=float(np.linalg.norm(Ns@s));assert nu>0
    partner=Ns@s/nu
    assert abs(np.vdot(s,partner))<1e-14
    assert np.linalg.norm(Ns@partner-nu*s)<1e-13
    anti=max(float(np.max(abs(g@Ns+Ns@g))) for g in original.GAMMA)
    assert anti<1e-13
    x,weights=np.polynomial.legendre.leggauss(192)
    chi=[]
    for center in (-.55,.55):
        v,_=packets.compact(x,center,.3);v/=np.sqrt(np.dot(weights,v*v));chi.append(v)
    q=1+x+.4*x*x
    # A trace-free metric-symbol calibration retaining original spin matrices.
    M=-.05*original.GAMMA[2]
    rows=[]
    for theta in (0.,.02,.05,.1):
        t=np.cos(theta)*s+np.sin(theta)*partner;ct=C@t.conj()
        isotropy=float(abs(np.vdot(ct,t)))
        ns=float(np.vdot(t,Ns@t).real)
        energy=float(np.vdot(t,original.GAMMA[2]@t).real)
        assert isotropy<1e-14 and abs(np.linalg.norm(t)-1)<1e-14
        assert abs(ns-nu*np.sin(2*theta))<1e-13
        assert abs(energy-np.cos(2*theta))<1e-13
        mode=np.stack([chi[0][:,None]*t,chi[1][:,None]*t,
                       chi[0][:,None]*ct,chi[1][:,None]*ct],axis=-1)
        gram=np.einsum('nik,n,nij->kj',mode.conj(),weights,mode)
        assert np.max(abs(gram-np.eye(4)))<1e-13
        action=np.einsum('ij,njk->nik',M,mode)*q[:,None,None]
        inside=np.einsum('nik,n,nij->kj',mode.conj(),weights,action)
        outside=action-np.einsum('nik,kj->nij',mode,inside)
        outgram=np.einsum('nik,n,nij->kj',outside.conj(),weights,outside)
        eig=np.linalg.eigvalsh(outgram);assert eig.min()>1e-8
        Dz=np.outer(t,t.conj())-np.outer(ct,ct.conj())
        Q_contrast=float(-.5*np.trace(Dz@Ns).real)
        assert abs(Q_contrast+ns)<1e-13
        # With weight sqrt(gamma)/sqrt(F), the acceleration contrast is -Delta Q.
        weighted_force_contrast=-Q_contrast
        rows.append(dict(theta=theta,CAR_isotropy_error=isotropy,
            singlet_mass_numerator_expectation=ns,principal_energy_polarization=energy,
            original_coefficient_leakage_Gram_min=float(eig.min()),
            Z_plus_minus_Q_Ns_contrast=Q_contrast,
            weighted_singlet_acceleration_contrast=weighted_force_contrast))
    assert rows[0]['singlet_mass_numerator_expectation']==0
    assert all(row['weighted_singlet_acceleration_contrast']>0 for row in rows[1:])
    return dict(round=824,all_checks_passed=True,
        all_64_original_mass_matrices_retained=True,
        target_metric_mass_derivative_identity_error=error,
        omitting_rank_one_kinverse_term_error=missing_metric,
        original_sterile_Ns_norm=nu,Clifford_Ns_anticommutator_error=anti,
        rows=rows,
        analytic_same_family_scalar_source_and_noise_existence=True,
        old_pure_particle_source_nonzero_claimed=False,
        selected_slice_is_Sigma_star_not_original_Sigma_zero=True,
        original_curved_propagation_or_source_at_Sigma_zero_computed=False,
        autonomous_preparation_or_all_geometry_no_go_proven=False,
        formal_test_groups_added=1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
