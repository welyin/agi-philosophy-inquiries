"""719: original sterile modes as system/reference and node-local records.

General full-H statements are analytic. Finite matrices retain original
coefficients but freeze bosonic configurations; no full Gibbs/field simulation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_record_history_transport as old
shared=old.shared
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_local_relational_record_results.json'


def annihilation(n,j):
    out=np.zeros((2**n,2**n),complex)
    for col in range(2**n):
        if col & (1<<j):
            sign=(-1)**((col & ((1<<j)-1)).bit_count())
            out[col^(1<<j),col]=sign
    return out


def projectors(X):
    eye=np.eye(len(X));xx=X@X
    ps=[(xx-X)/2,eye-xx,(xx+X)/2]
    assert np.linalg.norm(sum(ps)-eye)<1e-12
    for p in ps:assert np.linalg.norm(p@p-p)<1e-12
    return ps


def instrument(operators):
    joint=[np.eye(len(operators[0]))]
    for x in operators:joint=[p@k for p in projectors(x) for k in joint]
    return joint


def dephase(kraus,O):
    return sum(k@O@k.conj().T for k in kraus)


def pair_state(creators,sv,rv):
    vacuum=np.zeros(creators[0].shape[0],complex);vacuum[0]=1
    S=sum(v*c for v,c in zip(sv,creators))
    R=sum(v*c for v,c in zip(rv,creators))
    psi=S@R@vacuum
    assert abs(np.linalg.norm(psi)-1)<1e-12
    return psi


def setup(n,ids):
    cs=[annihilation(n,j) for j in range(n)]
    a,b,r,s=[cs[j] for j in ids]
    XA=a.conj().T@r+r.conj().T@a
    XB=b.conj().T@s+s.conj().T@b
    K=instrument([XA,XB])
    assert np.linalg.norm(XA@XB-XB@XA)<1e-12
    return cs,(a,b,r,s),XA,XB,K


def reference_record_check():
    # Order a,b,r,s; operators use CAR signs, not ordinary local tensor products.
    cs,(a,b,r,s),XA,XB,K=setup(4,[0,1,2,3])
    creators=[x.conj().T for x in cs]
    eye=np.eye(16);p=.75;q=1-p
    nf=(np.sqrt(p)*a+np.sqrt(q)*b).conj().T@(np.sqrt(p)*a+np.sqrt(q)*b)
    zr=r.conj().T@s
    effects=[k@k for k in K]
    RA=np.eye(16)-2*a.conj().T@a
    sumB=[sum(effects[3*j:3*j+3]) for j in range(3)]
    # instrument order: outer B result, inner A result.
    assert max(np.linalg.norm(RA@e@RA-e) for e in sumB)<1e-12
    rows=[]
    for eta in (1.,.4,0.):
        data=[]
        for sign in (1,-1):
            sv=np.array([1,sign,0,0])/np.sqrt(2)
            rays=[pair_state(creators,sv,np.array([0,0,1,z])/np.sqrt(2)) for z in (1,-1)]
            rho=sum(weight*np.outer(v,v.conj()) for weight,v in zip([(1+eta)/2,(1-eta)/2],rays))
            prob=np.array([np.trace(k@rho@k).real for k in K]).reshape(3,3)
            corr=float(np.trace(rho@XA@XB).real)
            original=float(np.trace(rho@nf).real)
            reconstructed=.5-2*np.sqrt(p*q)*corr/eta if eta else None
            after=dephase(K,rho)
            record_ref=max(abs(np.trace(zr@k@rho@k)) for k in K)
            assert abs(corr+sign*eta/2)<1e-12
            assert abs(np.trace(zr@after))<1e-12 and record_ref<1e-12
            if eta:assert abs(reconstructed-original)<1e-12
            data.append(dict(sign=sign,correlation=corr,ideal_mode_probability=original,
                             reconstructed_probability=reconstructed,
                             joint_probabilities=prob.tolist(),
                             B_marginal=prob.sum(axis=1).tolist(),
                             reference_coherence_after=float(abs(np.trace(zr@after)))))
        assert np.max(abs(np.array(data[0]['B_marginal'])-data[1]['B_marginal']))<1e-12
        joint_distance=float(np.sum(abs(np.array(data[0]['joint_probabilities'])-
                                        data[1]['joint_probabilities']))/2)
        if eta==0:assert joint_distance<1e-12
        rows.append(dict(reference_visibility=eta,states=data,joint_total_variation=joint_distance))
    kill=float(np.linalg.norm(dephase(K,zr)))
    assert kill<1e-12
    return dict(rows=rows,reference_coherence_operator_annihilation_error=kill,
                classical_statistics_not_same_Luders_instrument=True)


def spin_rotation(n,node,theta):
    u=np.eye(n,dtype=complex);sl=slice(32*node+30,32*node+32)
    u[sl,sl]=np.cos(theta)*np.eye(2)+1j*np.sin(theta)*np.array([[0,1],[1,0]])
    return u


def full_coefficient_check():
    rng=np.random.default_rng(71904);rows=[]
    sterile=np.array([30,31,62,63])
    for k in range(6):
        points=rng.normal(size=(2,5))*.25
        mass,dd,hop=shared.coefficients(points,[shared.group.sample(rng) for _ in range(2)],
                                       np.array([.21,.16]))
        h=mass+hop
        hh=np.zeros_like(h);ds=np.zeros_like(dd)
        for a in range(3):
            for b in range(3):
                u=spin_rotation(64,0,2*np.pi*a/3)@spin_rotation(64,1,2*np.pi*b/3)
                hh+=u@h@u.conj().T/9;ds+=u@dd@u.T/9
        expected=h.copy();expected[sterile,:]=0;expected[:,sterile]=0
        error=float(np.linalg.norm(hh-expected))
        pairerror=float(np.linalg.norm(ds-dd))
        assert max(error,pairerror)<1e-12
        rows.append(dict(sample=k,Dirac_and_hopping_removal_error=error,
                         original_Majorana_preservation_error=pairerror))
    # Exact 256-Fock source/energy calibration of the original neutral sector.
    data=old.fixture()
    H,G=data['H'],data['G']
    _,_,XA,XB,K=setup(8,[2,6,3,7])
    D=dephase(K,H)-H;DG=dephase(K,G)-G
    sp=old.spectrum(H);rho,rhod=old.thermal(sp,.73,G)
    mean=float(np.trace(rho@D).real)
    derivative=float(np.trace(rhod@D+rho@DG).real)
    finite=[]
    for eps in (2e-3,5e-4,1.25e-4):
        values=[]
        for lam in (.17+eps,.17-eps):
            new=old.fixture(lam)
            rr,_=old.thermal(old.spectrum(new['H']),.73)
            values.append(float(np.trace(rr@(dephase(K,new['H'])-new['H'])).real))
        fd=(values[0]-values[1])/(2*eps)
        finite.append(dict(step=eps,derivative=fd,error=abs(fd-derivative)))
    assert mean>0 and finite[-1]['error']<1e-8
    assert np.linalg.norm(dephase(K,XA@XB)-XA@XB)<1e-12
    return dict(original_64_mode_samples=rows,conditional_Gibbs_energy_increase=mean,
                conditional_Gibbs_total_source_derivative=derivative,
                finite_differences=finite,readout_product_preserved=True)


def shared_dynamics_check():
    data=old.fixture()
    cs,(a,b,r,s),XA,XB,K=setup(8,[2,6,3,7])
    sv=np.zeros(8,complex);sv[[2,6]]=[1,1j*.6];sv/=np.linalg.norm(sv)
    rv=np.zeros(8,complex);rv[[3,7]]=[1,np.exp(.3j)];rv/=np.linalg.norm(rv)
    psi=pair_state([c.conj().T for c in cs],sv,rv)
    zS=a.conj().T@b;zR=r.conj().T@s
    fS=a@b;fR=r@s
    sp=old.spectrum(data['H'])
    rows=[]
    for t in (0.,.6,1.7,3.1):
        v=old.evolve_columns(sp,t,psi[:,None])[:,0]
        expectation=lambda O:np.vdot(v,O@v)
        zs,zr,fs,fr=[expectation(O) for O in (zS,zR,fS,fR)]
        actual=float(expectation(XA@XB).real)
        factorized=float(2*np.real(np.conj(fs)*fr-zs*np.conj(zr)))
        correction=actual-factorized
        rows.append(dict(time=t,actual_local_record_correlation=actual,
                         product_of_current_marginals=factorized,
                         connected_system_reference_correction=correction,
                         reference_coherence_abs=float(abs(zr)),
                         system_reference_anomalous_product_real=float(np.real(np.conj(fs)*fr))))
    assert abs(rows[0]['connected_system_reference_correction'])<1e-12
    assert max(abs(r['connected_system_reference_correction']) for r in rows[1:])>1e-5
    return dict(rows=rows,original_Dirac_Majorana_hopping_all_retained=True,
                no_reset_or_new_autonomous_apparatus=True)


def run():
    dependencies=('research_note_524.md','research_note_598.md','research_note_623.md',
                  'research_note_633.md','research_note_634.md','research_note_706.md',
                  'research_note_709.md','research_note_718.md',
                  'joint_record_history_transport.py','joint_smooth_mode_contract.py')
    return dict(round=719,tests_run=3,failures=0,errors=0,
                relational_records=reference_record_check(),
                original_mass_and_source=full_coefficient_check(),
                actual_reference_dynamics=shared_dynamics_check(),
                scope=dict(local_Kraus_and_reference_preparation_are_explicit_inputs=True,
                           original_sterile_spin_modes_used_without_new_species=True,
                           full_H_domain_energy_and_coefficient_claims_analytic=True,
                           numerical_evolution_is_conditional_not_full_boson_Gibbs=True,
                           statistical_reconstruction_not_Luders_implementation=True,
                           no_free_reusable_reference_or_continuum_causal_device_claim=True,
                           old_633_counterexample_and_709_tools_reused=True,
                           goal_not_completed=True),
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest()
                                   for n in dependencies})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
