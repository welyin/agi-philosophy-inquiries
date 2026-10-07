"""835: pure internal compensation of the first two energy moments.

Finite spectral calibration of a proof for the original scalar Hamiltonian.
It is not a computation of that Hamiltonian's actual spectrum or controls.
"""
from pathlib import Path
import argparse,json
import numpy as np
HERE=Path(__file__).resolve().parent;TARGET=HERE/'pure_boson_moment_compensation_results.json'


def null_vector(psi,h,ids):
    columns=np.eye(len(psi),dtype=complex)[:,ids]
    constraints=np.vstack([psi.conj()@columns,(h@psi).conj()@columns])
    _,_,vh=np.linalg.svd(constraints,full_matrices=True)
    v=columns@vh[-1].conj()
    v/=np.linalg.norm(v)
    assert abs(np.vdot(psi,v))<2e-13 and abs(np.vdot(h@psi,v))<2e-12
    return v


def root_psd(k):
    eig,vec=np.linalg.eigh(k);assert min(eig)>-1e-12
    return (vec*np.sqrt(np.maximum(eig,0)))@vec.conj().T


def run():
    old=json.loads((HERE/'endpoint_energy_budget_results.json').read_text('utf-8'))
    delta=-old['organized_input_endpoint_packet']['full_H_squared_endpoint_difference']
    rng=np.random.default_rng(835)
    energies=np.array([0,1,2,3,5,8,13,21,35,100,300,900.],float)
    h=np.diag(energies)
    psi=(rng.normal(size=12)+1j*rng.normal(size=12))/(1+energies)
    psi[3:8]*=8
    psi/=np.linalg.norm(psi)
    random=rng.normal(size=(12,12))+1j*rng.normal(size=(12,12))
    k=.03*(random.conj().T@random)+.2*np.eye(12)
    assert np.linalg.norm(h@k-k@h)>1
    b=h@h+k
    mean=float(np.vdot(psi,h@psi).real)
    base=float(np.vdot(psi,b@psi).real)
    desired=base+delta
    assert mean>energies[2]
    low=null_vector(psi,h,[0,1,2]);high=null_vector(psi,h,[9,10,11])
    elo=float(np.vdot(low,h@low).real);ehi=float(np.vdot(high,h@high).real)
    p=(mean-elo)/(ehi-elo);assert 0<p<1
    chi=np.sqrt(1-p)*low+np.sqrt(p)*high
    upper=float(np.vdot(chi,b@chi).real)
    assert upper>=p*ehi**2-1e-9 and upper>desired
    assert abs(np.vdot(psi,chi))<2e-13
    assert abs(np.vdot(psi,h@chi))<2e-12
    a=0.;c=np.pi/2
    for _ in range(70):
        theta=(a+c)/2;eta=np.cos(theta)*psi+np.sin(theta)*chi
        value=float(np.vdot(eta,b@eta).real)
        if value<desired:a=theta
        else:c=theta
    eta=np.cos(theta)*psi+np.sin(theta)*chi
    errors=dict(norm=abs(np.vdot(eta,eta)-1),
                mean=abs(np.vdot(eta,h@eta)-mean),
                second=abs(np.vdot(eta,b@eta)-desired))
    assert max(errors.values())<2e-10
    # Embed both moment contracts into one finite Hamiltonian. Its two carrier
    # labels couple to different error sectors, retaining the positive K term.
    n=len(psi);full=np.zeros((4*n,4*n),complex)
    full[:n,:n]=h;full[n:2*n,n:2*n]=h
    rin=root_psd(k+delta*np.eye(n));rout=root_psd(k)
    full[:n,2*n:3*n]=rin;full[2*n:3*n,:n]=rin
    full[n:2*n,3*n:]=rout;full[3*n:,n:2*n]=rout
    before=np.zeros(4*n,complex);after=before.copy()
    before[:n]=psi;after[n:2*n]=eta
    moments=[]
    for power in (1,2,3):
        matrix=np.linalg.matrix_power(full,power)
        vi=float(np.vdot(before,matrix@before).real)
        vo=float(np.vdot(after,matrix@after).real)
        moments.append(dict(power=power,input=vi,output=vo,difference=vo-vi))
    assert abs(moments[0]['difference'])<1e-10
    assert abs(moments[1]['difference'])<1e-9
    assert abs(moments[2]['difference'])>1
    eig,vec=np.linalg.eigh(full)
    assert min(np.diff(eig))>1e-9
    pin=abs(vec.conj().T@before)**2;pout=abs(vec.conj().T@after)**2
    spectral_tv=float(.5*np.sum(abs(pin-pout)));assert spectral_tv>.1
    return dict(round=835,all_checks_passed=True,fresh_test_groups=1,
        compensation_is_one_pure_vector_not_a_discarded_mixture=True,
        native_endpoint_variance_deficit_used=delta,calibration_scalar_dimension=n,
        scalar_mean_energy=mean,low_energy=elo,high_energy=ehi,
        high_weight=p,initial_output_second=base,desired_second=desired,
        path_high_endpoint_second=upper,angle=theta,
        residuals={key:float(value) for key,value in errors.items()},
        full_single_H_moments=moments,full_spectral_weight_total_variation=spectral_tv,
        positive_error_term_not_dropped=True,
        original_scalar_spectrum_or_compensating_wavefunction_computed=False,
        original_other_source_statistics_preserved=False,
        native_autonomous_compensation_proven=False,
        scope='Pure fixed-mean-energy path with a noncommuting positive second-moment term. Two energy moments match in one diagnostic Hamiltonian while higher spectral data do not. The original-model existence proof uses its confining scalar spectrum; these finite eigenvalues are not the original physical spectrum.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();r=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
