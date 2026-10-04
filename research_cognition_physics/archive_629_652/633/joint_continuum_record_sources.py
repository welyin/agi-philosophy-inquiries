"""633: one original sterile mode ties record updates to continuum sources.

Free constant-background branch only; smooth compact Cauchy smearing.
Algebraic locality is explicitly separated from causal implementability.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_background_contact_matching as bg
import joint_quantum_response_matching as finite
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_continuum_record_sources_results.json'
IDS=[24,25,30,31]
SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1.,-1.])
PAULI=np.array([SX,SY,SZ])
PHI=np.array([0.,np.sqrt(bg.U0[0]),0.,0.,np.sqrt(bg.U0[1])])
hh,dd=bg.matter.mass_matrices(PHI)
H0=hh[np.ix_(IDS,IDS)]; DELTA=dd[np.ix_(IDS,IDS)]
B0=np.block([[H0,DELTA],[DELTA.conj().T,-H0.T]])
E=np.eye(8,dtype=complex)[:,2]
JCH=np.block([[np.zeros((4,4)),np.eye(4)],[np.eye(4),np.zeros((4,4))]])
Q=np.outer(E,E.conj())+np.outer(JCH@E.conj(),(JCH@E.conj()).conj())


def bdg(k):
    sk=np.einsum('a,aij->ij',k,PAULI)
    kin=np.block([[-sk,np.zeros((2,2))],[np.zeros((2,2)),sk]])
    return B0+np.block([[kin,np.zeros((4,4))],[np.zeros((4,4)),kin.T]])


def spectral(B):
    e,v=np.linalg.eigh(B)
    P=(v[:,e<0])@v[:,e<0].conj().T
    absolute=(v*abs(e))@v.conj().T
    return e,P,absolute


def mass_data():
    f=float(bg.matter.original.F(PHI))
    D=abs(bg.matter.Y['nu'])*PHI[1]/np.sqrt(f)
    M=abs(bg.matter.Y['s'])*PHI[4]/np.sqrt(f)
    root=np.sqrt(M*M+4*D*D)
    masses=np.array([(root-M)/2,(root+M)/2])
    weights=np.array([(1-M/root)/2,(1+M/root)/2])
    b=1/np.sqrt(6)/np.sqrt(np.sum(bg.U0)/(6*bg.matter.original.M))
    return masses,weights,float(D),float(M),float(b)


def original_mapping_check():
    masses,weights,D,M,b=mass_data()
    oldm,_=bg.canonical_mass_data(bg.U0)
    errs=[np.max(abs(oldm[-2:]-masses))]
    rng=np.random.default_rng(633)
    for k in rng.normal(size=(24,3)):
        B=bdg(k);kin=B-B0
        e,P,ab=spectral(B)
        expected=np.repeat(np.sqrt(k@k+masses*masses),2)
        errs.extend([np.max(abs(B-B.conj().T)),
                     np.max(abs(JCH@B.conj()@JCH+bdg(-k))),
                     np.max(abs(kin@B0+B0@kin)),
                     np.max(abs(e-np.r_[-expected[::-1],expected]))])
        errs.append(abs(np.vdot(E,ab@E).real-np.dot(weights,np.sqrt(k@k+masses*masses))))
    shell=[]
    for k0 in (0.,.3,1.1,3.):
        p=energy=source=trace_pressure=0.
        for axis in np.eye(3):
            for sign in (-1,1):
                k=sign*k0*axis;B=bdg(k);_,P,ab=spectral(B)
                p+=np.vdot(E,P@E).real/6
                energy+=np.vdot(E,ab@E).real/6
                # Reflection formula, integrated over parity-symmetric momenta.
                source+=(-2*np.vdot(E,P@B0@E).real)/6
                trace_pressure+=(-2*np.vdot(E,P@(B-B0)@E).real)/6
        expected_E=np.dot(weights,np.sqrt(k0*k0+masses*masses))
        expected_O=np.dot(weights,masses*masses/np.sqrt(k0*k0+masses*masses))
        errs.extend([abs(p-.5),abs(energy-expected_E),abs(source-expected_O),
                     abs(energy-source-trace_pressure)])
        shell.append(dict(k=k0,p=p,energy=energy,mass_source=source,
                          pressure_trace=trace_pressure))
    assert max(errs)<4e-14
    return dict(original_phi=PHI.tolist(),neutral_CAR_indices=IDS,
                neutral_masses=masses.tolist(),sterile_weights=weights.tolist(),
                Dirac_mass=D,Majorana_mass=M,radial_source_factor=b,
                max_mapping_error=float(max(errs)),shell_rows=shell,
                charged_spectators_retained=True,no_new_independent_massless_probe=True)


def exact_fock_check():
    aa=finite.annihilators(4);c=aa[2];n=c.conj().T@c;I=np.eye(16)
    H=finite.fock((H0,DELTA));e,v=np.linalg.eigh(H);vac=v[:,0]
    rho=np.outer(vac,vac.conj());Hrel=H-e[0]*I
    p=float(np.vdot(vac,n@vac).real)
    out=[(I-n)@vac,n@vac];probs=[1-p,p]
    post=[np.outer(x,x.conj())/pr for x,pr in zip(out,probs)]
    avg=sum(pr*r for pr,r in zip(probs,post))
    _,P,ab=spectral(B0)
    dP=-Q@P-P@Q+2*Q@P@Q
    psi=aa+[a.conj().T for a in aa]
    def covariance(state):
        return np.array([[np.trace(state@psi[j].conj().T@psi[i])
                          for j in range(8)] for i in range(8)])
    err=[np.max(abs(covariance(rho)-P)),
         np.max(abs(covariance(avg)-P-dP)),abs(p-.5),
         np.max(abs(n@n-n)),np.max(abs(sum(pr*r for pr,r in zip(probs,post))-avg))]
    energies=[float(np.trace(r@Hrel).real) for r in post]
    average_energy=float(np.trace(avg@Hrel).real)
    predicted=float(.5*np.trace(B0@dP).real)
    expected=float(np.vdot(E,ab@E).real)
    err.extend([abs(average_energy-predicted),abs(predicted-expected),
                max(abs(x-expected) for x in energies)])
    _,pairs=finite.matrices(finite.S0,IDS)
    J=finite.fock(pairs[1])
    JB=finite.matrices(finite.S0,IDS)[0][1]
    actual_source=float(np.trace((avg-rho)@J).real)
    spectral_source=float(.5*np.trace(JB@dP).real)
    err.append(abs(actual_source-spectral_source))
    assert min(np.linalg.eigvalsh(avg))>-2e-15
    assert max(err)<2e-14 and average_energy>.1 and abs(actual_source)>1e-4
    return dict(Fock_dimension=16,p=p,conditional_energy=energies,
                average_energy=average_energy,predicted_energy=expected,
                original_s_source_increment=actual_source,
                covariance_source_error=float(max(err)),
                nonselective_covariance_change_rank=int(np.linalg.matrix_rank(dP,tol=1e-10)),
                full_continuum_source_not_inferred_from_zero_mode_only=True)


def compact_bump(nr=192,nk=560,qmax=180.):
    r,wr=np.polynomial.legendre.leggauss(nr);r=(r+1)/2;wr=wr/2
    t=1-r*r;raw=np.exp(-1/t);N=1/np.sqrt(4*np.pi*np.dot(wr,r*r*raw*raw))
    f=N*raw
    grad=-2*r/t**2*f
    lap=(-6/t**2-8*r*r/t**3+4*r*r/t**4)*f
    K2=float(4*np.pi*np.dot(wr,r*r*grad*grad))
    K4=float(4*np.pi*np.dot(wr,r*r*lap*lap))
    q,wq=np.polynomial.legendre.leggauss(nk);q=(q+1)*qmax/2;wq=wq*qmax/2
    hat=np.sqrt(2/np.pi)*(np.sinc(np.outer(q,r)/np.pi)@(wr*r*r*f))
    probw=4*np.pi*wq*q*q*hat*hat
    return q,probw,N,K2,K4


def compact_continuum_check():
    masses,weights,D,M,b=mass_data()
    q,pw,N,K2,K4=compact_bump()
    q2,pw2,_,_,_=compact_bump(nr=232,nk=660)
    rows=[];errs=[]
    for ell in (1.,2.,4.):
        def moments(q,pw):
            k=q/ell;eng=np.sqrt(k[:,None]**2+masses**2)
            energy=float(np.dot(pw,eng@weights))
            source=float(np.dot(pw,(masses**2/eng)@weights))
            pressure=float(np.dot(pw,((k[:,None]**2/(3*eng))@weights)))
            return np.array([energy,source,pressure])
        energy,source,pressure=moments(q,pw)
        conv=float(np.max(abs(moments(q2,pw2)-[energy,source,pressure])))
        err=abs(energy-3*pressure-source);errs.append(err)
        upper=float(np.sqrt(K2/ell**2+D*D+M*M))
        lower=float(np.dot(weights,masses))
        # Analytic tail inequality, coefficient evaluated by quadrature.
        tail=float(np.sqrt(1+(masses.max()*ell/180)**2)*K4/(ell*180**3))
        assert lower<energy<upper and source>0 and pressure>0
        assert conv<2e-9 and err<2e-14
        rows.append(dict(support_radius=ell,record_probability=.5,
                         conditional_and_average_energy=energy,
                         mass_source_increment=source,radial_source_increment=b*source,
                         mean_integrated_pressure=pressure,energy_trace_error=err,
                         lower_energy_bound=lower,upper_energy_bound=upper,
                         tail_bound_coefficient_evaluated_numerically=tail,
                         quadrature_change=conv))
    norm=float(np.sum(pw))
    assert abs(norm-1)<2e-8
    # No numerical UV cutoff is used in the analytic existence/regularity proof.
    return dict(packet="N exp[-1/(1-r^2)] for r<1, zero otherwise",
                normalization=N,gradient_norm_squared=K2,laplacian_norm_squared=K4,
                momentum_cutoff_for_quadrature_only=180.,truncated_Fourier_norm=norm,
                norm_tail_upper_coefficient_evaluated_numerically=K4/180**4,
                rows=rows,Hadamard_and_all_energy_moments_proved_analytically=True,
                quadrature_not_interval_certified=True)


def causal_completion_check():
    a,b=finite.annihilators(2);theta=np.pi/6
    c=np.cos(theta)*a+np.sin(theta)*b
    n=c.conj().T@c;I=np.eye(4);nb=b.conj().T@b
    vacuum=np.array([1.,0.,0.,0.])
    plus=(a.conj().T+b.conj().T)@vacuum/np.sqrt(2)
    minus=(a.conj().T-b.conj().T)@vacuum/np.sqrt(2)
    U=I-2*a.conj().T@a
    err=float(np.linalg.norm(U@plus+minus))
    records=[];before=[];after=[]
    for v in (plus,minus):
        rho=np.outer(v,v.conj())
        avg=n@rho@n+(I-n)@rho@(I-n)
        before.append(float(np.trace(rho@nb).real))
        after.append(float(np.trace(avg@nb).real))
        records.append(float(np.trace(rho@n).real))
    signal=abs(after[1]-after[0])
    assert err<1e-15 and max(abs(x-.5) for x in before)<1e-15
    assert abs(signal-np.sqrt(3)/4)<1e-15
    assert np.max(abs((I-2*n)@(I-2*n)-I))<1e-15
    return dict(disjoint_smooth_sterile_modes=True,theta=float(theta),
                remote_initial_occupation=before,remote_after_instantaneous_readout=after,
                signal_difference=signal,local_phase_relation_error=err,
                record_probabilities=records,
                excluded="instantaneous distributed implementation of this whole-mode Lüders map",
                no_exclusion_of_finite_duration_causal_implementations=True,
                witness_state_differs_from_continuum_vacuum=True)


def run():
    mapping=original_mapping_check();fock=exact_fock_check()
    packet=compact_continuum_check();causal=causal_completion_check()
    deps=('joint_fermion_gauss_completion.py','joint_quantum_response_matching.py',
          'joint_background_contact_matching.py','joint_continuum_source_spectrum.py',
          'research_note_524.md','research_note_592.md','research_note_624.md','research_note_632.md')
    return dict(round=633,tests_run=4,failures=0,errors=0,mapping=mapping,
                exact_Fock=fock,compact_continuum=packet,causal_completion=causal,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope=dict(original_mass_mixing_retained=True,
                           free_constant_background_continuum_branch=True,
                           conditional_Hadamard_and_finite_energy_proved=True,
                           same_state_energy_mass_and_stress_sources=True,
                           local_counterterms_held_fixed_after_record=True,
                           finite_graph_to_continuum_mapping_open=True,
                           physical_record_device_and_internal_work_supply_open=True,
                           instantaneous_extended_instrument_rejected=True,
                           no_GR_or_dimension_or_SM_derivation_claim=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run();payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(payload)
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=633,tests=4,all_passed=True,
                         covariance_error=result['exact_Fock']['covariance_source_error'],
                         injection=result['compact_continuum']['rows'][0]['conditional_and_average_energy'],
                         instantaneous_signal=result['causal_completion']['signal_difference'])))

