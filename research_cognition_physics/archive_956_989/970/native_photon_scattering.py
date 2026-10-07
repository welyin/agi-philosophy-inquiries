"""970: exact one-excitation, positive-band scattering using the native charge.
Rotating-wave and fixed optical support are declared effective-model inputs.
Finite wavepacket endpoint errors use Cook/Duhamel tails, not a finite lattice.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
OUT=HERE/"native_photon_scattering_results.json"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text("utf-8-sig"))
def load(name,p):
    spec=importlib.util.spec_from_file_location(name,p)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def run():
    old=load("native956",STAGE/"956/native_material_interface.py")
    m=old.material(1.,.1);J=m["J"];d=m["occupancy"];Delta=1+J
    cs=[old.annihilate(i) for i in range(4)];ns=[c.T@c for c in cs]
    sec=np.eye(16)[:,[s for s in range(16) if s.bit_count()==2]]
    Q=sec.T@(ns[0]+ns[1]-ns[2]-ns[3])@sec/2
    ground=math.sqrt(1-d)*m["s"]+math.sqrt(d)*m["d"]
    excited=Q@m["d"]
    upper=-math.sqrt(d)*m["s"]+math.sqrt(1-d)*m["d"]
    assert np.linalg.norm(Q@ground-math.sqrt(d)*excited)<1e-14
    assert np.linalg.norm(m["H"]@excited-excited)<1e-14
    assert np.linalg.norm(m["H"]@upper-(1+J)*upper)<1e-14
    assert np.linalg.norm(Q@excited-math.sqrt(d)*ground-math.sqrt(1-d)*upper)<1e-14
    triplets=m["single"]@np.array([[1,0,0],[0,1/math.sqrt(2),0],
                                 [0,1/math.sqrt(2),0],[0,0,1.]])
    assert np.linalg.norm(Q@triplets)<1e-14
    # In the old total-spin-zero code, the spectator pair is unchanged.
    xi=.1;lam=.05;G=lam*math.sqrt(d);G2=G*G
    omega_c=Delta;delta=G2/(20*xi);eps=.001
    assert 0<delta<2*xi and omega_c-2*xi-G>0
    def coeff(E):
        vg=np.sqrt(4*xi*xi-(E-omega_c)**2)
        z=1j*vg*(E-Delta)
        r=G2/(z-G2);t=z/(z-G2)
        return r,t
    # Independently solve the three real-space matching equations.
    match=0.;unitarity=0.;record_identity=0.;samples=[]
    for off in (-delta,-delta/2,0.,delta/2,delta):
        E=Delta+off;k=math.acos((omega_c-E)/(2*xi))
        z=np.exp(1j*k)
        A=np.array([[-1,1,0],[xi*z,E-omega_c+xi*z,-G],[0,-G,E-Delta]],complex)
        b=np.array([1,-xi/z,0],complex)
        r0,t0,e=np.linalg.solve(A,b);r,t=coeff(E)
        match=max(match,abs(r-r0),abs(t-t0))
        SM=np.array([[t,r],[r,t]])
        unitarity=max(unitarity,float(np.linalg.norm(SM.conj().T@SM-np.eye(2),2)))
        se=t+r;so=t-r
        record_identity=max(record_identity,abs(so-1),abs((1+se.real)/2-abs(t)**2))
        samples.append(dict(energy=E,reflection_probability=float(abs(r)**2),
            transmission_probability=float(abs(t)**2),
            local_excited_stationary_amplitude_abs=float(abs(e)),
            even_plus_record_probability=float((1+se.real)/2)))
    assert max(match,unitarity,record_identity)<1e-12
    Rlower=100/101
    assert min(x["reflection_probability"] for x in samples)>Rlower-1e-12
    # Fixed compact spectral envelope phi(E)=A/sqrt(delta)*(1-x^2)^2.
    # Normalization is exact: int_-1^1 (1-x^2)^4 dx=256/315.
    Aphi=math.sqrt(315)/16
    quad=[]
    for n in (96,192):
        x,w=np.polynomial.legendre.leggauss(n)
        density=Aphi*Aphi*(1-x*x)**4
        E=Delta+delta*x;r,t=coeff(E)
        norm=float(w@density);ref=float(w@(density*abs(r)**2))
        p_even=float(w@(density*abs(t)**2))
        meanE=float(w@(density*E))
        se=t+r
        # With material initially |+>, reduced coherence is <se>/2.
        coherence=complex(w@(density*se))
        quad.append(dict(points=n,normalization=norm,mean_photon_energy=meanE,
            reflection_probability=ref,even_plus_probability=p_even,
            odd_plus_probability=1.,record_contrast=1-p_even,
            material_coherence_multiplier=[coherence.real,coherence.imag]))
        assert abs(norm-1)<1e-12 and abs(meanE-Delta)<1e-12
        assert ref>=Rlower and abs(p_even+ref-1)<1e-12
    assert abs(quad[0]["reflection_probability"]-quad[1]["reflection_probability"])<1e-12
    # Explicit finite-time Cook tails. For w(E)=phi(E)/sqrt(v(E)):
    # |<0|free(t)|phi_even>| <= ||w''||_1/(sqrt(pi)*t^2).
    # Outgoing even amplitude is se(E)*w(E).
    vmin=math.sqrt(4*xi*xi-delta*delta)
    q0=(4*xi*xi-delta*delta)**(-.25)
    q1=delta/2*(4*xi*xi-delta*delta)**(-1.25)
    q2=.5*(4*xi*xi-delta*delta)**(-1.25)+1.25*delta*delta*(4*xi*xi-delta*delta)**(-2.25)
    phi0=Aphi*math.sqrt(delta)*16/15
    phi1=2*Aphi/math.sqrt(delta)
    phi2=Aphi/delta**1.5*32/(3*math.sqrt(3))
    W0=q0*phi0;W1=q0*phi1+q1*phi0;W2=q0*phi2+2*q1*phi1+q2*phi0
    x1=(2*xi+delta*delta/vmin)/G2
    x2=(3*delta/vmin+delta**3/vmin**3)/G2
    s1=2*x1;s2=4*x1*x1+2*x2
    Wout2=W2+2*s1*W1+s2*W0
    # Round tau upward with a further 1% slack. No time-window optimization.
    tau=math.ceil(1.01*G*(W2+Wout2)/(math.sqrt(math.pi)*eps))
    err=G*(W2+Wout2)/(math.sqrt(math.pi)*tau)
    finite_lower=Rlower-4*eps
    assert err<eps and finite_lower>.986
    # Coupling energy and exchange in the same one-excitation H, finite diagnostic.
    # No claim that this small open chain approximates the infinite scattering.
    n=7;F=omega_c*np.eye(n)-xi*(np.eye(n,k=1)+np.eye(n,k=-1))
    Hf=np.zeros((n+1,n+1));Hf[:n,:n]=F
    Hm=np.zeros_like(Hf);Hm[-1,-1]=Delta
    V=np.zeros_like(Hf);V[n//2,-1]=V[-1,n//2]=G
    H=Hf+Hm+V
    currents=[1j*(H@p-p@H) for p in (Hf,Hm,V)]
    current_balance=float(np.linalg.norm(sum(currents),2))
    assert current_balance<1e-14
    assert np.linalg.norm(currents[0],2)>0
    # Same parity inputs have the same FREE photon energy distribution.
    # Exact incoming scattering states inherit total spectral equality.
    # Finite bare preparations can differ in higher interacting-energy moments.
    cosmos=read(STAGE/"969/local_field_cosmology_results.json")
    par=cosmos["shared_background_parameters"];mass=cosmos["local_sources"]["complete_mass"]
    hub=math.sqrt((mass+par["R"])/par["mu"])
    mismatch=dict(cook_duration_over_969_T=2*tau/par["T"],
        initial_969_hubble_times_cook_duration=2*tau*hub,
        interpretation="this sufficient stationary window does not certify transport to the fixed 969 history")
    files=[STAGE/"956/native_material_interface.py",STAGE/"research_note_960.md",
           STAGE/"969/local_field_cosmology_results.json",STAGE/"research_note_969.md",
           HERE/"drafts/scattering_adoption_decision.md"]
    return dict(round=970,all_scientific_checks_passed=True,
        model="positive-band one-excitation rotating-wave native-charge scattering",
        parameters=dict(U=1.,v=.1,J=J,d=d,Delta=Delta,xi=xi,lambda_charge=lam,
            G=G,band_center=omega_c,band_min=omega_c-2*xi,band_max=omega_c+2*xi,
            half_packet_energy_width=delta),
        real_space_matching_error=float(match),two_port_unitarity_error=float(unitarity),
        record_probability_identity_error=float(record_identity),stationary_samples=samples,
        packet_rows=quad,
        analytic_bounds=dict(reflection_and_asymptotic_record_lower=Rlower,
            fixed_packet_phi_normalization="sqrt(315)/16",
            w_L1=W0,w_prime_L1=W1,w_second_L1=W2,
            outgoing_w_second_L1=Wout2,
            start_time=-tau,end_time=tau,
            endpoint_state_error_bound=err,reported_endpoint_error=eps,
            finite_record_contrast_lower=finite_lower,
            ideal_scattering_operator_error=2/math.sqrt(101),
            one_excitation_positive_energy_lower=omega_c-2*xi-G),
        resources=dict(incoming_parities_same_free_photon_energy_distribution=True,
            asymptotic_scattering_states_same_total_energy_distribution=True,
            finite_bare_preparations_same_total_mean_energy=True,
            finite_bare_preparations_full_spectral_equality_claimed=False,
            stationary_emitter_returns_to_its_initial_energy=True,
            original_relational_code_used=True,
            total_energy_current_balance=current_balance,
            interaction_energy_current_norm=float(np.linalg.norm(currents[2],2)),
            finite_packet_energy=Delta,
            phase_reference_and_waveguide_support_are_inputs=True,
            mechanical_recoil_dynamics_included=False),
        overlap_with_969=mismatch,
        scope=dict(actual_spatial_channel_input=True,spatial_dimension_derived=False,
            standard_bounded_task_effect=True,all_task_inputs_and_passive_references_tail_bound=True,
            same_full_965_Hamiltonian=False,Pauli_Fierz_matching_error_proved=False,
            same_969_cosmology_window_proved=False,irreversible_arrow_proved=False,
            full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files},
        references=["https://arxiv.org/abs/0802.4204"])
def compare(a,b,path=""):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+"/"+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for k,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f"/{k}")
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-12),(path,a,b)
    else:assert a==b,(path,a,b)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--write",action="store_true");a=p.parse_args()
    result=run()
    if a.write:
        with OUT.open("x",encoding="utf-8") as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
    else:compare(result,read(OUT))
    print(json.dumps({k:result[k] for k in ("round","all_scientific_checks_passed","parameters",
        "packet_rows","analytic_bounds","overlap_with_969")},ensure_ascii=False,indent=2))
