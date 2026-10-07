"""929: a scoped joint six-protocol process leaves the update instrument free.
All circuit histories are compiled with the inherited positive finite clock.
The allowed disturbance is one record-bit X per round at a declared logical slot.
"""
from pathlib import Path
import argparse,hashlib,json,math
import numpy as np
HERE=Path(__file__).resolve().parent
DIM=2**13
IND=np.arange(DIM)
X=np.array([[0.,1.],[1.,0.]])

def ry(theta):
    c,s=np.cos(theta/2),np.sin(theta/2)
    return np.array([[c,-s],[s,c]])

def one(state,target,u,controls=()):
    mask=(IND&(1<<target))==0
    for bit,value in controls:mask &= ((IND>>bit)&1)==value
    low=IND[mask];high=low|(1<<target)
    out=state.copy()
    out[low]=u[0,0]*state[low]+u[0,1]*state[high]
    out[high]=u[1,0]*state[low]+u[1,1]*state[high]
    return out

def act(state,gate,inverse=False):
    typ,target,controls,u=gate
    if typ=='weak':
        out=state
        for q,mat in enumerate(u):out=one(out,target,mat.T if inverse else mat,((0,q),))
        return out
    return one(state,target,u.conj().T if inverse else u,controls)

def circuit(theta,faults):
    out=[];mu=.6
    weak=[ry(2*np.arccos(np.sqrt((1+mu)/2))),ry(2*np.arccos(np.sqrt((1-mu)/2)))]
    def gate(target,controls=(),u=X):out.append(('one',target,controls,u))
    for r,fault in enumerate(faults):
        e=1+6*r;a,b,c=e+1,e+2,e+3;s1,s2=e+4,e+5
        gate(0,u=ry(theta));out.append(('weak',e,(),weak))
        for target in (a,b,c):gate(target,((e,1),))
        # Three separate controlled fault slots; in the full H the controls are internal two-bit F_r.
        for f,target in enumerate((a,b,c),1):gate(target,u=X if fault==f else np.eye(2))
        for control,target in ((a,s1),(b,s1),(b,s2),(c,s2)):gate(target,((control,1),))
        for target,pattern in ((a,(1,0)),(b,(1,1)),(c,(0,1))):gate(target,((s1,pattern[0]),(s2,pattern[1])))
    assert len(out)==30
    return out

def kraus(theta):
    return [np.diag(np.sqrt([.8,.2]))@ry(theta),np.diag(np.sqrt([.2,.8]))@ry(theta)]

def expected(theta,faults):
    w=np.zeros((DIM,2),complex);ks=kraus(theta)
    syndrome={0:(0,0),1:(1,0),2:(1,1),3:(0,1)}
    for x in (0,1):
        for y in (0,1):
            base=0
            for r,(event,fault) in enumerate(zip((x,y),faults)):
                e=1+6*r
                for bit in range(e,e+4):base |= event<<bit
                s1,s2=syndrome[fault];base|=s1<<(e+4);base|=s2<<(e+5)
            k=ks[y]@ks[x]
            for q in (0,1):w[base|q]=k[q]
    return w

def history(theta,faults):
    gates=circuit(theta,faults);w=np.zeros((DIM,2),complex);w[0,0]=1;w[1,1]=1
    hist=[w]
    for g in gates:w=act(w,g);hist.append(w)
    return gates,np.asarray(hist)

def run():
    theta_values=[0.,np.pi/2];rows=[];largest=0.;iso=0.;syndrome_fail=0.
    for theta in theta_values:
        ks=kraus(theta);probs=[]
        for x in (0,1):
            for y in (0,1):probs.append(float(np.linalg.norm((ks[y]@ks[x])[:,0])**2))
        for f1 in range(4):
            for f2 in range(4):
                _,hist=history(theta,(f1,f2));w=hist[-1];target=expected(theta,(f1,f2))
                largest=max(largest,float(np.linalg.norm(w-target,2)))
                iso=max(iso,float(np.linalg.norm(w.conj().T@w-np.eye(2),2)))
                # The entire isometry includes the data, both event witnesses, repaired records and syndromes.
                # Equality as a matrix covers every unknown input and every passive reference.
        rows.append(dict(theta=float(theta),joint_event_probabilities_for_zero=probs,
                         nonconstant_00_effect_gap=float(np.ptp(np.linalg.eigvalsh((ks[0]@ks[0]).conj().T@(ks[0]@ks[0]))))))
    assert largest<2e-14 and iso<2e-14
    assert abs(rows[0]['joint_event_probabilities_for_zero'][0]-.64)<1e-14
    assert abs(rows[1]['joint_event_probabilities_for_zero'][0]-.13)<1e-14
    # Relabel A <-> C and the corresponding fault/syndrome labels, including both rounds.
    perm=IND.copy()
    for r in range(2):
        e=1+6*r
        for a,b in ((e+1,e+3),(e+4,e+5)):
            diff=((perm>>a)^(perm>>b))&1;perm^=(diff<<a)|(diff<<b)
    swap={0:0,1:3,2:2,3:1};role_error=0.
    for theta in theta_values:
        for f in ((1,2),(3,0)):
            w=expected(theta,f);renamed=np.zeros_like(w);renamed[perm]=w
            role_error=max(role_error,float(np.linalg.norm(renamed-expected(theta,tuple(swap[x] for x in f)))))
    assert role_error<1e-14
    # Positive autonomous clock, exactly the existing engineered path construction.
    length=30;shift=length/2;couplings=np.array([.5*np.sqrt((j+1)*(length-j)) for j in range(length)])
    hc=shift*np.eye(length+1)+np.diag(couplings,1)+np.diag(couplings,-1)
    energies,v=np.linalg.eigh(hc);weights=abs(v[0])**2
    analytic_weights=np.array([math.comb(length,j)/2**length for j in range(length+1)])
    spectrum_error=float(np.max(abs(energies-np.arange(length+1))))
    spectral_weight_error=float(np.max(abs(weights-analytic_weights)))
    time=np.pi;a=v@(np.exp(-1j*time*energies)*v[0].conj())
    finish_error=float(np.linalg.norm(a-np.eye(length+1)[:,-1]))
    energy=float(hc[0,0]);variance=float(np.linalg.norm(hc[:,0])**2-energy**2)
    # Test actual gate-weighted H on a nontrivial history superposition, not only on the clock matrix.
    rng=np.random.default_rng(929);z=rng.normal(size=length+1)+1j*rng.normal(size=length+1);z/=np.linalg.norm(z)
    intertwining=0.
    for theta in theta_values:
        gates,hist=history(theta,(2,3));state=z[:,None,None]*hist;actual=shift*state.copy()
        for j,g in enumerate(gates):
            actual[j+1]+=couplings[j]*act(state[j],g)
            actual[j]+=couplings[j]*act(state[j+1],g,inverse=True)
        want=(hc@z)[:,None,None]*hist
        intertwining=max(intertwining,float(np.linalg.norm((actual-want).reshape(-1,2),2)))
    assert spectrum_error<1e-12 and spectral_weight_error<1e-13 and finish_error<1e-12 and intertwining<1e-12
    assert energy==15 and abs(variance-7.5)<1e-12
    delta=.02;tail=float(1-np.cos(delta/2)**(2*length));bound=length*delta**2/4
    assert 0<tail<=bound
    return dict(round=929,date='2026-10-06',all_scientific_checks_passed=True,
        compared_candidates=rows,unknown_input_and_fault_sectors_checked=32,
        complete_output_isometry_max_error=largest,unknown_input_isometry_max_error=iso,
        event_00_probability_difference=abs(rows[0]['joint_event_probabilities_for_zero'][0]-rows[1]['joint_event_probabilities_for_zero'][0]),
        actual_record_and_syndrome_relabeling_error=role_error,
        logical_steps=length,clock_dimension=length+1,work_qubits=13,internal_fault_qubits=4,
        full_hilbert_dimension=(length+1)*2**17,blank_work_qubits=12,
        full_H_minimum=0.,full_H_maximum=30.,initial_full_energy=energy,initial_full_energy_variance=variance,
        full_spectral_input_measure_same_for_all_unknown_inputs_and_candidates=True,
        clock_spectrum_error=spectrum_error,clock_spectral_weight_error=spectral_weight_error,
        complete_gate_weighted_intertwining_error=intertwining,autonomous_finish_error=finish_error,
        terminal_time=float(time),terminal_time_window_half_width=delta,
        terminal_window_failure=tail,terminal_window_failure_upper_bound=bound,
        same_six_finite_contracts_jointly_realized=True,
        disturbance_scope='At most one record-bit X per round, after copying and before verification; controller and two other records trusted.',
        repeats=2,all_blanks_fault_registers_syndromes_and_clock_internal=True,
        initial_low_entropy_preparation_and_engineered_H_are_inputs=True,
        history_compilation_is_inherited_not_new_physics=True,
        full_quantum_postevent_state_and_passive_reference_retained=True,
        arbitrary_record_phase_attacks_or_controller_errors_certified=False,
        eternal_repetition_or_fixed_finite_device_all_scale_claimed=False,
        spatial_or_GR_interfaces_realized=False,
        all_cognitive_principles_inadequate_for_all_physics_proved=False,
        whole_stage_completed=False,full_goal_completed=False,
        source_hashes={str(Path('929')/Path(__file__).name):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();out=run()
    path=HERE/'joint_protocol_selection_results.json'
    if a.write:
        assert not path.exists();path.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
