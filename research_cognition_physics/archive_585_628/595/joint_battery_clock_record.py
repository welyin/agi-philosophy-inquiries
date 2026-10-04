"""595: conserving battery dilation connected to the existing internal clock.

Original-H application is analytic. Finite matrices below are explicitly
labelled algebra fixtures for the clock/free-energy/source connection, never
claimed to be a spectrum or simulation of the original full research H0.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import persistent_prefix_history_audit as clock

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_battery_clock_record_results.json'
spec=importlib.util.spec_from_file_location('entry595',HERE/'round595_drafts/positive_battery_entry.py')
entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)


def matrix_fixture():
    # System, positive finite diagnostic battery, pointer, memory.
    dims=(2,3,2,2);size=int(np.prod(dims))
    def index(s,b,p,m):return np.ravel_multi_index((s,b,p,m),dims)
    energies=np.array([s+b for s in range(2) for b in range(3) for p in range(2) for m in range(2)],float)
    system=np.array([s for s in range(2) for b in range(3) for p in range(2) for m in range(2)],float)
    A=np.array([[.72,.11],[.11,.90]])
    small=clock.unitary(np.kron(A,clock.Y),1.)
    gate=np.eye(size,dtype=complex)
    for total in (1,2):
        for m in range(2):
            ix=[index(s,total-s,p,m) for s in range(2) for p in range(2)]
            gate[np.ix_(ix,ix)]=small
    copy=np.zeros((size,size),complex)
    for s in range(2):
        for b in range(3):
            for p in range(2):
                for m in range(2):copy[index(s,b,p,m^p),index(s,b,p,m)]=1
    gate=copy@gate
    Hf=np.diag(energies)
    assert np.linalg.norm(gate.conj().T@gate-np.eye(size))<4e-15
    assert np.linalg.norm(Hf@gate-gate@Hf)<1e-14
    rng=np.random.default_rng(595)
    payload=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2));payload/=np.linalg.norm(payload)
    battery=np.array([0.,.8,.6])
    psi=np.zeros((size,2),complex)
    for s in range(2):
        for b in range(3):psi[index(s,b,0,0)]=battery[b]*payload[s]
    # A surrogate noncommuting observable, not the physical geometry momentum.
    Pi=np.kron(clock.X,np.eye(12))*.31
    return Hf,system,gate,psi,Pi


def reduced_record(v,length):
    arr=v.reshape(length,2,3,2,2,2).transpose(1,4,5,0,2,3).reshape(8,-1)
    return arr@arr.conj().T


def joint_clock_checks():
    Hf,Es,W,psi,Pi=matrix_fixture();dim=len(Hf);length=18;J=.35;hbar=.7
    hc,D,prefixes=clock.history_matrices([W]+[np.eye(dim)]*(length-2))
    Hfree=np.kron(np.eye(length),Hf)
    H=Hfree+J*hc
    bareclock=2*np.eye(length)+np.diag(np.ones(length-1),1)+np.diag(np.ones(length-1),-1)
    conjugated=D@(Hfree+J*np.kron(bareclock,np.eye(dim)))@D.conj().T
    dictionary_error=float(np.linalg.norm(H-conjugated))
    energy_commutator=float(np.linalg.norm(H@Hfree-Hfree@H))
    assert dictionary_error<2e-13 and energy_commutator<2e-13
    eig,vec=np.linalg.eigh(H);assert eig[0]>0
    initial=np.zeros((length*dim,2),complex);initial[:dim]=psi
    coefficients=vec.conj().T@initial
    Pi_full=np.kron(np.eye(length),Pi)
    Fapp=1j/hbar*(J*hc@Pi_full-Pi_full@ (J*hc))
    deltaPi=W.conj().T@Pi@W-Pi
    input_free=float(np.vdot(initial,Hfree@initial).real)
    input_full=float(np.vdot(initial,H@initial).real)
    input_system=float(np.sum(Es[:,None]*abs(psi)**2))
    gate_system=float(np.sum(Es[:,None]*abs(W@psi)**2))
    rows=[];forces=[]
    for t in (.5,1.5,3.):
        evolved=vec@(np.exp(-1j*t*eig/hbar)[:,None]*coefficients)
        u=J*t/hbar
        amplitudes=clock.path_amplitudes(length,u)*np.exp(-2j*u)
        free=np.exp(-1j*t*np.diag(Hf)/hbar)[:,None]*psi
        expected=np.array([amplitudes[n]*(free if n==0 else W@free) for n in range(length)]).reshape(length*dim,2)
        history_error=float(np.linalg.norm(evolved-expected))
        rho=reduced_record(evolved,length)
        undone=reduced_record(free,1);done=reduced_record(W@free,1)
        p0=float(abs(amplitudes[0])**2)
        mixture=p0*undone+(1-p0)*done
        mixture_error=float(np.linalg.norm(rho-mixture))
        record_distance=clock.trace_distance(rho,done)
        output_system=float(np.sum(abs(evolved.reshape(length,dim,2))**2*Es[None,:,None]))
        balance=output_system-input_system-(1-p0)*(gate_system-input_system)
        free_error=float(np.vdot(evolved,Hfree@evolved).real-input_free)
        full_error=float(np.vdot(evolved,H@evolved).real-input_full)
        assert history_error<4e-14 and mixture_error<4e-14
        assert record_distance<=p0+2e-14
        assert abs(balance)<5e-14 and abs(free_error)<5e-14 and abs(full_error)<5e-14
        rows.append(dict(time=t,history_error=history_error,record_mixture_error=mixture_error,
                         p_unfinished=p0,reference_retaining_record_distance=record_distance,
                         system_energy_gain=output_system-input_system,
                         completed_weight_energy_residual=balance,
                         free_total_energy_error=free_error,full_autonomous_energy_error=full_error,
                         clock_truncation_bound=clock.truncation_bound(length,u)))
        direct=float(np.vdot(evolved,Fapp@evolved).real)
        pdot=-2*J/hbar*np.imag(amplitudes[0].conjugate()*amplitudes[1])
        delta=float(np.vdot(free,deltaPi@free).real)
        predicted=float(pdot*delta)
        force_error=abs(direct-predicted)
        assert force_error<3e-14
        forces.append(dict(time=t,apparatus_force=direct,completion_current_force=predicted,
                           completion_probability_derivative=float(pdot),momentum_change_in_completed_branch=delta,
                           identity_error=force_error))
    assert max(abs(r['apparatus_force']) for r in forces)>1e-5
    late=[]
    comm_norm=float(np.linalg.norm(W@Pi-Pi@W,2))
    for u in (2.,5.,10.,20.):
        a=clock.halfline_amplitudes(2,u)
        coefficient=2*J/hbar*abs(a[0]*a[1])*comm_norm
        envelope=4*J/(hbar*u*u)*comm_norm
        assert coefficient<=envelope+1e-14
        late.append(dict(scaled_time=u,force_norm_expectation_bound=float(coefficient),late_envelope=envelope))
    return (dict(dictionary_error=dictionary_error,free_energy_commutator=energy_commutator,
                 minimum_finite_fixture_eigenvalue=float(eig[0]),rows=rows,
                 original_full_H0_not_numerically_propagated=True),
            dict(rows=forces,late_bounds=late,nonzero_force_despite_total_energy_conservation=True,
                 fixture_observable_is_not_original_geometry_momentum=True))


def run():
    frozen=entry.run();assert frozen==json.loads(entry.TARGET.read_text('utf8'))
    joint,force=joint_clock_checks()
    evidence=dict(positive_energy_profile=frozen['sine_profile'],
                  energy_conserving_reference_dilation=frozen['energy_fibre'],
                  joint_clock_energy_and_record=joint,apparatus_force=force)
    names=('research_note_594.md','research_note_592.md','research_note_401.md','persistent_prefix_history_audit.py',
           'round595_drafts/positive_battery_entry.md','round595_drafts/positive_battery_entry.py',
           'round595_drafts/positive_battery_entry_results.json','round595_drafts/positive_battery_entry_checks.json')
    return dict(round=595,tests_run=4,failures=0,errors=0,evidence=evidence,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
                scope='Conditional original-H0 dilation with a positive continuous battery, explicit coherence resources, neutral record registers and a half-infinite history clock. One fixed positive autonomous H per accuracy gives late record/old-source approximation and exact total energy conservation. Gate geometry force is included and bounded, not a total gravitational stress tensor. Spectral/global couplings, extra hardware and preparation are inputs. Matrix fixtures are not original H0 spectra. No local matter realization, GR constraints or continuum completion.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==r
    print(json.dumps(r,ensure_ascii=False,indent=2))
