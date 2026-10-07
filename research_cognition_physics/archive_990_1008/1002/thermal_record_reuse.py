"""Reuse the actual round-980 thermal state for a bounded population record.

The extra finite auxiliary supplies energy AND takes records/entropy.
Its pure initial energy state and the engineered interaction are inputs.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
TARGET=HERE/'thermal_record_reuse_results.json'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text('utf-8-sig'))
def entropy(p):
    p=np.asarray(p);p=p[p>1e-14]
    return float(-np.dot(p,np.log(p)))
def partial(rho):
    a=rho.reshape(2,4,5,2,4,5)
    return (np.einsum('ijkljk->il',a),np.einsum('ijkimk->jm',a),
            np.einsum('ijkijn->kn',a))


def calculate():
    old=read(STAGE/'980/finite_thermal_records_results.json')
    e=np.array(old['parameters']['energies']);p=np.array(old['parameters']['thermal_populations'])
    beta=old['parameters']['beta'];rev=np.arange(3,-1,-1)
    rebuilt=np.exp(-beta*(e-e[0]));rebuilt/=sum(rebuilt)
    assert np.max(abs(p-rebuilt))<1e-14
    base=2.;aux=np.r_[base,base+e-e[rev]];g=.01
    assert min(aux)>0 and len(np.unique(aux))==5
    label=np.array([0.,.2])
    energies=(label[:,None,None]+e[None,:,None]+aux[None,None,:]).ravel()
    H0=np.diag(energies);V=np.zeros((40,40))
    def ix(x,i,b):return (x*4+i)*5+b
    for i in range(4):
        left=ix(1,i,0);right=ix(1,rev[i],i+1)
        V[left,right]=V[right,left]=1
        assert abs(energies[left]-energies[right])<1e-15
    H=H0+g*V
    evals,vecs=np.linalg.eigh(H)
    commute=float(np.linalg.norm(H0@V-V@H0,2));assert commute<1e-14
    assert evals[0]>0
    Plow=np.diag(np.tile(np.repeat([1,1,0,0],5),2))
    active_comm=float(np.linalg.norm(Plow@H-H@Plow,2));assert active_comm>0
    t0=math.pi/(2*g);halfwidth=.1/g
    def evolve(rho,t):
        U=(vecs*np.exp(-1j*evals*t))@vecs.conj().T
        return U@rho@U.conj().T
    def initial(x):
        rho=np.zeros((40,40),complex)
        for i in range(4):rho[ix(x,i,0),ix(x,i,0)]=p[i]
        return rho
    low=float(p[0]+p[1]);q=1-low;contrast=2*low-1
    samples=[];endpoint=[]
    for t in (t0-halfwidth,t0,t0+halfwidth):
        probabilities=[]
        for x in (0,1):
            state=evolve(initial(x),t)
            probabilities.append(float(np.trace(Plow@state).real))
            assert abs(np.trace(H@state)-np.trace(H@initial(x)))<1e-12
            assert abs(np.trace(g*V@state))<1e-12
            if t==t0:endpoint.append(state)
        theory_p1=low-contrast*math.sin(g*t)**2
        assert abs(probabilities[0]-low)<2e-12
        assert abs(probabilities[1]-theory_p1)<2e-12
        avg=(1-probabilities[0]+probabilities[1])/2
        samples.append(dict(time=t,p_low_given_0=probabilities[0],p_low_given_1=probabilities[1],
                            equal_prior_error=avg))
    source,memory,battery=partial(endpoint[1])
    assert np.max(abs(np.diag(memory).real-p[rev]))<2e-12
    assert np.max(abs(np.diag(battery).real-np.r_[0.,p]))<2e-12
    transfer=float(np.dot(e,p[rev]-p));db=float(np.dot(aux,np.diag(battery).real)-base)
    relative=float(np.dot(p[rev],np.log(p[rev]/p)))
    assert transfer>0 and abs(transfer+db)<1e-12
    assert abs(relative-beta*transfer)<1e-13
    # Auxiliary alone has disjoint label supports at the selected endpoint.
    b0=partial(endpoint[0])[2]
    td_B=float(np.sum(abs(np.linalg.eigvalsh(battery-b0)))/2)
    assert abs(td_B-1)<1e-12
    # Source superposition is NOT preserved locally. Its global coherence
    # stays in the joint process; no assumption of classical input collapse.
    pure_source=np.zeros((40,40),complex)
    for i in range(4):
        for x in (0,1):
            for y in (0,1):pure_source[ix(x,i,0),ix(y,i,0)]=p[i]/2
    S_out=partial(evolve(pure_source,t0))[0]
    assert np.max(abs(S_out-np.eye(2)/2))<1e-12
    # Reuse only the existing rigorous reset bound, not a new bath simulation.
    reset=float(13/14000+1e-4)
    assert any(isinstance(value,dict) and value.get('full_charge_reset_bound')==reset for value in old.values())
    window=q+contrast*math.sin(.1)**2/2
    propagated=window+reset
    assert propagated<.193032
    # Inherit the certified thermal population intervals, not rounded spectrum
    # decimals, for the reported finite-window bound. sin(.1)^2 <= .01.
    intervals=old['eigensystem_certificate']['coarse_analytic_certificate']['thermal_intervals']
    assert sum(F(str(intervals[i][1])) for i in (2,3))<F(188903,1000000)
    certified=F(99,100)*F(188903,1000000)+F(1,200)+F(13,14000)+F(1,10000)
    assert certified<F(193043,1000000)
    refs=[Path(__file__),HERE/'drafts/STATUS.md',HERE/'drafts/adoption_decision.md',
          STAGE/'980/finite_thermal_records_results.json',STAGE/'research_note_980.md',
          STAGE/'1001/research_round_1001_checks.json']
    return dict(round=1002,date='2026-10-07',all_scientific_checks_passed=True,
        kind='finite_thermal_population_rewrite_with_energy_and_record_auxiliary',
        parameters=dict(energies=e.tolist(),thermal_populations=p.tolist(),beta=beta,
            eta=.05,auxiliary_energies=aux.tolist(),g=g,source_energies=label.tolist()),
        finite_dimension=40,H_min=float(evals[0]),energy_commutator_norm=commute,
        active_read_commutator_norm=active_comm,peak_time=t0,window_halfwidth=halfwidth,
        samples=samples,population_contrast=contrast,peak_equal_prior_error=q,
        window_equal_prior_error_bound=window,old_reset_error_bound=reset,
        after_old_reset_equal_prior_error_bound=propagated,
        certified_after_reset_error_bound=dict(exact=str(certified),value=float(certified)),
        conditional_1_memory_energy_gain=transfer,conditional_1_auxiliary_energy_change=db,
        conditional_1_memory_relative_entropy=relative,
        conditional_1_auxiliary_entropy=entropy(np.diag(battery).real),
        conditional_1_memory_auxiliary_mutual_information=entropy(p),
        auxiliary_label_trace_distance=td_B,unknown_source_dephased_locally=True,
        round980_state_reused=True,auxiliary_initial_pure_energy_state_is_input=True,
        engineered_energy_matched_interaction_is_input=True,
        original_charge_interaction_implements_writer=False,
        all_resources_prepared_from_thermal_states=False,
        population_effect_is_qnd_under_full_active_H=False,
        infinite_retention_certified=False,cyclic_auxiliary_restoration_certified=False,
        physical_readout_instrument_constructed=False,joint_clock_communication_lifecycle_certified=False,
        full_goal_completed=False,
        source_hashes={str(pth.relative_to(ROOT)):sha(pth) for pth in refs})


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-10,abs_tol=2e-12),(a,b)
    else:assert a==b,(a,b)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();out=calculate()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as dest:
            json.dump(out,dest,ensure_ascii=False,indent=2);dest.write('\n')
    else:compare(out,read(TARGET))
    print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
