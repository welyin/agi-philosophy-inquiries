"""741: complete same-record source hierarchy and reference-family lift.

Continuum first-order construction is analytic in note741. Finite full-species
matrices test the object mapping, not continuum stress or Einstein evolution.
"""
import argparse
import hashlib
import json
from functools import lru_cache
from pathlib import Path
import sys
import numpy as np
import joint_dynamic_continuum_reference as ref
import joint_relative_source_development as prior
import joint_local_source_normalization as local
import joint_retarded_reference_response as response
HERE=Path(__file__).resolve().parent;TARGET=HERE/'joint_absolute_source_development_results.json'
sys.path.insert(0,str(HERE/'round741_drafts'))
import absolute_source_order_entry as entry

def matrices(t,k,parameter,choice):
    base,_=ref.matrices(t,.08,k)
    _,phi,_,_,_=ref.collar(t,.08)
    M=ref.old.bdg(*ref.old.matter.mass_matrices(phi))
    extra=choice*response.entry.pulse(t)
    ramp=ref.smoothstep((t+.16)/.08)
    h=ramp+extra
    K=np.exp(-parameter*h)*(base-M)
    phip=np.exp(parameter*extra)*phi
    MM=ref.old.bdg(*ref.old.matter.mass_matrices(phip))
    D=-h*K+extra*2/ref.old.matter.original.F(phip)*MM
    return K+MM,D,-K,-h*(-K)

@lru_cache(maxsize=18)
def lift(parameter,choice,steps=96):
    start=-.24;dt=-start/steps;k=np.array([.31,-.27,.19])
    ps=[];dps=[];ends=[];source=[];dsource=[]
    for sign in (1,-1):
        B,*_=matrices(start,sign*k,parameter,choice)
        P=ref.projector(B);dP=np.zeros_like(P)
        for j in range(steps):
            B,D,*_=matrices(start+(j+.5)*dt,sign*k,parameter,choice)
            U,dU=ref.step(B,D,dt);before=P
            dP=dU@before@U.conj().T+U@dP@U.conj().T+U@before@dU.conj().T
            P=U@before@U.conj().T
        B,D,G,dG=matrices(0.,sign*k,parameter,choice)
        ps.append(P);dps.append(dP);ends.append(B);source.append(G);dsource.append(dG)
    P,dP,B,G,dG=(prior.assemble(xs) for xs in (ps,dps,ends,source,dsource))
    post=P+response.family.record_change(P)
    dpost=dP+response.family.record_change(dP)
    _,phi,_,_,_=ref.collar(0.,.08)
    J=local.jacobian(phi)
    N=[prior.assemble([n,n]) for n in local.basis()]
    Q=np.r_[prior.source(post,G),J.T@np.array([prior.source(post,n) for n in N])]
    dQ=np.r_[prior.source(dpost,G)+prior.source(post,dG),
             J.T@np.array([prior.source(dpost,n) for n in N])]
    C=np.block([[np.zeros((64,64)),np.eye(64)],[np.eye(64),np.zeros((64,64))]])
    ev=np.linalg.eigvalsh(post)
    reality=float(np.max(abs(C@post.conj()@C-(np.eye(128)-post))))
    assert reality<2e-12 and ev.min()>-2e-12 and ev.max()<1+2e-12
    return Q,dQ,B,dict(post_record_min_eigenvalue=float(ev.min()),
        post_record_max_eigenvalue=float(ev.max()),post_record_reality_error=reality)

def reference_extension_check():
    a,da,B0,info=lift(0.,0);b,db,BB,_=lift(0.,1)
    assert np.max(abs(a-b))<1e-13 and np.max(abs(B0-BB))==0
    inherited=entry.actual(0.,0.)[0]
    assert np.max(abs(a-inherited))<3e-12
    delta=db-da
    assert np.max(abs(delta))>1e-3
    rows=[]
    for eps in (.02,.01,.005):
        ap,_,Ba,ia=lift(eps,0);bp,_,Bb,ib=lift(eps,1)
        am,*_=lift(-eps,0);bm,*_=lift(-eps,1)
        error=float(np.max(abs(((bp-ap)-(bm-am))/(2*eps)-delta)))
        weighted=eps*(bp-ap)
        assert np.max(abs(Ba-Bb))==0
        assert np.max(abs(Ba-B0))>1e-3
        rows.append(dict(parameter=eps,endpoint_background_difference_between_lifts=0.,
            common_endpoint_change_from_original=float(np.max(abs(Ba-B0))),
            complete_source_history_difference=(bp-ap).tolist(),
            loop_weighted_history_difference=weighted.tolist(),
            normalized_second_order_coefficient=(weighted/eps**2).tolist(),
            centered_derivative_error=error,reference_lift_a=ia,reference_lift_b=ib))
    ratios=[rows[i]['centered_derivative_error']/rows[i+1]['centered_derivative_error'] for i in range(2)]
    assert all(3.6<r<4.4 for r in ratios),ratios
    return dict(original_dimension=128,zero_order_source=a.tolist(),
        first_order_source_derivative_a=da.tolist(),first_order_source_derivative_b=db.tolist(),
        second_order_backreaction_extension_ambiguity=delta.tolist(),rows=rows,
        derivative_error_ratios=ratios,common_preparation_and_record=True,
        same_entire_endpoint_neighborhood_not_only_same_value=True,
        no_claim_that_two_finite_parameter_states_are_identical=True,
        no_claim_of_continuum_stress_or_a_constraint_solved_numerical_metric=True)

def run():
    first=entry.run();assert first==json.loads(entry.TARGET.read_text('utf8'))
    second=reference_extension_check()
    names=('research_note_573.md','research_note_731.md','research_note_732.md','research_note_734.md',
        'research_note_735.md','research_note_740.md','joint_source_constraint_response.py',
        'joint_retarded_reference_response.py','joint_local_source_normalization.py',
        'round741_drafts/absolute_source_order_entry.py','round741_drafts/absolute_source_order_entry_results.json')
    return dict(round=741,tests_run=2,failures=0,errors=0,
        checks=['same_record_absolute_source_order_check','common_initial_reference_extension_check'],
        results=dict(source_order=first,reference_extension=second),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
        scope='Under the declared one-fermion-loop conserving normalization, original smooth classical development and chosen reference/record, a full absolute first-order backreaction and smooth compatible state family give an O(lambda^2) complete-equation and constraint residual. Fixed zeroth-order preparation makes the first coefficient independent of first-order auxiliary extensions, while second order can differ. No exact semiclassical solution, physical lambda=1 error bound, autonomous instrument, full quantum Gauss or UV completion is claimed.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:r[k] for k in ('round','tests_run','failures','errors')}))
