"""734 entry: common-past metric pulse versus endpoint record-preserving reset.

Original complete local matter matrices; no continuum stress prediction.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE));sys.path.insert(0,str(ARCHIVE/'round732_drafts'))
import joint_dynamic_continuum_reference as reference
import joint_relative_source_development as prior
import source_feedback_entry as family
TARGET=HERE/'causal_reference_entry_results.json'

def pulse(t):
    s=(t+.11)/.08
    return float(np.exp(4-1/(s*(1-s)))) if 0<s<1 else 0.

def matrices(t,k,gamma):
    B,_=reference.matrices(t,.08,k)
    _,phi,_,_,_=reference.collar(t,.08)
    M=reference.old.bdg(*reference.old.matter.mass_matrices(phi))
    w=pulse(t);K=np.exp(-gamma*w)*(B-M)
    return K+M,-w*K

def flow(gamma,derivative=True):
    start=-.24;steps=96;dt=-start/steps;k=np.array([.31,-.27,.19])
    ps=[];ds=[];us=[];bs=[];pre=0.
    for sign in (1,-1):
        B,_=matrices(start,sign*k,gamma);P=reference.projector(B)
        dP=np.zeros_like(P);Ut=np.eye(64,dtype=complex)
        for j in range(steps):
            t=start+(j+.5)*dt;B,D=matrices(t,sign*k,gamma)
            U,dU=reference.step(B,D,dt);before=P
            if derivative:
                dP=dU@before@U.conj().T+U@dP@U.conj().T+U@before@dU.conj().T
            P=U@before@U.conj().T;Ut=U@Ut
            if t<-.11:pre=max(pre,float(np.max(abs(dP))))
        ps.append(P);ds.append(dP);us.append(Ut);bs.append(matrices(0,sign*k,gamma)[0])
    return *(prior.assemble(x) for x in (ps,ds,us,bs)),pre

def block(P0,C):
    f=np.zeros(128,complex);f[30]=f[62]=1/np.sqrt(2);cf=C@f.conj()
    U,s,_=np.linalg.svd(np.column_stack((f,cf,P0@f,P0@cf)),full_matrices=False)
    basis=U[:,s>1e-12]
    return basis@basis.conj().T

def run():
    P0,dP0,U0,B0,pre=flow(0.)
    original,_,_,_,C=prior.initial_modes();assert np.max(abs(P0-original))<1e-12
    K=block(P0,C);Z=np.eye(128)-K;delta0=family.record_change(P0)
    gamma=.2;P,dP,U,B,_=flow(gamma)
    patched=K@P0@K+Z@P@Z
    assert np.max(abs(B-B0))<1e-14 and pre==0.
    old_delta=family.record_change(P);new_delta=family.record_change(patched)
    record_derivative=family.record_change(dP0)
    assert np.linalg.norm(old_delta-delta0)>1e-6
    assert np.linalg.norm(new_delta-delta0)<1e-12
    # Pure common-past evolution remains pure. The reset generally does not.
    purity=float(np.trace(patched@(np.eye(128)-patched)).real)
    exact_purity=float(np.linalg.norm(K@P@Z,'fro')**2)
    assert purity>1e-9 and abs(purity-exact_purity)<2e-11
    reset=patched-P;past_reset=U.conj().T@reset@U
    assert np.linalg.norm(past_reset-reset)>1e-9
    assert abs(np.linalg.norm(past_reset,'fro')-np.linalg.norm(reset,'fro'))<1e-12
    _,_,_,Ds=family.flow(0.,0.,derivatives=False);G=Ds[0]
    direct_response=prior.source(dP0,G)
    reset_response=prior.source(Z@dP0@Z,G)
    state_correction=prior.source(Z@dP0@Z-dP0,G)
    h=2e-5;Pp,_,_,_,_=flow(h,False);Pm,_,_,_,_=flow(-h,False)
    fd=(Pp-Pm)/(2*h);fd_error=float(np.max(abs(fd-dP0)))
    assert fd_error<3e-9
    assert abs(reset_response-direct_response-state_correction)<1e-14
    assert np.linalg.norm(record_derivative)>1e-6
    assert abs(state_correction)>1e-9
    deps=('research_note_733.md','joint_reference_polarization_boundary.py',
          'joint_reference_polarization_boundary_results.json','joint_dynamic_continuum_reference.py',
          'joint_relative_source_development.py','round732_drafts/source_feedback_entry.py')
    return dict(entry_round=734,new_formal_round=False,physical_modes=64,Nambu_dimension=128,
        pulse_support=[-.11,-.03],amplitude=gamma,endpoint_H_change=float(np.max(abs(B-B0))),
        before_pulse_response=pre,base_original_reference_error=float(np.max(abs(P0-original))),
        common_past_record_change=float(np.linalg.norm(old_delta-delta0,'fro')),
        pinned_record_error=float(np.linalg.norm(new_delta-delta0,'fro')),
        common_past_purity_error=float(np.linalg.norm(P@P-P,'fro')),
        reset_mixedness=purity,exact_cross_block_mixedness=exact_purity,
        required_past_state_correction_norm=float(np.linalg.norm(past_reset,'fro')),
        common_past_geometric_source_response=direct_response,
        pinned_geometric_source_response=reset_response,
        additional_preparation_source_response=state_correction,
        actual_record_response_norm=float(np.linalg.norm(record_derivative,'fro')),
        finite_difference_response_error=fd_error,
        dependencies={name:hashlib.sha256((ARCHIVE/name).read_bytes()).hexdigest() for name in deps},
        scope='Original common-past pure reference under a smooth past metric pulse has a record and source response even at unchanged endpoint Hamiltonian. Endpoint pinning removes record response but changes the preparation and generally mixes the state. Not a no-go for apparatus/environment-assisted preparation, not a renormalized continuum stress or full semiclassical solution.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
