"""734: original joint-background retarded response and contact checks.

Finite local complete-species calibration only. Continuum point-splitting,
its differentiability scope and Ward conditions are stated in note734.
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
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'round734_drafts'))
import causal_reference_entry as entry
sys.path.insert(0,str(HERE/'round732_drafts'))
import source_feedback_entry as family
TARGET=HERE/'joint_retarded_reference_response_results.json'

def maxabs(x):return float(np.max(abs(x)))

def local(t,k,gamma=0.):
    g,phi,a,a0,_=ref.collar(t,.08);w=entry.pulse(t)
    base,_=ref.matrices(t,.08,k)
    M0=ref.old.bdg(*ref.old.matter.mass_matrices(phi))
    K=np.exp(-gamma*w)*(base-M0)
    phi=np.exp(gamma*w)*phi
    M=ref.old.bdg(*ref.old.matter.mass_matrices(phi))
    F=ref.old.matter.original.F(phi)
    return K+M,-w*K+w*2/F*M,-K,-w*(-K)

def matrices(t,gamma=0.):
    k=np.array([.31,-.27,.19]);xs=[local(t,s*k,gamma) for s in (1,-1)]
    return tuple(prior.assemble([x[j] for x in xs]) for j in range(4))

@lru_cache(maxsize=12)
def flow(end=0.,gamma=0.,steps=96,history=False):
    start=-.24;dt=(end-start)/steps
    P=ref.projector(matrices(start,gamma)[0]);dP=np.zeros_like(P);saved=[]
    for j in range(steps):
        t=start+(j+.5)*dt;B,D,_,_=matrices(t,gamma)
        U,dU=ref.step(B,D,dt);old=P
        dP=dU@old@U.conj().T+U@dP@U.conj().T+U@old@dU.conj().T
        P=U@old@U.conj().T
        if history:saved.append((t,U,dU))
    return P,dP,saved

def retarded_kernel_check():
    P,dP,saved=flow(history=True)
    A=np.zeros_like(P);tail=np.eye(128,dtype=complex);before=0.
    for t,U,dU in reversed(saved):
        localA=1j*dU@U.conj().T
        assert maxabs(localA-localA.conj().T)<2e-14
        if t<-.11:before=max(before,maxabs(localA))
        A+=tail@localA@tail.conj().T
        tail=tail@U
    response=-1j*(A@P-P@A)
    _,_,G,_=matrices(0.)
    source=prior.source(dP,G)
    commutator=prior.source(P,-1j*(G@A-A@G))
    record=prior.source(family.record_change(dP),G)
    record_commutator=prior.source(P,-1j*(family.record_change(G)@A-A@family.record_change(G)))
    err=max(maxabs(response-dP),abs(source-commutator),abs(record-record_commutator))
    assert err<2e-12 and before==0 and abs(source)>1e-4 and abs(record)>1e-6
    return dict(original_dimension=128,propagated_response_vs_retarded_kernel_error=err,
        pre_support_kernel=before,geometric_reference_response=source,
        actual_record_geometric_response=record,
        independently_assembled_commutator_response=commutator,
        kernel_is_discrete_exact_Frechet_response_not_a_continuum_stress=True)

def gauge_sources(t,gamma,kind,index):
    g,phi,a,a0,_=ref.collar(t,.08);w=entry.pulse(t)
    g=np.exp(2*gamma*w)*g;phi=np.exp(gamma*w)*phi
    dphi=np.zeros(5);da=np.zeros_like(a)
    X=phi[:2]+1j*phi[2:4]
    if kind=='weak':
        change=1j*(ref.SIG[index]/2)@X;dphi[:4]=np.r_[change.real,change.imag]
        da=np.cross(a,np.eye(3)[index])
    elif kind=='circle':
        change=3j*X;dphi[:4]=np.r_[change.real,change.imag]
    bg=(g,phi,a,a0);zero=[np.zeros_like(x) for x in bg]
    _,S=prior.matrix_and_source(*bg,zero[0],dphi,zero[2],zero[3])
    _,V=prior.matrix_and_source(*bg,zero[0],zero[1],da,zero[3])
    return S,V

def linearized_exchange_check():
    t=-.07;P,dP,_=flow(t);B,D,_,_=matrices(t)
    dotdP=-1j*(D@P-P@D+B@dP-dP@B)
    rows=[];maximum=0.;operator_error=0.;h=2e-5
    for kind,index,H in prior.physical_generators():
        Q=prior.assemble([ref.old.bdg(H,np.zeros_like(H))]*2)
        S,V=gauge_sources(t,0.,kind,index);K=S+V
        Sp,Vp=gauge_sources(t,h,kind,index);Sm,Vm=gauge_sources(t,-h,kind,index)
        dK=(Sp+Vp-Sm-Vm)/(2*h)
        operator_error=max(operator_error,maxabs(K-1j*(Q@B-B@Q)),maxabs(dK-1j*(Q@D-D@Q)))
        rate=prior.source(dotdP,Q)
        state=prior.source(dP,K);contact=prior.source(P,dK)
        maximum=max(maximum,abs(rate+state+contact))
        rows.append(dict(kind=kind,index=index,charge_rate_response=rate,
                         state_response=state,background_source_contact=contact,
                         missing_contact_defect=rate+state))
    assert operator_error<3e-8 and maximum<3e-9
    assert max(abs(r['background_source_contact']) for r in rows)>1e-5
    return dict(generators=12,full_matrix_linearized_identity_error=operator_error,
        complete_exchange_balance_error=maximum,rows=rows,
        numerical_check_not_an_anomaly_cancellation_proof=True)

def same_record_source_contact_check():
    t=-.07;P,dP,_=flow(t);B,D,G,dG=matrices(t)
    post=P+family.record_change(P);dpost=dP+family.record_change(dP)
    state=prior.source(dpost,G);contact=prior.source(post,dG);total=state+contact
    h=2e-5;forces=[]
    for gamma in (h,-h):
        Q,_,_=flow(t,gamma);_,_,GG,_=matrices(t,gamma)
        forces.append(prior.source(Q+family.record_change(Q),GG))
    fd=(forces[0]-forces[1])/(2*h)
    assert abs(fd-total)<2e-7 and abs(contact)>1
    relative_state=prior.source(family.record_change(dP),G)
    relative_contact=prior.source(family.record_change(P),dG)
    assert abs(relative_state)>1e-6 and abs(relative_contact)>1e-3
    return dict(actual_post_record_state_response=state,actual_operator_contact=contact,
        total_source_response=total,independent_difference=fd,error=abs(fd-total),
        record_only_state_response=relative_state,record_only_operator_contact=relative_contact,
        finite_matrix_contact_not_Hadamard_counterterm=True)

def run():
    checks=('retarded_kernel_check','linearized_exchange_check','same_record_source_contact_check')
    results={name:globals()[name]() for name in checks}
    deps=('research_note_620.md','research_note_628.md','research_note_632.md','research_note_733.md',
          'joint_dynamic_continuum_reference.py','joint_relative_source_development.py',
          'joint_reference_polarization_boundary.py','joint_reference_polarization_boundary_results.json',
          'round734_drafts/causal_reference_entry.py','round734_drafts/causal_reference_entry_results.json')
    return dict(round=734,tests_run=3,failures=0,errors=0,checks=list(checks),results=results,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Same-past retarded reference response on the declared smooth continuum background family requires simultaneous state, source-operator, parametrix and local-term variation. The joint differentiated Ward identity has background contacts; zero anomaly requires the declared compatible renormalization contract. Original complete finite matrices verify causal kernel and nonzero contacts, not renormalized continuum stress or nonlinear semiclassical existence.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:r[k] for k in ('round','tests_run','failures','errors')}))
