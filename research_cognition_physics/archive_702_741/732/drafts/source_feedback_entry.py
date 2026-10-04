"""732 entry: actual source response versus a prepared-energy Hessian.

Original730 complete local matrix calibration; not a continuum stress solver.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_dynamic_continuum_reference as old
TARGET=HERE/'source_feedback_entry_results.json'


def matrices(t,k,gamma,beta):
    g,phi,a,a0,eta=old.collar(t,.08)
    original_mass=old.old.bdg(*old.old.matter.mass_matrices(phi))
    base,_=old.matrices(t,.08,k)
    kinetic=np.exp(-gamma*eta)*(base-original_mass)
    phi=np.exp(beta*eta)*phi
    F=old.old.matter.original.F(phi)
    mass=old.old.bdg(*old.old.matter.mass_matrices(phi))
    return kinetic+mass,(-eta*kinetic,eta*2/F*mass)


def flow(gamma=0.,beta=0.,steps=96,derivatives=True):
    start=-.24;dt=-start/steps;k=np.array([.31,-.27,.19])
    ps=[];dps=[[],[]];bs=[];sources=[[],[]]
    for sign in (1,-1):
        B,_=matrices(start,sign*k,gamma,beta);P=old.projector(B)
        dP=[np.zeros_like(P),np.zeros_like(P)]
        for j in range(steps):
            B,Ds=matrices(start+(j+.5)*dt,sign*k,gamma,beta)
            U,dU=old.step(B,Ds[0],dt);before=P
            if derivatives:
                _,dV=old.step(B,Ds[1],dt)
                dP=[dQ@before@U.conj().T+U@der@U.conj().T+U@before@dQ.conj().T
                    for dQ,der in zip((dU,dV),dP)]
            P=U@before@U.conj().T
        B,Ds=matrices(0,sign*k,gamma,beta)
        ps.append(P);bs.append(B)
        for i in (0,1):dps[i].append(dP[i]);sources[i].append(Ds[i])
    order=np.r_[np.arange(32),np.arange(64,96),np.arange(96,128),np.arange(32,64)]
    def assemble(xs):return old.old.block(xs)[np.ix_(order,order)]
    return assemble(ps),[assemble(xs) for xs in dps],assemble(bs),[assemble(xs) for xs in sources]


def record_change(P):
    e=np.zeros(128,complex);e[30]=e[62]=1/np.sqrt(2)
    conjugation=np.block([[np.zeros((64,64)),np.eye(64)],[np.eye(64),np.zeros((64,64))]])
    ce=conjugation@e.conj();Q=np.outer(e,e.conj())+np.outer(ce,ce.conj())
    R=np.eye(128)-2*Q
    return (R@P@R-P)/2


def evaluate(gamma=0.,beta=0.,steps=96,derivatives=True):
    P,dPs,B,Ds=flow(gamma,beta,steps,derivatives)
    change=record_change(P)
    source=np.array([np.trace(D@change).real/2 for D in Ds])
    energy=float(np.trace(B@change).real/2)
    # Cross operator derivatives vanish: metric kinetic scaling and radial mass
    # scaling are independent coordinates of the stated original local family.
    cross=np.array([[np.trace(D@record_change(dp)).real/2 for dp in dPs] for D in Ds])
    state_response=np.array([np.trace(B@record_change(dp)).real/2 for dp in dPs])
    return source,energy,cross,state_response


def run():
    source,energy,cross,state=evaluate()
    reference=old.record_source_check()
    assert abs(source[0]-reference['direct_fixed_state_source'])<1e-12
    assert abs(state[0]-reference['preparation_response'])<1e-12
    h=2e-5;fd=np.zeros((2,2));energy_fd=np.zeros(2)
    for j in (0,1):
        plus=np.zeros(2);minus=np.zeros(2);plus[j]=h;minus[j]=-h
        sp,ep,_,_=evaluate(*plus,derivatives=False)
        sm,em,_,_=evaluate(*minus,derivatives=False)
        fd[:,j]=(sp-sm)/(2*h);energy_fd[j]=(ep-em)/(2*h)
    cross_error=max(abs(fd[0,1]-cross[0,1]),abs(fd[1,0]-cross[1,0]))
    total_error=float(np.max(abs(energy_fd-source-state)))
    fine,_,fine_cross,_=evaluate(steps=192)
    curl=float(cross[0,1]-cross[1,0]);fine_curl=float(fine_cross[0,1]-fine_cross[1,0])
    assert cross_error<2e-8 and total_error<2e-8
    assert abs(curl)>1e-5 and abs(curl-fine_curl)<abs(curl)*.01
    deps=('research_note_601.md','research_note_730.md','research_note_731.md',
          'joint_dynamic_continuum_reference.py','joint_dynamic_continuum_reference_results.json')
    return dict(entry_round=732,new_formal_round=False,physical_modes=64,Nambu_dimension=128,
        nonselective_record_energy=energy,direct_sources=source.tolist(),preparation_response=state.tolist(),
        offdiag_source_response=dict(metric_to_radial=float(cross[1,0]),radial_to_metric=float(cross[0,1])),
        antisymmetric_source_response=curl,antisymmetric_source_response_192_steps=fine_curl,
        offdiag_finite_difference_error=float(cross_error),prepared_energy_total_derivative_error=total_error,
        source_96_vs_192_error=float(np.max(abs(source-fine))),
        continuum_source_regularization_and_self_consistency_not_checked=True,
        dependencies={name:hashlib.sha256((ARCHIVE/name).read_bytes()).hexdigest() for name in deps},
        scope='Original730 transported full local matter and same record. Metric kinetic/radial mass direct forces have nonsymmetric cross response; they cannot be replaced by the Hessian of total prepared record energy. Finite local symbol only, not a no-go for causal feedback or a semiclassical spacetime solution.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
