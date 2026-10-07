"""811: exact smooth comparison connection, source family, and state transport.

The ten-dimensional BV model calibrates identities, not the original continuous
QME coefficients, physical elapsed-time dynamics, or a quantum apparatus.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'800'))
import physical_time_slice_transport as old
p=old.p
C=p.C;ZI=p.Z;I=p.I
TARGET=HERE/'smooth_slice_connection_results.json'
ZERO=(0,0,0,0)  # epsilon,a,b,source
ONE={ZERO:I}
def clean(s):return {k:v for k,v in s.items() if any(v.flat)}
def add(*ss):
    out={}
    for s in ss:
        for k,v in s.items():out[k]=out.get(k,ZI)+v
    return clean(out)
def sc(s,c):return clean({k:p.scale(v,c) for k,v in s.items()})
def mul(s,t):
    out={}
    for k,v in s.items():
        for l,w in t.items():
            n=tuple(x+y for x,y in zip(k,l))
            if n[0]<=3 and n[3]<=2:out[n]=out.get(n,ZI)+v@w
    return clean(out)
def adj(s):return clean({k:p.sharp(v) for k,v in s.items()})
def diff(s,n):
    out={}
    for k,v in s.items():
        if k[n]:
            l=list(k);l[n]-=1;l=tuple(l)
            out[l]=out.get(l,ZI)+p.scale(v,k[n])
    return clean(out)
def exp(s):
    out=ONE;power=ONE
    for n in range(1,7):
        power=mul(power,s)
        out=add(out,sc(power,Q(1,factorial(n))))
    return out
def equal(s,t):return not add(s,sc(t,-1))
def comm(s,t):return add(mul(s,t),sc(mul(t,s),-1))
def ad(u,s):return mul(mul(u,s),adj(u))
def at(s,a,b):
    out={}
    for k,v in s.items():
        l=(k[0],0,0,k[3])
        out[l]=out.get(l,ZI)+p.scale(v,a**k[1]*b**k[2])
    return clean(out)
def expectation(rho,s):
    out={}
    for k,v in s.items():
        x=np.trace(rho@old.physical(v))
        if x:out[k]=x
    return out
def dump(s):return {str(k):str(v) for k,v in sorted(s.items())}
def first(s):
    k=min(s);i,j=next((i,j) for i in range(10) for j in range(10) if s[k][i,j])
    return dict(powers=list(k),row=i,column=j,value=str(s[k][i,j]))

def run():
    hx=old.representative(p.base.X,Q(1,5))
    hz=old.representative(p.base.Z,Q(1,7))
    u=mul(exp({(1,1,0,0):p.scale(hx,C(0,1))}),
          exp({(1,0,1,0):p.scale(hz,C(0,1))}))
    assert equal(mul(u,adj(u)),ONE)
    ga=sc(mul(diff(u,1),adj(u)),C(0,-1))
    gb=sc(mul(diff(u,2),adj(u)),C(0,-1))
    curl=add(diff(gb,1),sc(diff(ga,2),-1))
    curvature=add(curl,sc(comm(ga,gb),C(0,-1)))
    assert not curvature and curl
    for s in (u,ga,gb):
        assert all(not any(p.delta(v).flat) for v in s.values())
    assert equal(ga,adj(ga)) and equal(gb,adj(gb))
    a0={ZERO:old.lifted(p.base.Z)}
    sources=add(a0,{(0,0,0,1):old.lifted(p.base.X),
                    (2,0,0,2):p.scale(old.lifted(p.base.Y),Q(1,13))})
    carried=ad(u,sources)
    assert equal(diff(carried,3),ad(u,diff(sources,3)))
    assert equal(diff(carried,1),sc(comm(ga,carried),C(0,1)))
    assert equal(diff(carried,2),sc(comm(gb,carried),C(0,1)))
    # Removing old quantum source contact is detected by its second derivative.
    bare=add(a0,{(0,0,0,1):old.lifted(p.base.X)})
    lost=add(diff(diff(carried,3),3),
             sc(diff(diff(ad(u,bare),3),3),-1))
    assert lost
    frames=[at(u,a,b) for a,b in ((Q(1,3),Q(2,5)),(Q(2,3),Q(-1,5)),(Q(4,5),Q(1,7)))]
    v10=mul(frames[1],adj(frames[0]))
    v21=mul(frames[2],adj(frames[1]))
    v20=mul(frames[2],adj(frames[0]))
    assert equal(mul(v21,v10),v20)
    rho=p.scale(p.base.I+p.scale(p.base.X,Q(1,3))+p.scale(p.base.Y,Q(1,7))
                +p.scale(p.base.Z,Q(1,5)),Q(1,2))
    assert Q(1,9)+Q(1,49)+Q(1,25)<1
    localrho=ad(u,{ZERO:old.lifted(rho)})
    transformed=ad(u,a0)
    scalar=lambda s:{k:np.trace(old.physical(v)) for k,v in s.items()
                     if np.trace(old.physical(v))}
    joint=scalar(mul(localrho,transformed))
    original=expectation(rho,a0)
    assert joint==original
    cancellation=add(mul(diff(localrho,1),transformed),
                     mul(localrho,diff(transformed,1)))
    assert not scalar(cancellation)
    frozen_change=expectation(rho,diff(transformed,1))
    assert frozen_change
    # A real change in the original query is retained by a moving dictionary.
    dynamics=exp({(1,1,0,0):p.scale(old.lifted(p.base.Y),C(0,1))})
    actual=ad(dynamics,a0)
    full=ad(u,actual)
    transported_derivative=add(diff(full,1),sc(comm(ga,full),C(0,-1)))
    assert equal(transported_derivative,ad(u,diff(actual,1)))
    assert scalar(mul(localrho,full))==expectation(rho,actual)
    actual_change=expectation(rho,diff(actual,1))
    assert actual_change
    return dict(round=811,all_checks_passed=True,epsilon_order=3,source_order=2,
        BV_matrix_dimension=10,smooth_parameter_count=2,
        closed_selfadjoint_generators=True,flat_comparison_connection=True,
        three_frame_composition_exact=True,full_source_derivatives_commute=True,
        state_and_fixed_observable_changes_cancel=True,
        actual_query_change_retained=True,
        omitted_commutator_curvature_defect=first(curl),
        lost_original_second_source_contact=first(lost),
        falsely_frozen_state_record_change=dump(frozen_change),
        retained_actual_query_change=dump(actual_change),
        original_continuous_smooth_family_proof_in_note=True,
        comparison_generator_identified_with_physical_time_Hamiltonian=False,
        original_sharp_time_operator_nonexistence_proven=False,
        autonomous_apparatus_or_finite_coupling_proven=False,
        new_numbered_test_groups=1)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();result=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite saved evidence.'
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))

