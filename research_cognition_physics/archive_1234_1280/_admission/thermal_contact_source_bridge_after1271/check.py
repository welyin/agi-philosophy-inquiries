"""Read-only finite source/quantum comparison diagnostics, not a device model."""
from pathlib import Path
from fractions import Fraction as F
import json
import sys
import numpy as np

D = Path(__file__).resolve().parent


def add(a, b):
    out = [F(0)] * max(len(a), len(b))
    for i, x in enumerate(a): out[i] += x
    for i, x in enumerate(b): out[i] += x
    return out


def scale(a, c): return [c*x for x in a]


def mul(a, b):
    out = [F(0)] * (len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i+j] += x*y
    return out


def der(a): return [F(i)*a[i] for i in range(1, len(a))] or [F(0)]
def integ(a): return sum((x/F(i+1) for i,x in enumerate(a)), F(0))
def end(a): return sum(a, F(0))-a[0]


def source_identity(f, q, a, b):
    # Integrated stress P=(Q''-B-A')/2, no stationarity assumption.
    p = scale(add(add(der(der(q)), scale(b, -1)), scale(der(a), -1)), F(1,2))
    lhs = integ(mul(f,p))
    boundary = F(1,2)*end(add(add(mul(f,der(q)),scale(mul(der(f),q),-1)),scale(mul(f,a),-1)))
    rhs = boundary + F(1,2)*(integ(mul(der(der(f)),q))-integ(mul(f,b))+integ(mul(der(f),a)))
    assert lhs == rhs
    return {'lhs':str(lhs), 'rhs':str(rhs), 'endpoint':str(boundary),
            'incorrect_if_drop_endpoint':str(rhs-boundary)}


def record_prob(theta, kept):
    # Full A4 B4 plus M,E, with each thermal preparation label retained.
    S=np.zeros((16,16),complex)
    for a in range(4):
        for b in range(4): S[b*4+a,a*4+b]=1
    U=lambda x:np.cos(x)*np.eye(16)-1j*np.sin(x)*S
    prob_m=prob_a_minus=prob_joint=0.0
    for n,qn in enumerate([.75,.25]):
        for m,qm in enumerate([.6,.4]):
            psi=np.zeros((16,2,2),complex)
            psi[(2*n)*4+(2*m+1),0,0]=1
            psi=np.einsum('ab,bij->aij',U(theta),psi)
            if kept:
                nxt=np.zeros_like(psi)
                for k in range(16):
                    r=(k//4)%2
                    for i in range(2):
                        for j in range(2): nxt[k,i^r,j^r]=psi[k,i,j]
                psi=nxt
            psi=np.einsum('ab,bij->aij',U(np.pi/4),psi)
            assert abs(np.vdot(psi,psi).real-1)<1e-12
            weight=qn*qm
            prob_m+=weight*np.sum(abs(psi[:,1,:])**2)
            for k in range(16):
                aa=k//4
                if aa%2:
                    p=weight*np.sum(abs(psi[k])**2)
                    prob_a_minus+=p
                    if aa//2: prob_joint+=p
    return {'M_1':float(prob_m),'A_minus':float(prob_a_minus),'J1_Aminus':float(prob_joint)}


def calculate():
    total_p=F(0); mass=F(0); pressure=F(0); energy=F(0)
    # Targets+stocks (each side mass7) and one seed mass1 on each side.
    for m,r in [(7,1),(7,-1),(1,1),(1,-1)]:
        total_p+=F(4,3)*m*r
        mass+=m; pressure+=F(16,15)*m; energy+=F(5,3)*m
    assert (total_p,mass,pressure,energy)==(0,16,F(256,15),F(80,3))
    assert energy+4==F(92,3)
    pressure_check={'total_spatial_P':str(total_p),'moving_rest_mass':str(mass),
                    'free_bundle_pressure_xx':str(pressure),'new_energy':str(energy+4)}
    target=np.array([np.cos(np.pi/4),-1j*np.sin(np.pi/4)])
    identity=np.array([1,0],complex)
    delta=np.outer(target,target.conj())-np.outer(identity,identity.conj())
    trace_half=float(np.sum(abs(np.linalg.eigvalsh(delta)))/2)
    assert abs(trace_half-1/np.sqrt(2))<1e-14
    # Dimensionless polynomial tests; temporal rescaling preserves identities.
    q=list(map(F,[1,2,3,-4,2])); a=list(map(F,[2,-1,3])); b=list(map(F,[-1,4,-2]))
    interior=source_identity(list(map(F,[0,0,1,-2,1])),q,a,b)
    endpoints=source_identity(list(map(F,[1,1])),q,a,b)
    assert endpoints['endpoint'] != '0'
    flat=source_identity(list(map(F,[0,0,1,-2,1])),q,[F(0)],[F(0)])
    lam=F(1,100); W=F(1,50); L=Lbar=F(1); remainder=F(1,10)
    bound=abs(lam)*W/2+lam*lam*(L*L+Lbar*Lbar)/2+lam*lam*remainder
    assert bound==F(21,100000)
    records=[]
    for off in [0,np.pi/12,-np.pi/12,.01]:
        kept=record_prob(np.pi/4+off,True)
        bare=record_prob(np.pi/4+off,False)
        assert abs(kept['M_1']-np.sin(np.pi/4+off)**2)<1e-13
        assert abs(kept['A_minus']-.5)<1e-13
        assert abs(bare['A_minus']-np.cos(off)**2)<1e-13
        assert abs(kept['J1_Aminus']-.2)<1e-13
        records.append({'angle_offset':float(off),'kept':kept,'unrecorded':bare})
    return {'status':'PASS_CONDITIONAL_ARITHMETIC_ONLY','new_actual_contact_model':False,
            'new_rounds':0,'new_axioms':0,'pressure':pressure_check,
            'free_wait_vs_gate_half_trace':trace_half,
            'finite_source_identities':{'flat_interior':flat,'forced_interior':interior,'forced_endpoints':endpoints},
            'illustrative_not_certified_channel_budget':str(bound),
            'original_record_calibration':records}


def main():
    data=calculate()
    p=D/'results.json'
    if '--write' in sys.argv:
        with p.open('x',encoding='utf-8') as f: json.dump(data,f,ensure_ascii=False,indent=2); f.write('\n')
        print('SAVED_FINITE_SOURCE_AND_ORIGINAL_RECORD_DIAGNOSTIC')
    else:
        old=json.loads(p.read_text(encoding='utf-8'))
        assert old==data
        print('READ_ONLY_MATCH_CONDITIONAL_SOURCE_BRIDGE_DIAGNOSTIC')


if __name__=='__main__': main()
