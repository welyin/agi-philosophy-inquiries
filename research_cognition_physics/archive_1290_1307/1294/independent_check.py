#!/usr/bin/env python3
"""1294 independent finite-ring verification. No author implementation imported.
Default/--check: read-only comparison. --write: exclusive first creation.
Exact rational certificates and integer trace moments carry the strict bounds;
finite matrices diagnose the stated finite graph, not a continuum apparatus.
"""
import argparse, hashlib, itertools, json, math
from fractions import Fraction as F
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if (p/'research_cognition_physics').is_dir())
OUT=HERE/'independent_results.json'

def rat(x):
    x=F(x)
    return {'fraction':str(x),'decimal':float(x)}

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def err(a,b=0): return float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
def small(x): return float(f'{x:.10g}')
def mm(a,b): return a@b

def exact_block(n,q):
    # Real multiplication representation of Q or Q(omega), omega^2+omega+1=0.
    d=2 if n==3 else 1
    I=np.eye(d,dtype=object)
    if q==0: R=I.copy(); RI=I.copy()
    elif n==2: R=-I; RI=-I
    else:
        R=np.array([[0,-1],[1,-1]],dtype=object)
        RI=R@R
    h=np.zeros((8*d,8*d),dtype=object)
    def add(i,j,m): h[i*d:(i+1)*d,j*d:(j+1)*d]+=m
    for x in range(4):
        add(2*x+1,2*x+1,I)
        y=(x+1)%4
        for r in range(2):
            add(2*y+r,2*x+r,-RI if x==3 else -I)
            add(2*x+r,2*y+r,-R if x==3 else -I)
    add(4,5,I); add(5,4,I)
    return h,d

def exact_certificates():
    tmax=F(1,512); tmin=F(1,1024); beta=F(1,1024)
    R4=(4*tmax)**4/F(24)/(1-4*tmax/5)
    assert R4/tmax**3==F(40,1917)<F(1,40)
    amp=F(1,6)-F(1,40)
    assert amp**2>F(1,50)
    thermal_tail=16*(4*beta)**5/F(120)/(1-4*beta/6)
    assert thermal_tail/beta**4==F(1024,7675)<F(1,4)
    traces={}
    for n in (2,3):
        moments=[]
        for q in (0,1):
            h,d=exact_block(n,q); power=np.eye(len(h),dtype=object)
            m=[]
            for k in range(5):
                tr=F(int(np.trace(power)),d)
                assert tr.denominator==1
                m.append(int(tr)); power=power@h
            moments.append(m)
        assert moments[0][:4]==moments[1][:4]
        diff=moments[0][4]-moments[1][4]
        assert diff=={2:32,3:24}[n]
        traces[str(n)]={'q0':moments[0],'q1':moments[1],
            'leading_trace_difference_coefficient':rat(F(diff,24))}
    h,_=exact_block(2,1)
    e=np.zeros(8,dtype=object); e[0]=1
    assert np.array_equal(h@h@e,2*e)
    assert all((h@e)[j]==0 for j in (1,3,5,7))
    zsum_bound={str(n):rat(F(8*n)/(1-4*beta)) for n in (2,3)}
    for n in (2,3): assert F(8*n)/(1-4*beta)<9*n
    return {'record_time_interval':[rat(tmin),rat(tmax)],
        'record_relative_tail':rat(R4/tmax**3),
        'record_amplitude_coefficient_lower':rat(amp),
        'record_probability_coefficient_lower':rat(amp**2),
        'uniform_record_gap_strict_lower':rat(tmin**6/50),
        'endpoint_record_gap_strict_lower':rat(tmax**6/50),
        'thermal_beta':rat(beta),'thermal_difference_tail_over_beta4':rat(thermal_tail/beta**4),
        'thermal_trace_moments':traces,'thermal_partition_sum_upper':zsum_bound,
        'thermal_sector_bias_common_strict_lower':rat(beta**4/54),
        'dark_subspace_exact_relation':'h^2 e0=2 e0; e0,h e0 have record zero',
        'hamiltonian_lower_bound':1,'hamiltonian_upper_bound':8,
        'h_norm_upper_bound':4}

def car_checks():
    dim=16; cs=[]
    for j in range(4):
        c=np.zeros((dim,dim))
        for b in range(dim):
            if b&(1<<j): c[b^(1<<j),b]=(-1)**((b&((1<<j)-1)).bit_count())
        cs.append(c)
    one=[1<<j for j in range(4)]
    residual=0.
    for i in range(4):
        for j in range(4):
            residual=max(residual,err(cs[i]@cs[j].T+cs[j].T@cs[i],np.eye(dim)*(i==j)))
            actual=(cs[i].T@cs[j])[np.ix_(one,one)]
            expected=np.zeros((4,4)); expected[i,j]=1
            residual=max(residual,err(actual,expected))
    assert residual==0
    return {'positive_charge_car_modes':4,'full_positive_car_dimension':16,
        'fixed_one_particle_dimension':4,'exact_car_and_hopping_residual':residual}

def physical_ring(n):
    # Independent enumeration of all electric configurations and the Gauss kernel.
    states=[]
    for x in range(4):
        for es in itertools.product(range(n),repeat=4):
            if all((es[j]-es[(j-1)%4]-(int(j==x)-int(j==0)))%n==0 for j in range(4)):
                for r in (0,1): states.append((x,es,r))
    assert len(states)==8*n
    lookup={s:i for i,s in enumerate(states)}; D=len(states)
    T=np.zeros((D,D)); nr=np.zeros((D,D)); V=np.zeros((D,D))
    local=[]
    for edge in range(4):
        te=np.zeros((D,D))
        for j,(x,es,r) in enumerate(states):
            if x==edge:
                out=list(es); out[x]=(out[x]-1)%n
                te[lookup[((x+1)%4,tuple(out),r)],j]=1
        local.append(te)
        T+=te
    for j,(x,es,r) in enumerate(states):
        nr[j,j]=r
        if x==2: V[lookup[(x,es,1-r)],j]=1
    K=T+T.T; h=-K+nr+V; H=4*np.eye(D)+h
    J=np.zeros((D,D),complex)
    for q in range(n):
        for x in range(4):
            for r in (0,1):
                c=8*q+2*x+r
                for l in range(n):
                    es=tuple((l-int(j<x))%n for j in range(4))
                    J[lookup[(x,es,r)],c]=np.exp(-2j*np.pi*q*l/n)/np.sqrt(n)
    return states,local,K,nr,V,h,H,J

def taylor_unitary(h,t,order=30):
    out=np.eye(len(h),dtype=complex); term=out.copy()
    for k in range(1,order+1):
        term=term@(-1j*t*h)/k; out+=term
    return out

def finite_diagnostics(n):
    states,local,K,nr,V,h,H,J=physical_ring(n); D=len(states)
    resid=err(J.conj().T@J,np.eye(D))
    hf=J.conj().T@h@J
    off=hf.copy()
    for q in range(n): off[8*q:8*q+8,8*q:8*q+8]=0
    resid=max(resid,err(off))
    thirds=[]
    for q in range(n):
        phi=2*np.pi*q/n
        block=hf[8*q:8*q+8,8*q:8*q+8]
        e=np.zeros(8,complex); e[0]=1
        for k in (0,1,2): resid=max(resid,err((np.linalg.matrix_power(block,k)@e)[1::2]))
        expected=np.zeros(4,complex); expected[2]=1+np.exp(1j*phi)
        third=(block@block@block@e)[1::2]
        resid=max(resid,err(third,expected)); thirds.append([small(abs(z)) for z in third])
    code=J[:,[0,8]]
    resid=max(resid,err(code.conj().T@H@code,4*np.eye(2)),
        err(code.conj().T@H@H@code,18*np.eye(2)))
    P=[np.diag([int(s[2]==r) for s in states]) for r in (0,1)]
    smallU=taylor_unitary(h,1/512)
    effects=[code.conj().T@smallU.conj().T@p@smallU@code for p in P]
    probs=[float(effects[1][i,i].real) for i in range(2)]
    resid=max(resid,err(effects[0]+effects[1],np.eye(2)))
    if n==3: assert probs[1]>float(F(1,512)**6/F(50))
    else: assert probs[1]<1e-30
    # An order-one diagnostic time checks the complete poststate and work ledger.
    ev,Ue=np.linalg.eigh(H); U=(Ue*np.exp(-1j*ev/3))@Ue.conj().T
    W=U@code; Ks=[p@W for p in P]
    resid=max(resid,err(sum(k.conj().T@k for k in Ks),np.eye(2)))
    bell=np.array([1,0,0,1],complex)/np.sqrt(2)
    joint=np.kron(W,np.eye(2))@bell
    rho=np.outer(joint,joint.conj()).reshape(D,2,D,2)
    rr=np.einsum('iaib->ab',rho)
    resid=max(resid,err(rr,np.eye(2)/2))
    rin=np.array([[.37,.12+.09j],[.12-.09j,.63]],complex)
    out=W@rin@W.conj().T
    deph=sum(p@out@p for p in P)
    work=float(np.trace(H@(deph-out)).real)
    prescribed=-float(np.trace(V@out).real)
    resid=max(resid,abs(work-prescribed),abs(float(np.trace(H@out).real)-4))
    # Both finite masses and all three couplings are retained as separate sources.
    directions=[np.eye(D),np.eye(D),-K,nr,V]
    pars=np.array([2.,2.,1.,1.,1.]); eps=1e-5
    for j,s in enumerate(directions):
        plus=sum((pars[k]+(eps if k==j else 0))*directions[k] for k in range(5))
        minus=sum((pars[k]-(eps if k==j else 0))*directions[k] for k in range(5))
        assert err((plus-minus)/(2*eps),s)<1e-9
    residual_sources=err(sum(pars[k]*(code.conj().T@U.conj().T@directions[k]@U@code) for k in range(5)),4*np.eye(2))
    resid=max(resid,residual_sources)
    pieces=[pars[k]*directions[k] for k in range(5)]
    resid=max(resid,err(sum(1j*(H@piece-piece@H) for piece in pieces)))
    # Removing the local coupling keeps the record identically zero.
    U0=taylor_unitary(-K+nr,1/3)
    resid=max(resid,err(P[1]@U0@code))
    # Disable the return hopping while keeping the full ring Hilbert/Gauss/code data.
    hopen=h+(local[3]+local[3].T)
    hopf=J.conj().T@hopen@J
    resid=max(resid,err(hopf[:8,:8],hopf[8:16,8:16]))
    beta=1/1024
    expH=(Ue*np.exp(-beta*ev))@Ue.conj().T
    thermal=expH/np.trace(expH)
    weights=[]
    for q in range(n):
        jq=J[:,8*q:8*q+8]
        weights.append(float(np.trace(jq.conj().T@thermal@jq).real))
    assert weights[0]>1/n
    assert resid<1e-12
    return {'group_order':n,'physical_dimension':D,'kinematic_fixed_charge_dimension':8*n**4,
        'max_algebraic_residual':small(resid),'H_spectrum_extrema':[small(ev[0]),small(ev[-1])],
        'third_order_record_coefficients_abs':thirds,
        'endpoint_record_probabilities_q0_q1':[small(p) for p in probs],
        'Gibbs_sector_weights_diagnostic_only':[float(f'{x:.16g}') for x in weights],
        'Gibbs_q0_bias_diagnostic_only':small(weights[0]-1/n),
        'source_parameters':['m_c','m_d','g','a_p','r'],
        'initial_energy':4,'initial_second_energy_moment':18,
        'dephasing_work_t_one_third':small(work),'minus_interaction_mean':small(prescribed),
        'arbitrary_reference_certificate':'Kraus completeness as operator equality; Bell is a diagnostic',
        'deletions':['local record coupling zero: no record','return hopping disabled on retained ring: q0 and q1 blocks equal']}

def calculate():
    return {'status':'PASS_WITH_SCOPE','group_count':1,
        'exact_certificates':exact_certificates(),'CAR':car_checks(),
        'finite_diagnostics':[finite_diagnostics(n) for n in (2,3)],
        'scope':['fixed one-positive/one-negative charge sector; negative position fixed at zero',
            'finite graph sources are parameter derivatives, not continuum Hilbert stress',
            'local record formed autonomously; final ideal readout is separately licensed',
            'full H Gibbs on 8N physical dimensions, not a uniform topological code reference',
            'no claim of all FUCP controls, four-dimensional covariance, or HG completion']}

def bindings():
    paths=[HERE/'proof.md',HERE/'check.py',HERE/'results.json',Path(__file__).resolve(),
        HERE.parent/'_shared/round1294_admission.md']
    return {str(p.relative_to(ROOT)).replace('\\','/'):digest(p) for p in paths}

def main():
    ap=argparse.ArgumentParser(); group=ap.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true'); group.add_argument('--check',action='store_true')
    group.add_argument('--diagnose',action='store_true',help='read-only unbound calculation before author freeze')
    args=ap.parse_args(); science=calculate()
    if args.diagnose:
        print(json.dumps(science,ensure_ascii=False,indent=2)); return
    data={'schema':'1294-independent-v1','scientific_payload':science,'files_sha256':bindings()}
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as f: json.dump(data,f,ensure_ascii=False,indent=2); f.write('\n')
        print('CREATED PASS_WITH_SCOPE')
    else:
        old=json.loads(OUT.read_text(encoding='utf-8'))
        if old!=data: raise AssertionError('saved independent result/bindings differ')
        print('READ_ONLY_MATCH PASS_WITH_SCOPE')
    print('results_sha256='+digest(OUT))

if __name__=='__main__': main()
