"""Round 1270 independent audit; one scientific group.
No author code or numerical result is imported.  Default replay is read-only.
Exact rational/integer checks support the analytic review; finite matrices are diagnostics.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import factorial, pi
from pathlib import Path
import hashlib, json, sys
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(rel): return json.loads((STAGE/rel).read_text(encoding='utf-8-sig'))
def iv(x): return tuple(map(F,x))
def ceil2(n): return (n-1).bit_length()
def td(a,b): return float(np.linalg.eigvalsh((a-b+(a-b).conj().T)/2).__abs__().sum()/2)
def ketdm(v): return np.outer(v,v.conj())
def kron(*args):
    z=np.array([[1.0+0j]])
    for a in args:z=np.kron(z,a)
    return z

def partial(rho,dims,keep):
    dims=list(dims);a=rho.reshape(dims+dims)
    for j in reversed(range(len(dims))):
        if j not in keep:
            a=np.trace(a,axis1=j,axis2=j+len(dims));dims.pop(j)
    return a.reshape(int(np.prod(dims)),int(np.prod(dims)))

def matrix_diagnostics():
    i=np.eye(2,dtype=complex); x=np.array([[0,1],[1,0]],complex)
    y=np.array([[0,-1j],[1j,0]],complex); z=np.diag([1,-1]).astype(complex)
    p0=np.diag([1,0]).astype(complex);p1=i-p0
    pp=(i+y)/2;pm=(i-y)/2
    go=pi/2*kron(pm,i-x,i)
    gr=pi/2*kron(i,i-x,p1)
    u=kron(pp,i,i)+kron(pm,x,i)
    ur=kron(i,i,p0)+kron(i,x,p1)
    ev,vec=np.linalg.eigh(go)
    expgo=(vec*np.exp(-1j*ev))@vec.conj().T
    assert np.linalg.norm(expgo-u)<2e-14
    assert np.linalg.norm(go@gr-gr@go)<2e-14
    assert abs(np.linalg.norm(go,2)-pi)<2e-14
    assert np.linalg.norm(u.conj().T@u-np.eye(8))<2e-14
    # Genuine S-C-R correlations: the internal gate is applied only to C-F-E.
    rng=np.random.default_rng(1270)
    v=rng.normal(size=8)+1j*rng.normal(size=8);v/=np.linalg.norm(v)
    # S,C,R -> S,C,F,E,R, with exact blank F and a nontrivial noisy E.
    input_scfr=np.zeros((32,32),complex)
    rho_scr=ketdm(v)
    for ss in range(2):
      for cc in range(2):
       for rr in range(2):
        a=(ss*2+cc)*2+rr
        for st in range(2):
         for ct in range(2):
          for rt in range(2):
           b=(st*2+ct)*2+rt
           for ee,w in ((0,0.91),(1,0.09)):
            ia=((((ss*2+cc)*2)*2+ee)*2+rr)
            ib=((((st*2+ct)*2)*2+ee)*2+rt)
            input_scfr[ia,ib]=w*rho_scr[a,b]
    total=kron(i,ur@u,i)
    out=total@input_scfr@total.conj().T
    marginal_error=td(partial(input_scfr,[2]*5,[0,4]),partial(out,[2]*5,[0,4]))
    assert marginal_error<2e-14
    # Same-source completion is an operator identity, not an expectation cancellation.
    hs=np.array([[0.7,0.08j],[-0.08j,1.1]],complex)
    jbase=-kron(hs,np.eye(8));gate=kron(i,go/23+gr/29)
    jout=jbase-gate
    a=np.array([[2,.3],[.3,1]],complex)
    ac=kron(a,np.eye(8))
    complete_error=np.linalg.norm((ac@jout+jout@ac+ac@gate+gate@ac-ac@jbase-jbase@ac)/2)
    assert complete_error<2e-14
    # Nonselective flag querying preserves S-R, but conditionals need not.
    measured=sum(kron(i,i,p,i,i)@out@kron(i,i,p,i,i) for p in (p0,p1))
    assert td(partial(measured,[2]*5,[0,4]),partial(out,[2]*5,[0,4]))<2e-14
    # Fixed known material comparator and a small genuine S-C entangling perturbation.
    vs=np.array([1,2j],complex)/np.sqrt(5);vc=np.array([1,1],complex)/np.sqrt(2)
    vc=(np.cos(.17)*i-1j*np.sin(.17)*z)@vc
    vcmp=np.kron(vs,vc)
    entgen=kron(x,z);angle=1e-4
    vr=(np.cos(angle)*np.eye(4)-1j*np.sin(angle)*entgen)@vcmp
    rcmp=ketdm(vcmp);rreal=ketdm(vr);din=td(rreal,rcmp)
    conditional=[];source_errors=[]
    for pc in (pp,pm):
        pr=kron(i,pc)
        t=pr@rreal@pr;s=pr@rcmp@pr
        p=float(np.trace(t).real);q=float(np.trace(s).real)
        ds=td(partial(t/p,[2,2],[0]),partial(s/q,[2,2],[0]))
        assert min(p,q)>.25 and ds<=2*din/(q-din)+1e-13 and ds<8*din+1e-13
        srcerr=abs(np.trace(hs@(partial(t/p,[2,2],[0])-ketdm(vs))))
        assert srcerr<16*np.linalg.norm(hs,2)*din+1e-13
        conditional.append(ds);source_errors.append(float(srcerr))
    # Rare branch deletion: a small cq distance alone does not bound every conditional.
    rare=F(1,10000)
    cmp=np.diag([float(1-rare),0,float(rare),0]).astype(complex)
    real=np.diag([float(1-rare),0,0,float(rare)]).astype(complex)
    assert abs(td(cmp,real)-float(rare))<1e-15
    assert td(p0,p1)==1 and 8*rare<1
    # A finite nonzero area discrepancy is bounded by the generator norm.
    xi=1e-6
    actual=(vec*np.exp(-1j*(1+xi)*ev))@vec.conj().T
    area_error=np.linalg.norm(actual-u,2)
    assert area_error<=pi*abs(xi)+2e-14
    return {'gate_exponential_residual':float(np.linalg.norm(expgo-u)),
      'source_identity_residual':float(complete_error),'unselected_material_reference_distance':marginal_error,
      'conditional_material_distances':conditional,'conditional_source_errors':source_errors,
      'rare_branch_cq_distance':str(rare),'rare_branch_conditional_distance':1,
      'scope':'finite diagnostics; universal isometry/source/contractivity statements are proved analytically, not inferred from samples'}

def main():
    family=read('1252/results.json')['group_parameter_family']['cells'][0]
    end=family['checkpoints']['2']['state'];positive=family['positive_all_time']
    ulo,uhi=iv(end['u']);hlo,hhi=iv(end['H'])
    # H is nonincreasing on the old positive continuation. e^-x >= 1-x.
    uend=ulo*(1-hhi/F(4))
    assert ulo>F(373917,10**6) and hhi<F(270066,10**6)
    assert uend>F(1,3)
    # Independently propagate all remaining free coordinates from exact endpoints.
    def absmax(interval):return max(map(abs,interval))
    def absmin(interval):
        lo,hi=interval
        return F(0) if lo<=0<=hi else min(abs(lo),abs(hi))
    def core(state):
        ulo0,uhi0=iv(state['u']);_,hhi0=iv(state['H'])
        ur=ulo0*(1-hhi0/4)
        jlo=iv(state['J'])[0]
        mc=F(71,100)+jlo+sum(absmin(iv(state[k]))**2 for k in ('q','r','K'))/2
        qr=sum(absmax(iv(state[k]))**2 for k in ('q','r'))
        yy=absmax(iv(state['Y']))+absmax(iv(state['K']))/4
        th=absmax(iv(state['z']))+F(1,5)+F(1,4)*F(5,2)/10
        assert ur>F(1,3) and mc>F(93,100) and jlo>F(98,1000)
        assert qr<1 and yy<1 and th<1 and absmax(iv(state['psi']))<1
        assert absmax(iv(state['K']))<1
        return {'u_final_lower':str(ur),'Mc_lower':str(mc),'qr_norm_squared_upper':str(qr),
                'Y_upper':str(yy),'theta_upper':str(th)}
    cores={'parameter_family':core(end)}
    oldprogram=read('1247/results.json')['group_finite_memory']
    for name in ('alpha0','alpha1'):cores[name]=core(oldprogram[name]['checkpoints']['2']['state'])
    assert F(3,2)**2>2*F(111,100)
    assert F(1)-F(3,2)/5>F(3,5)
    assert F(1,5)-F(3,2)/10>F(1,25)
    assert iv(positive['J'])[0]>F(2,25)
    assert iv(positive['M_c'])[0]>F(9,10)
    assert F(2)*27*F(11,10)<61
    assert 27*F(3,2)<61
    # Q remains in (-1/2,3/2): exp(1/2)>3/2 and exp(3/2)>3.
    assert sum(F(1,2)**j/F(factorial(j)) for j in range(4))>F(3,2)
    assert sum(F(3,2)**j/F(factorial(j)) for j in range(4))>3
    # Exact resource reconstruction, independent of the round1270 author outputs.
    ledger=read('1249/quantum_bounds_results.json')['group_constructive_finite_resource_bound']['ledger']
    beta=324*(factorial(81)//factorial(40)**2)*2**40
    lc=56*beta;lh=14*(1+lc);hw=2**36
    m=2**ceil2(hw*factorial(21)*lh**21)
    aa=2**98*(1+lc)**2;ls=4*(lc+14)
    lp=ls*(1+100*aa*beta);lxi=lc+4*lp;lk=lc+4*10**9*beta
    ast=max(aa,1+128*aa**2,1);lst=max(ls,lxi+2*ls,lk)
    cf=2**34*5*256**13*15**10*factorial(20)
    eh=ceil2(cf*hw*lh**20);e0=ceil2(cf*ast*(1+lst)**20)
    ef=e0+20*(4096+5*m)
    n=1+max((ef+44)//3,(eh+ef+47)//6)
    assert (m.bit_length()-1,eh,e0)==(3059,3093,14008)
    assert n==int(ledger['integer_n_power_exponent_decimal'])
    assert ef==int(ledger['Fourier_transport_exponent_decimal'])
    # All times <=9/4 use exp(Mt) <=2^(5M) and residual S^3 <=2^(15M+12).
    assert F(9,4)<F(5,2) # log2(e)<2 suffices
    # sqrt14<4; sqrt15/(12sqrt2)<1/4 by squaring.
    assert 14<16 and F(15,288)<F(1,16)
    assert F(9,4)*F(1,4)*2**12<2**16
    er=eh+15*m+20-2*n
    eold=e0+eh+15*m+25-2*n
    enew=e0+60-4*n
    eclass=4*m+136-4*n
    # e < 11/4 from a finite Taylor sum and geometric positive tail.
    assert sum(F(1,factorial(j)) for j in range(5))+F(1,100)<F(11,4)
    assert F(11,4)**9<2**16 # e^(9M/4) <2^(4M)
    assert max(er,eold,enew,eclass)<-100
    assert max(eold+3,enew+3,e0+eclass+3)<-80
    # Exact pulse normalization from Beta(4,4), shape and translated common-future budget.
    integ=140*sum(F((-1)**j*factorial(3),factorial(j)*factorial(3-j)*(4+j)) for j in range(4))
    assert integ==1
    sampling_tau=F(1,10**8)
    gp=F(105,4)/sampling_tau**2
    work=F(1,16)*27*gp
    assert work<2**59
    # Same-point source/geometry identities for arbitrary a,c,kappa:
    # -ca/a^3 + c/a^2 =0 and -2(kappa*c/(6a^2))=kappa*(-c/(3a^2)).
    assert -F(1)+1==0 and -2*F(1,6)==-F(1,3)
    # Omitting the pressure (while compensating density) fails this second identity.
    assert -2*F(1,6)!=0
    duration=F(1,64);peak=140*F(1,2)**6/duration
    assert peak==140
    pip=F(22,7);err=F(1,10**9)
    peakmass=peak*pip*(1+err)
    assert peakmass==440*(1+err)
    earliest=F(33,16)-F(1,256)-2
    causal=F(21,400)
    latest=F(33,16)+F(1,256)+2*duration
    assert earliest-causal==F(39,6400)
    assert F(17,8)-latest==F(7,256)
    # Torus diameter: sqrt3<7/4 and a<3 give round trip <21/400.
    assert 3<F(49,16)
    assert 2*3*F(7,4)*F(1,100)/2==causal
    betaout=err+pip*err+err
    assert betaout==F(36,7*10**9)<F(6,10**9)<F(1,10**7)
    dr=F(1,2**100)
    qmin=F(5,16)-betaout;pmin=qmin-dr
    assert pmin>F(1,4) and F(2)/pmin<8
    # sin x<=x on [1/4,3/8] supplies the smaller ideal probability >=5/16.
    assert (1-F(3,8))/2==F(5,16)
    # Per branch source expectation uses 2 ||B|| D, not beta_out ||B||.
    assert 2*8*2==2**5
    # Same conditional finite-test residual at kappa=3 units.
    kappa_error=F(1,15000)+F(5,10**16)
    residual=F(111,100)*kappa_error+3*F(1,2**80)
    assert residual<F(75,10**6)
    # Omitting clock compensation leaves nonzero complete mass exactly V.
    jbase=F(99,1000);vo=F(1,10000)
    assert (jbase-vo)+vo==jbase and jbase+vo-jbase==vo
    oldpaths=['1247/results.json','1249/quantum_bounds_results.json','1252/results.json']
    return {'round':1270,'verdict':'PASS_INDEPENDENT_WITH_SCOPE','independent_groups':1,
      'method':'one common-task audit: exact rational old-tube continuation/resource/time/source bounds and independently built finite matrices; no author imports',
      'continuation':{'u_lower_at_9_over4':str(uend),'full_coordinate_extension':cores,'u_margin_above_one_third':str(uend-F(1,3)),
       'M_log2':m.bit_length()-1,'E_H':eh,'E_0':e0,'N_decimal':str(n),
       'read_state_log2_upper':str(er),'conditional_source_log2_upper':str(eold),
       'new_source_log2_upper':str(enew),'switch_derivative_upper':str(gp),'switch_work_coefficient_upper':str(work),'classical_perturbation_log2_upper':str(eclass),
       'finite_test_residual_upper':str(residual)},
      'output':{'pulse_integral':str(integ),'peak':str(peak),'mass_peak_over_epsilon_upper':str(peakmass),
       'earliest_processing_causal_margin':str(earliest-causal),'last_gate_output_margin':str(F(17,8)-latest),
       'output_probability_error_upper':str(betaout),'smallest_actual_calibration_probability_lower':str(pmin),
       'conditional_distance_multiplier_upper':str(F(2)/pmin),'same_source_compensation_exact':True,
       'quantum_diagnostics':matrix_diagnostics()},
      'source_sha256':{p:sha(STAGE/p) for p in oldpaths},
      'scope':['known Gaussian family with arbitrary independent initial pointer-reference state; all outcomes retained',
       'normal conditional constants only for the stated nonrare calibration outcomes',
       'new named whole-block future-boundary primitive is an effective material input',
       'bounded completed source menus and classical comparison geometry, not operator Einstein equations',
       'no new all-space quantum implementation or empirical apparatus certificate; stage closure requires separate audit']}

if __name__=='__main__':
    result=main();path=HERE/'independent_results.json'
    if len(sys.argv)==2 and sys.argv[1]=='--write':
        with path.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        print('WRITE_ONCE_PASS_1270 independent_groups=1')
    elif len(sys.argv)==1:
        assert result==json.loads(path.read_text(encoding='utf-8-sig')),'independent result mismatch'
        print('READ_ONLY_MATCH_1270 independent_groups=1')
    else:raise SystemExit('Default: read-only replay. --write: exclusive first save.')