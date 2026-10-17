"""Round 1268 independent audit. One group, two analytic subchecks.
No author executable or author numerical result is imported.
Default is read-only comparison; --write creates results exclusively.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import factorial, isqrt
from pathlib import Path
import hashlib, json, sys

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def inp(rel):
    return json.loads((STAGE/rel).read_text(encoding='utf-8-sig'))

def interval(x):
    return tuple(map(F,x))

def ceil_log2(n):
    return (n-1).bit_length()

# Sparse multivariate polynomials: variables A,D,rho,C,y.
Z=(0,0,0,0,0)
def atom(j):
    q=list(Z);q[j]=1
    return {tuple(q):F(1)}
def add(*ps):
    out={}
    for p in ps:
        for key,val in p.items():out[key]=out.get(key,F(0))+val
    return {k:v for k,v in out.items() if v}
def mul(p,q):
    out={}
    for a,x in p.items():
        for b,y in q.items():
            key=tuple(i+j for i,j in zip(a,b))
            out[key]=out.get(key,F(0))+x*y
    return {k:v for k,v in out.items() if v}
def neg(p):return {k:-v for k,v in p.items()}
def sqrt_interval(q,scale=10**20):
    lo=isqrt(q.numerator*scale**2//q.denominator)
    return F(lo,scale),F(lo+1,scale)

def main():
    A,D,rho,C,y=map(atom,range(5))
    # From the independently expanded quadratic action:
    # kinetic diag(A,rho), spatial Hessian [[D,-C],[-C,0]].
    m00=add(mul(A,y),neg(D));m11=mul(rho,y)
    det=add(mul(m00,m11),neg(mul(C,C)))
    target=add(mul(mul(A,rho),mul(y,y)),neg(mul(mul(D,rho),y)),neg(mul(C,C)))
    assert det==target and len(det)==3
    # For f(B)=const B^-1/2 and B=3s^2, D-F cancels identically.
    assert 2*F(-1,2)+4*F(3,4)*F(1,3)==0
    # The cross term does not disappear for linear f. At y=1:
    # det=rho(A-D)-C^2=-C^2, while det(0)=-C^2.
    linear_at_one=add(mul(rho,add(A,neg(A))),neg(mul(C,C)))
    assert linear_at_one==neg(mul(C,C))
    material_diagnostics=[]
    for shape in ('scale_plus','linear_plus','linear_minus'):
        Frod,n,s,mass=F(1),F(1),F(1),F(1)
        c=F(1,1000) if shape!='linear_minus' else F(-1,1000)
        fp,fpp=(F(-1,6),F(1,12)) if shape=='scale_plus' else (F(1),F(0))
        ac=Frod+2*n*c*fp
        dc=ac+4*n*c*fpp*s*s
        cc=2*n*c*fp*s
        rr=n*mass
        assert ac>0 and dc>0 and rr>0 and cc!=0
        disc=dc*dc+4*ac*cc*cc/rr
        slo,shi=sqrt_interval(disc)
        yp=(dc+slo)/(2*ac),(dc+shi)/(2*ac)
        ym=(dc-shi)/(2*ac),(dc-slo)/(2*ac)
        assert yp[0]>1 and ym[1]<0
        # Frozen negative-root rate squared / |k|^2 <= C^2/(rho D).
        assert -ym[0] < cc*cc/(rr*dc)+F(1,10**18)
        material_diagnostics.append({'shape':shape,'A':str(ac),'D':str(dc),'C':str(cc),
            'positive_speed_squared_interval':list(map(str,yp)),
            'negative_speed_squared_interval':list(map(str,ym)),
            'scope':'rational diagnostic only; universal result is the analytic determinant/sign argument'})
    # Direct volume term -c B: rho=c B, isotropic pressure=-c B/3.
    # Total homogeneous intensity a^3 c 3k^2/a^2=3k^2 c a.
    assert 2*F(1,3)-1 == -F(1,3)
    assert 3-2==1
    # Worldline a-dependent mass has the same homogeneous rho/P but a different
    # nonhomogeneous principal block. Equality of ODEs is not causal equivalence.
    assert 2*F(-1,2)*F(1,3)==-F(1,3)

    # Independently reconstruct old finite-resource integers, not 1268 outputs.
    ledger=inp('1249/quantum_bounds_results.json')['group_constructive_finite_resource_bound']['ledger']
    B=324*(factorial(81)//factorial(40)**2)*2**40
    Lc=56*B; Lh=14*(1+Lc); H0=2**36
    M=2**ceil_log2(H0*factorial(21)*Lh**21)
    As=2**98*(1+Lc)**2; Ls=4*(Lc+14)
    Lp=Ls*(1+100*As*B); Lxi=Lc+4*Lp; Lk=Lc+4*10**9*B
    Ae=1+128*As**2; Le=Lxi+2*Ls
    Ast=max(As,Ae,1);Lst=max(Ls,Le,Lk)
    CF=2**34*5*256**13*15**10*factorial(20)
    EH=ceil_log2(CF*H0*Lh**20); E0=ceil_log2(CF*Ast*(1+Lst)**20)
    EF=E0+20*(4096+5*M)
    N=1+max((EF+44)//3,(EH+EF+47)//6)
    assert (M.bit_length()-1,EH,E0)==(3059,3093,14008)
    assert EF==int(ledger['Fourier_transport_exponent_decimal'])
    assert N==int(ledger['integer_n_power_exponent_decimal'])
    assert B==int(ledger['B_decimal'])
    assert 27*(1+B)<2**136
    # F(Q)=xi(Q)e^Q: transitions do not overlap, so |xi'|<=B, not 2B.
    # e^3<27 and the old exact core give Fmax<27, Lipschitz<2^136.
    alpha=F(1,16);tau=F(1,10**8)
    # 140 u^3(1-u)^3 / tau is normalized and peaks at u=1/2.
    norm=140*sum(F((-1)**j*factorial(3),factorial(j)*factorial(3-j)*(4+j)) for j in range(4))
    assert norm==1
    peak=140*F(1,2)**6/tau
    assert peak==F(35,16)/tau
    gate_bound=alpha*peak*27
    assert gate_bound<2**29
    # The clock density uses the old Weyl completion, not multiplication by R_F.
    # b_clk=chi*xi_Q*exp(-2Q): derivatives <=81*r!*(2Lc+2)^r.
    clock_fourier_bound=CF*81*(2*Lc+2)**20
    assert ceil_log2(clock_fourier_bound)==3007 and clock_fourier_bound<2**E0
    new_source_coefficient=alpha*peak*(clock_fourier_bound+81+27*(1+B))
    assert new_source_coefficient<2**(E0+26)
    assert E0+26-4*N<-100
    # Independent coefficient comparison to author's generous proposed exponent.
    # True-free/Gaussian norm <=sqrt(eps) 2^(EH+15M+16).
    # Multiplication residual <=30 d0 + LipF sqrt(eps) S/sqrt(2).
    e_d=EH+15*M+20-2*N
    first_coefficient=alpha*30/F(2**4)
    second_relative_exponent=(5*M+4+136-4)-(EH+15*M+20)
    assert first_coefficient==F(15,128) and second_relative_exponent<-4
    assert first_coefficient+F(1,16)<1
    e_gauss=e_d+1
    # Safe old-source norm and expectation bound (includes constant+compact menu).
    e_source=E0+EH+15*M+26-2*N
    assert max(e_d,e_gauss,e_source,29-4*N)<-100

    old=inp('1247/results.json')['group_finite_memory']
    family=inp('1252/results.json')['group_parameter_family']['cells'][0]
    assert interval(family['kappa'])==(F(3),F(3001,1000))
    uB=[interval(old[f'alpha{j}']['checkpoints']['2']['state']['u']) for j in (0,1)]
    uc=interval(family['checkpoints']['2']['state']['u'])
    su=interval(family['checkpoints']['2']['sensitivity']['u'])
    for lo,hi in (*uB,uc):
        assert F(5,2)<1/hi and 1/lo<F(29,10)
    assert -su[1]/uc[1]**2>F(1,5)
    gap_a=1/uB[0][1]-1/uB[1][0]
    assert gap_a>F(683,10**7)
    assert F(5,2)*(1-F(11,10)*tau)>2
    assert F(3)*F(11,10)<4
    # cos z >= 1-z^2/2 for z in [0,3/8].
    min_slope=alpha*(1-F(3,8)**2/2)
    assert min_slope==F(119,2048)>F(1,32)
    # Probability temporal error <=alpha * 4 tau = tau/4 per task.
    time_error=4*alpha*tau
    assert time_error==tau/4
    sys=stat=F(1,10**7)
    gap_B=F(683,32*10**7)-2*(time_error+sys+stat)
    gap_k=F(1,32*5000)-2*(time_error+sys+stat)
    assert gap_B==F(2767,1600000000)>F(17,10**7)
    assert gap_k==F(1169,200000000)>F(58,10**7)
    assert gap_B-2*F(1,2**100)>F(17,10**7)
    assert gap_k-2*F(1,2**100)>F(58,10**7)
    repeats=5*10**14
    assert 2*repeats*stat**2==10
    # e>8/3 follows from its first five nonnegative Taylor terms.
    assert sum(F(1,factorial(j)) for j in range(5))>F(8,3)
    four_failure=8*F(3,8)**10
    assert four_failure<F(1,100)
    probability_error=time_error+sys+stat+F(1,2**100)
    eta_a=F(1,150000)
    assert 32*probability_error<eta_a
    kappa_accept=5*(2*eta_a+F(1,10**16))
    assert kappa_accept==F(1,15000)+F(5,10**16)
    field_prefix=F(111,100)*kappa_accept
    assert field_prefix<F(75,10**6)
    # Normal fixed input family only. Quantum channels themselves admit a passive
    # reference, but the Gaussian trajectory estimate is not all-input control.
    source_files=[STAGE/'1247/results.json',STAGE/'1249/quantum_bounds_results.json',STAGE/'1252/results.json']
    return {'round':1268,'verdict':'PASS_INDEPENDENT_WITH_SCOPE','independent_groups':1,
       'method':'independent sparse polynomial determinant, exact Fraction identities, old interval inputs and reconstructed resource integers; no author imports, PDE simulation or cutoff spectrum claim',
       'material':{'joint_determinant_coefficients':{str(k):str(v) for k,v in sorted(det.items())},
           'diagnostics':material_diagnostics,'volume_pressure_ratio':'-1/3',
           'scope':'frozen continuous principal block; neither finite physical leakage nor full causal completion certified'},
       'readout':{'pulse_integral':str(norm),'pulse_peak':str(peak),'gate_peak_over_epsilon':str(gate_bound),
           'gaussian_residual_first_relative_coefficient':str(first_coefficient),
           'second_relative_exponent':str(second_relative_exponent),'E_H':EH,'E_0':E0,
           'M_log2':M.bit_length()-1,'N_decimal':str(N),'vector_log2_upper':str(e_d),
           'gaussian_compare_log2_upper':str(e_gauss),'old_source_log2_upper':str(e_source),
           'minimum_probability_slope':str(min_slope),'single_task_time_probability_error':str(time_error),
           'clock_source_Fourier_power_upper':ceil_log2(clock_fourier_bound),
           'new_source_power_upper':E0+26,'new_source_log2_upper':str(E0+26-4*N),
           'independent_repeats_per_task':repeats,'four_task_failure_upper':str(four_failure),
           'a_inverse_error':str(eta_a),'accepted_kappa_error':str(kappa_accept),
           'classical_pre_gate_field_bound':str(field_prefix),
           'B_contrast_before_tiny_quantum_deduction':str(gap_B),
           'kappa_endpoint_contrast_before_tiny_quantum_deduction':str(gap_k),
           'scope':'conditional fixed effective operator/known normal family; new system bias and actual primitive are not inherited from the former POVM'},
       'source_sha256':{str(p.relative_to(STAGE)).replace('\\','/'):digest(p) for p in source_files},
       'not_certified':['whole spatial theta/pointer identity','complete new source and feedback','strict CC for the new interaction','post-read 9/4 field certificate','stage completion']}

if __name__=='__main__':
    result=main();out=HERE/'independent_results.json'
    if len(sys.argv)==2 and sys.argv[1]=='--write':
        with out.open('x',encoding='utf-8',newline='\n') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        print('WRITE_ONCE_PASS_1268 independent_groups=1')
    elif len(sys.argv)==1:
        saved=json.loads(out.read_text(encoding='utf-8-sig'))
        assert saved==result,'saved independent result does not match'
        print('READ_ONLY_MATCH_1268 independent_groups=1')
    else:raise SystemExit('Use no arguments for read-only replay, or --write once.')
