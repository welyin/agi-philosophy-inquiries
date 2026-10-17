"""1283: exact Gaussian moments and finite phase-cost selection.
No old checker import. Default read-only; --write exclusively creates results.
"""
import argparse
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
DIM=8
ZERO=(0,)*DIM


def add(*ps):
    out={}
    for p in ps:
        for k,v in p.items():out[k]=out.get(k,F(0))+v
    return {k:v for k,v in out.items() if v}


def scale(p,c):return {k:v*c for k,v in p.items() if v*c}


def mul(p,q):
    out={}
    for k,v in p.items():
        for l,w in q.items():
            key=tuple(a+b for a,b in zip(k,l))
            out[key]=out.get(key,F(0))+v*w
    return {k:v for k,v in out.items() if v}


def variable(i):
    key=[0]*DIM;key[i]=1
    return {tuple(key):F(1)}


def moment(poly,s,t):
    total=F(0)
    variances=[s]*4+[t]*4
    for powers,co in poly.items():
        value=co
        for n,v in zip(powers,variances):
            if n%2:value=0;break
            value*=math.prod(range(1,n,2))*v**(n//2)
        total+=value
    return total


def fs(x):return str(x.numerator)+'/'+str(x.denominator)


def quat(q):
    a,b,c,d=q
    return np.array([[a,-b,-c,-d],[b,a,-d,c],[c,d,a,-b],[d,-c,b,a]],float)


def numerical_group_check():
    # Fixed points and amplitudes, not a scan or a finite-group truncation.
    x1=np.array([.7,-.4,.2,1.1]);x2=np.array([-.3,.8,.5,-.6])
    uq=np.array([1.,2.,-1.,3.]);uq/=np.linalg.norm(uq)
    rx=quat(uq)@x2;s=2.;t=1.;kap=.8;q=(s-t)/(4*s*t)
    zeta=float(x1@rx);D=float((x1@x1)*(x2@x2)-zeta*zeta)
    ell=3*q*zeta/4-q*q*D/4
    G=kap*kap*D/4;d=kap*(3*zeta/4-q*D/2)
    h=F(1,2000);hf=float(h)
    rows=[]
    for theta in (0.,.5,1.):
        summ=0j
        for a in range(1,4):
            ratios=[]
            for sign in (-1,1):
                qa=np.zeros(4);qa[0]=math.cos(hf/2);qa[a]=sign*math.sin(hf/2)
                zetah=float(x1@(quat(qa)@rx))
                dz=zetah-zeta
                ratios.append(np.exp(q*dz-1j*theta*kap*dz))
            summ-=(sum(ratios)-2)/(hf*hf)
        predicted=ell+theta*theta*G-1j*theta*d
        residual=abs(summ-predicted)
        assert residual<2e-6,(theta,residual)
        rows.append({'theta':theta,'second_difference_residual':float(residual)})
    # Physical Gauss invariance of both the real preparation and phase function.
    g1=np.array([2.,1.,0.,-1.]);g1/=np.linalg.norm(g1)
    g2=np.array([1.,-1.,2.,1.]);g2/=np.linalg.norm(g2)
    r1=quat(g1);r2=quat(g2)
    transformed=float((r1@x1)@(r1@quat(uq)@r2.T)@(r2@x2))
    assert abs(transformed-zeta)<1e-13
    return rows


def calculate():
    # Fixed before the first computation in proof section 5.
    s,t,kap=F(2),F(1),F(4,5)
    akin,b,theta=F(1,2),F(3,10),F(1,2)
    v=[variable(i) for i in range(8)]
    Y=add(*(mul(x,x) for x in v[:4]));Z=add(*(mul(x,x) for x in v[4:]))
    P=add(*(mul(v[i],v[i+4]) for i in range(4)))
    zeta=scale(add(Y,scale(Z,-1)),F(1,2))
    D=add(mul(Y,Z),scale(mul(P,P),-1))
    q=(s-t)/(4*s*t)
    ell=add(scale(zeta,3*q/4),scale(D,-q*q/4))
    G=scale(D,kap*kap/4)
    d=add(scale(zeta,3*kap/4),scale(D,-kap*q/2))
    exact={
        'D':moment(D,s,t),'D2':moment(mul(D,D),s,t),
        'zeta2':moment(mul(zeta,zeta),s,t),
        'zeta_D':moment(mul(zeta,D),s,t),
        'L':moment(ell,s,t),'L2':moment(mul(ell,ell),s,t)}
    expected={'D':12*s*t,'D2':360*s*s*t*t,
              'zeta2':6*(s*s+t*t)-8*s*t,'zeta_D':36*s*t*(s-t),
              'L':3*(s-t)**2/(16*s*t),
              'L2':9*(s-t)**2/(64*s*t)+45*(s-t)**4/(512*s*s*t*t)}
    assert exact==expected,(exact,expected)
    # Direct Wick evaluation of the squared covariant differential norm.
    quartic=moment(mul(G,G),s,t)
    quadratic_L2=moment(add(scale(mul(ell,G),2),mul(d,d)),s,t)
    assert quartic==45*kap**4*s*s*t*t/2
    assert quadratic_L2==9*kap*kap*(5*(s*s+t*t)-6*s*t)/16
    Ebase_without_W=akin*(1/s+1/t+F(1,2))+4*kap*t+b*exact['L']
    candidates=[]
    for c in (F(0),F(1,10)):
        A2=16*akin*kap*kap*t+b*moment(G,s,t)+c*quadratic_L2
        E0=Ebase_without_W+c*exact['L2']
        energies=[]
        for th in (F(0),theta,2*theta):
            direct_L2=moment(add(mul(add(ell,scale(G,th*th)),add(ell,scale(G,th*th))),
                                 scale(mul(d,d),th*th)),s,t)
            e=Ebase_without_W+16*akin*kap*kap*t*th*th+b*moment(G,s,t)*th*th+c*direct_L2
            assert e==E0+A2*th*th+c*quartic*th**4
            energies.append(e)
        residual=energies[2]-4*energies[1]+3*energies[0]
        assert residual==12*c*quartic*theta**4
        candidates.append({'c':fs(c),'A2':fs(A2),
                           'energies_minus_common_W':[fs(x) for x in energies],
                           'doubling_residual':fs(residual)})
    denominator=12*quartic*theta**4
    eta=F(1,10000);eta_R=8*eta
    assert quartic==F(4608,125)
    assert denominator==F(3456,125)
    assert eta_R/denominator==F(1,34560)
    # General positive highest moment: recurrence of the two chi-square factors.
    n=3
    rising=lambda a,n:math.prod(a+j for j in range(n))
    expected_top=(kap*kap*s*t)**n*rising(F(2),n)*rising(F(3,2),n)
    direct_top=moment(mul(mul(G,G),G),s,t)
    assert direct_top==expected_top>0
    sources=[ROOT/'archive_554_584/research_note_557.md',HERE.parent/'1282/proof.md']
    return {
        'round':1283,'verdict':'PASS_WITH_SCOPE','author_groups':1,'new_axioms':0,
        'method':'exact eight-coordinate Wick moments; fixed group-path second difference',
        'parameters':{'s':fs(s),'t':fs(t),'kappa':fs(kap),'a_kin':fs(akin),'b':fs(b),'theta':fs(theta)},
        'source_sha256':{str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        'exact_moments':{k:fs(v) for k,v in exact.items()},
        'L2_quadratic_coefficient':fs(quadratic_L2),'C4':fs(quartic),
        'candidates':candidates,
        'finite_error':{'each_energy_error':fs(eta),'residual_error':fs(eta_R),
                        'denominator':fs(denominator),'c_upper_if_observed_residual_zero':fs(eta_R/denominator)},
        'polynomial_example_highest_moment_n3':fs(direct_top),
        'group_derivative_diagnostic':numerical_group_check(),
        'scope':{'old_W_kept_as_common_finite_term':True,'normal_Gauss_preparation':True,
                 'controller_energy_or_actual_energy_readout_certified':False,
                 'new_phase_permission_is_explicit_effective_input':True,
                 'leading_quadratic_Taylor_alone_selects_c_zero':False,
                 'exact_open_window_selects_finite_polynomial_degree_at_most_one':True,
                 'finite_error_classification_for_arbitrary_polynomial':False,
                 'condition_applies_to_all_physical_phase_operations':False,
                 'gauge_group_or_b_or_W_selected':False}}


def equivalent(a,b):
    if isinstance(a,dict):return isinstance(b,dict) and a.keys()==b.keys() and all(equivalent(a[k],b[k]) for k in a)
    if isinstance(a,list):return isinstance(b,list) and len(a)==len(b) and all(equivalent(x,y) for x,y in zip(a,b))
    if isinstance(a,float):return isinstance(b,(int,float)) and abs(a-b)<=1e-11*max(1,abs(a),abs(b))
    return a==b


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=calculate();path=HERE/'results.json'
    if args.write:
        with path.open('x',encoding='utf-8',newline='\n') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2);stream.write('\n')
        print('PASS: exclusive first write')
    else:
        assert equivalent(result,json.loads(path.read_text('utf-8'))),'Saved result mismatch'
        print('PASS: READ_ONLY_MATCH')
    print('result_sha256='+hashlib.sha256(path.read_bytes()).hexdigest())
    print('C4='+result['C4']+'; residual(c=1/10)='+result['candidates'][1]['doubling_residual'])
    print('max_group_diagnostic_residual='+str(max(x['second_difference_residual'] for x in result['group_derivative_diagnostic'])))


if __name__=='__main__':main()
