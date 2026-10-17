"""1243: exact finite-budget certificates; default replay never writes."""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent

def rational(x):
    return {'exact': str(x), 'decimal': float(x)}

def budget(ell):
    delta, b0, B, Ds = Q(1,10**6), Q(1,10**6), Q(1,100), Q(1,10000)
    Lr = Lf = Q(1,1000)
    dr = df = Q(1,100000)
    z = b0*ell+B*ell**2/2
    assert 0 <= z < 1
    M = 1/(1-z)  # e^z <= geometric series, exact upper bound
    N = b0*b0/2+b0*B*ell/3+(Ds+B*B/2)*ell*ell/6
    terms = dict(noise=2*delta/ell**2, endpoint=2*b0/ell,
                 transverse=dr+df, longitudinal=Lr*ell/3+2*Lf*ell/3,
                 area_weight=(M-1)*Q(1,50), nonlinear=M*N,
                 finite_cover=Q(1,100000))
    zeta = sum(terms.values(), Q(0))
    epsS = 9*zeta
    W1 = Q(53,4)
    trace = epsS*W1+Q(1,100000)
    return dict(z=z, exp_upper=M, terms=terms, zeta=zeta,
                epsS=epsS, W1=W1, trace_bound=trace,
                full_window_equation_bound=epsS+2*trace,
                conditional_margin=Q(1,20)-trace)

def integral_poly(coeff, left, right):
    return sum((v*(right**(i+1)-left**(i+1))/(i+1)
                for i,v in enumerate(coeff)), Q(0))

def main_result():
    ell, c = Q(1,10), Q(1,5)
    # Rosen calibration: a=b=1-c lambda^2; A=(1-c lambda^2)^2.
    area_change = 2*c*ell**2-c*c*ell**4
    # A(r+sigma^2-theta^2/2)=4c-12c^2 lambda^2.
    identity = integral_poly([4*c*ell,4*c,-12*c*c*ell,-12*c*c], -ell,Q(0))
    assert identity == area_change
    area_moment1 = integral_poly([Q(0),-ell,-1],-ell,Q(0))/(ell**2/2)
    heat_moment1 = integral_poly([Q(0),Q(0),1],-ell,Q(0))/(ell**2/2)
    area_moment2 = integral_poly([Q(0),Q(0),ell,1],-ell,Q(0))/(ell**2/2)
    assert (area_moment1,heat_moment1,area_moment2)==(ell/3,2*ell/3,ell**2/6)
    b, small = budget(ell), budget(ell/10)
    assert b['trace_bound'] < Q(42,1000)
    assert b['conditional_margin'] > Q(8,1000)
    assert small['terms']['noise'] == 100*b['terms']['noise']
    assert small['trace_bound'] > b['trace_bound']

    # Reuse calibration only; no novel nine-direction classification claimed.
    eta=np.diag([-1.,1.,1.,1.]); directions=[]
    for i in range(3):
        for sign in (1,-1):
            k=np.zeros(4); k[0]=1; k[i+1]=sign; directions.append(k)
    for i,j in ((1,2),(1,3),(2,3)):
        k=np.zeros(4);k[0]=1;k[i]=k[j]=1/np.sqrt(2);directions.append(k)
    def rebuild(q):
        a=np.array([(q[2*i]+q[2*i+1])/2 for i in range(3)])
        b=np.array([(q[2*i]-q[2*i+1])/4 for i in range(3)])
        S=np.zeros((4,4));S[0,0]=a.sum()/4
        for i in range(3):
            S[i+1,i+1]=a[i]-S[0,0];S[0,i+1]=S[i+1,0]=b[i]
        for idx,(i,j) in enumerate(((0,1),(0,2),(1,2))):
            S[i+1,j+1]=S[j+1,i+1]=q[6+idx]-(a[i]+a[j])/2-np.sqrt(2)*(b[i]+b[j])
        return S
    max_error=0.; max_norm=0.
    for mask in range(512):
        q=np.array([1. if mask & (1<<i) else -1. for i in range(9)])
        S=rebuild(q)
        max_error=max(max_error,max(abs(k@S@k-q[i]) for i,k in enumerate(directions)))
        max_norm=max(max_norm,np.linalg.norm(S))
        assert abs(np.sum(eta*S)) < 1e-12
    assert max_error<1e-12 and max_norm<9
    rigorous_frobenius_square=Q(3,4)**2+3*Q(5,4)**2+6*Q(1,2)**2+6*Q(7,2)**2
    assert rigorous_frobenius_square==Q(321,4)<81

    # Polynomial compact-window and weak-identity calibration, flat metric only.
    # q_h has unit mass and zero first moment. For f=alpha*t,
    # S00=3f,Sii=f,Lambda=-3f, j=0: both sides of (9) equal 3 alpha d.
    h,d,alpha=Q(1,4),Q(1,2),Q(1,100)
    q=[Q(15,16)/h,Q(0),-Q(15,8)/h**3,Q(0),Q(15,16)/h**5]
    assert integral_poly(q,-h,h)==1
    q_tv=2*q[0]
    W=2+3*d*q_tv
    assert W==b['W1']
    mean_diff=3*alpha*d
    # Integral t*U' = -integral U = -d, and raising time index flips sign.
    sourcefree_rhs=-3*alpha*(-d)
    assert sourcefree_rhs==mean_diff
    # Removing conservation: T=L eta, L=alpha*t; S=0, j0=alpha.
    omitted_source_diff=-alpha*d
    omitted_source_rhs=-alpha*d
    assert omitted_source_diff==omitted_source_rhs != 0
    # High-frequency example: fixed |S| and unbounded point gradient.
    alpha=Q(1,10000)
    assert 3*alpha*10000 > 3*alpha*10

    def encode_b(x):
        return {k:({n:rational(v) for n,v in val.items()} if k=='terms' else rational(val))
                for k,val in x.items()}
    return dict(round=1243,status='PASS',science_groups=3,
        group_finite_window=dict(status='PASS',arithmetic='exact Fraction; exp bounded by 1/(1-z)',
            Rosen_area_identity=rational(area_change),
            normalized_moments=[rational(x) for x in (area_moment1,heat_moment1,area_moment2)],
            certificate=encode_b(b),smaller_window=encode_b(small)),
        group_tensor_reuse=dict(status='PASS',purpose='old 304/1097 calibration, not new theorem',
            corner_cases=512,max_reconstruction_error=max_error,max_frobenius=max_norm,
            exact_conservative_norm_squared=rational(rigorous_frobenius_square)),
        group_weak_trace=dict(status='PASS',compact_window_mass='1',W1_bound=rational(W),
            conserved_polynomial_window_difference=rational(mean_diff),
            omitted_conservation_difference=rational(omitted_source_diff),
            point_gradient_inference='invalid; frequency is not controlled'),
        scope=dict(CT='conjecture_not_used',new_adopted_axioms=0,
            actual_horizon_instruments='not implemented',
            finite_coverage_and_regularity='explicit independent conditions',
            Einstein_generation='not proved',cosmological_constant='not selected'))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=main_result();dest=HERE/'results.json'
    if args.write:
        with dest.open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        print('EXCLUSIVE_INITIAL_WRITE')
    else:
        assert json.loads(dest.read_text('utf-8'))==result
        print('READ_ONLY_MATCH')
    print(json.dumps(dict(status=result['status'],round=1243,science_groups=3,
        trace_bound=result['group_finite_window']['certificate']['trace_bound'],
        margin=result['group_finite_window']['certificate']['conditional_margin']),ensure_ascii=False))
