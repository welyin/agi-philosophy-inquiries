"""Outward decimal quadrature plus analytic complex-disk and tail bounds.

This certifies an unnormalized sixth-derivative integral for the exact-decimal
parameter interpretation of exact_history_density.py, not a finite-time value.
No external numerical library is required.
"""
from decimal import Decimal as D, getcontext
from fractions import Fraction as Q
from math import comb,factorial
from pathlib import Path
import json
getcontext().prec=60
HERE=Path(__file__).resolve().parent

class Interval:
    __slots__=('lo','hi')
    def __init__(self,lo,hi=None):
        self.lo=D(lo);self.hi=D(lo if hi is None else hi)
    @staticmethod
    def rational(q):
        q=Q(q);v=D(q.numerator)/D(q.denominator)
        return Interval(v.next_minus(),v.next_plus())
    def __add__(self,o):
        if not isinstance(o,Interval):o=Interval.rational(o)
        return Interval((self.lo+o.lo).next_minus(),(self.hi+o.hi).next_plus())
    __radd__=__add__
    def __neg__(self):return Interval(-self.hi,-self.lo)
    def __sub__(self,o):return self+-o if isinstance(o,Interval) else self+(-Q(o))
    def __mul__(self,o):
        if not isinstance(o,Interval):o=Interval.rational(o)
        v=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi]
        return Interval(min(v).next_minus(),max(v).next_plus())
    __rmul__=__mul__
    def reciprocal(self):
        assert self.lo>0 or self.hi<0
        return Interval((1/self.hi).next_minus(),(1/self.lo).next_plus())
    def __truediv__(self,o):
        if not isinstance(o,Interval):o=Interval.rational(o)
        return self*o.reciprocal()
    def sqrt(self):
        assert self.lo>=0
        return Interval(self.lo.sqrt().next_minus(),self.hi.sqrt().next_plus())
    def exp(self):return Interval(self.lo.exp().next_minus(),self.hi.exp().next_plus())

def newton_cotes(n):
    # Exact rational integrals of Lagrange basis on [-1,1].
    x=[Q(2*j,n)-1 for j in range(n+1)];weights=[]
    for j,xj in enumerate(x):
        p=[Q(1)];den=Q(1)
        for k,xk in enumerate(x):
            if j==k:continue
            out=[Q(0)]*(len(p)+1)
            for m,a in enumerate(p):out[m]-=xk*a;out[m+1]+=a
            p=out;den*=xj-xk
        weights.append(sum(2*a/Q(m+1) for m,a in enumerate(p) if m%2==0)/den)
    for k in range(n+1):
        assert sum(w*z**k for w,z in zip(weights,x))==(Q(2,k+1) if k%2==0 else 0)
    return weights

def prepare(K=20):
    data=json.loads((HERE/'exact_history_density_results.json').read_text('utf8'))
    coef={(v['r']//2,v['s']//2):Q(int(v['numerator']),int(v['denominator'])) for v in data['coefficients']['6']}
    assert max(a+b for a,b in coef)<=6
    total=sum(abs(v) for v in coef.values());assert total<2000
    poly=[[Q(0)]*(K+1) for _ in range(7)]
    for k in range(1,K+1):
        sine=Q((-1)**(k+1)*2**(2*k-3),factorial(2*k)) if k>=2 else Q(1,4)
        for (a,b),value in coef.items():
            moment=sum(Q(2*(-1)**j*comb(a+1,j),2*(k+b+j)+1) for j in range(a+2))
            assert moment>0
            poly[a+b][k]+=value*sine*moment
    return [[Interval.rational(v) for v in row] for row in poly],total

def evaluate(x,poly):
    if not x:return Interval(0)
    r=Interval.rational(x);v=r*r;B=12*v/(6+v)
    g=Interval(0)
    for row in reversed(poly):
        f=Interval(0)
        for c in reversed(row):f=f*B+c
        g=g*v+f
    return v*v*(-v).exp()/(1+v/6).sqrt()*g

def run():
    K=20;n=20;cells=500;R=Q(10);h=R/(2*cells);delta=Q(1,2)
    poly,total=prepare(K);weights=newton_cotes(n)
    # Adjacent cells share endpoint values, accumulated weights remain rational.
    grid={}
    for cell in range(cells):
        for j,w in enumerate(weights):
            index=cell*n+j;grid[index]=grid.get(index,Q(0))+h*w
    assert sum(grid.values())==R
    answer=Interval(0)
    for index,w in grid.items():
        if w:answer=answer+evaluate(R*index/(cells*n),poly)*Interval.rational(w)
    # On complex disks |z-x|<=1/2, x in [0,10]:
    # |z|<11, |6+z²|>=23/4, |B|<25, inverse square-root<2,
    # |exp(-z²)|<2, angular sine polynomial <=exp(10)/8<3^10/8.
    M=Q(11**4)*2*2*(2*total*11**12)*Q(3**10,8)
    assert M<10**26
    ratio=h/delta
    quad_error=cells*h*(2+sum(abs(w) for w in weights))*M*ratio**(n+1)/(1-ratio)
    # Entire real sine remainder. 2 int r^(4+2d) exp(-r²) <= (d+2)!,
    # since Gamma(d+5/2) <= Gamma(d+3) for d>=0.
    sine_error=total*factorial(8)*Q(48**(K+1),8*factorial(2*K+2))
    # True integrand tail: |A|<=1/4, angular absolute integral<=2.
    # exp(-100)<2^-100, and int_10∞ r^(4+2d)e^-r²
    # <= exp(-100)*(d+2)! sum(100^j/j!)/(20).
    max_tail=factorial(8)*sum(Q(100**j,factorial(j)) for j in range(9))
    tail_error=total*Q(1,40)*Q(1,2**100)*max_tail
    analytic_error=quad_error+sine_error+tail_error
    error=Interval.rational(analytic_error)
    lower=(answer.lo-error.hi).next_minus();upper=(answer.hi+error.hi).next_plus()
    assert upper<0
    return dict(exact_parameter_interpretation='Frozen decimal table and Yukawa values as real rationals; exact u solve.',
        derivative_order=6,unnormalized_integral_interval=[str(lower),str(upper)],
        quadrature_interval=[str(answer.lo),str(answer.hi)],quadrature_decimal_precision=60,
        quadrature_nodes=len(grid),interpolation_degree=n,sine_series_terms=K,
        complex_disk_radius=str(delta),complex_disk_bound=str(M),
        quadrature_analytic_error_upper=str(Interval.rational(quad_error).hi),
        sine_series_error_upper=str(Interval.rational(sine_error).hi),
        tail_error_upper=str(Interval.rational(tail_error).hi),
        newton_cotes_absolute_weight_sum=str(sum(abs(w) for w in weights)),
        coefficient_l1_norm=str(total),strict_negative=True,
        normalized_derivative_has_same_sign=True,finite_time_interval_not_computed=True,
        connected_graph_not_certified=True,all_checks_passed=True)

if __name__=='__main__':
    r=run()
    with Path(__file__).with_name('certify_history_integral_results.json').open('x',encoding='utf8') as f:
        json.dump(r,f,ensure_ascii=False,indent=2)
    print(json.dumps({k:r[k] for k in ('unnormalized_integral_interval','quadrature_analytic_error_upper','sine_series_error_upper','strict_negative')}))
