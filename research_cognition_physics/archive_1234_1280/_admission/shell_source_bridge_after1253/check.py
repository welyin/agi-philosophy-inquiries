"""Exact-rational scalar admission diagnostic; not a shell quantum model.
Default mode compares saved results without writing. No existing research imported.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial,isqrt
import hashlib,json,sys

class I:
    def __init__(self,a,b=None):
        self.lo=F(a); self.hi=F(a if b is None else b)
        assert self.lo<=self.hi
    def __add__(self,o):
        o=cv(o);return I(self.lo+o.lo,self.hi+o.hi)
    __radd__=__add__
    def __neg__(self): return I(-self.hi,-self.lo)
    def __sub__(self,o):return self+-cv(o)
    def __rsub__(self,o):return cv(o)+-self
    def __mul__(self,o):
        o=cv(o);v=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi]
        return I(min(v),max(v))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=cv(o);assert not o.lo<=0<=o.hi
        return self*I(1/o.hi,1/o.lo)
    def __rtruediv__(self,o):return cv(o)/self
    def __pow__(self,n):
        assert n>=0
        r=I(1)
        for _ in range(n):r=r*self
        return r
def cv(x):return x if isinstance(x,I) else I(x)
def pseries(u,kind,N=12):
    """Positive series on 0<=u<=1/10 with first-omitted geometric tail."""
    u=F(u);assert 0<=u<=F(1,10)
    if kind=='F':
        term=lambda k: F(2*(2*k-1),factorial(2*k))*u**(2*k)
        s=sum((term(k) for k in range(2,N+1)),F(0))
        r=u*u/F((2*N+1)*(2*N+4))
    elif kind=='cosh':
        term=lambda k:u**(2*k)/factorial(2*k)
        s=sum((term(k) for k in range(N+1)),F(0))
        r=u*u/F((2*N+3)*(2*N+4))
    else:
        assert kind=='sinh'
        term=lambda k:u**(2*k+1)/factorial(2*k+1)
        s=sum((term(k) for k in range(N+1)),F(0))
        r=u*u/F((2*N+4)*(2*N+5))
    return I(s,s+term(N+1)/(1-r))
def fun(u,kind):return I(pseries(u.lo,kind).lo,pseries(u.hi,kind).hi)
def root(target):
    a,b=F(0),F(1,10)
    assert pseries(b,'F').lo>target>0
    for _ in range(90):
        mid=(a+b)/2; f=pseries(mid,'F')
        if f.hi<target:a=mid
        elif f.lo>target:b=mid
        else:raise AssertionError('Increase tail precision')
    return I(a,b)
def out(x):
    x=cv(x);scale=10**20
    a=x.lo*scale;b=x.hi*scale
    low=a.numerator//a.denominator
    high=-((-b.numerator)//b.denominator)
    def fmt(v):
        sign='-' if v<0 else '';v=abs(v)
        return sign+str(v//scale)+'.'+str(v%scale).zfill(20)
    return [fmt(low),fmt(high)]
def compute():
    scale=10**40;k=isqrt(19*scale**2)
    assert k*k<=19*scale**2<(k+1)**2
    sqrt19=I(F(k,scale),F(k+1,scale))
    X=F(10000); rows=[]; hs=[]; vs=[]
    for n in range(6):
        mu=18+n+2*sqrt19
        target=mu**2/(X*X)
        u=I(root(target.lo).lo,root(target.hi).hi)
        c=fun(u,'cosh');s=fun(u,'sinh'); d=c-1
        e=X*u*u/2
        Fu=2*u*d
        v=(s-u*c)/d
        hm=mu/(X*d)
        hmm=1/(X*d)-mu**2*s/(X**3*u*d**3)
        assert Fu.lo>0 and v.hi<0
        assert F(99,100)<hm.lo<hm.hi<1 and hmm.hi<0
        assert e.lo>0 and e.hi<mu.lo
        hs.append(e);vs.append(v)
        # Source relation: p_unreduced from B=exp(-u)-1+u.
        exp_minus=c-s;B=exp_minus-1+u;w=mu/X
        assert B.lo>0
        P=X*(w*w/B-B)/2
        V=X*(w*w/B+B)/2
        constraint=e-(V+mu*mu/(2*X)-P*u)
        assert constraint.lo<=0<=constraint.hi
        assert (1-u).lo>0
        rows.append(dict(n=n,mass=out(mu),u=out(u),ADM_energy=out(e),
                         radial_velocity_at_Pc0=out(v),energy_mass_derivative=out(hm),
                         energy_mass_second_derivative=out(hmm),implicit_jacobian=out(Fu),
                         original_constraint_residual=out(constraint)))
    second=[hs[i+2]-2*hs[i+1]+hs[i] for i in range(4)]
    assert all(x.hi<0 for x in second)
    diff=[vs[i+1]-vs[i] for i in range(5)]
    assert all(x.hi<0 for x in diff)
    # Positive unital point evaluation Q[f]=f(z0) I has zero commutator.
    # Bounded a(X)=1/2+(1/4)tanh((X-X0)/L), a'(X0)=1/(4L).
    L=F(1000); residual=-vs[0]/(4*L)
    assert residual.lo>F(1,100000)
    return dict(status='PASS',kind='EXACT_RATIONAL_ADMISSION_DIAGNOSTIC',
                new_scientific_rounds=0,source_scope='Classical reduced symbol at Pc=0 only; no finite quantum propagation certified.',
                X=str(X),rational_bisections_per_endpoint=90,positive_series_last_index=12,
                rows=rows,energy_second_differences=[out(x) for x in second],
                neighboring_velocity_differences=[out(x) for x in diff],
                positive_unital_point_map=dict(commutator=0,transport_residual_absolute=out(residual),
                    conclusion='Positive unital quantization alone does not match even this one area-derived readout derivative.'),
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
def main():
    result=compute();path=Path(__file__).with_name('results.json')
    if '--write-once' in sys.argv:
        with path.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        print('FIRST_WRITE_PASS_SHELL_ADMISSION')
    else:
        assert json.loads(path.read_text('utf-8'))==result
        print('READ_ONLY_MATCH_SHELL_ADMISSION')
if __name__=='__main__':main()
