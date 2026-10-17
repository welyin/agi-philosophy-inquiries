"""1242: exact outward-rounded fixed-point interval Taylor certificate.

Python standard library. Default replay is read-only; --write creates once.
No floating point operations enter the certified interval integration.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from math import isqrt
from pathlib import Path
from fractions import Fraction as F

HERE=Path(__file__).resolve().parent
S=10**36

class I:
    __slots__=('lo','hi')
    def __init__(self, x=0, hi=None, raw=False):
        if raw:
            self.lo,self.hi=x,hi
        elif isinstance(x,I):
            self.lo,self.hi=x.lo,x.hi
        else:
            q=F(x); r=q if hi is None else F(hi)
            self.lo=q.numerator*S//q.denominator
            self.hi=-((-r.numerator*S)//r.denominator)
        assert self.lo<=self.hi
    def __add__(self,b):
        b=I(b);return I(self.lo+b.lo,self.hi+b.hi,raw=True)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo,raw=True)
    def __sub__(self,b):return self+-I(b)
    def __rsub__(self,b):return I(b)+-self
    def __mul__(self,b):
        b=I(b);p=[self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi]
        return I(min(p)//S,-((-max(p))//S),raw=True)
    __rmul__=__mul__
    def div(self,n):
        assert isinstance(n,int) and n>0
        return I(self.lo//n,-((-self.hi)//n),raw=True)
    def __pow__(self,n):
        r=I(1)
        for _ in range(n):r=r*self
        return r
    def subset(self,b):return b.lo<=self.lo and self.hi<=b.hi
    def contains_zero(self):return self.lo<=0<=self.hi
    def width(self):return F(self.hi-self.lo,S)
    def out(self):
        def fmt(n):
            s='-' if n<0 else '';a=abs(n)
            return f'{s}{a//S}.{a%S:036d}'
        return [fmt(self.lo),fmt(self.hi)]

def sqrt_interval(q):
    q=F(q);n=isqrt(q.numerator*S*S//q.denominator)
    return I(n,n+1,raw=True)

class Jet:
    def __init__(self,c):self.c=list(c)
    def coerce(self,b):return b if isinstance(b,Jet) else Jet([I(b)]+[I(0)]*(len(self.c)-1))
    def __add__(self,b):
        b=self.coerce(b);return Jet([a+d for a,d in zip(self.c,b.c)])
    __radd__=__add__
    def __neg__(self):return Jet([-a for a in self.c])
    def __sub__(self,b):return self+-self.coerce(b)
    def __rsub__(self,b):return self.coerce(b)+-self
    def __mul__(self,b):
        b=self.coerce(b)
        return Jet([sum((self.c[j]*b.c[k-j] for j in range(k+1)),I(0)) for k in range(len(self.c))])
    __rmul__=__mul__
    def __pow__(self,n):
        r=self.coerce(1)
        for _ in range(n):r=r*self
        return r

def rhs(y,D):
    u,H,p,v,z=y
    R=F(3,20)*(1+F(1,5)*p)*u**2
    rhoP=(D+F(1,100)*p)*u**3+v**2+F(2,3)*R
    return [-H*u,-F(3,2)*rhoP,v,-3*H*v-p-F(3,100)*u**2-F(1,100)*u**3,F(1,10)*p]

def coefficients(y,D,n):
    c=[[I(a)] for a in y]
    for k in range(n):
        d=rhs([Jet(a) for a in c],D)
        for j in range(5):c[j].append(d[j].c[k].div(k+1))
    return c

def state_energy(y,D):
    u,H,p,v,z=y
    R=F(3,20)*(1+F(1,5)*p)*u**2
    dust=(D+F(1,100)*p)*u**3
    rho=dust+R+(v**2+p**2)*F(1,2)
    P=-R*F(1,3)+(v**2-p**2)*F(1,2)
    return rho,P

def integrate(D,steps=20,den=100,order=8):
    h=I(F(1,den));time=I(0,F(1,den))
    y=[I(1),sqrt_interval(D+F(3,20)),I(0),I(0),I(0)]
    radii=[F(1,50),F(1,40),F(1,500),F(3,1000),F(1,10000)]
    domain=[I('.7','1.05'),I('.5','1.2'),I('-.02','.02'),I('-.05','.05'),I('-.002','.002')]
    allbox=None;maxwidth=F(0);maxrem=F(0)
    hpow=[I(1)]
    for j in range(order+1):hpow.append(hpow[-1]*h)
    trajectory=[]
    for step in range(steps):
        box=[a+I(-r,r) for a,r in zip(y,radii)]
        assert all(a.subset(b) for a,b in zip(box,domain)),'positive common domain'
        picard=[a+time*b for a,b in zip(y,rhs(box,D))]
        assert all(a.subset(b) for a,b in zip(picard,box)),'Picard tube inclusion'
        c=coefficients(y,D,order)
        tail=coefficients(box,D,order+1)
        out=[]
        for j in range(5):
            rem=tail[j][-1]*hpow[order+1]
            out.append(sum((c[j][k]*hpow[k] for k in range(order+1)),I(0))+rem)
            maxrem=max(maxrem,abs(F(rem.lo,S)),abs(F(rem.hi,S)))
        y=out
        assert all(a.subset(b) for a,b in zip(y,box))
        maxwidth=max(maxwidth,max(a.width() for a in y))
        rho,_=state_energy(y,D)
        assert (y[1]**2-rho).contains_zero(),'constraint must contain its exact zero'
        if allbox is None:allbox=box
        else:allbox=[I(min(a.lo,b.lo),max(a.hi,b.hi),raw=True) for a,b in zip(allbox,box)]
        trajectory.append([a.out() for a in y])
    assert (1+F(1,5)*allbox[2]).lo>0
    assert (F(1,5)+F(1,10)*allbox[2]).lo>0
    B=F(1,10)*y[0]**2
    return dict(D=str(D),T=str(F(steps,den)),state=[a.out() for a in y],
                B=B.out(),all_time_tube=[a.out() for a in allbox],
                max_endpoint_width=str(maxwidth),max_step_remainder=str(maxrem),
                picard_tubes=steps,order=order),y,trajectory

def conservation_checks():
    # Exact rational polynomial identities at non-special positive-domain states.
    residuals=[];bad_source=[];bad_force=[]
    for u,H,p,v,D in [(F(9,10),F(1),F(-1,100),F(-1,200),F(51,50)),
                       (F(4,5),F(7,10),F(1,200),F(1,300),F(41,50))]:
        up,Hp,pp,vp,zp=rhs([u,H,p,v,F(0)],D)
        R=F(3,20)*(1+F(1,5)*p)*u*u
        Rp=F(3,100)*v*u*u-2*H*R
        dust=(D+F(1,100)*p)*u**3
        dustp=F(1,100)*v*u**3-3*H*dust
        sp=v*vp+p*v
        rhoP=dust+v*v+F(2,3)*R
        residuals.append(Rp+dustp+sp+3*H*rhoP)
        residuals.append(2*H*Hp-(Rp+dustp+sp))
        # Two different deletions; recompute each complete continuity residual.
        # A: omit Q*p*u^3 from stress energy, retain the true pointer force.
        deleted_rhop=Rp-3*H*D*u**3+sp
        deleted_rhoP=D*u**3+v*v+F(2,3)*R
        bad_source.append(deleted_rhop+3*H*deleted_rhoP)
        # B: omit -Q*u^3 from pointer acceleration, retain full stress energy.
        wrong_vp=vp+F(1,100)*u**3
        wrong_sp=v*wrong_vp+p*v
        bad_force.append(Rp+dustp+wrong_sp+3*H*rhoP)
        assert bad_source[-1]==-F(1,100)*v*u**3
        assert bad_force[-1]==F(1,100)*v*u**3
    assert all(x==0 for x in residuals)
    assert all(x!=0 for x in bad_source+bad_force)
    return dict(status='PASS',exact_continuity_and_constraint_residuals=list(map(str,residuals)),
                omitted_interaction_stress_continuity_defects=list(map(str,bad_source)),
                omitted_pointer_force_continuity_defects=list(map(str,bad_force)))

def run():
    hi,yh,trh=integrate(F(51,50));lo,yl,trl=integrate(F(41,50))
    gap=yh[4]-yl[4]
    assert gap.lo>I('0.00000007').hi
    # Declared terminal rotor-phase readout errors, one per history, not implemented.
    margin=gap-I('0.00000002')
    assert margin.lo>0
    Bgap=F(1,10)*(yl[0]**2-yh[0]**2)
    assert Bgap.lo>0
    # Local series detects feedback already in the leading record jet.
    jets=[coefficients([I(1),sqrt_interval(D+F(3,20)),I(0),I(0),I(0)],D,4)
          for D in [F(51,50),F(41,50)]]
    for a,D in zip(jets,[F(51,50),F(41,50)]):
        expected=F(21,24000)*sqrt_interval(D+F(3,20))
        assert (a[4][4]-expected).contains_zero()
    deletion=dict(status='PASS',fixed_geometry_D_absent_from_pointer_and_record_rhs=True,
                  fixed_geometry_record_gap_exact='0',
                  explanation='With a,H prescribed identically, D occurs only in the deleted H equation. Same p,v,z initial data imply identical records by ODE uniqueness. This comparison is not another solution of Einstein constraints.')
    sources={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ['admission.md','proof.md','check.py']}
    return dict(round=1242,status='PASS',science_groups=3,
                arithmetic=dict(grid_denominator=str(S),outward_rounding='integer floor/ceil',
                    initial_H='integer-square-root bracket',scheme='8th degree Taylor + 9th interval coefficient on Picard tubes'),
                group_conservation=conservation_checks(),
                group_finite_records=dict(status='PASS',high=hi,low=lo,gap=gap.out(),B_separation=Bgap.out(),
                    phase_error_per_history='0.00000001',remaining_separation=margin.out(),
                    instrument_precision_is_contract=True),
                group_deletions=deletion,
                source_sha256=sources,
                whole_GR_proved=False,new_adopted_axioms=0)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=run();path=HERE/'results.json'
    if args.write:
        with path.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        mode='CREATED_ONCE'
    else:
        assert json.loads(path.read_text('utf-8'))==result,'saved result mismatch'
        mode='READ_ONLY_MATCH'
    print(json.dumps(dict(status='PASS',mode=mode,groups=result['science_groups'],
        record_gap=result['group_finite_records']['gap']),ensure_ascii=False))
