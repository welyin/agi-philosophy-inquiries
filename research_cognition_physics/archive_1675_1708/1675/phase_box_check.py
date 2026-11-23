"""Certify that the full two-width input box misses gamma mod 180 in [60,66.5].

Use the positive squared-minor expression of round 1658 and interval first
derivatives. If a partial has a fixed sign, move that coordinate to its
adverse face before mean-value bounding/subdivision. If there is a zero,
the adverse extremum cannot have a strict outward-certified sign.
"""
import argparse
import itertools
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OLD = runpy.run_path(str(HERE/'check.py'))
N = 10


def down(x):
    return np.nextafter(x, -np.inf)


def up(x):
    return np.nextafter(x, np.inf)


def prod(x,y):
    a,b=x
    c,d=y
    t1,t2,t3,t4=a*c,a*d,b*c,b*d
    return down(np.minimum(np.minimum(t1,t2),np.minimum(t3,t4))), up(np.maximum(np.maximum(t1,t2),np.maximum(t3,t4)))


def add(x,y):
    return down(x[0]+y[0]),up(x[1]+y[1])


class Jet:
    """Separate primal and ten partial derivatives, all enclosed by intervals."""
    def __init__(self,v,d):
        self.v,self.d=v,d

    def __add__(self,z):
        if isinstance(z,Jet):
            return Jet(add(self.v,z.v),add(self.d,z.d))
        return Jet(add(self.v,(z,z)),self.d)

    __radd__=__add__

    def __neg__(self):
        return Jet((-self.v[1],-self.v[0]),(-self.d[1],-self.d[0]))

    def __sub__(self,z):
        return self+-z

    def __rsub__(self,z):
        return (-self)+z

    def __mul__(self,z):
        if isinstance(z,Jet):
            return Jet(prod(self.v,z.v),add(prod(self.d, z.v),prod(self.v,z.d)))
        return Jet(prod(self.v,(z,z)),prod(self.d,(z,z)))

    __rmul__=__mul__

    def reciprocal(self):
        a,b=self.v
        assert np.all((a>0)|(b<0))
        v=(down(1/b),up(1/a))
        der=prod(v,v)
        return Jet(v,prod((-der[1],-der[0]),self.d))

    def __truediv__(self,z):
        return self*(z.reciprocal() if isinstance(z,Jet) else 1/z)

    def sq(self):
        a,b=self.v
        low=np.where((a<=0)&(b>=0),0,np.minimum(a*a,b*b))
        val=(down(low),up(np.maximum(a*a,b*b)))
        return Jet(val,prod((down(2*a),up(2*b)),self.d))

    def sqrt(self):
        a,b=self.v
        assert np.all(a>0)
        val=down(np.sqrt(a)),up(np.sqrt(b))
        return Jet(val,prod((down(.5/val[1]),up(.5/val[0])),self.d))


def roots(alpha, rising):
    a,b=alpha.v
    assert np.all(a>0) and np.all(b<4/27)
    def solve(target):
        lo=np.zeros_like(target)+(0 if rising else 1/3)
        hi=np.zeros_like(target)+(1/3 if rising else 1)
        for _ in range(52):
            x=(lo+hi)/2
            take=((x*(1-x)**2<target)==rising)
            lo,hi=np.where(take,x,lo),np.where(take,hi,x)
        return (hi+lo)/2
    lo=down(solve(a if rising else b)-1e-11)
    hi=up(solve(b if rising else a)+1e-11)
    assert np.all(lo>0) and np.all(hi<1)
    assert np.all(hi < 1/3) if rising else np.all(lo > 1/3)
    # Verify the root enclosure by outward evaluation, independently of bisection precision.
    for end,target,lower in ((lo,a if rising else b,rising),(hi,b if rising else a,not rising)):
        om=(down(1-end),up(1-end))
        value=prod((end,end),prod(om,om))
        assert np.all(value[1]<target) if lower else np.all(value[0]>target)
    w=Jet((lo,hi),alpha.d)
    slope=(1-w)*(1-3*w)
    return Jet(w.v,prod(slope.reciprocal().v,alpha.d))


def residual(box, rising, sign):
    coords=[]
    for i in range(N):
        d=np.zeros((len(box),N))
        d[:,i]=1
        coords.append(Jet((box[:,i,0,None],box[:,i,1,None]),(d,d)))
    downmass,strange,bottom,p,u,t,ele,mu,tau,cd=coords
    aa,ss=(downmass/bottom).sq(),(strange/bottom).sq()
    cp,cu,ct,sd=[(1-x.sq()).sqrt() for x in (p,u,t,cd)]
    base=aa*cp.sq()+ss*p.sq()
    mix=(ss-aa)*p*cp
    A=base*cu.sq()+u.sq()
    h12=(cu*(ct*mix+t*u*(1-base)*cd),-cu*t*u*(1-base)*sd)
    h13=(cu*(-t*mix+ct*u*(1-base)*cd),-cu*ct*u*(1-base)*sd)
    def sqabs(z):
        return z[0].sq()+z[1].sq()
    def conjprod(v,w):
        return v[0]*w[0]+v[1]*w[1],v[1]*w[0]-v[0]*w[1]
    v2=[(-p*ct-cp*t*u*cd,-cp*t*u*sd),(cp*ct-p*t*u*cd,-p*t*u*sd),(t*cu,t*0)]
    v3=[(p*t-cp*ct*u*cd,-cp*ct*u*sd),(-cp*t-p*ct*u*cd,-p*ct*u*sd),(ct*cu,t*0)]
    weights=(ss,aa,aa*ss)
    minor=sum(v*sqabs(z) for v,z in zip(weights,v2))
    C,E=sqabs(h12)/A,sqabs(h13)/A
    B,F=aa*ss/minor,minor/A
    ap=[conjprod(v,k) for v,k in zip(v2,v3)]
    adj=[sum(w*z[j] for w,z in zip(weights,ap)) for j in range(2)]
    Z=sqabs(adj)/(A*minor)
    y=conjprod(h13,h12)
    signed=-6*(y[0]*adj[0]+y[1]*adj[1])/A.sq()
    K=A+B+9*C
    h33=aa*sqabs(v3[0])+ss*sqabs(v3[1])+sqabs(v3[2])
    L=Z+h33
    h=(downmass/bottom)*(strange/bottom)
    tr=ele.sq()+mu.sq()+tau.sq()
    D=ele*mu*tau/(tr*tr.sqrt())
    P=(ele.sq()*mu.sq()+ele.sq()*tau.sq()+mu.sq()*tau.sq())/tr.sq()
    ww=roots((D*K*L.sqrt()/h).sq(),rising)
    om=1-ww
    # Positive squared-minor form of 1658, with a signed interference term.
    r=D*ww.sqrt()/(h*L.sqrt())
    v=ww/L
    cf=(C*F).sqrt()
    jh=cp*p*ct*t*cu.sq()*u*sd*(ss-aa)*(1-ss)*(1-aa)
    er=signed/(6*cf)
    ei=jh/(A*cf)
    unit=A*B/K.sq()*om.sq()+(A*(Z+F)+B*h33)*r+(ei*v).sq()
    return unit+((9*C*F*r).sqrt()-sign*er*v).sq()-P


def domain():
    iv=OLD['IV']
    def marg(row):
        x,e=row
        return iv['I'].point(x)+iv['bounds'](-2,2)*iv['I'].point(e)
    cols=[marg(row) for row in OLD['qrows']]
    ls=[marg(row) for row in OLD['lrows']]
    vals=cols[:3]+[iv['trig'](x,True) for x in cols[3:]]+ls
    common=np.array([[down(v.values()[0]),up(v.values()[1])] for v in vals])
    # delta = gamma + asin(k sin gamma) on (0,pi). For gamma or pi-gamma
    # in the enlarged measured arc the angle correction is smaller than .001001
    # radians since k<.001. Cosine is one-Lipschitz and decreasing on (0,pi).
    k=iv['trig'](cols[3],False)*iv['trig'](cols[5],True)*iv['trig'](cols[4],True)/(iv['trig'](cols[3],True)*iv['trig'](cols[5],False))
    assert k.values()[1]<.001
    from fractions import Fraction
    pi=iv['bounds']('3.14159265358979323846','3.14159265358979323847')
    cosmin=iv['trig'](pi*iv['I'].point(Fraction(133,360)),False)
    cosmax=iv['trig'](pi*iv['I'].point(Fraction(1,3)),False)
    assert cosmin.values()[0]-.001001>.397 and cosmax.values()[1]<.501
    assert -cosmax.values()[1]-.001001>-.502 and -cosmin.values()[0]<-.398
    return [np.vstack([common, end]) for end in ((.397,.501),(-.502,-.398))]


def classify(dom, rising, sign):
    work=dom[None,:,:].copy()
    center=dom.mean(axis=1)
    value=residual(np.repeat(center[None,:,None],2,axis=2),rising,sign)
    positive=bool(value.v[0][0,0]>0)
    assert positive or value.v[1][0,0]<0
    terminal=0
    mingap=1.
    for depth in range(22):
        if not len(work):
            return [terminal,depth,mingap]
        batch=[]
        for box in np.array_split(work, math.ceil(len(work)/256)):
            v=residual(box,rising,sign)
            # Prove only the sign at the adverse extremum. A derivative of
            # definite sign pushes the extremum to one face of this box.
            for _ in range(3):
                fixed=box[:,:,0]==box[:,:,1]
                increasing=(v.d[0]>=0)&~fixed
                decreasing=(v.d[1]<=0)&~fixed&~increasing
                if not (increasing|decreasing).any():
                    break
                minimize=increasing if positive else decreasing
                maximize=decreasing if positive else increasing
                box[:,:,1]=np.where(minimize,box[:,:,0],box[:,:,1])
                box[:,:,0]=np.where(maximize,box[:,:,1],box[:,:,0])
                v=residual(box,rising,sign)
            # Enclose the value at the center plus interval gradient times displacement.
            mid=box.mean(axis=2)
            vc=residual(np.repeat(mid[:,:,None],2,axis=2),rising,sign)
            d=prod(v.d,(down(box[:,:,0]-mid),up(box[:,:,1]-mid)))
            lo,hi=vc.v[0][:,0],vc.v[1][:,0]
            for i in range(N):
                lo,hi=add((lo,hi),(d[0][:,i],d[1][:,i]))
            lo=np.maximum(lo,v.v[0][:,0])
            hi=np.minimum(hi,v.v[1][:,0])
            assert np.all(lo<=hi)
            good=lo>0 if positive else hi<0
            if np.any(good):
                terminal+=int(good.sum())
                mingap=min(mingap,float(np.maximum(lo[good],-hi[good]).min()))
            if np.all(good):
                continue
            todo=box[~good]
            # Split the coordinate dominating the certified mean-value width.
            index=np.argmax((d[1]-d[0])[~good],axis=1)
            j=np.arange(len(todo))
            left,right=todo.copy(),todo.copy()
            left[j,index,1]=right[j,index,0]=mid[~good][j,index]
            batch.extend([left,right])
        work=np.concatenate(batch) if batch else []
        if len(work)>32768:
            raise ValueError('unclosed cover',depth,len(work))
    raise ValueError('depth cap',len(work))


def main():
    cases={}
    checks=0
    for orientation,dom in zip(('positive','conjugate'),domain()):
        for frac in (.2,.5,.8):
            x=dom[:,0]+frac*(dom[:,1]-dom[:,0])
            quarks=np.r_[x[:3],np.arcsin(x[3:6])]
            actual=OLD['ratio_shapes'](x[-1],quarks,x[6:9])
            point=np.repeat(x[None,:,None],2,axis=2)
            for j,(rising,sign) in enumerate(itertools.product((True,False),(-1,1))):
                val=residual(point,rising,sign).v
                assert val[0][0,0]<=actual[j][0]<=val[1][0,0]
                checks+=1
        for rising,sign in itertools.product((True,False),(-1,1)):
            label=f'{orientation}_{"small" if rising else "large"}_{"minus" if sign<0 else "plus"}'
            cases[label]=classify(dom,rising,sign)
    return dict(excluded_gamma_mod180_deg=[60.,66.5], experimental_one_width_deg=[61.,65.5],
                maximum_transfer_displacement_for_one_width_deg=1.,
                cosdelta_enclosures=[[.397,.501],[-.502,-.398]],
                extra_central_and_noncentral_singular_value_checks=checks,
                terminal_boxes=sum(c[0] for c in cases.values()),
                branch_columns=['terminal_boxes','depth_counter','minimum_certified_absolute_shape_gap'],
                branches=cases)


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--check',action='store_true')
    args=cli.parse_args()
    out=main()
    if args.check:
        old=json.loads((HERE/'phase_box_results.json').read_bytes())
        for k,v in out.items():
            if k=='branches':
                for branch,certificate in v.items():
                    np.testing.assert_allclose(certificate,old[k][branch],rtol=1e-10,atol=0)
            else:
                assert v==old[k]
        print('PASS: outward phase cover, both CP conjugates and all four determinant branches')
    else:
        with (HERE/'phase_box_results.json').open('x',encoding='utf8') as f:
            f.write(json.dumps(out,indent=2)+'\n')
