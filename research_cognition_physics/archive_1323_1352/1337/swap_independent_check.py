"""Independent finite partial-SWAP construction for round 1337. NumPy only."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import math
import numpy as np


def prepare(r, theta):
    n=len(r)
    ident=np.eye(n)
    involution=np.transpose(np.eye(n*n).reshape(n,n,n,n),(1,0,2,3)).reshape(n*n,n*n)
    w=(math.cos(theta)*np.eye(n*n)-1j*math.sin(theta)*involution)@np.kron(r[:,None],ident)
    q,_=np.linalg.qr(w,mode="complete")
    return w,q[:,n:]


def spectrum(r,chi,theta,t,M):
    w,complement=prepare(r,theta)
    p=np.kron(chi[:,None],np.eye(len(r)))
    a=p.conj().T@w
    leading=np.linalg.svd(a,compute_uv=False)[::-1]
    mass=np.column_stack([t*p,M*complement])
    full=np.linalg.svd(mass,compute_uv=False)[::-1]
    return leading,full[:len(r)]


def ratio(y):
    return float(y[0]*y[-1]/(y[1]**2))


def historic():
    data={"u":[[.00085,.00127,.00177],[.535,.619,.703],[168.7,171.7,174.7]],
          "d":[[.00171,.00290,.00414],[.040,.055,.071],[2.8,2.89,2.98]]}
    bounds={}
    for f,(m1,m2,m3) in data.items():
        bounds[f]=[m1[0]*m3[0]/m2[2]**2,m1[1]*m3[1]/m2[1]**2,m1[2]*m3[2]/m2[0]**2]
    return data,bounds


def interior_example():
    target=1.04
    theta=math.acos(1/target)
    r=np.eye(3)[:,0]
    w,_=prepare(r,theta)
    def one(overlap):
        chi=np.array([overlap, math.sqrt(1-overlap**2),0.])
        a=np.kron(chi[:,None],np.eye(3)).conj().T@w
        return chi,np.linalg.svd(a,compute_uv=False)[::-1]
    examples={}
    for f,low,high in [('u',.00175,173.),('d',.00180,2.83)]:
        middle=math.sqrt(low*high/target)
        left,right=0.,1.
        for _ in range(64):
            z=(left+right)/2
            _,y=one(z)
            if y[1]/y[2] < middle/high:
                left=z
            else:
                right=z
        z=(left+right)/2
        _,y=one(z)
        t=high/y[2]
        assert np.allclose(t*y,[low,middle,high],rtol=1e-10,atol=1e-12)
        examples[f]=dict(z=z,t=float(t),leading_GeV=(t*y).tolist())
    M=20*max(d['t'] for d in examples.values())
    data,_=historic()
    for f,d in examples.items():
        chi,_=one(d['z'])
        _,actual=spectrum(r,chi,theta,d['t'],M)
        assert all(lo<x<hi for (lo,_,hi),x in zip(data[f],actual))
        d['complete_GeV']=actual.tolist()
        d['complete_curvature']=ratio(actual)
    return dict(cos_theta=math.cos(theta),M_GeV=M,departments=examples)


def compute():
    rng=np.random.default_rng(133703)
    def random_vec(n):
        v=rng.normal(size=n)+1j*rng.normal(size=n)
        return v/np.linalg.norm(v)
    asym_error=0.; spectator_error=0.; finite_spectator_error=0.; excess=0.; number=0
    for n in [3,4,6]:
        for theta in [0.,.33,.81,1.2,2.3,-.94]:
            r=random_vec(n);chi=random_vec(n)
            c=abs(math.cos(theta)); z=abs(np.vdot(chi,r)); b=c*z
            for t,M in [(.08,3.),(.9,.9),(2.,.5)]:
                leading,full=spectrum(r,chi,theta,t,M)
                asym_error=max(asym_error,abs(ratio(leading)-1/c))
                spectator_error=max(spectator_error,np.max(abs(leading[1:-1]-b)))
                finite_spectator_error=max(finite_spectator_error,float(np.ptp(full[1:-1])))
                k=min(t/M,M/t)
                R=ratio(full)
                excess=max(excess,1/c-R,R-(1+k*k)/c)
                number+=1
    assert max(asym_error,spectator_error,finite_spectator_error,excess)<1e-9
    data,bounds=historic()
    example=interior_example()
    return dict(finite_cases=number, max_leading_curvature_error=asym_error,
                max_leading_middle_mass_error=spectator_error,max_complete_middle_mass_spread=finite_spectator_error,
                max_finite_two_sided_bound_excess=excess, historical_GeV=data,
                historical_curvature=[dict(sector=f,low=x[0],central=x[1],high=x[2]) for f,x in bounds.items()],
                down_over_up_central=bounds['d'][1]/bounds['u'][1],rectangle_witness=example)


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,str):
        assert a==b
    else:
        assert np.isclose(a,b,rtol=3e-8,atol=3e-9),(a,b)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--check',action='store_true')
    args=p.parse_args()
    result=compute()
    destination=Path(__file__).with_name('swap_independent_results.json')
    if args.check:
        compare(result,json.loads(destination.read_text(encoding='utf8')))
    else:
        with destination.open('x',encoding='utf8') as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write('\n')
    print(json.dumps(result,allow_nan=False))
