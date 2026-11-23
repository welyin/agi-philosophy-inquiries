"""Finite physical-representation and rational-window checks for round 1348."""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path

import numpy as np


def adj(a): return a.conj().T


def compare(a,b):
    np.testing.assert_allclose(a,b,atol=4e-11,rtol=4e-11)


def rep(cases):
    # (twice-spin, twice-T3 at the reference electric charge). Hypercharge
    # is chosen separately in each summand to put that magnetic weight at b charge.
    dims=[j+1 for j,m in cases]
    n=sum(dims)
    plus=np.zeros((n,n),complex)
    t3=np.zeros(n)
    charge=np.zeros(n)
    origin=0
    for (j,m),d in zip(cases,dims):
        for k in range(d):
            eig=-j/2+k
            t3[origin+k]=eig
            charge[origin+k]=eig-m/2
            if k<d-1: plus[origin+k+1,origin+k]=math.sqrt((d-1-k)*(k+1))
        origin+=d
    compare(plus@adj(plus)-adj(plus)@plus,2*np.diag(t3))
    return plus,np.diag(t3),charge


def frame(rng,n,k):
    u,_=np.linalg.qr(rng.normal(size=(n,k))+1j*rng.normal(size=(n,k)))
    return u[:,:k]


def bounds(h,e,b):
    gh=adj(h)@h
    ge=adj(e)@e
    compare(gh-ge,b)
    ev=np.linalg.eigvalsh(b)
    tol=1e-9
    rp=int(np.sum(ev>tol))
    rn=int(np.sum(ev< -tol))
    assert np.linalg.matrix_rank(h,tol)>=rp
    assert np.linalg.matrix_rank(e,tol)>=rn
    assert np.trace(gh).real+tol>=np.maximum(ev,0).sum()
    assert np.trace(ge).real+tol>=np.maximum(-ev,0).sum()
    return rp,rn


def tests():
    rng=np.random.default_rng(1348)
    q,t,charge=rep([(1,-1),(1,-1),(1,1),(2,0),(3,-1),(0,0)])
    eye=np.eye(len(q))
    base=eye[:,charge==0]; upper=eye[:,charge==1]
    hist={}
    for _ in range(100):
        x=base@frame(rng,base.shape[1],3)
        u=upper@frame(rng,upper.shape[1],2)
        a=adj(u)@q@x
        h=(eye-u@adj(u))@q@x
        e=adj(q)@x
        b=-2*adj(x)@t@x-adj(a)@a
        key=str(bounds(h,e,b))
        hist[key]=hist.get(key,0)+1
    return hist


def historical():
    i3_upper=F(-490+15,1000)
    i3_lower=F(-490-12,1000)
    form=F(80,100)
    vcb=(F(351,10000)+F(19+18+8,10000))/form
    r2=F(5,100)
    obs_sq=vcb**2*(1+r2)
    low=-2*i3_upper-obs_sq
    return dict(vcb_upper=float(vcb),open_squared_upper=float(obs_sq),
                unseen_lower_fraction=str(low),unseen_lower=float(low),
                single_partner_amplitude_lower=math.sqrt(low),
                single_no_exotic_amplitude_upper=math.sqrt(-2*i3_lower),
                independent_fLV_Vtb_box=[.96,1.08])


def wrong_sign():
    q,t,charge=rep([(1,-1),(1,1),(0,0)])
    # First doublet can be resolved into u and c by adding an inert up charge;
    # its open norm below is unchanged by that basis embedding.
    native=-.490
    sin2=.2327
    shift=-sin2/3
    gl=native-shift
    mirror=-gl+shift
    g_right=-.028-shift
    gram=lambda g:(g*g+g_right*g_right,g*g-g_right*g_right)
    compare(gram(gl),gram(-gl))
    w_minus=.002
    w_plus=2*mirror+w_minus
    w_zero=1-w_plus-w_minus
    x=np.eye(len(t))[:,charge==0]@np.sqrt([w_minus,w_plus,w_zero])
    compare(np.vdot(x,t@x),mirror)
    compare(np.linalg.norm(q@x)**2,w_minus)
    compare(np.linalg.norm(adj(q)@x)**2,w_plus)
    return dict(z_only_reflected_left=mirror,open_up_squared=w_minus,
                extra_down_squared=w_plus,additional_up_squared=0.)


def run():
    historical_bound=historical()
    mixed=tests()
    assert any(a>0 and b>0 for a,b in (map(int,k[1:-1].split(", ")) for k in mixed))
    # Realizable rank two with two input lower doublets and one observed
    # coherent raising direction contained partly in a third doublet.
    q,t,ch=rep([(1,-1)]*3)
    eye=np.eye(6)
    x=eye[:,ch==0][:,:2]
    u=eye[:,ch==1]@np.array([[1.],[0.],[1.]])/np.sqrt(2)
    a=adj(u)@q@x; h=(eye-u@adj(u))@q@x
    b=-2*adj(x)@t@x-adj(a)@a
    compare(np.linalg.eigvalsh(b),[.5,1.])
    assert bounds(h,adj(q)@x,b)==(2,0)
    hgram=np.array([[2,math.sqrt(2)],[math.sqrt(2),1]])
    b0=np.diag([1,-1])
    assert min(np.linalg.eigvalsh(hgram-b0))>=-1e-12
    assert min(np.linalg.eigvalsh(hgram-np.diag([1,0])))<0
    return dict(round=1348,historical=historical_bound,mirror=wrong_sign(),
                general_representations=mixed,two_hidden_positive_eigenvalues=[.5,1.],
                pole_mass_is_input=False)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--write",action="store_true")
    args=p.parse_args()
    result=run()
    target=Path(__file__).with_name("results.json")
    if args.write:
        with target.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2); f.write("\n")
    else:
        saved=json.loads(target.read_text(encoding="utf-8"))
        assert saved==result
    print(json.dumps(result["historical"]))

