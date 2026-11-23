"""Inverse mass branches and continuous-box certificates for a geometric singular spectrum."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


class Box:
    """Interval arithmetic; widen elementary calculations beyond float round-off."""
    def __init__(self,lo,hi=None):
        self.lo,self.hi=lo,lo if hi is None else hi

    @staticmethod
    def cast(x):
        return x if isinstance(x,Box) else Box(x)

    @staticmethod
    def padded(lo,hi):
        pad=2e-14*max(abs(lo),abs(hi),1e-200)
        return Box(lo-pad,hi+pad)

    def __add__(self,b):
        b=self.cast(b)
        return self.padded(self.lo+b.lo,self.hi+b.hi)
    __radd__=__add__

    def __neg__(self):
        return Box(-self.hi,-self.lo)

    def __sub__(self,b):
        return self+(-self.cast(b))

    def __rsub__(self,b):
        return self.cast(b)+(-self)

    def __mul__(self,b):
        b=self.cast(b)
        edges=[u*v for u in (self.lo,self.hi) for v in (b.lo,b.hi)]
        return self.padded(min(edges),max(edges))
    __rmul__=__mul__

    def __truediv__(self,b):
        b=self.cast(b)
        assert b.lo>0
        return self*self.padded(1/b.hi,1/b.lo)

    def __rtruediv__(self,b):
        return self.cast(b)/self

    def __pow__(self,exponent):
        assert self.lo>=0 and exponent>0
        return self.padded(self.lo**exponent,self.hi**exponent)

    def center(self):
        return (self.lo+self.hi)/2


def intervals(items,n):
    return [Box(x-n*err,x+n*err) for x,err in items]


def calibration(cfg,scale,n):
    if scale=="MZ":
        source=read(HERE/cfg["low"])
        u,d=source["yukawa"]["up"],source["yukawa"]["down"]
        angles=intervals(source["quark_sines_and_phase"],n)
    else:
        source=read(HERE/cfg["high"])["high_scale"]
        u,d=source["up"],source["down"]
        angles=intervals(source["angles_theta12_theta13_theta23_delta"],n)
        for i in range(3):
            angles[i]=Box.padded(math.sin(angles[i].lo),math.sin(angles[i].hi))
    return intervals(u,n),intervals(d,n),angles


def ckm_probabilities(angles):
    a,b,c,phase=angles
    A,B,C=[(1-k**2)**.5 for k in (a,b,c)]
    cosdelta=Box.padded(math.cos(phase.hi),math.cos(phase.lo))
    interference=2*a*A*c*C*b*cosdelta
    return [[(A*B)**2,(a*B)**2,b**2],
            [(a*C)**2+(A*c*b)**2+interference,(A*C)**2+(a*c*b)**2-interference,(c*B)**2],
            [(a*c)**2+(A*C*b)**2-interference,(A*c)**2+(a*C*b)**2+interference,(C*B)**2]]


def invariants(c,t,d,s,b,p):
    small_weights=[(d/j)**2 for j in (d,s,b)]
    big_weights=[(j/b)**2 for j in (d,s,b)]
    small_weights[0]=big_weights[2]=Box(1)
    small=[sum(x*y for x,y in zip(row,small_weights)) for row in p]
    big=[sum(x*y for x,y in zip(row,big_weights)) for row in p]
    G=big[2]+(c/t)**2*big[1]
    T=t**2*b**2*G
    det=(c*t*d*s*b)**2
    Pcoeff=det*(small[1]/c**2+small[2]/t**2)/d**2
    quartic=(b**2*big[0])/T
    cubic=Pcoeff/(det**(1/3)*T)
    h0=(c**2*s**2/(t*b*d))**(2/3)*small[0]/G
    return quartic,cubic,h0,small,G


def roots_certificate(coeff):
    quartic,cubic,h0=coeff[:3]
    def sign(x):
        h=x**(2/3)
        return h-h0+quartic*h**4-cubic*h**3
    # Exactly three positive roots by these three sign changes and Descartes' rule.
    bins=[(1e-7,.0001,1),(.0001,10.,-1),(10.,1e8,1)]
    enclosures=[]
    for left,right,direction in bins:
        assert (direction*sign(left)).hi<0 and (direction*sign(right)).lo>0
        edges=[]
        for test in ("all_negative","all_positive"):
            lo,hi=left,right
            for _ in range(70):
                mid=math.sqrt(lo*hi)
                box=direction*sign(mid)
                if (box.hi<0 if test=="all_negative" else box.lo<=0):
                    lo=mid
                else:
                    hi=mid
            edges.append(lo if test=="all_negative" else hi)
        enclosures.append(edges)
    return enclosures


def central_svd(up,down,angles):
    a,b,c,delta=[z.center() for z in angles]
    A,B,C=np.sqrt(1-np.array([a,b,c])**2)
    ph=np.exp(1j*delta)
    v=np.array([[A*B,a*B,b/ph],[-a*C-A*c*b*ph,A*C-a*c*b*ph,c*B],
                [a*c-A*C*b*ph,-A*c-a*C*b*ph,C*B]])
    vals=np.linalg.svd(np.array(up)[:,None]*v*np.array(down),compute_uv=False)
    # Use the exact determinant to avoid losing the tiny last singular value.
    return float(vals[1]**3/(np.prod(up)*np.prod(down)))


def evaluate_scale(cfg,scale):
    outputs={}
    for n in (0,1,2):
        (x,c,t),(d,s,b),angles=calibration(cfg,scale,n)
        p=ckm_probabilities(angles)
        coeff=invariants(c,t,d,s,b,p)
        roots=roots_certificate(coeff)
        h=coeff[3]
        rho_upper=(c**2*s**2/(t*b*d*x)*((h[0]+x**2*(h[1]/c**2+h[2]/t**2))/p[2][2])**1.5).hi
        assert roots[0][1]<x.lo and roots[1][0]>c.hi and rho_upper<1
        out={"ordered_positive_root_ranges":roots,"held_out_yu":[x.lo,x.hi],
             "minimum_fractional_reduction_of_yu_if_only_it_is_revised":1-roots[0][1]/x.lo,
             "upper_s2_squared_over_s1_s3_at_measured_mass":rho_upper}
        if n==0:
            up=[z.center() for z in (x,c,t)]
            down=[z.center() for z in (d,s,b)]
            out["direct_geometric_ratio"]=central_svd(up,down,angles)
            for window in roots:
                up[0]=sum(window)/2
                assert abs(central_svd(up,down,angles)-1)<2e-7
        outputs[str(n)]=out
    return outputs


def eta_budget(cfg,relative_half_width):
    (unused,c,t),(_,s,b),angles=calibration(cfg,"MZ",2)
    p=ckm_probabilities(angles)
    sc,stat,sys=cfg["ms_over_average"]
    S=Box(sc-2*math.hypot(stat,sys),sc+2*math.hypot(stat,sys))
    transfer=Box(1-relative_half_width,1+relative_half_width)
    SH=S*transfer
    # Positive masses with r_QCD=u/d<1, hence r_H<=1+residual.
    d=2*s/(SH*Box(1,1+transfer.hi))
    xmax=2*s.hi/SH.lo
    _,cubic,_,small,G=invariants(c,t,d,s,b,p)
    slack=cubic.hi*xmax**(4/3)
    assert 0<slack<.05
    upper=(SH**2*c**2/(4*t*b)*(small[0]/G)**1.5).hi/(1-slack)**1.5
    assert 0<upper<.25
    r_limit=2*upper/(1-2*upper+math.sqrt(1-4*upper))
    assert 1/r_limit>transfer.hi
    rq=r_limit/transfer.lo
    assert rq<1
    q_limit=.5*math.sqrt((S.hi**2-1)*(1+rq)/(1-rq))*(1+1e-12)
    return {"matching_half_width":relative_half_width,"rH_over_one_plus_rH_squared_max":upper,
            "rQCD_upper":rq,"Q_upper":q_limit,"small_branch_slack":slack}


def run():
    cfg=read(HERE/"inputs.json")
    variants={scale:evaluate_scale(cfg,scale) for scale in cfg["scales"]}
    eta=eta_budget(cfg,cfg["QCD_to_MZ_relative_matching"][1]-1)
    qc,err=cfg["eta_q"]
    assert eta["Q_upper"]<qc-2*err
    return {"per_scale":variants,"independent_eta_holdout":eta,
            "registered_eta_working_interval":[qc-2*err,qc+2*err]}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    data=run()
    destination=HERE/"results.json"
    if args.write_new:
        with destination.open("x",encoding="utf-8") as stream:
            json.dump(data,stream,indent=2)
            stream.write("\n")
    elif args.check:
        def equal(a,b):
            if isinstance(a,dict):
                return a.keys()==b.keys() and all(equal(a[k],b[k]) for k in a)
            return np.allclose(a,b,rtol=5e-9,atol=5e-16)
        assert equal(data,read(destination))
        print("1440 three inverse branches and independent eta holdout verified")
    else:
        print(json.dumps(data,indent=2))
