"""Whole five-angle-free bounds on the two-plane commutator norm (extra input: theta23, theta13)."""
import argparse
import json
from fractions import Fraction as F
from math import sqrt
from pathlib import Path

import numpy as np


def centered(h):
    return h-np.trace(h)*np.eye(3)/3


def area_invariant(h,k):
    h,k=centered(h),centered(k)
    c=1j*(h@k-k@h)
    gram=np.vdot(h,h).real*np.vdot(k,k).real-np.vdot(h,k).real**2
    return np.vdot(c,c).real/gram,c


def leptons():
    # Exact rational interval on tau overlap x, shape parameters h and r. No phi or theta12.
    h,r=F(4,1000),F(4,100)
    den1=F(4,9)*(1-h+h*h)*(1-r+r*r)
    candidates=[]
    for i in range(250,610):
        lo,hi=F(i,1000),F(i+1,1000)
        xc=max(lo,min(hi,F(1,2)))
        tmax=max(abs(lo-(1+h)*(1+r)/3),abs(hi+r*(1-hi)+h-F(1,3)))
        gram=den1-tmax*tmax
        assert gram>0
        candidates.append(2*(xc*(1-xc)+h*h/4)/gram)
    upper=max(candidates)
    assert upper<F(142,100)
    return float(upper)


def projections():
    # Both widths together; inequalities are monotone and use raw endpoints.
    def normalized(ylo,ymid,yhi):
        return (ymid**2-ylo**2)/(yhi**2-ylo**2)
    ru=normalized(0,.00356+2*.00006,.967-2*.004)
    rd=normalized(0,.000306+2*.000004,.01630-2*.00009)
    re=normalized(0,.000585042+2*.000000075,.0099378-2*.0000014)
    q=ru+rd+ru*rd
    assert q<.000405 and re<.004
    kq=(.998,.9984)
    for sign in [-1,1]:
        w=(1-(.00370+sign*2*.00008)**2)*(1-(.04193+sign*2*.00041)**2)
        assert kq[0]<w<kq[1]
    lepton_boxes={}
    for name,solar,ds,atm,da,mix,ms in [
            ("normal",.0000750,.0000012,.00251,.00010,.514,.056),
            ("inverted",.0000753,.0000018,.00248,.00010,.511,.055)]:
        gap=solar+2*ds
        r=gap/(atm-2*da+gap)
        assert r<.04
        pmin=(1-sqrt(1-(.0851-2*.0024)))/2
        pmax=(1-sqrt(1-(.0851+2*.0024)))/2
        lo=(1-mix-2*ms)*(1-pmax)
        hi=(1-mix+2*ms)*(1-pmin)
        assert .25<lo<=hi<.61
        lepton_boxes[name]=[lo,hi,r]
    # Relative deviations of commutator and wedge, both positive at x=0.9984.
    a,b=kq
    cc=sqrt(2*b*(1-b))
    ww=sqrt((1-b)*(b+1/3))
    q=.000405
    lower=2*a/(a+1/3)*((1-q/sqrt(2)/cc)/(1+2*q/(3*ww)))**2
    assert lower>1.466
    return {"quark_B_lower":lower,"lepton_B_upper":leptons(),
        "two_width_lepton_bounds_x_and_rmax":lepton_boxes,"two_width_quark_h":ru,"two_width_quark_r":rd}


def competition():
    a=np.diag([0.,1.,3.])
    b=np.array([[0,1,-1j],[1,0,1],[1j,1,0]],dtype=complex)
    v,c=area_invariant(a,b)
    w,d=area_invariant(a,b+np.diag([0.,5.,0.]))
    assert np.allclose(c,d) and abs(np.linalg.det(c))>1 and abs(v-w)>.3
    changed,_=area_invariant(4*a+b+2*np.eye(3),2*a-3*b+7*np.eye(3))
    assert abs(changed-v)<1e-14
    return {"commutator_spectrum":np.linalg.eigvalsh(c).tolist(),
            "nonzero_determinant":np.linalg.det(c).real,"area_invariants":[v,w]}


def standard(s12,s13,s23,phase):
    c12,c13,c23=np.sqrt(1-np.array([s12,s13,s23])**2)
    e=np.exp(1j*phase)
    return np.array([[c12*c13,s12*c13,s13/e],[-s12*c23-c12*s23*s13*e,c12*c23-s12*s23*s13*e,s23*c13],
        [s12*s23-c12*c23*s13*e,-c12*s23-s12*c23*s13*e,c23*c13]])


def examples():
    v=standard(.2251,.0037,.04193,1.139)
    h=np.diag(np.array([.00000704,.00356,.967])**2)
    k=v@np.diag(np.array([.0000154,.000306,.01630])**2)@v.conj().T
    bq,_=area_invariant(h,k)
    assert bq>1.466
    out=[]
    for phase in [0.,np.pi/2,np.pi,3*np.pi/2]:
        p13=(1-sqrt(1-.0851))/2
        u=standard(sqrt(.3092),sqrt(p13),sqrt(.514),phase)
        charged=np.diag(np.array([.00000277713,.000585042,.0099378])**2)
        nu=u@np.diag([0,.000075,.00251+.000075])@u.conj().T
        bl,_=area_invariant(charged,nu)
        assert bl<leptons()
        out.append(bl)
    return {"central_quark_B":bq,"central_NO_B_phases_0_pi2_pi_3pi2":out}


def rounded(x):
    if isinstance(x,dict):
        return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,list):
        return [rounded(v) for v in x]
    return float(format(x,".10g")) if isinstance(x,(float,np.floating)) else x


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--write",action="store_true")
    ap.add_argument("--check",action="store_true")
    out=Path(__file__).with_name("independent_results.json")
    result=rounded({"certified_domain":projections(),"counterexample":competition(),"examples":examples()})
    if ap.parse_args().write:
        with out.open("x",encoding="utf-8") as fp:
            json.dump(result,fp,indent=2,ensure_ascii=False)
            fp.write("\n")
    else:
        assert result==json.loads(out.read_text(encoding="utf-8"))
    print("1386 independent area obstruction, full CP range and explicit competition verified.")
