"""Audit the one-loop threshold polytope with exact representation weights."""
import argparse
from fractions import Fraction as F
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
SM=(F(41,10),F(-19,6),F(-7))
T=(F(3,5),F(2,3),F(0))
L=(F(1,30),F(1,2),F(1,3))
S=(F(8,15),F(0),F(5,6))
W=(F(0),F(1,3),F(0))
assert [sum(m[i] for m in [T,L,S]) for i in range(3)]==[F(7,6)]*3


def beta(h):
    return tuple(SM[i]+sum(b[i]*s for b,s in zip([T,L,S,W],h)) for i in range(3))


def interp(data):
    lo,hi=[],[]
    for h in itertools.product((0,1),repeat=3):
        v=beta([*h,0]); ratio=(v[1]-v[2])/(v[0]-v[1])
        lo.append((ratio,h)); hi.append((ratio,h))
    q0=min(lo)[0]; q=max(hi)[0]
    assert (q0,q)==(F(5,13),F(70,101))
    x,y=data[0]-data[1],data[1]-data[2]; q=float(q);q0=float(q0)
    r=y-q*x
    budget=25*x-32*y
    return {"alpha3_observed":1/data[2],"alpha3_holdout_low":1/(data[1]-q0*x),"alpha3_holdout_high":1/(data[1]-q*x),
       "color_only_extra_inverse_shift_needed":r,"min_span_of_matching_offsets":r/(1+q),
       "extra_real_triplet_min_active_log_over_2pi":3*r/(1+q),"augmented_weight_budget":budget}


def run():
    inp=json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    orig=json.loads((HERE/inp["sm_inputs"]).read_text(encoding="utf-8"))
    gv=np.array([orig[k] for k in ("g1_GUT","g2","g3")])
    low=orig["mu_ref_GeV"];twopi=2*math.pi
    centers=4*math.pi/gv[:,0]**2
    val=interp(centers)
    box=[interp(4*math.pi/(gv[:,0]+inp["two_width"]*gv[:,1]*np.array(c))**2) for c in itertools.product([-1,1],repeat=3)]
    ranges={k:[min(x[k] for x in box),max(x[k] for x in box)] for k in val}
    # Type-II threshold ceilings with one additional real adjoint triplet, closed at MZ.
    maxbud=ranges["augmented_weight_budget"][1]
    def bounds(limit):
        return {"T_GeV":low*math.exp(twopi*limit/23),"L_GeV":low*math.exp(twopi*limit/17),
              "one_of_T_and_L_GeV":low*math.exp(twopi*limit/40),"min_M6_over_MX":math.exp(-twopi*limit/40)}
    # Construct an exact common threshold witness (charge-triplet and L at a fixed mass).
    mass=inp["example_triplet_leptoquark_mass_GeV"]; l=math.log(mass/low)/twopi
    x=centers[0]-centers[1]; y=centers[1]-centers[2]
    tau=(x+y+.3*l)/11.4
    lW=(25*x-32*y-40*l)/19
    massW=low*math.exp(twopi*lW); mx=low*math.exp(twopi*tau)
    vecs=[np.asarray(v,dtype=float) for v in [SM,T,L,S,W]]
    inv=centers.copy()
    spec=[(mass,vecs[1]+vecs[2]),(massW,vecs[4]),(mx,vecs[3])]
    slope=vecs[0].copy(); e0=low
    for scale,change in sorted(spec,key=lambda p:p[0]):
        inv-=slope*math.log(scale/e0)/twopi
        slope+=change; e0=scale
    assert np.ptp(inv)<2e-12 and min(inv)>0 and low<massW<mx
    # No special choice of splitting can overcome the scalar-15 certificate.
    for h in itertools.product([0,.7,1],repeat=3):
        tt=3.5;tr,lq,sex=h
        b=beta([F(str(i)) for i in (*h,0)])
        aa=[40+float(i)*tt for i in b]; xx=aa[0]-aa[1];yy=aa[1]-aa[2]
        assert abs(140*xx-202*yy-tt*(144*(1-tr)+99*(1-lq)+243*sex))<1e-10
    return {"ratio_rational_ends":["5/13","70/101"],"central":val,"two_width":ranges,
       "adjoint_central_threshold_limits":bounds(val["augmented_weight_budget"]),
       "adjoint_two_width_threshold_limits":bounds(maxbud),
       "one_real_adjoint_explicit_central_solution":{"fixed_M_T_M_L":mass,"M_Wtriplet":massW,"M_X_and_M6":mx,
             "common_inverse_alpha":float(np.mean(inv))}}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    a=p.add_mutually_exclusive_group();a.add_argument("--write",action="store_true");a.add_argument("--check",action="store_true")
    opt=p.parse_args(); out=json.dumps(run(),indent=2,allow_nan=False)+"\n";name=HERE/"results.json"
    if opt.write:
        with name.open("x",encoding="utf-8") as f: f.write(out)
    elif out!=name.read_text(encoding="utf-8"): raise SystemExit("1458 mismatch")
    else: print("1458 OK")
