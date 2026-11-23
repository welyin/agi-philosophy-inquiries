"""Piecewise QCD validity and EW inverse-coupling gap on the same 1754 particles."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent
TAU=2*math.pi
SOFT=30000.
WMAX=1000.
SCMAX=31000.


def read(path):
    return json.loads((HERE/path).read_bytes())


def make(case,new_soft=False,ordinary_soft=True):
    low=read("../1369/inputs.json")
    base=read("../1698/results.json")["shared_witness"]
    mc,mw=read("../1744/inputs.json")["pole_mass_GeV"]
    f=read("../1754/results.json")["family"]
    mZ=low["mu_ref_GeV"]
    ranges=np.array([low[k] for k in ("g1_GUT","g2","g3")])
    couplings=ranges[:,0]+np.array([0,0,-.67])*ranges[:,1]
    inputs=4*math.pi/couplings**2
    sm=np.array([4.1,-19/6,-7.])
    susy=np.array([6.6,1,-3.])
    events=[]
    if ordinary_soft:
        events.append((SOFT,susy-sm,-np.array([0,2,3])/(12*np.pi)))
    for name,inp,charge in (("colored",mc,np.array([.4,0,1.])),("weak",mw,np.array([.6,1.,0]))):
        for mode,frac in (("fermions",2/3),("scalars",1/3)):
            if mode=="scalars" and not ordinary_soft:
                continue
            ext=inp if mode=="fermions" else math.hypot(SOFT,inp)
            events.extend((m,frac*charge,np.zeros(3)) for m in [*base[name+"_"+mode+"_GeV"],ext,ext])
    beta=np.zeros(3)
    if case:
        masses=f["masses_chiral_UDEN"][:2]
        # Each chiral family contributes 4/3 to each gauge beta coefficient.
        beta=np.array([4/3,4/3,4/3])+np.array([2/15,2,4/3])+np.array([16/15,0,2/3])
        if case=="single":
            masses=masses+f["vector_poles_including_inactive_weak_partner"]
        else:
            pair=f["extra_down_singlet_precision_holdout"]
            masses=masses+[f["hard_vector_mass"]*(1+sg*pair[key])
                           for key in ("q_up","q_down") for sg in (-1,1)]
            beta+=np.array([4/15,0,2/3])
        assert abs(beta[2]-(2/3)*len(masses))<1e-12
        assert min(masses)>mZ and max(masses)<WMAX
        events.extend((m,np.array([0.,0.,2/3]),np.zeros(3)) for m in masses)
        events.append((WMAX,np.r_[beta[:2],0.],np.zeros(3)))
        if new_soft:
            for m in masses:
                sc=math.hypot(SOFT,m)
                assert SOFT<=sc<=SCMAX
                events.append((sc,np.array([0.,0.,1/3]),np.zeros(3)))
            events.append((SCMAX,np.r_[beta[:2]/2,0.],np.zeros(3)))
    # EW hard masses are a reference envelope; colour thresholds are the actual masses.
    return mZ,inputs,sm,sorted(events,key=lambda e:e[0]),beta,ranges


def point(lo,a,b,events,scale):
    plus=a.copy()
    for mass,db,match in events:
        if mass>scale: break
        plus=plus+match-db*math.log(scale/mass)/TAU
    return plus-b*math.log(scale/lo)/TAU


def check_case(case,new_soft,ordinary_soft):
    lo,a,b,events,beta,ranges=make(case,new_soft,ordinary_soft)
    U=read("../1744/independent_results.json")["new_spectrum"]["cross_GeV"]
    at_U=point(lo,a,b,events,U)
    high=b+sum((e[1] for e in events),np.zeros(3))
    max_m=max(m for m,_,_ in events)
    alast=point(lo,a,b,events,max_m)
    def first_cross(inverse_limit):
        # Find the first crossing, including jumps and non-monotone low-scale flow.
        scale=lo
        val=a[2]
        rate=b[2]
        for m,db,offset in events+[(np.inf,np.zeros(3),np.zeros(3))]:
            if rate>0:
                ln=TAU*(val-inverse_limit)/rate
                if ln>=0 and ln<math.log(m/scale):
                    return scale*math.exp(ln)
            if not np.isfinite(m): break
            val-=rate*math.log(m/scale)/TAU
            val+=offset[2]
            if val<=inverse_limit:return m
            rate+=db[2]
            scale=m
        return None
    stop=first_cross(1.)
    if stop:
        np.testing.assert_allclose(point(lo,a,b,events,stop)[2],1.,atol=2e-13)
    # Select each matching mass independently in [MZ,1 TeV] and [30,31 TeV].
    # A positive SU(2) shift improves the gap. Displace every U(1) logarithm earliest.
    cost=beta[0]*math.log(WMAX/lo)/TAU
    if new_soft:
        cost+=beta[0]*math.log(SCMAX/SOFT)/(2*TAU)
    end=stop or U
    probes=[lo,end]
    probes.extend(s for m,_,_ in events if m<end for s in [m*(1-1e-12),m])
    amp=min((point(lo,a,b,events,s)[0]-point(lo,a,b,events,s)[1],s) for s in probes)
    # If high b1-b2>0 this formal EW meeting is unique after the last threshold.
    meet=max_m*math.exp(TAU*(alast[0]-alast[1])/(high[0]-high[1]))
    if stop and case: assert stop<meet and amp[0]-cost>0
    return dict(high_b=high,extra_IR_EW_tilt_budget=cost,
                alpha3_inverse_at_10PeV=point(lo,a,b,events,1e7)[2],
                alpha3_inverse_at_old_U_formal=at_U[2],
                alpha3_half_scale=first_cross(2.),alpha3_one_scale=stop,
                one_loop_formal_QCD_pole=first_cross(0.),
                reference_EW_meeting_formal=meet,
                EW_a1_minus_a2_minimum_up_to_stop_or_old_U=amp[0],
                scale_at_minimum=amp[1],uniform_gap_lower_bound=amp[0]-cost)


def native(x):
    if isinstance(x,dict):return {k:native(v) for k,v in x.items()}
    if isinstance(x,(tuple,list,np.ndarray)):return [native(v) for v in x]
    if isinstance(x,(float,np.floating)):return float(f"{x:.10g}")
    return x


def calculate():
    baseline=check_case(None,False,True)
    assert abs(baseline["alpha3_inverse_at_old_U_formal"]-10.4151944)<6e-8
    rows={}
    for case in ("single","double"):
        rows[case]={}
        for ns,os,label in ((False,True,"old_soft_new_fermions"),(True,True,"old_soft_new_completed"),
                            (False,False,"all_fermions_without_soft_sectors")):
            row=check_case(case,ns,os)
            if os:
                assert row["alpha3_one_scale"]<1e13
            else:
                assert row["alpha3_one_scale"] is None and row["high_b"][2]<0
            rows[case][label]=row
    return dict(baseline_old_U_inverse_alpha3=baseline["alpha3_inverse_at_old_U_formal"],
                independently_recovered_old_EW_crossing_GeV=baseline["reference_EW_meeting_formal"],models=rows)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    check=cli.parse_args().check
    result=native(calculate())
    path=HERE/"independent_results.json"
    if check:
        assert json.loads(path.read_bytes())==result
        print("1756 independent: colour stopping scales, weak gaps and asymptotically free controls verified.")
    else:
        with path.open("x",encoding="utf8") as out:
            json.dump(result,out,indent=2,allow_nan=False)
            out.write("\n")

