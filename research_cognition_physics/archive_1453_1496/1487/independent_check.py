"""Independently convert masses and profile physical two-body Higgs widths."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

BASE=Path(__file__).resolve().parent


def read(name):
    return json.loads((BASE/name).read_text(encoding="utf-8"))


def rank_ratio(q,s):
    b=math.hypot(q*math.sqrt(1-s*s),s)
    other=(1-q*q)*s*math.sqrt(1-s*s)/b
    charged=np.array([[0., q/(2*b)],[6*b,6*other]])
    values=np.linalg.svd(charged,compute_uv=False)
    return float(values[0]/values[1])


def rate_test(t):
    # Distances |3-k| for real k with measured |k²-1|<=x, all four signs.
    def possible(err):
        l,h=math.sqrt(1-err),math.sqrt(1+err)
        intervals=[(3-h,3-l),(3+l,3+h)]
        return any(a/d<=t<=b/c for (a,b),(c,d) in itertools.product(intervals,repeat=2))
    lo,hi=0.,.7
    assert possible(hi)
    for _ in range(75):
        mid=(lo+hi)/2
        if possible(mid): hi=mid
        else: lo=mid
    # Two exactly SM partial widths with complex k, electron pure scalar.
    # k_S=cos(phi), k_P=beta sin(phi) retains the exact physical two-body rate.
    beta=math.sqrt(1-(2*.1056583755/125.0)**2)
    low,high=0.,math.pi
    for _ in range(75):
        phase=(low+high)/2
        k=complex(math.cos(phase),beta*math.sin(phase))
        if abs(3-k)<2/t: low=phase
        else: high=phase
    phase=(low+high)/2
    k=complex(math.cos(phase),beta*math.sin(phase))
    assert np.allclose([t*abs(3-k),k.real**2+(k.imag/beta)**2],[2.,1.])
    def cp_intervals(cap, vitesse):
        edge=math.sqrt(max(0.,1-(cap/vitesse)**2))
        def distance(c):
            return math.hypot(3-c,vitesse*math.sqrt(max(0.,1-c*c)))
        return [(2.,distance(edge)),(distance(-edge),4.)]
    beta_e=math.sqrt(1-(.001022/125.)**2)
    cp_lo,cp_hi=0.,beta
    for _ in range(75):
        mid=(cp_lo+cp_hi)/2
        feasible=any(a/dd<=t<=bb/cc for a,bb in cp_intervals(mid,beta_e)
                     for cc,dd in cp_intervals(mid,beta))
        if feasible: cp_hi=mid
        else: cp_lo=mid
    assert np.isclose(cp_hi,abs(k.imag),atol=1e-12)
    return {"two_real_absolute_widths_minimax":hi,
            "all_scalar_signs_minimax_absolute_pseudoscalar":cp_hi,
            "exact_SM_width_CP_phase_parameter_deg":math.degrees(phase),
            "actual_mu_chiral_phase_deg":math.degrees(math.atan2(k.imag,k.real)),
            "imaginary_kappa_in_saturating_example":k.imag}


def audit():
    cfg=read("inputs.json")
    m=read(cfg["quark_masses_and_muonium"])
    d,ed=m["lattice"]["md"]
    strange,es=m["lattice"]["ms"]
    pole,ep=m["muonium"]["pole_ratio"]
    _,s,errs=next(row for row in read(cfg["charm"])["edges"] if row[0]=="cd")
    a=3/(2*math.pi*cfg["inverse_alpha_on_shell"])
    common=lambda ratio: ratio*(1+a*math.log(ratio))
    normalization=rank_ratio(d/strange,s)
    # Near mass MZ numerical expansion from each individual species.
    mu=.1056583755
    electron=mu/pole
    corr=lambda m: -(1+.75*math.log((91.1876/m)**2))/(math.pi*cfg["inverse_alpha_on_shell"])
    assert np.isclose(pole*(1+corr(mu)-corr(electron)),common(pole))
    qhi=(d+2*ed)/(strange-2*es)
    clo=(s-2*sum(errs))/1.02
    t=common(pole)/normalization
    hi=common(pole+2*ep)*cfg["differential_transport_remainder"][1]/rank_ratio(qhi,clo)
    return {"same_scale_central_t":t,"widest_t_upper":hi,
            "central":rate_test(t),"weakest_requirement_over_box":rate_test(hi)}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    m = p.add_mutually_exclusive_group()
    m.add_argument("--write-new", action="store_true")
    m.add_argument("--check", action="store_true")
    args = p.parse_args()
    path = BASE/"independent_results.json"
    out = json.loads(json.dumps(audit(), default=float))
    if args.write_new:
        with path.open("x", encoding="utf-8") as f:
            f.write(json.dumps(out, indent=2)+"\n")
    elif args.check:
        assert out == json.loads(path.read_text(encoding="utf-8"))
        print("Independent tests passed.")
    else:
        print(json.dumps(out, indent=2))
