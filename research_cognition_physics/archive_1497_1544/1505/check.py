"""Cross-channel oriented zeros and the scope of a folded two-bin dip."""

import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def kinematics(energy, mw, mz):
    rw,rz = (mw/energy)**2,(mz/energy)**2
    alpha = 1-rw-rz
    beta = math.sqrt(alpha*alpha-4*rw*rz)
    return alpha,beta,beta/(1+rw-rz),beta/(1-rw+rz)


def edge(bound, vw, vz):
    if math.atanh(vw)+math.atanh(vz) <= bound:
        return 1.0
    lo,hi = 0.0,1.0
    for _ in range(60):
        mid = (lo+hi)/2
        if math.atanh(vw*mid)+math.atanh(vz*mid)<bound:
            lo=mid
        else:
            hi=mid
    return (lo+hi)/2


def ideal_depth(q, mw, mz, r, bins, n):
    a,b,vw,vz = kinematics(q,mw,mz)
    x0 = r*a/b
    edges = [edge(y,vw,vz) for y in bins]
    nodes,weights = np.polynomial.legendre.leggauss(n)
    counts = []
    for lo,hi in zip(edges,edges[1:]):
        x = (hi-lo)/2*nodes+(hi+lo)/2
        # Fold the two incoming directions and sum the opposite-helicity channels.
        shape = (1-x**4)*(x*x+x0*x0)/(a*a-b*b*x*x)**2
        counts.append(float((hi-lo)/2*(weights@shape)))
    return 1-2*counts[0]/counts[1]


def evaluate():
    cfg = read(HERE/"inputs.json")
    cal = read(HERE/cfg["calibration"])
    sin2, ds = cal["SLD_sin2_eff_e"]
    mw, dw = cal["CMS_MW_GeV"]
    mz, dz = cal["LEP_MZ_GeV"]
    half = ds+cfg["additional_quark_weak_angle_transfer"]
    xgamma=1/3
    r=-xgamma*sin2/(1-sin2)
    threshold=math.sqrt(mw*mw+mz*mz+2*mw*mz/math.sqrt(1-r*r))
    rows={}
    for energy in cfg["test_invariant_masses_GeV"]:
        a,b,vw,vz = kinematics(energy,mw,mz)
        x = r*a/b
        g1,g2=.5-2*sin2/3,-.5+sin2/3
        assert abs(g1/(a+b*x)+g2/(a-b*x))<1e-12
        def zero(s2, wmass, zmass):
            alpha,beta,*_ = kinematics(energy,wmass,zmass)
            return -alpha/beta*xgamma*s2/(1-s2)
        limits=[zero(k,wi,zi) for k in [sin2-half,sin2+half]
                for wi in [mw-dw,mw+dw] for zi in [mz-dz,mz+dz]]
        depth=ideal_depth(energy,mw,mz,r,cfg["absolute_delta_y_bins"],32)
        assert math.isclose(depth,ideal_depth(energy,mw,mz,r,cfg["absolute_delta_y_bins"],64),abs_tol=1e-11)
        rapidity=math.atanh(vw*x)+math.atanh(vz*x)
        rows[str(energy)]={"W_angle_cos":x,"W_angle_box":[min(limits),max(limits)],
                           "signed_yW_minus_yZ_at_zero":rapidity,"ideal_opposite_TT_depth":depth}
    error=cfg["additional_quark_weak_angle_transfer"]
    # Same isoscalar vertex error h in both quark Z charges.
    tests={}
    for h in [0,0.02]:
        gu,gd=.5-2*sin2/3+h,-.5+sin2/3+h
        zz=(gu+gd)/(gu-gd)
        recovered=((1-sin2)*zz+sin2*xgamma)/2
        assert math.isclose(h,recovered,abs_tol=1e-15)
        tests[str(h)]=zz
    measured,err=cfg["atlas_TT_depth"]
    d0=cfg["D0"]
    for xg in [-.27,1/3,.59]:
        up=(1+xg)/2
        down=(xg-1)/2
        gu,gd=.5-sin2*up,-.5-sin2*down
        zz=(gu+gd)/(gu-gd)
        assert math.isclose(-zz/(xg-zz),sin2,abs_tol=1e-13)
    return {"weak_angle_central":sin2,"quark_angle_extra_budget":error,
            "gamma_zero_cos_Wplus":xgamma,"z_reduced_zero":r,
            "first_physical_z_zero_GeV":threshold,"rows":rows,
            "reduced_z_if_common_quark_vertex_shift":tests,
            "D0_MUH_lower_tail":d0["MUH_at_or_below"]/d0["pseudo_trials"],
            "ATLAS_no_dip_component_max_sides_one_two_widths":
                [1-measured+k*err for k in [1,2]],
            "at_central_depth_minimum_central_suppression_vs_flat":measured}


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            compare(a[k],b[k])
    else:
        np.testing.assert_allclose(a,b,rtol=2e-9,atol=1e-11)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    mode=p.add_mutually_exclusive_group()
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write",action="store_true")
    args=p.parse_args()
    output=HERE/"results.json"
    data=evaluate()
    if args.write:
        serialized=json.dumps(data,indent=2,allow_nan=False)
        with output.open("x",encoding="utf-8") as stream:
            stream.write(serialized+"\n")
    else:
        compare(data,read(output))
    print("1505 verified",data["z_reduced_zero"])
