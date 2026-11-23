"""Independent 1539 audit: inverse triangle, monotone m_bb and finite tau->mu gamma floor."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parent


def load(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def compute():
    cfg=load(ROOT/"inputs.json")
    data=load(ROOT/cfg["oscillation"])
    main=load(ROOT/"results.json")
    n=cfg["widths"]
    pm=lambda pair:(pair[0]-n*pair[1],pair[0]+n*pair[-1])
    smin,smax=pm(data["IO"]["solar_sine_squared"])
    clo,chi=pm(data["IO"]["atmospheric_sine_squared"])
    tmin,tmax=[(1-math.sqrt(1-s))/2 for s in pm(data["DayaBay_amplitude"])]
    ratio=math.sqrt(tmax/(1-tmax))/(1-2*smax)
    qlow=[clo*(1-tmax),(1-chi)*(1-tmax)]
    def component(q):
        return max(0,math.sqrt(q)-ratio*math.sqrt(1-q))**2/(1+ratio**2)
    low_mu,low_tau=map(component,qlow)
    product=min(x*(1-x) for x in (low_mu,1-low_tau))
    assert 0<product<.25
    br_opal_min=.1734-2*math.hypot(.0009,.0006)
    gf=1.1663787e-5
    alpha_lower=1/138
    effective_transfer=.90
    min_br_at_one_TeV=effective_transfer*br_opal_min*alpha_lower*product/(48*math.pi*gf**2)*1e-12
    bound=(4.2e-8/min_br_at_one_TeV)**.25
    # Verify the parent's finite witness by reconstructing the triangle directly.
    z=cfg["m3_witness_eV"]
    a=data["IO"]["atmospheric_gap"][0]
    d=data["IO"]["solar_gap"][0]
    s=data["IO"]["solar_sine_squared"][0]
    t=(1-math.sqrt(1-data["DayaBay_amplitude"][0]))/2
    p=np.array([(1-t)*(1-s),(1-t)*s,t])
    def square_terms(z):
        return np.array([a+z*z,a+d+z*z,z*z])
    def mbb_sq(z):
        terms=p**2/square_terms(z)
        return a*(a+d)*terms[2]-a*d*terms[0]+(a+d)*d*terms[1]
    masses=np.sqrt(square_terms(z))
    legs=p/masses
    phi=math.acos((legs[2]**2-legs[0]**2-legs[1]**2)/(2*legs[0]*legs[1]))
    phases=np.array([1,np.exp(1j*phi),-(legs[0]+np.exp(1j*phi)*legs[1])/legs[2]])
    assert abs(legs@phases)<1e-13
    assert math.isclose(abs((p*masses)@phases),math.sqrt(mbb_sq(z)),rel_tol=1e-12)
    assert math.isclose(math.sqrt(mbb_sq(z)),main["witness"]["mbb_eV"])
    # Dominant negative derivative with y=m3² throughout the parent's certified outer window.
    amin,amax=pm(data["IO"]["atmospheric_gap"])
    dmin,dmax=pm(data["IO"]["solar_gap"])
    ymin,ymax=[v*v for v in main["IO_m3_outer_eV"]]
    neg_min=amin*(amin+dmin)*tmin**2/ymax**2
    pos_max=amax*dmax*((1-tmin)*(1-smin))**2/(amin+ymin)**2
    zmin,zmax=main["IO_m3_central_interval_eV"]
    ends=[math.sqrt(mbb_sq(zmax)),math.sqrt(mbb_sq(zmin))]
    assert neg_min>pos_max and zmin < z < zmax
    return {"max_tan_angle_n_vs_nu3":ratio,"min_x_mu":low_mu,"min_x_tau":low_tau,
            "min_product_x_mu_x_tau":product,"min_norm_Br_tau_mu_nu":br_opal_min,
            "tau_mu_gamma_floor_at_f0_over_M_one_per_TeV":min_br_at_one_TeV,
            "bound_f0_over_M_inverse_TeV":bound, "soft_matching_rate_floor":effective_transfer,
            "sharp_mbb_central_eV":ends,"dominant_slope_min":neg_min,
            "opposing_slope_max":pos_max}


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    opt=cli.add_mutually_exclusive_group(required=True)
    opt.add_argument("--write-new",action="store_true")
    opt.add_argument("--check",action="store_true")
    args=cli.parse_args()
    res=compute()
    path=ROOT/"independent_results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as dst:
            json.dump(res,dst,indent=2,allow_nan=False)
            dst.write("\n")
    else:
        saved=load(path)
        assert res.keys()==saved.keys()
        for key in res:
            np.testing.assert_allclose(res[key],saved[key],rtol=1e-10,atol=1e-16)
        print("1539 inverse triangle, m_bb endpoints, flavor floor and Belle readout verified")
