"""Independent 1533 audit. Use charged-Higgs rather than neutral-production templates."""
import argparse
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT/"inputs.json").read_text(encoding="utf-8-sig"))


def widths(mplus, xi_d):
    # Two-body tree widths for t->Wb, t->H+b in common G_F m_t^3 |Vtb|²/(8*pi*sqrt2) units.
    mt, mb, mw = 172.5, 4.6, 80.38
    rb, rw, rh = [(m/mt)**2 for m in (mb, mw, mplus)]
    def ps(r):
        return math.sqrt((1-r-rb)**2-4*r*rb)
    w = ps(rw) * ((1-rb)**2+rw*(1+rb)-2*rw**2)
    h = ps(rh) * ((1+rb-rh)*(1+rb*xi_d**2)+4*rb*xi_d)
    br = h/(h+w)
    # Independent stress: H width relative to W reduced by 15%, other widths <=5% W,
    # only 80% of predicted charged search signal reaches published template.
    br_floor = .85*h/(.85*h+1.05*w)
    gf, mtau = 1.1663787e-5, DATA["tau_mass_GeV"]
    tau_width = gf*mplus*mtau**2*(1-(mtau/mplus)**2)**2/(4*math.pi*math.sqrt(2))
    tau_max_nominal = .0027/br
    tau_max_stress = .0027/(.8*br_floor)
    return {"br_t_Hb":br, "stress_floor_br_t_Hb":br_floor,
            "max_br_H_taunu_using_weakest_limit":tau_max_nominal,
            "stress_max_br_H_taunu":tau_max_stress,
            "tau_partial_GeV":tau_width,
            "min_other_partial_GeV_nominal":tau_width*(1/tau_max_nominal-1)}


def calculate():
    mh, v, mtau = [DATA[k] for k in ("higgs_mass_GeV","vev_GeV","tau_mass_GeV")]
    source = json.loads((ROOT/"results.json").read_text(encoding="utf-8"))
    mwitness = math.hypot(mh, math.sqrt(DATA["witness_gap_squared_GeV2"][0]))
    # Linear reconstruction from the measured spin ratio rather than the author's polar solver.
    angle = math.radians(DATA["cms_alpha_tau_deg"][0])
    rate = DATA["cms_tau_over_z_sm"][0]
    result = {}
    for s, label in ((1,"positive"),(-1,"negative")):
        t = math.tan(angle)
        u = s*math.sqrt(rate/(1+t*t/(1-4*(mtau/mh)**2)))
        w = u*t
        n = math.hypot(1,1-u,w)
        assert math.isclose(source[label]["Stau_a"],u)
        assert math.isclose(source[label]["Ptau_a"],w)
        assert math.isclose(source[label]["scalar_ttH_floor"],(2-u)**2)
        assert math.isclose(source[label]["h"][0],1/n)
        soft = source[label]["soft_Z2_stationary_witness"]
        assert math.isclose(soft["charged_mass_GeV"],mwitness)
        # Raise charged mass to the experimental 150 GeV hypothesis leaving all neutral data;
        # lambda3 -> lambda3+shift, lambda4 -> lambda4-shift, lambda345 unchanged.
        shift = 2*(150**2-mwitness**2)/v**2
        assert min(soft["BFB_margins"])>0 and shift>0
        result[label] = {"a":1/n, "tau_even_over_a":u, "tau_odd_over_a":w,
                         "scalar_top_over_a":2-u,"shift_lambda3":shift}
    for width in (1,2):
        row = source["windows"][str(width)]
        assert row["transferred_mirror_min"] > row["ATLAS_range"][1]
    result["charged"] = {str(m):{label:widths(m,xi) for label,xi in
                                 (("II",1),("X",-1))}
                         for m in (mwitness,150.)}
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser()
    arg=p.add_mutually_exclusive_group(required=True)
    arg.add_argument("--check",action="store_true")
    arg.add_argument("--write-new",action="store_true")
    args=p.parse_args()
    file=ROOT/"independent_results.json"
    out=calculate()
    if args.write_new:
        with file.open("x",encoding="utf-8") as f:
            json.dump(out, f, indent=2, allow_nan=False)
            f.write("\n")
    else:
        assert out==json.loads(file.read_text(encoding="utf-8"))
        print("1533: neutral amplitudes, charged-mass shift and independent ATLAS top-decay bounds pass")
