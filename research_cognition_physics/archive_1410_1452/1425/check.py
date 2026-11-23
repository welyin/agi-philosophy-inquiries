"""Common-pole Higgs flavour rates, explicit monotone endpoints, round 1425."""
import argparse
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def compute():
    cfg = read("inputs.json")
    ref = read(cfg["angles_and_matching"])
    def width_values(data, n):
        mid, down, *up = data
        return [mid] if n == 0 else [mid - n*down, mid + n*(up[0] if up else down)]
    def rates(n):
        coords = [ref["DayaBay_amplitude"], ref["NO"]["atmospheric_sine_squared"],
                  ref["quark_sines_and_phase"][0], *ref["yukawa"]["charged_leptons"],
                  *ref["yukawa"]["down"][:2]]
        axes = [width_values(v,n) for v in coords]
        if min(axes[1]) < .5 < max(axes[1]):
            axes[1].append(.5)  # stationary point of y(1-y) in the exchanged row
        rec = []
        for amp, y, cq, e, mu, tau, down, strange in itertools.product(*axes):
            x = (1-math.sqrt(1-amp))/2
            p = (1-x,x*y,x*(1-y))
            ds_proj = cq*cq*(1-cq*cq)
            ds = 3*(down**2+strange**2)*ds_proj
            em = p[0]*p[1]*(e**2+mu**2)
            et = p[0]*p[2]*(e**2+tau**2)
            mt = p[1]*p[2]*(mu**2+tau**2)
            exchanged_mt = y*(1-y)*(mu**2+tau**2)
            rec.append(dict(mt_et_original=mt/et, em_et_original=em/et,
                            ds_mt_exchanged=ds/exchanged_mt, em_ds_original=em/ds,
                            et_ds_original=et/ds, chiral_em_over_ds=math.sqrt(p[0]*p[1]/ds_proj)))
        return {key:[min(row[key] for row in rec),max(row[key] for row in rec)] for key in rec[0]}
    center = {k:v[0] for k,v in rates(0).items()}
    fit = rates(cfg["width_multipliers"])
    b = cfg["leptonic_relative_rate_factor_half_width"]
    # Independent per-channel transfer in [1-b, 1+b]. Hadronic QCD factor kept symbolic.
    r_max = fit["mt_et_original"][1]*(1+b)/(1-b)
    assert r_max < 1
    atlas, cms = cfg["atlas"], cfg["cms"]
    lower = atlas["symmetry_difference"] - cfg["width_multipliers"]*atlas["symmetry_width"]
    et_lim = atlas["MC_marginal95"]["et"]
    assert lower > 0 and lower < atlas["MC_marginal95"]["mt"]
    orig_em_ceil = et_lim * fit["em_et_original"][1]*(1+b)/(1-b)
    # ds below is the pre-hadronic matching branch B_ds/T_ds.
    ds_threshold_em = cms["em_upper95"]/(fit["em_ds_original"][0]*(1-b))
    ds_threshold_et = et_lim/(fit["et_ds_original"][0]*(1-b))
    return {"round":1425,"NO_central":center,"NO_two_width_tree":fit,
            "with_lepton_factor_budget":{"mt_et_original_max":r_max,
            "electron_relative_efficiency_loss_needed_for_fake_positive_difference":1-r_max,
            "original_em_upper_from_et":orig_em_ceil,
            "improvement_over_direct_em":cms["em_upper95"]/orig_em_ceil,
            "ds_over_hadronic_multiplier_lower_required_to_exclude_original_with_em":ds_threshold_em,
            "ds_over_hadronic_multiplier_lower_required_to_exclude_original_with_et":ds_threshold_et,
            "two_width_signed_difference_lower":lower,
            "two_width_signed_difference_response_error_threshold":lower,
            "exchanged_ds_lower_over_hadronic_multiplier_if_positive_lower":lower*fit["ds_mt_exchanged"][0]/(1+b),
            "three_width_signed_difference_lower":atlas["symmetry_difference"]-3*atlas["symmetry_width"]},
            "working_exchanged_witness":{"et":0,"em":0,"mt":0.001,
              "passes_listed_marginal_bounds_and_two_width_difference":bool(
                  lower<.001<min(atlas["MC_marginal95"]["mt"],cms["tau_upper95"]["mt"]))}}


def same(x,y):
    if isinstance(x,dict):
        return isinstance(y,dict) and x.keys()==y.keys() and all(same(v,y[k]) for k,v in x.items())
    if isinstance(x,list):
        return isinstance(y,list) and len(x)==len(y) and all(same(a,b) for a,b in zip(x,y))
    if isinstance(x,bool):
        return x is y
    return math.isclose(x,y,rel_tol=2e-10,abs_tol=1e-16)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    ans=compute()
    out=HERE/"results.json"
    if args.write:
        with out.open("x",encoding="utf-8") as stream:
            json.dump(ans,stream,indent=2)
            stream.write("\n")
    else:
        assert same(ans,read("results.json"))
    print(json.dumps(ans))
