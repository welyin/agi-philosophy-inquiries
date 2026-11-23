"""Certified triangle and interval bounds on a triplet partial-width combination."""
import argparse
import json
import math
from pathlib import Path

BASE = Path(__file__).resolve().parent


def read(rel):
    return json.loads((BASE / rel).read_text(encoding="utf-8"))


def interval(row, widths):
    return row[0] - widths * row[1], row[0] + widths * row[-1]


def flavor_box(order, widths, cfg):
    prior = read("../1386/inputs.json")
    juno = read("../1446/inputs.json")["juno_207_days"]
    mixing = interval(prior["DayaBay_amplitude"], widths)
    reactor = tuple((1 - math.sqrt(1 - a)) / 2 for a in mixing)
    solar = interval(cfg["solar_angle"], widths + cfg["solar_center_padding"])
    small = interval(juno["solar_21"], widths + cfg["solar_center_padding"])
    key = "normal_atmospheric_31" if order == "NO" else "inverted_abs_31"
    atmospheric = interval(juno[key], widths)
    # For IO the primary result is for *signed* Δ31. v1 quotes -2.482^(+.025)_(-.026).
    # The archived absolute-value box happens to have the right lower and upper sides.
    return solar, reactor, interval(prior[order]["atmospheric_sine_squared"], widths), small, atmospheric


def bounds(order, box):
    x, y, z, small, atm = box
    # x(1-x) is increasing throughout the adopted solar rectangle.
    background_max = small[1] * math.sqrt(x[1]*(1-x[1])*(1-z[0]))
    if order == "NO":
        scale_min = atm[0] - small[1]*x[1]
    else:
        scale_min = atm[0] + small[0]*x[0]
    p_lower = math.sqrt(1-y[1])*(math.sqrt(y[0]*z[0])*scale_min - background_max)
    solar_g = (small[0]*(1-3*x[1]*(1-y[0])), small[1]*(1-3*x[0]*(1-y[1])))
    atm_g = (atm[0]*(1-3*y[1]), atm[1]*(1-3*y[0]))
    assert 0 < x[0] <= x[1] < .5 and solar_g[0] > 0 and atm_g[0] > solar_g[1]
    g_range = ((solar_g[0]+atm_g[0], solar_g[1]+atm_g[1]) if order == "NO"
               else (solar_g[0]-atm_g[1], solar_g[1]-atm_g[0]))
    assert p_lower > 0
    return {"P_lower_eV2": p_lower, "G_interval_eV2": list(g_range),
            "R_lower": p_lower / max(map(abs, g_range))}


def calculate():
    cfg = read("inputs.json")
    gf = read("../1449/inputs.json")["charged_current_check"]["GF_GeVminus2"]
    inverse_alpha = read("../1392/inputs.json")["qcd"]["inverse_alpha_zero_rounded"]
    m = cfg["mass_GeV"]
    norm = math.sqrt(cfg["meg_limit"]*192*math.pi*gf**2*inverse_alpha)
    orders = {}
    for order in ("NO","IO"):
        result = bounds(order, flavor_box(order, cfg["widths"], cfg))
        limits = {}
        for name, weight in (("equal_charged_masses",9),("any_single_charged_mass",8)):
            w = norm*m**2/(weight*8*math.pi*result["R_lower"])
            limits[name] = {"max_abs_W": w, "max_abs_signed_width_GeV": m*w,
                            "at_90pct_amplitude_retention_GeV": m*w/.9}
        result["benchmarks"] = limits
        orders[order] = result
    return {"orders":orders, "benchmark_Mpp_GeV":m,
            "diagonal_competitor_W":.3**2/(8*math.pi),"diagonal_competitor_dipole":0}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args = ap.parse_args()
    result = calculate()
    if args.write_new:
        with (BASE/"results.json").open("x",encoding="utf-8") as out:
            json.dump(result, out, indent=2)
            out.write("\n")
    else:
        assert result == read("results.json")
    print(json.dumps(result["orders"], ensure_ascii=True))
