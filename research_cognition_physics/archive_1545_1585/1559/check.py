"""Optimise the common positive mass, inverse-mass and logarithmic moments."""
import argparse
import json
from math import cosh, exp, log, pi, sinh, sqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sharp(budget, eta, cfg):
    """budget is the sum of positive quadratic threshold terms, in GeV^2."""
    b = 8 * pi**2 * budget
    c = 2 * eta / cfg["v_over_sqrt2_GeV"]**2
    mass = (b / c)**0.25
    scale = cfg["cutoff_GeV"]
    ell = log(scale / mass)
    assert cfg["weak_scale_GeV"] < mass < scale and ell > 0.5
    weight = sqrt(b * c)
    # Two degenerate Majoranas with a relative i cancel each other's dim-5 term.
    assert abs(weight * mass**2 / b - 1) < 1e-12
    assert abs(weight / mass**2 / c - 1) < 1e-12
    for shift in [-5, -1, -0.01, 0, 0.01, 1, 5]:
        majorant = ell * cosh(2 * shift) - 0.5 * sinh(2 * shift)
        assert majorant >= ell - shift - 1e-12
    return {"mass_at_equality_GeV": mass, "total_weight_at_equality": weight,
            "log_moment_upper": weight * ell,
            "ratio_gain_from_singlets_upper": exp(3 * weight * ell / (32 * pi**2))}


def minimum_budget(required, eta, cfg, baseline):
    low = high = baseline
    if sharp(baseline, eta, cfg)["log_moment_upper"] >= required:
        raise ValueError("Inverse calculation needs a deficit at the nominal budget")
    while sharp(high, eta, cfg)["log_moment_upper"] < required:
        high *= 2
    for _ in range(90):
        mid = (low + high) / 2
        if sharp(mid, eta, cfg)["log_moment_upper"] < required:
            low = mid
        else:
            high = mid
    return high


def evaluate():
    cfg = read(HERE / "inputs.json")
    low = read(HERE / cfg["weak_scale_charged"])["yukawa"]
    points = read(HERE / cfg["quark_lepton_points"])
    labels = points["labels"]
    d, s = low["down"][:2]
    e, mu = low["charged_leptons"][:2]
    low_ratio = (d[0]-2*d[1])*(mu[0]-2*mu[1])/((s[0]+2*s[1])*(e[0]+2*e[1]))
    cap = cfg["quark_trajectory_caps"]
    quarks = log(cfg["cutoff_GeV"]/cfg["weak_scale_GeV"]) * (
        cap["ys"]**2 + cap["yc"]**2 + (cap["yt"]*cap["Vtd"])**2)
    mh, width = cfg["higgs_mass_GeV"]
    budget = (mh+2*width)**2 / cfg["higgs_transfer_lower"]
    bounds = {}
    for channel in ["mu_mu", "trace"]:
        eta = cfg["eta_"+channel+"_upper"]
        bound = sharp(budget, eta, cfg)
        holdouts = {}
        for key in ["release_down_A", "release_down_B"]:
            v = dict(zip(labels, points["candidates"][key]["values"]))
            ratio = (v["yd"] / v["ys"]) * (v["ymu"] / v["ye"])
            cost = (32*pi**2/3) * log(low_ratio/ratio) - quarks
            threshold_needed = minimum_budget(cost, eta, cfg, budget)
            holdouts[key] = {"high_scale_double_ratio":ratio,"required_log_moment":cost,
                             "additional_integrated_flavour_term_needed":cost-bound["log_moment_upper"],
                             "min_positive_threshold_GeV2_with_same_weak_limit":threshold_needed,
                             "threshold_to_admitted_budget":threshold_needed/budget}
        low_high_ratio = exp(3*(bound["log_moment_upper"]+quarks)/(32*pi**2))
        bounds[channel] = {"upper":bound,"points":holdouts,
                           "necessary_high_double_ratio_lower":low_ratio/low_high_ratio}
    return {"admitted_positive_threshold_GeV2":budget,"weak_double_ratio_lower":low_ratio,
            "quark_integral_allowance":quarks,"bounds":bounds}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = evaluate()
    path = HERE / "results.json"
    if args.check:
        data = read(path)
        def check_values(actual, saved):
            if isinstance(actual, dict):
                assert actual.keys() == saved.keys()
                for key in actual:
                    check_values(actual[key], saved[key])
            else:
                assert abs(actual-saved) <= 2e-12*max(1., abs(actual))
        check_values(result, data)
        print("1559 shared threshold, current and logarithmic budgets verified")
    else:
        path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
