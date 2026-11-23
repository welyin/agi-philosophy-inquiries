"""Two charge levels: exact branch dimensions and a three-flavour Gram enclosure."""
import argparse
from fractions import Fraction
import json
from math import ceil, floor, pi, sin, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def enclosure(value, direction, digits=7):
    scale = 10**digits
    return Fraction((floor if direction < 0 else ceil)(value*scale), scale)


def positive_shift(diagonal, edges, shift):
    # Principal minors of diag(diagonal) - edges - shift I (real symmetric).
    x, y, z = [d-shift for d in diagonal]
    a, b, c = edges
    return x > 0 and x*y-a*a > 0 and x*y*z-x*c*c-y*b*b-z*a*a-2*a*b*c > 0


def lower_eigenvalue(diagonal, edges):
    left, right = Fraction(-3), Fraction(3)
    for _ in range(52):
        middle = (left+right)/2
        if positive_shift(diagonal, edges, middle):
            left = middle
        else:
            right = middle
    return float(left)


def run():
    cfg = read("inputs.json")
    w, z = cfg["W_ratio"], cfg["Z"]
    wm, wt = w["muon_over_e"], w["tau_over_e"]
    n = w["width_multiplicity"]
    lows, highs = [Fraction(1)], [Fraction(1)]
    transfer_low, transfer_high = w["relative_rate_to_gram"]
    for center, error in [(wm["value"],wm["work_halfwidth_per_unit"]),
                          (wt["value"],wt["stat"]+wt["syst"])]:
        lows.append(enclosure((center-n*error)/transfer_high,-1))
        highs.append(enclosure((center+n*error)/transfer_low,1))

    coefficient = z["min_rate_response"]*z["Gmu"]*z["min_mZ_GeV"]**3/(6*sqrt(2)*pi*z["max_width_GeV"])
    coefficient = float(enclosure(coefficient,-1))
    limits = [z[key] for key in ("branch_limit_em","branch_limit_et","branch_limit_mt")]
    entries = [enclosure(sqrt(limit/coefficient),1) for limit in limits]
    lower = lower_eigenvalue(lows,entries)
    upper = -lower_eigenvalue([-x for x in highs],entries)
    distance = max(1-sqrt(lower),sqrt(upper)-1)
    defect = max(sqrt(lower)*(1-lower),sqrt(upper)*(upper-1))
    assert lower > 1/3 and upper >= 1
    result = {
        "diagonal_intervals":list(zip(map(float,lows),map(float,highs))),
        "branch_coefficient_lower":coefficient, "offdiagonal_upper":list(map(float,entries)),
        "squared_singular_value_enclosure":[lower,upper],
        "all_three_nonzero_singular_values":[sqrt(lower),sqrt(upper)],
        "distance_upper_to_one_isometry":distance,"cubic_residual_upper":defect,
    }
    delta = cfg["test"]["delta"]
    strengths = np.array([1,1+delta,1-delta])
    A = np.diag(strengths)
    T = np.zeros((6,6))
    T[:3,3:] = A
    H = T@T.T-T.T@T
    residual = (H@T-T@H)/2-T
    assert np.allclose(residual[:3,3:],A@A.T@A-A)
    assert np.all((strengths**2 >= np.array(lows,dtype=float)) &
                  (strengths**2 <= np.array(highs,dtype=float)))
    N = cfg["test"]["pulses_N"]
    coupled, uncoupled = sin((1+delta)*N*pi)**2, sin(N*pi)**2
    result["distinct_blocks_witness"] = {
        "algebra_dimension":3*len(set(strengths)), "distance":delta,
        "operator_residual":float(np.linalg.norm(residual,2)),
        "common_pulse_area":N*pi,
        "excitation_common_vs_split":[uncoupled,coupled],
    }
    h = cfg["test"]["hidden_loss"]
    heavy_row = np.array([sqrt(h),-sqrt(h)])
    G = np.eye(2)-np.outer(heavy_row,heavy_row)
    eig, vectors = np.linalg.eigh(G)
    light = (vectors*np.sqrt(eig))@vectors.T
    full = np.vstack([light,heavy_row])
    assert np.allclose(full.T@full,np.eye(2))
    result["missing_neutral_partner"] = {
        "inclusive_light_mu_over_e":G[1,1]/G[0,0],
        "normalized_projected_overlap":G[0,1]/G[0,0],
        "full_neutral_LFV_charge":0,
    }
    return json.loads(json.dumps(result),parse_float=lambda x:float(f"{float(x):.11g}"))


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    cli.add_argument("--write-new",action="store_true")
    args = cli.parse_args()
    results = run()
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as out:
            json.dump(results,out,indent=2)
            out.write("\n")
    else:
        assert results == read("results.json")
        print("1473: three-flavour spectral bounds, algebra defects, and missing partner checked")
