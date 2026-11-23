"""Exact extremes for a common two-family Yukawa block with one Clebsch."""
import argparse
import json
from itertools import product
from math import acosh, cosh, exp, hypot, log, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def maximum_double_ratio(ratio, clebsch):
    gap = clebsch * (1 - ratio)
    heavy = (gap + hypot(gap, 2 * sqrt(ratio))) / 2
    return heavy**2


def attainable_trace(r, gain_squared, w):
    """Maximize P_e at fixed P_d=1/r+r over the imbalance of the two down edges."""
    pd = 1/r + r
    q_bound = acosh(pd/2)
    left = gain_squared - w*w
    right = 1/gain_squared - w*w
    candidates = [0, q_bound]
    if left < 0 and right < 0:
        q_star = .5*log(right/left)
        if 0 < q_star < q_bound:
            candidates.append(q_star)
    return max(2*cosh(q+log(gain_squared)) + w*w*(pd-2*cosh(q)) for q in candidates)


def minimum_diagonal_weight(r, charged_ratio, gain_squared):
    target = charged_ratio + 1/charged_ratio
    lo, hi = 1., 20.
    for _ in range(70):
        mid = (lo+hi)/2
        if attainable_trace(r,gain_squared,mid) < target:
            lo=mid
        else:
            hi=mid
    return hi


def run():
    cfg = json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))
    old = json.loads((HERE / cfg["prior"]).read_text(encoding="utf-8"))
    rows = old["scales"][cfg["scale_key"]]
    d, s, e, mu = [rows[i] for i in [0, 1, 3, 4]]
    clebsch = abs(cfg["lepton_clebsch"])
    def estimates(values):
        down, strange, electron, muon = values
        r = down / strange
        observed = r * muon / electron
        cap = maximum_double_ratio(r, clebsch)
        c_min = (sqrt(observed) - r / sqrt(observed)) / (1 - r)
        return {"down_ratio": r, "double_ratio": observed, "upper_at_three": cap,
                "electron_prediction_floor": muon * r / cap, "electron_holdout": electron,
                "observed_over_upper": observed / cap, "necessary_clebsch": c_min,
                "diagonal_singlet_fraction_floor": (c_min-clebsch)/(1+clebsch)}
    inputs = [d, s, e, mu]
    central = estimates([pair[0] for pair in inputs])
    width = cfg["marginal_widths"]
    corners = [estimates(z) for z in product(*[(v-width*w, v+width*w) for v,w in inputs])]
    bounds = {k:[min(v[k] for v in corners),max(v[k] for v in corners)] for k in central}
    # Construct all the way from the balanced to an off-diagonal transpose block.
    r = central["down_ratio"]
    count = 0
    for off_sq in np.geomspace(r*r, r, 23):
        a, b = sqrt(off_sq), r/sqrt(off_sq)
        c = sqrt(max(0, 1+r*r-a*a-b*b))
        for phase in [0.0, 0.71, 2.0]:
            dmat = np.array([[0,a*np.exp(1j*phase)],[b,c*1j]])
            emat = dmat.copy()
            emat[1,1] *= -clebsch
            sd = np.linalg.svd(dmat,compute_uv=False)
            se = np.linalg.svd(emat,compute_uv=False)
            assert np.allclose(sd,[1,r]) and np.isclose(np.prod(se),r)
            assert 1-1e-11 <= r*se[0]/se[1] <= maximum_double_ratio(r,clebsch)+1e-10
            count += 1
    assert np.isclose(r*se[0]/se[1],maximum_double_ratio(r,clebsch))
    assert bounds["electron_prediction_floor"][0] > bounds["electron_holdout"][1]
    assert bounds["necessary_clebsch"][0] > 3
    tradeoff = {}
    for k in [1.,2.,4.,8.,9.,10.]:
        required=[minimum_diagonal_weight(row[0]/row[1],row[3]/row[2],k)
                  for row in product(*[(v-width*w, v+width*w) for v,w in inputs])]
        wc=minimum_diagonal_weight(r,mu[0]/e[0],k)
        assert abs(attainable_trace(r,k,wc)-(mu[0]/e[0]+e[0]/mu[0])) < 2e-10
        tradeoff[str(k)]={"center":wc,"easiest_two_width_corner":min(required)}
    return {"round":1429, "central":central,"two_width_rectangle":bounds,"svd_examples":count,
            "joint_edge_gain_squared_vs_diagonal_weight":tradeoff}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args = parser.parse_args()
    output = run()
    path = HERE / "results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(output,stream,indent=2)
            stream.write("\n")
    else:
        assert output == json.loads(path.read_text(encoding="utf-8"))
        print("1429 joint gain and diagonal-source tradeoff verified")
