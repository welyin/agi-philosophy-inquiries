"""Affine inverse top Yukawa evolution with an analytic gauge kernel."""
import argparse
import json
from math import log, pi, sqrt
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
INPUTS = json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))
MATCH = json.loads((HERE / INPUTS["matching"]).read_text(encoding="utf-8"))
BETA = np.array(INPUTS["beta_g"])
COEFF = np.array(INPUTS["C_gut_norm"])
G0 = np.array([MATCH[key] for key in ("g1_GUT", "g2", "g3")])
K = 8 * pi * pi


def affine(scale, gauge, order=48):
    upper = log(scale / MATCH["mu_ref_GeV"])
    def kernel(times):
        times = np.atleast_1d(times)
        args = -np.outer(times, BETA * gauge**2 / K)
        assert np.min(args) > -1
        return np.exp(np.log1p(args) @ (COEFF / BETA))
    nodes, weights = np.polynomial.legendre.leggauss(order)
    integral = upper / 2 * np.dot(weights, kernel((nodes + 1) * upper / 2))
    return float(kernel(upper)[0]), float(integral * INPUTS["a"] / K)


def infer(y, slope, intercept):
    remainder = 1/y**2 - intercept
    assert remainder > 0
    return sqrt(slope/remainder)


def calculate():
    y, dy = MATCH["yt"]
    report = []
    for scale in INPUTS["scale_GeV"]:
        slope, intercept = affine(scale, G0[:, 0])
        check_slope, check_intercept = affine(scale, G0[:, 0], 96)
        assert abs(intercept-check_intercept) < 1e-12 and slope == check_slope
        _, high_intercept = affine(scale, G0[:, 0] - G0[:, 1])
        _, low_intercept = affine(scale, G0[:, 0] + G0[:, 1])
        eta_min = 1-(y+dy)**2 * high_intercept
        eta_max = 1-(y-dy)**2 * low_intercept
        threshold_min = sqrt((1-INPUTS["sensitivity"])/high_intercept)
        outputs = [1/sqrt(intercept + slope/x**2) for x in INPUTS["UV_tests"]]
        assert eta_min > INPUTS["sensitivity"] and threshold_min > y+dy
        report.append({"scale_GeV": scale, "slope_in_inverse_square": slope,
                       "intercept_in_inverse_square": intercept,
                       "mathematical_upper_IR": 1/sqrt(intercept),
                       "inferred_UV_central": infer(y, slope, intercept),
                       "observed_sensitivity_central": 1-y*y*intercept,
                       "observed_sensitivity_box": [eta_min, eta_max],
                       "min_IR_for_one_percent_focus": threshold_min,
                       "min_additive_IR_matching_change": threshold_min-(y+dy),
                       "finite_boundary_IR_outputs": outputs})
    return {"scheme": "GUT-normalized g1, top-only one loop, dimensionless Yukawa",
            "reference_GeV": MATCH["mu_ref_GeV"], "working_box_yt": [y-dy, y+dy],
            "rounds": report}


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            compare(a[key], b[key])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            compare(x, y)
    elif isinstance(a, str):
        assert a == b
    else:
        assert np.isclose(a, b, rtol=5e-11, atol=5e-13), (a, b)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    answer = calculate()
    if args.check:
        compare(answer, json.loads((HERE / "results.json").read_text(encoding="utf-8")))
        print("1395 analytic affine-map reconstruction verified")
    else:
        print(json.dumps(answer, indent=2))
