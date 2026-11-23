"""Verify the two-family shared Yukawa texture with NumPy and published input rows."""
import argparse
import json
from math import sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def trace_ratio(ratio):
    return ratio + 1 / ratio


def ratio_from_trace(trace):
    return 0.5 * (trace + sqrt((trace - 2) * (trace + 2)))


def minimum_weight(down, lepton):
    return (sqrt(lepton) - 1 / sqrt(lepton)) / (sqrt(down) - 1 / sqrt(down))


def ceiling(down, weight):
    return ratio_from_trace(2 + weight**2 * (trace_ratio(down) - 2))


def min_common_distortion(down, lepton, weight):
    return max(1.0, (sqrt(lepton) + weight / sqrt(down)) /
               (weight * sqrt(down) + 1 / sqrt(lepton)))


def matrix_ratio(matrix):
    masses = np.linalg.svd(matrix, compute_uv=False)
    return masses[0] / masses[1]


def witness(down, lepton, weight):
    # Normalize |a*b| to one. Solve separately for off-diagonal imbalance and |c|.
    strength = (trace_ratio(lepton) - trace_ratio(down)) / (weight**2 - 1)
    imbalance = trace_ratio(down) - strength
    assert strength >= 0 and imbalance >= 2
    a = sqrt(ratio_from_trace(imbalance))
    b = 1 / a
    c = sqrt(strength)
    first = np.array([[0, a*1j], [b*np.exp(0.31j), c]], dtype=complex)
    second = first.copy()
    second[1, 1] *= -weight
    np.testing.assert_allclose([matrix_ratio(first), matrix_ratio(2.3*second.T)],
                               [down, lepton], rtol=1e-12)
    return [a, b, c]


def evaluate():
    data = json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))
    units = np.array(data["units"])
    results = {}
    for label, point in data["points"].items():
        middle = np.array(point["central"]) * units
        delta = np.array(point["half_width"]) * units
        low, high = middle - delta, middle + delta
        down, lepton = middle[1]/middle[0], middle[3]/middle[2]
        down_lo, down_hi = low[1]/high[0], high[1]/low[0]
        lepton_lo, lepton_hi = low[3]/high[2], high[3]/low[2]
        k = minimum_weight(down, lepton)
        k_bounds = [minimum_weight(down_hi, lepton_lo), minimum_weight(down_lo, lepton_hi)]
        ceiling_3 = ceiling(down, 3)
        assert ceiling(down_hi, 3) < lepton_lo
        distortion = min_common_distortion(down, lepton, 3)
        np.testing.assert_allclose(ceiling(distortion * down, 3), lepton/distortion)
        assert 3 < k_bounds[0] <= k <= k_bounds[1] < 4
        results[label] = {
            "down_ratio": down, "lepton_ratio": lepton,
            "down_ratio_box": [down_lo, down_hi], "lepton_ratio_box": [lepton_lo, lepton_hi],
            "allowed_light_lepton_ratio_at_central_down": {
                str(w): [down, ceiling(down, w)] for w in [1, 3, 4]},
            "min_abs_weight": k, "min_abs_weight_box": k_bounds,
            "max_lepton_ratio_over_input_box_C3": ceiling(down_hi, 3),
            "relative_to_C3_ceiling": lepton/ceiling_3,
            "min_common_condition_distortion_C3": distortion,
            "distortion_lower_on_box_C3": min_common_distortion(down_hi, lepton_lo, 3),
            "min_abs_identity_admixture": (k-3)/(k+1),
            "C4_witness_abs_a_b_c": witness(down, lepton, 4)}
    return results


def compare(actual, stored):
    if isinstance(actual, dict):
        assert actual.keys() == stored.keys()
        for key, value in actual.items():
            compare(value, stored[key])
    else:
        np.testing.assert_allclose(actual, stored, atol=1e-12, rtol=1e-12)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    results = evaluate()
    target = HERE / "results.json"
    if args.write_new:
        with target.open("x", encoding="utf-8") as stream:
            json.dump(results, stream, indent=2, allow_nan=False)
            stream.write("\n")
    else:
        compare(results, json.loads(target.read_text(encoding="utf-8")))
    for scale, values in results.items():
        print(scale, values["min_abs_weight"], values["min_common_condition_distortion_C3"])
