"""Independent one-loop matching and thresholds for the |beta|=sqrt(3) embedding."""

import argparse
import itertools
import json
import math
from pathlib import Path


BASE = Path(__file__).resolve().parent
LOOP = 8 * math.pi**2


def threshold_at(scale, mz, g1, g2):
    log_span = math.log(scale / mz)
    inv_y = 1 / (0.6 * g1**2) - (41 / 6) * log_span / LOOP
    inv_w = 1 / g2**2 + (19 / 6) * log_span / LOOP
    deficit = -LOOP * (inv_y - 3 * inv_w)
    assert deficit > 0
    # One real X=0 adjoint adds log(m_D/m_T) to the signed budget in loop units.
    min_adjoints = math.floor(deficit / log_span) + 1
    perturbative_deficit = deficit + LOOP / (4 * math.pi)
    min_for_alpha_one = math.ceil(perturbative_deficit / log_span)
    # Scale connection coordinates so that the neutral kinetic metric has unit diagonal.
    min_kinetic_mixing = math.sqrt(1 - inv_y / (3 * inv_w))
    g_l = 1 / math.sqrt(inv_w)
    g_x = g_l / (math.sqrt(3) * min_kinetic_mixing)
    beta = -math.sqrt(3)
    matching = beta**2 * inv_w + 1 / g_x**2 + 2 * beta * min_kinetic_mixing / (g_x * g_l)
    assert math.isclose(matching, inv_y) and min_kinetic_mixing < 1
    return {
        "required_sum_log_split": deficit, "required_split_product": math.exp(deficit),
        "maximum_log_split_per_octet": log_span,
        "minimum_number_for_positive_canonical_matching": min_adjoints,
        "minimum_number_for_alphaX_at_most_one": min_for_alpha_one,
        "minimum_abs_kinetic_mixing": min_kinetic_mixing, "gX_at_equality": g_x,
    }


def certificate():
    inputs = json.loads((BASE.parent / "1369" / "inputs.json").read_text(encoding="utf-8"))
    mz = inputs["mu_ref_GeV"]
    g1, err1 = inputs["g1_GUT"]
    g2, err2 = inputs["g2"]
    pole_scale = mz * math.exp(LOOP * (1 / (0.6 * g1**2) - 3 / g2**2) / (49 / 3))
    scale = 10000
    center = threshold_at(scale, mz, g1, g2)
    corners = [
        threshold_at(scale, mz, g1 + a * err1, g2 + b * err2)
        for a, b in itertools.product([-2, 2], repeat=2)
    ]
    # The threshold exponents and mixing bounds are monotone in each matched coupling.
    rectangles = {
        key: [min(row[key] for row in corners), max(row[key] for row in corners)]
        for key in center
        if key != "gX_at_equality"
    }
    # Direct SU(3)-weight traces including a symmetric six representation and its conjugate.
    fundamental = [(0.5, 1/(2*math.sqrt(3))), (-0.5, 1/(2*math.sqrt(3))),
                   (0.0, -1/math.sqrt(3))]
    six = [(i[0]+j[0], i[1]+j[1])
           for idx, i in enumerate(fundamental) for j in fundamental[idx:]]
    weight_checks = {}
    for label, weights in (("fundamental", fundamental), ("sextet", six),
                           ("conjugate_sextet", [(-t, -h) for t, h in six])):
        for charge in [0.0, 0.5, -1/3]:
            surplus = sum((charge - math.sqrt(3)*h)**2 - 3*t**2 for t, h in weights)
            assert math.isclose(surplus, len(weights) * charge**2, abs_tol=1e-13)
            weight_checks[label + "_" + str(charge)] = round(surplus, 12)
    return {"canonical_SM_face_GeV": pole_scale,
            "adjoint_real_triplet_delta_slope": -1, "adjoint_complex_doublet_delta_slope": 1,
            "central_at_10_TeV": center, "independent_two_width_corners": rectangles,
            "explicit_weight_traces_before_spin_color_factor": weight_checks}


def same(a, b):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(same(value, b[key]) for key, value in a.items())
    if isinstance(a, list):
        return len(a) == len(b) and all(same(i, j) for i, j in zip(a, b))
    return math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-12)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    path = BASE / "independent_results.json"
    result = certificate()
    if args.write:
        with path.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")
    elif not same(result, json.loads(path.read_text(encoding="utf-8"))):
        raise ValueError("Independent result mismatch")
    print(json.dumps(result, indent=2))
