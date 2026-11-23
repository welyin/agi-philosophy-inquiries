"""Reproduce the pole-residue boxes, chiral moments, and isospin witnesses."""
import argparse
import json
from math import hypot, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE = json.loads((HERE / "inputs.json").read_text(encoding="utf-8-sig"))


def interval(name, widths):
    center, *errors = SOURCE[name]
    return center - widths * sum(errors), center + widths * sum(errors)


def polar_tests(n):
    electron = interval("SLD_Ae_1996_98", n)
    tau = interval("OPAL_Atau", n)
    forward = interval("DELPHI_AFB0tau_9parameter", n)
    transverse = interval("ALEPH_CTT", n)
    normal = interval("ALEPH_CTN", n)
    expected_forward = [0.75 * electron[0] * tau[0], 0.75 * electron[1] * tau[1]]
    near = [sqrt(1 - tau[1]**2), sqrt(1 - tau[0]**2)]
    aligned_tau = max(tau[0], 4 * forward[0] / (3 * electron[1]))
    normalized_grade = hypot(transverse[0], tau[0])
    combined = None
    if aligned_tau <= tau[1]:
        combined = (1 + hypot(aligned_tau, transverse[0])) / 2
    # With a single complex residue all three normalized moments lie on a sphere.
    min_x_pure = max(transverse[0], sqrt(1 - tau[1]**2 - max(map(abs, normal))**2))
    return {"forward_predicted": expected_forward, "forward_measured": forward,
            "minimum_forward_transfer_gap": max(0, forward[0] - expected_forward[1]),
            "CTT_if_axial": near, "CTT_if_vector": [-near[1], -near[0]],
            "ALEPH_CTT": transverse, "ALEPH_CTN": normal,
            "mirror_minimum_additive_shift": transverse[0] + near[0],
            "unresolved_labels_axial_weight_lower": (1+transverse[0])/2,
            "residue_purity_lower": (1 + normalized_grade**2)/2,
            "largest_eigenweight_lower": (1 + normalized_grade)/2,
            "largest_eigenweight_with_forward_if_matched": combined,
            "pure_abs_v_over_a_upper": sqrt((1 - min_x_pure)/(1 + min_x_pure))}


def moments(a, v):
    joined = np.array([a, v], dtype=complex)
    gram = joined @ joined.conj().T
    gram /= np.trace(gram)
    return (float(2 * gram[1,0].real), float(2 * gram[0,0].real - 1),
            float(-2 * gram[1,0].imag))


def representations():
    # Basis D-, S-, T- with distinct orthonormal neutral partners D0, T0.
    third = np.diag([-0.5, 0., -1.])
    up_from_minus = np.array([[1., 0., 0.], [0., 0., sqrt(2)]])
    up_from_zero = np.array([[0., sqrt(2)]])
    points = []
    for x in [0., 0.125, 0.5]:
        light = np.sqrt([1 - 2*x, x, x])
        diagonal = float(light @ third @ light)
        off = third @ light - diagonal * light
        first = up_from_minus @ light
        second = up_from_zero @ first
        vals = [diagonal, float(off @ off), float(first @ first), float(second @ second)]
        assert np.allclose(vals, [-.5, x/2, 1., 4*x], atol=1e-15)
        points.append({"x": x, "T3_diag_offnorm2_raise1norm2_raise2norm2": vals})
    return points


def compute():
    a_tau = SOURCE["OPAL_Atau"][0]
    axial_ratio = a_tau / (1 + sqrt(1 - a_tau*a_tau))
    assert np.allclose(moments([1], [axial_ratio])[:2], [a_tau, sqrt(1 - a_tau**2)])
    assert np.allclose(moments([axial_ratio], [1])[:2], [a_tau, -sqrt(1 - a_tau**2)])
    cfg = SOURCE["two_width_witness"]
    fixed_a = cfg["Atau"]
    shapes = []
    for x in [sqrt(1-fixed_a**2), cfg["hidden_normalized_CTT"]]:
        matrix = .5*np.array([[1+x, fixed_a], [fixed_a,1-x]])
        vals, vecs = np.linalg.eigh(matrix)
        factors = vecs @ np.diag(np.sqrt(np.maximum(vals,0)))
        measured = moments(factors[0], factors[1])
        assert np.allclose(measured, [fixed_a,x,0.])
        afb = .75 * cfg["Ae"] * measured[0]
        for name, value in [("DELPHI_AFB0tau_9parameter", afb),
                            ("OPAL_Atau", measured[0]), ("ALEPH_CTT", measured[1]),
                            ("ALEPH_CTN", measured[2]), ("SLD_Ae_1996_98", cfg["Ae"])]:
            lo, hi = interval(name,2)
            assert lo <= value <= hi
        shapes.append({"moments_A_X_Y": measured, "normalized_eigenvalues": vals.tolist(),
                       "AFB": afb})
    return {"central_real_v_over_a": [axial_ratio,1/axial_ratio],
            "central_predicted_forward": .75 * SOURCE["SLD_Ae_1996_98"][0] * a_tau,
            "working_boxes": {str(n): polar_tests(n) for n in SOURCE["widths"]},
            "single_and_multiple_label_witnesses": shapes, "same_diagonal_representations": representations()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    out = HERE/"results.json"
    result = compute()
    if args.check:
        assert json.loads(out.read_text(encoding="utf-8-sig")) == json.loads(json.dumps(result))
        print("1506: four independent detector samples, residue and representation witnesses checked.")
    else:
        out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
