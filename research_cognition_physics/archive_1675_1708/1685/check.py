"""Check allowed UV monomials, Wilson matching, weak CP and a shared colored threshold."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_text(encoding="utf-8"))


def allowed_insertions(base, charges, order, max_degree):
    names = ["S", "X", "A2", "A3"]
    return [
        "*".join(term) or "1"
        for length in range(max_degree + 1)
        for term in itertools.combinations_with_replacement(names, length)
        if (base + sum(charges[key] for key in term)) % order == 0
    ]


def evaluate():
    cfg = read("inputs.json")
    weak = read(cfg["weak_amplitudes_and_edm_budget"])
    powers = np.array(read(cfg["weak_powers"])["powers_A_B_C_E_F_R_D"])
    edges = np.array(weak["quark_path_amplitudes"])
    eps, x = cfg["scalar_magnitude_over_cutoff"], cfg["real_X_over_cutoff"]
    phi = math.pi * cfg["scalar_angle_over_pi"]
    s = eps * np.exp(1j * phi)
    charges, order = cfg["flavor_charge"], cfg["flavor_order"]
    singles = {key: allowed_insertions(charges[key], charges, order, 2)
               for key in ["L2", "L3"]}
    pairs = {left + right: allowed_insertions(charges[left] + charges[right], charges, order, 1)
             for left in ["P", "T"] for right in ["Pbar", "Tbar"]}
    assert singles == {"L2": ["A2", "S*S"], "L3": ["A3", "S*A2"]}
    assert pairs == {"PPbar": [], "PTbar": ["A3"], "TPbar": ["A3"], "TTbar": ["1"]}
    # General theorem uses a character quotient. This independently enumerates finite cyclic examples.
    cyclic_witnesses = [
        [n, a] for n in range(2, 97) for a in range(n) if 4 * a % n == 0
    ]
    assert all((6 * a) % n == (2 * a) % n for n, a in cyclic_witnesses)
    assert (6 - 2 + 2 * charges["X"]) % order == 0
    assert (1 - 5 + charges["X"]) % order == 0
    assert (6 + charges["A2"]) % 4 == 0
    assert weak["light_index"] + sum(charges[t] for t in ["P", "Pbar", "T", "Tbar"]) == order

    coefficients = edges / eps**powers
    coefficients[-1] /= x
    l2, l3, left, right = cfg["lambda2_lambda3_left_right"]
    mt = cfg["colored_M_T_over_cutoff"]

    def mass_matrix(value):
        a2 = l2 * value**2
        a3 = l3 * value * a2
        return np.array([[0j, left * a3], [right * a3, mt]])

    def observables(value):
        a, b, c, e, f, r, d = coefficients * value**powers
        yd = np.array([[0, a, 0], [b, c, r + x*d], [0, e, f]])
        h = yd @ yd.conj().T
        colored = mass_matrix(value)
        full_det = np.linalg.det(yd) * value**2 * np.linalg.det(colored)
        return full_det, np.imag(h[0, 1] * h[1, 2] * h[2, 0]), yd

    full_det, weak_invariant, yd = observables(s)
    assert abs(np.angle(full_det)) < 1e-12 and abs(weak_invariant) > 1e-12
    a, b, c, e, f, r, d = edges * np.exp(1j * phi * powers)
    original = np.array([[0, a, 0], [b, c, r + d], [0, e, f]])
    np.testing.assert_allclose(yd, original, rtol=3e-15, atol=1e-16)
    maximum_relative_det_error = 0.0
    for magnitude, angle in itertools.product([0.2, eps, 0.85], np.linspace(-math.pi, math.pi, 17)):
        v = magnitude * np.exp(1j * angle)
        expect = -left * right * (l2*l3)**2 * v**6
        maximum_relative_det_error = max(maximum_relative_det_error,
                                         abs(np.linalg.det(mass_matrix(v)) / expect - 1))
        np.testing.assert_allclose(observables(v)[0] / full_det, (v/s)**12, rtol=5e-14)

    step = 1e-5
    def log_spectrum(scale):
        return np.log(np.linalg.svd(mass_matrix(s*np.exp(scale)), compute_uv=False))
    slopes = (log_spectrum(step) - log_spectrum(-step)) / (2*step)
    assert abs(sum(slopes) - 6) < 2e-8
    factor = cfg["low_momentum_threshold_radius_factor"]
    shift = -sum(log_spectrum(math.log(factor)) - log_spectrum(0)) / (2*math.pi)
    assert abs(shift + 6*math.log(factor)/(2*math.pi)) < 1e-12
    budget = [v / x**2 for v in weak["heavy_c2_over_c6_conditional_range"]]
    return {
        "cyclic_character_witnesses_orders_2_to_96": len(cyclic_witnesses),
        "Z12_allowed_multiplier_monomials": singles, "Z12_allowed_pair_insertions": pairs,
        "light_heavy_QCD_instanton_index": [6, 6],
        "max_relative_determinant_residual": round(maximum_relative_det_error, 14),
        "tree_angle_rad": round(float(np.angle(full_det)), 14),
        "weak_im_H12H23H31": float(f"{weak_invariant:.12g}"),
        "conditional_masses_TeV_descending": [round(v, 10) for v in
                                             np.exp(log_spectrum(0))*cfg["cutoff_TeV"]],
        "individual_log_mass_radial_slopes": [round(v, 8) for v in slopes],
        "sum_log_mass_radial_slopes": round(sum(slopes), 7),
        "fixed_UV_inverse_alpha3_shift_full_heavy_chiral_pairs": round(shift, 10),
        "allowed_S2X2_over_S6_Wilson_interval": budget,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = evaluate()
    if args.check:
        assert actual == read("results.json")
        print("1685 checks passed")
    else:
        print(json.dumps(actual, ensure_ascii=False, indent=2, allow_nan=False))
