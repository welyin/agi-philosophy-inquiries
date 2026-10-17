"""Bounded same-action field-band screen; default is read-only verification.
No ODE integration, no complete constrained perturbation or thermal-state claim.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent.parent
SOURCES = {
    "1247/proof.md": "f3a79f8f43aa235d2fa13cc87a232e2b64643fc1bade375bfee96c204ab4bf3d",
    "1247/results.json": "dbcd590e7816bc423c49a4e1bf83bab6fa50a84dd5578249cc493a15cd72bf06",
    "1250/proof.md": "5e9f5e1495156f6bcd81c48a5046ac538ea72b2823a6cb7a9e3b7a59ff499d9f",
    "1250/results.json": "84fa683b1e343f94ba2f5a4f2e26653dfa49c7d6c04d3833a436daa8cb880284",
}


def sf(value: F) -> str:
    return str(value)


def calculate() -> dict:
    for relative, expected in SOURCES.items():
        actual = hashlib.sha256((ARCHIVE / relative).read_bytes()).hexdigest()
        assert actual == expected, (relative, actual, expected)
    original = json.loads((ARCHIVE / "1247/results.json").read_text(encoding="utf-8"))
    endpoint_upper = {}
    for branch in ("alpha0", "alpha1"):
        upper = F(original["group_finite_memory"][branch]["checkpoints"]["2"]["state"]["u"][1])
        assert upper < F(3, 8)
        endpoint_upper[branch] = sf(upper)

    # Exact coefficient audit of U_p after Hdot=-3(rho+P)/2 and H^2=rho.
    up_rho = F(3, 2) * F(-3, 2) + F(9, 4)
    up_pressure = F(3, 2) * F(-3, 2)
    assert up_rho == 0 and up_pressure == F(-9, 4)
    # U_r-U_p = 3H Fdot/(2F) + Fddot/(2F) - Fdot^2/(4F^2).
    # Fdot=v/5, Fddot=(-3Hv-p-3u^2/100-u^3/100)/5.
    ur_coeff = {
        "H_v_over_F": F(3, 2) * F(1, 5) + F(1, 2) * F(-3, 5),
        "p_over_F": F(1, 2) * F(-1, 5),
        "u_squared_over_F": F(1, 2) * F(-3, 500),
        "u_cubed_over_F": F(1, 2) * F(-1, 500),
        "v_squared_over_F_squared": F(-1, 4) * F(1, 25),
    }
    assert ur_coeff == {
        "H_v_over_F": F(0), "p_over_F": F(-1, 10),
        "u_squared_over_F": F(-3, 1000), "u_cubed_over_F": F(-1, 1000),
        "v_squared_over_F_squared": F(-1, 100),
    }

    # Pure rod expansion: -1/2 (Fbar + eps*s/5) (g + eps*d)^2.
    # Exact eps^2 coefficient is -Fbar*d^2/2 - s*g*d/5;
    # the second summand is the single cross term in README (2).
    Fbar, s, g, d = F(7, 5), F(2, 3), F(-3, 7), F(5, 11)
    eps2_by_product = F(-1, 2) * (Fbar * d * d + F(2, 5) * s * g * d)
    eps2_expected = -Fbar * d * d / 2 - s * g * d / 5
    assert eps2_by_product == eps2_expected
    assert F(2) * F(-1, 2) * F(1, 5) == F(-1, 5)

    rho_max, p_abs, f_min = F(111, 100), F(3, 2), F(7, 10)
    up_max = F(9, 4) * rho_max
    ur_max = up_max + p_abs / (10 * f_min)
    ur_simple = F(339, 125)
    gamma_sq_max = F(1, 10) / (25 * f_min)
    gamma_simple = F(2, 25)
    assert up_max == F(999, 400)
    assert ur_max < ur_simple
    assert gamma_sq_max == F(1, 175) < gamma_simple**2 < 1
    assert up_max - 1 < ur_simple

    q_wave = F(2)
    p_diag_min = q_wave**2 + 1 - up_max
    rod_diag_min = q_wave**2 - ur_max
    det_min = p_diag_min * rod_diag_min - gamma_sq_max * q_wave**2
    assert p_diag_min > 0 and rod_diag_min > 0 and det_min > 0
    frozen_det_at_two = q_wave**2 * (q_wave**2 + 1 - gamma_sq_max)
    assert frozen_det_at_two > 0
    # For q=2+y, y>=0: q^2-ur_simple-gamma_simple*q
    # = margin + (4-gamma_simple)*y + y^2.
    margin = q_wave**2 - ur_simple - gamma_simple * q_wave
    coeff_y = 2*q_wave - gamma_simple
    assert margin == F(141, 125) and coeff_y == F(98, 25)
    assert margin > 0 and coeff_y > 0

    u_lower = F(29, 40) * F(364929449, 10**9)
    pi_lower, pi_upper = F(157, 50), F(22, 7)
    q13_lower = 2*pi_lower*F(13, 10)*u_lower
    q13_upper = 2*pi_upper*F(13, 10)*F(3, 8)
    assert 2 < q13_lower < q13_upper == F(429, 140) < F(13, 4)

    return {
        "status": "PASS",
        "kind": "unnumbered_bounded_direction_screen",
        "source_sha256": SOURCES,
        "scope": {
            "fixed_original_metric_and_two_fluid_background": True,
            "field_field_quadratic_block_only": True,
            "full_constraint_perturbations_proved": False,
            "actual_quantum_preparation_proved": False,
            "KMS_Unruh_or_Jacobson_reference_proved": False,
            "compact_rod_zero_mode_claim": False,
            "ODE_reintegrated": False,
        },
        "exact_coefficient_audit": {
            "U_p_rho_coefficient": sf(up_rho),
            "U_p_pressure_coefficient": sf(up_pressure),
            "U_r_minus_U_p": {k: sf(v) for k, v in ur_coeff.items()},
            "rod_eps_squared_product": sf(eps2_by_product),
            "rod_eps_squared_expected": sf(eps2_expected),
            "rod_cross_coefficient": "-1/5",
            "rod_F_pp": "0",
        },
        "uniform_bounds_both_original_branches": {
            "interval": ["2", "9/4"], "U_p_upper": sf(up_max),
            "U_r_upper": sf(ur_max), "U_r_strict_simple_upper": sf(ur_simple),
            "gamma_squared_upper": sf(gamma_sq_max),
            "gamma_strict_simple_upper": sf(gamma_simple),
        },
        "fixed_q_two": {
            "p_diagonal_lower": sf(p_diag_min), "rod_diagonal_lower": sf(rod_diag_min),
            "determinant_lower": sf(det_min), "frozen_determinant_lower": sf(frozen_det_at_two),
        },
        "all_q_at_least_two": {
            "strict_eigenvalue_lower": sf(margin),
            "q_equals_two_plus_y_remainder_coefficients": {"y": sf(coeff_y), "y_squared": "1"},
        },
        "fixed_finite_mode_witness": {
            "real_sine_cosine_pair_n": [13, 0, 0], "torus_period": "10",
            "u_two_saved_upper": endpoint_upper, "u_hold_lower": sf(u_lower),
            "pi_rational_bounds_used": [sf(pi_lower), sf(pi_upper)],
            "physical_q_strict_lower": sf(q13_lower),
            "physical_q_strict_upper": sf(q13_upper),
            "included_in_fixed_band": ["2", "13/4"],
        },
        "limitation": "Positive field Hessian band is not a full constrained state, local apparatus, thermal reference, or independently closed quantum source certificate.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--save", action="store_true", help="first creation only; refuse overwrite")
    args = parser.parse_args()
    actual = calculate()
    result_path = HERE / "results.json"
    if args.save:
        with result_path.open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(actual, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        print("PASS: exact finite-band audit; results created once.")
    else:
        expected = json.loads(result_path.read_text(encoding="utf-8"))
        assert actual == expected, "saved result mismatch"
        print("PASS: read-only exact audit matches saved results; no ODE or full-model claim.")


if __name__ == "__main__":
    main()
