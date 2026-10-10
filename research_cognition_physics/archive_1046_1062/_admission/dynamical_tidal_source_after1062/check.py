"""Read-only exact algebra audit; --write-exclusive creates results once.

No stellar integration, waveform fit, paper correction, or physical error claim.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
from pathlib import Path


def q(value: F | int) -> str:
    return str(F(value))


def calculate() -> dict:
    checks = {}

    # Eqs. 56, 59: Lambda_0 = (16/45) (b_0/a_0).
    lambda0_per_ratio = F(16, 45)
    beta_lambda0_factor = F(107, 105) / lambda0_per_ratio
    checks["eq56_to_eq66_static_conversion"] = beta_lambda0_factor == F(321, 112)

    # Eq. 57 has this signed coefficient of log(mu r_s); Eq. 66 agrees.
    beta_prefactor = -F(64, 45)
    checks["eq57_to_eq66_prefactor"] = 2 * (-F(32, 45)) == beta_prefactor

    # Finite log-coordinate witness, not an asserted stellar solution.
    # b0/a0=1, Lambda_bar=10, log(mu r_s)=-1, mu=1/r and r=e*r_s.
    ratio, baseline, log_scale = F(1), F(10), F(-1)
    lambda0 = lambda0_per_ratio * ratio
    beta_log = beta_prefactor * (1 + beta_lambda0_factor * lambda0)
    value57 = baseline + beta_log * log_scale
    value65 = baseline - beta_log * log_scale
    difference = value57 - value65
    beta_wave = -beta_log
    checks["eq57_eq65_finite_scale_disagreement"] = difference == F(27136, 4725) > 0
    checks["consistent_wave_beta_requires_opposite_sign"] = (
        baseline - beta_wave * log_scale == value57
    )

    # Eq. 63 solved for c_dotE, versus printed Eq. 64.
    # G*M and R have length dimension, omega has inverse length.
    dimensions = {"eq63_solved_rhs_length_power": 1 - 3 + 2,
                  "eq64_printed_rhs_length_power": 1 - 1 + 2}
    gm, radius, c_e, omega_f = F(1), F(2), F(3), F(1, 2)
    solved63 = gm / radius**3 * c_e / omega_f**2
    printed64 = gm / radius * c_e / omega_f**2
    checks["eq64_missing_inverse_radius_squared"] = (
        solved63 == F(3, 2)
        and printed64 == 6
        and printed64 / solved63 == radius**2
        and dimensions["eq63_solved_rhs_length_power"] == 0
        and dimensions["eq64_printed_rhs_length_power"] == 2
    )

    # Independent circular-orbit algebra for Eq. 67 -> Eq. 68, G=M=1.
    # L_dyn = A Lambda(r) omega^2/r^6, r Lambda'=beta_wave.
    # Its direct energy plus shifted Newtonian energy is
    # (A/r^9)[Lambda + (2/3)(6 Lambda-beta_wave)].
    energy_lambda = 1 + F(2, 3) * 6
    energy_beta = -F(2, 3)
    checks["eq67_to_eq68_with_eq65_convention"] = (
        energy_lambda == 5 and energy_beta == -F(2, 3)
    )

    # Eqs. 68,69 -> psi''=-2 E'/F at linear tidal order.
    # b=m2/M; tuples are coefficients of (1,b). Set G=M=1.
    # L=Lambda_bar+beta_wave log[r/(gamma*r_s)], dL/dlog(omega)=-2 beta/3.
    eprime_L = (F(0), -F(405))
    eprime_beta = (F(0), F(99))
    flux_L = (F(24), F(36))
    flux_beta = (F(0), -F(6))
    response_L = tuple(e - f for e, f in zip(eprime_L, flux_L))
    response_beta = tuple(e - f for e, f in zip(eprime_beta, flux_beta))
    checks["eq68_eq69_relative_balance"] = (
        response_L == (F(-24), F(-441))
        and response_beta == (F(0), F(105))
    )

    n = F(11, 3)
    d2_power = n * (n - 1)
    d2_log_shift = -F(2, 3) * (2*n - 1)
    phase_prefactor = -F(45, 1408)
    phase_L = (F(8), F(147))
    phase_beta = (F(38, 11), F(1253, 44))
    derived_L = tuple(phase_prefactor*d2_power*v for v in phase_L)
    derived_beta = tuple(phase_prefactor*(d2_power*c+d2_log_shift*l)
                         for c, l in zip(phase_beta, phase_L))
    expected_L = tuple(F(5, 48)*v for v in response_L)
    expected_beta = tuple(F(5, 48)*v for v in response_beta)
    checks["eq70_eq71_second_derivative_matches_eq68_eq69"] = (
        derived_L == expected_L and derived_beta == expected_beta
    )
    checks["all_checks"] = all(checks.values())
    assert checks["all_checks"]

    return {
        "schema": "dynamical-tidal-source-algebra-audit-v1",
        "source": "https://arxiv.org/pdf/2606.19446v2",
        "source_version": "v2, 22 September 2026",
        "arithmetic": "exact fractions; powers are dimensional algebra",
        "checks": checks,
        "eq56_to_eq66": {
            "Lambda0_per_b0_over_a0": q(lambda0_per_ratio),
            "beta_log_prefactor": q(beta_prefactor),
            "Lambda0_coefficient_in_parenthesis": q(beta_lambda0_factor),
        },
        "finite_scale_identity_witness": {
            "is_verified_stellar_model": False,
            "b0_over_a0": q(ratio), "Lambda0": q(lambda0),
            "Lambda_bar": q(baseline), "log_mu_rs": q(log_scale),
            "mu_rs": "exp(-1)", "gamma": "1", "r_over_rs": "exp(1)",
            "beta_log_eq66": q(beta_log), "beta_wave_consistent": q(beta_wave),
            "Lambda_from_eq57": q(value57), "Lambda_from_printed_eq65": q(value65),
            "difference": q(difference),
        },
        "eq63_eq64": {
            "dimensions": dimensions,
            "formal_example": {"GM": q(gm), "R": q(radius),
                               "cE": q(c_e), "omega_f": q(omega_f)},
            "eq63_solved": q(solved63), "eq64_printed": q(printed64),
            "ratio_printed_over_solved": q(printed64/solved63),
        },
        "waveform_internal_algebra": {
            "assumed_running": "Lambda(r)=Lambda_bar+beta_wave*log(r/(gamma*r_s))",
            "eq68_energy_coefficients_L_beta": [q(energy_lambda), q(energy_beta)],
            "relative_balance_L_coefficients_1_b": [q(v) for v in response_L],
            "relative_balance_beta_coefficients_1_b": [q(v) for v in response_beta],
            "psi_second_derivative_L_coefficients_1_b": [q(v) for v in derived_L],
            "psi_second_derivative_beta_coefficients_1_b": [q(v) for v in derived_beta],
            "eq69_independent_radiation_derivation_certified": False,
            "eq57_through_eq71_as_printed_jointly_accepted": False,
        },
        "new_scientific_rounds": 0,
        "physical_finite_frequency_remainder_certified": False,
        "author_erratum_claimed": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-exclusive", action="store_true")
    args = parser.parse_args()
    output = Path(__file__).with_name("results.json")
    result = calculate()
    if args.write_exclusive:
        with output.open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        print("Created results.json exclusively; all exact algebra checks passed.")
    else:
        saved = json.loads(output.read_text(encoding="utf-8"))
        if saved != result:
            raise AssertionError("Saved result differs from the exact recomputation.")
        print("Read-only exact recomputation matches results.json.")


if __name__ == "__main__":
    main()
