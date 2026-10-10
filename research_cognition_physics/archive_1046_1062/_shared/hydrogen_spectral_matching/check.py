"""Count-0 algebra and published-budget audit; default execution is read-only.

No spectroscopy fit, QED calculation, confidence interval, or new empirical group.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
import json
from pathlib import Path


def frac(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def run() -> dict:
    A, B = F(3, 4), F(2, 9)
    a, b = F(-7, 8), F(-1, 8)
    r = B / A
    slope = b - r * a
    slope_13 = F(-26, 27) - F(8, 9) / A * a
    assert r == F(8, 27)
    assert slope == F(29, 216)
    assert slope_13 == F(2, 27)
    assert slope / slope_13 == F(29, 16)
    assert slope > 0 > b

    # A shared C is one nuisance, so this 2x2 covariance has rank one.
    covariance = ((a*a, a*b), (a*b, b*b))
    var_common = r*r*covariance[0][0] + covariance[1][1] - 2*r*covariance[0][1]
    var_wrong_independent = r*r*a*a + b*b
    assert covariance[0][0] * covariance[1][1] == covariance[0][1]**2
    assert var_common == slope**2
    assert var_wrong_independent > var_common

    # This bounded C is an algebra diagnostic, not a physical QED remainder.
    rho = F(1, 10**6)
    denom_fraction = F(7, 6) * rho
    uniform_remainder = abs(slope) * rho * denom_fraction / (1 - denom_fraction)
    finite_family = []
    for C in (-rho, F(0), rho):
        exact = A * (B + b*C) / (A + a*C) - B
        identity = slope*C / (1 - F(7, 6)*C)
        remainder = abs(exact - slope*C)
        assert exact == identity
        assert remainder <= uniform_remainder
        finite_family.append({"C": frac(C), "exact_normalized_shift": frac(exact),
                              "linearized_shift": frac(slope*C),
                              "absolute_remainder": frac(remainder)})

    # Independent finite numerator/denominator changes; set nu_A=A so s0=1.
    eta_A, eta_B = F(1, 10**5), F(2, 10**5)
    general_bound = eta_A * (eta_B + abs(r)*eta_A) / (A - eta_A)
    finite_corners = []
    for da in (-eta_A, eta_A):
        for db in (-eta_B, eta_B):
            actual = A*(B + db)/(A + da) - B
            linear = db - r*da
            identity = A*(db - r*da)/(A + da)
            assert actual == identity
            assert abs(actual - linear) <= general_bound
            finite_corners.append({"delta_A": frac(da), "delta_B": frac(db),
                                   "remainder": frac(abs(actual-linear))})

    with localcontext() as ctx:
        ctx.prec = 55
        D = Decimal
        components = {"QED": D(179), "muonic_radius": D(138),
                      "alpha": D('1.7'), "calibration_line": D(3),
                      "mass_ratio": D('0.0001')}
        variance_hz2 = sum(v*v for v in components.values())
        budget = variance_hz2.sqrt()
        combined = (variance_hz2 + D(480)**2).sqrt()
        frequency_hz = D('730690248610790')
        relative_ppt = combined/frequency_hz * D(10)**12
        assert abs(budget-D(226)) < D('0.1')
        assert abs(combined-D(530)) < D(1)

        # Muonic Table I: central values as printed, coefficient held fixed.
        eqed, ens, eexp, coeff = map(D, ('206.0344', '.0289', '202.3706', '5.2259'))
        radius_from_rounded_table = ((eqed+ens-eexp)/coeff).sqrt()
        assert abs(radius_from_rounded_table-D('.84060')) < D('.000005')
        covariance_ratio = (D(var_wrong_independent.numerator)/D(var_wrong_independent.denominator)
                            / (D(var_common.numerator)/D(var_common.denominator))).sqrt()
        budgets = {
            "inputs_standard_uncertainty_Hz": {k: str(v) for k, v in components.items()},
            "quadrature_variance_Hz2": str(variance_hz2),
            "quadrature_standard_uncertainty_Hz": str(budget),
            "combined_with_480Hz_Hz": str(combined),
            "combined_relative_ppt": str(relative_ppt),
            "printed_observation_minus_prediction_kHz": "0.00",
            "is_new_pvalue_or_hard_error_bound": False,
        }
        muonic = {"units_energy": "meV", "units_radius": "fm",
                   "radius_from_rounded_table_central_only": str(radius_from_rounded_table),
                   "published_radius": "0.84060(39)",
                   "recomputed_radius_uncertainty": False}
        ratio_std = str(covariance_ratio)

    return {
        "passed": True,
        "status": "mature adoption; algebra and published-budget audit only",
        "science_group_delta": 0,
        "empirical_group_delta": 0,
        "new_cognitive_axiom_delta": 0,
        "new_round": False,
        "exact_coefficients": {"A0": frac(A), "B0": frac(B), "B0_over_A0": frac(r),
                               "shared_C_calibrated_slope": frac(slope),
                               "fixed_scale_target_only_slope": frac(b),
                               "one_S_three_S_calibrated_slope": frac(slope_13),
                               "sensitivity_ratio": frac(slope/slope_13)},
        "shared_nuisance_covariance": {
            "rank_one": True, "variance_with_covariance": frac(var_common),
            "incorrect_independent_variance": frac(var_wrong_independent),
            "incorrect_to_correct_std_ratio": ratio_std},
        "algebra_C_diagnostic": {"abs_C_bound": frac(rho),
                                 "uniform_remainder_bound": frac(uniform_remainder),
                                 "is_physical_QED_bound": False, "samples": finite_family},
        "general_ratio_corner_diagnostic": {"bound": frac(general_bound),
                                              "samples": finite_corners},
        "published_budget": budgets,
        "muonic_theory_extraction_diagnostic": muonic,
        "not_performed": ["raw spectra or line-shape reanalysis", "full QED energies or covariance reconstruction",
                          "new empirical fit or p value", "62MB dataset download", "new-particle scan"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--save-exclusive', action='store_true',
                        help='Create results.json once; refuse to overwrite.')
    args = parser.parse_args()
    result = run()
    path = Path(__file__).with_name('results.json')
    if args.save_exclusive:
        with path.open('x', encoding='utf-8', newline='\n') as out:
            json.dump(result, out, ensure_ascii=False, indent=2)
            out.write('\n')
        print('Hydrogen spectral matching: saved exclusively; count 0.')
    else:
        saved = json.loads(path.read_text(encoding='utf-8'))
        if saved != result:
            raise AssertionError('Results differ from saved count-0 snapshot.')
        print('Hydrogen spectral matching: read-only checks passed; count 0.')


if __name__ == '__main__':
    main()
