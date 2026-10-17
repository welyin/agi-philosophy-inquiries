"""Independent exact-arithmetic audit of round 1258.

Two groups: continuum energy/scaling algebra, and the full double-Dyson
majorant with a finite budget. No author code is imported. Universal field
claims belong to the analytic proof; finite coefficient checks below are
diagnostics of stated identities, not replacements for that proof.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def coeff_add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def coeff_scale(c, a):
    return tuple(c*x for x in a)


def calculate():
    # Independent derivation in the nonnegative monomial basis
    # (|u_dot|^2, |u|^2, (k^2/a^2)|u|^2), after dividing by a^3.
    energy = (F(1), F(1), F(1))
    # KG friction contributes -6*A; volume derivative contributes 3*e;
    # frequency derivative contributes -2*y*B.
    rate_over_H = coeff_add(coeff_scale(3, energy), (F(-6), F(0), F(-2)))
    assert rate_over_H == (-3, 3, 1)
    positive_margin = coeff_add(coeff_scale(3, energy), rate_over_H)
    negative_margin = coeff_add(coeff_scale(3, energy), coeff_scale(-1, rate_over_H))
    assert positive_margin == (0, 6, 4)
    assert negative_margin == (6, 0, 2)
    assert all(x >= 0 for x in (*positive_margin, *negative_margin))
    # Thus |e_dot| <= 3 |H| e independently of the mode number.

    interval = (F(2), F(9, 4))
    base_time = F(17, 8)
    half_window = max(abs(t-base_time) for t in interval)
    assert half_window == F(1, 8)
    H_bound = F(11, 10)
    gronwall_exponent = 3*H_bound*half_window
    assert gronwall_exponent == F(33, 80)
    # For 0<x<1, the positive exponential series is strictly below
    # the geometric series; this avoids a numerical exp comparison.
    B_rational_upper = 1/(1-gronwall_exponent)
    assert B_rational_upper == F(80, 47) < 2

    # b1=C exp[-1/(1-s^2)] on |s|<1. On |s|<=1/2,
    # the integral is >= exp(-4/3), hence C<9 from e<3.
    # With z=(1-s^2)^-1, |b1''| <=
    # C exp(-z)*(4*z^4+8*z^3+2*z^2).
    # sup z^m exp(-z)=(m/e)^m<(m/2)^m gives a strict bound.
    derivative_point_upper = sum(c*F(m, 2)**m
                                 for c, m in [(4, 4), (8, 3), (2, 2)])
    assert derivative_point_upper == 93
    b1_second_L1_upper = 9*2*derivative_point_upper
    assert b1_second_L1_upper == 1674
    h_L1_upper = 2*b1_second_L1_upper
    C2 = F(1)
    assert h_L1_upper == 3348
    # h=(d1^2-d2^2)(b1*b1*b1); zero and first moments follow by
    # compact smooth integration by parts. |b_hat|<=||b||_1=1,
    # hence |h_hat(xi)| <= |xi|^2 as well as ||h||_1.
    #
    # No cutoff: min(C0^2/y, C2^2*l^4*y) <= C0*C2*l^2.
    # Either y is below their crossing and the increasing term suffices,
    # or above it and the decreasing term suffices.
    ell = F(1, 16)
    C0 = F(h_L1_upper)
    crossing_y = C0/(C2*ell**2)
    spectral_majorant = C0*C2*ell**2
    assert C0*C0/crossing_y == spectral_majorant
    assert C2*C2*ell**4*crossing_y == spectral_majorant
    spatial_dimension = 3
    smear_prefactor_power = -3
    dt_power = 1
    dx_power = spatial_dimension
    amplitude_power = smear_prefactor_power+dt_power+dx_power
    variance_power = 2*amplitude_power+2
    assert amplitude_power == 1 and variance_power == 4

    volume, amin, chi_L1, E_star = F(1000), F(1), F(1), F(1)
    coefficient_with_B2 = 2*2*chi_L1**2*C0*C2/(volume*amin)
    assert coefficient_with_B2 == F(1674, 125) < 14
    coefficient_with_rational_B = 2*B_rational_upper*C0/(volume*amin)
    assert coefficient_with_rational_B == F(13392, 1175)
    D_upper = F(14)*E_star*ell**4
    rho_upper = F(2)*E_star/(volume*amin**3)
    assert rho_upper == F(1, 500)

    # A time center at 17/8 and |supp chi|<=1 gives a true interior
    # support for ell=1/16. This is support of the field test function,
    # not permission to apply an instantaneous localized detector gate.
    smear_window = (base_time-ell, base_time+ell)
    assert interval[0] < smear_window[0] < smear_window[1] < interval[1]

    # Independently recombine the two Dyson branches by total degree r.
    # Sum_{m+n=r}1/(m!n!)=2^r/r!, and the field-product error <=r^2 D/2.
    coefficient_diagnostics = []
    for r in range(21):
        direct = F(r*r, 2)*sum(
            F(1, math.factorial(m)*math.factorial(r-m))
            for m in range(r+1))
        collapsed = F(r*r*2**r, 2*math.factorial(r))
        assert direct == collapsed
        coefficient_diagnostics.append(str(direct))
    # Full-series identity follows from applying (z*d/dz)^2 to e^z:
    # sum r^2 z^r/r! = (z+z^2)e^z, z=2L.
    L = F(1, 8)
    polynomial = L*(1+2*L)
    assert polynomial == F(5, 32)
    exp_2L_upper = 1/(1-2*L)
    assert exp_2L_upper == F(4, 3)
    all_order_factor_upper = polynomial*exp_2L_upper
    assert all_order_factor_upper == F(5, 24)
    probability_upper = D_upper*all_order_factor_upper
    assert probability_upper == F(35, 786432)
    advertised_upper = F(3, 65536)
    assert probability_upper < advertised_upper < F(1, 10000)
    stronger_probability_upper = coefficient_with_B2*ell**4*all_order_factor_upper
    assert stronger_probability_upper < probability_upper

    return {
        "status": "PASS",
        "independent_groups": 2,
        "independent_check_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "group_energy_and_scaling": {
            "energy_basis": ["u_dot_squared", "u_squared", "k_squared_over_a_squared_times_u_squared"],
            "energy_coefficients": list(map(str, energy)),
            "rate_over_H_coefficients": list(map(str, rate_over_H)),
            "three_energy_plus_rate_coefficients": list(map(str, positive_margin)),
            "three_energy_minus_rate_coefficients": list(map(str, negative_margin)),
            "reference_time": str(base_time),
            "maximum_time_displacement": str(half_window),
            "gronwall_exponent_upper": str(gronwall_exponent),
            "B_strict_rational_upper": str(B_rational_upper),
            "b1_second_derivative_L1_strict_upper": str(b1_second_L1_upper),
            "h_L1_strict_upper": str(h_L1_upper),
            "quadratic_fourier_constant": str(C2),
            "spectral_crossing_k_squared": str(crossing_y),
            "smear_amplitude_scale_power": amplitude_power,
            "extra_variance_scale_power": variance_power,
            "variance_coefficient_strict_upper_B2": str(coefficient_with_B2),
            "variance_coefficient_strict_upper_rational_B": str(coefficient_with_rational_B),
            "chosen_ell": str(ell),
            "chosen_smear_time_support": list(map(str, smear_window)),
            "additional_occupation_source_density_strict_upper": str(rho_upper),
            "scope": "All modes are covered by the analytic piecewise spectral bound; the script checks algebra and constants, not the Hadamard property or physical preparation."
        },
        "group_full_dyson_and_budget": {
            "double_series_coefficients_checked_through": 20,
            "coefficient_diagnostics": coefficient_diagnostics,
            "full_series_identity": "sum r^2 z^r/r! = (z+z^2)exp(z), z=2L",
            "integrated_interaction_norm_L": str(L),
            "all_order_factor_strict_upper": str(all_order_factor_upper),
            "variance_strict_upper_at_E1": str(D_upper),
            "probability_strict_upper": str(probability_upper),
            "stronger_probability_strict_upper_B2": str(stronger_probability_upper),
            "advertised_probability_strict_upper": str(advertised_upper),
            "margin_to_advertised_upper": str(advertised_upper-probability_upper),
            "margin_advertised_to_one_ten_thousandth": str(F(1, 10000)-advertised_upper),
            "scope": "Arbitrary probe plus spectator input and bounded joint final effect, provided the field is initially independent and the two Gaussian states share CCR and mean with positive additional covariance."
        },
        "not_certified": [
            "Hadamard regularity from finite excitation energy alone",
            "automatic thermality, KMS qualification, or an Unruh ratio for the baseline",
            "implementation of spacetime-smear cosine as an instantaneous causal gate",
            "reference-independent comparison for amplitude-renormalized menus",
            "free preparation, readout, control, or vanishing device source",
            "curvature transport error or self-consistent geometric backreaction",
            "initial field correlations with the arbitrary probe or spectator"
        ]
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true",
                        help="exclusive first creation of independent_results.json")
    args = parser.parse_args()
    actual = calculate()
    path = HERE/"independent_results.json"
    if args.write:
        with path.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(actual, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        action = "FIRST_WRITE"
    else:
        expected = json.loads(path.read_text(encoding="utf-8"))
        assert expected == actual, "saved independent result differs"
        action = "READ_ONLY_MATCH"
    print(json.dumps({"status": "PASS", "action": action,
                      "independent_groups": 2,
                      "probability_strict_upper": str(F(35, 786432))}))


if __name__ == "__main__":
    main()
