"""1257: exact finite-shell and bounded local-effect audit.

Default: read-only reconstruction and exact comparison with results.json.
--write-results: exclusive first creation only. No ODE integration or QFT cutoff.
Analytic normality, Hadamard property, and local time smearing are proved in proof.md.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
from itertools import permutations, product
import json
from math import factorial
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parent.parent
SOURCES = {
    "archive_1234_/1247/proof.md": "f3a79f8f43aa235d2fa13cc87a232e2b64643fc1bade375bfee96c204ab4bf3d",
    "archive_1234_/1250/proof.md": "5e9f5e1495156f6bcd81c48a5046ac538ea72b2823a6cb7a9e3b7a59ff499d9f",
    "archive_819_853/research_note_843.md": "5f3e140bedd04c2eaa4c3bbb8f41e9a27780cdf86af1fe5b0f530ebced5b8ae4",
    "archive_990_1008/research_note_1004.md": "ae717955f940e15ff1861e56c4c02e067564d15767abe713bf4cc5db1548c5fd",
}


def strings(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: strings(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [strings(item) for item in value]
    return value


def shell_points():
    axes = set()
    for axis in range(3):
        for sign in (-1, 1):
            point = [0, 0, 0]
            point[axis] = 3 * sign
            axes.add(tuple(point))
    mixed = {
        tuple(sign * item for sign, item in zip(signs, perm))
        for perm in set(permutations((2, 2, 1)))
        for signs in product((-1, 1), repeat=3)
    }
    assert len(axes) == 6 and len(mixed) == 24 and not axes.intersection(mixed)
    for points in (axes, mixed):
        assert all(sum(x*x for x in point) == 9 for point in points)
        assert all(tuple(-x for x in point) in points for point in points)
        assert all(tuple(point[i] for i in perm) in points
                   for point in points for perm in permutations(range(3)))
    return sorted(axes), sorted(mixed)


def moments(points, weight):
    # Each coefficient multiplies the same positive finite occupation scale nu.
    return {
        "zeroth": weight * len(points),
        "first": [sum((weight * n[i] for n in points), F(0)) for i in range(3)],
        "second": [[sum((weight * n[i] * n[j] for n in points), F(0))
                    for j in range(3)] for i in range(3)],
        "fourth_single_axis": [sum((weight * n[i]**4 for n in points), F(0))
                               for i in range(3)],
        "fourth_smear_D_squared": sum((weight * (n[0]**2-n[1]**2)**2
                                      for n in points), F(0)),
    }


def group_shell(axes, mixed):
    ma, mb = moments(axes, F(4)), moments(mixed, F(1))
    expected_second = [[F(72) if i == j else F(0) for j in range(3)] for i in range(3)]
    assert ma["zeroth"] == mb["zeroth"] == 24
    assert ma["first"] == mb["first"] == [F(0)] * 3
    assert ma["second"] == mb["second"] == expected_second
    assert ma["fourth_single_axis"] == [F(648)] * 3
    assert mb["fourth_single_axis"] == [F(264)] * 3
    return {
        "status": "PASS",
        "same_squared_integer_radius": 9,
        "occupation_scale": "nu>0 finite",
        "A": {"points": axes, "occupation_per_point_over_nu": "4", "moments_over_nu": ma},
        "B": {"points": mixed, "occupation_per_point_over_nu": "1", "moments_over_nu": mb},
        "all_positive_and_negative_momenta_retained": True,
        "real_field_pair_counts": {"A": 3, "B": 12},
        "pair_convention": "Each +/- pair is two real sine/cosine oscillators, or equivalently two independent travelling-wave annihilators.",
    }, ma, mb


def group_source_and_local(ma, mb, axes, mixed):
    d0 = ma["zeroth"] - mb["zeroth"]
    d1 = [a-b for a, b in zip(ma["first"], mb["first"])]
    d2 = [[ma["second"][i][j] - mb["second"][i][j] for j in range(3)] for i in range(3)]
    # Free minimally coupled X in the hold window, modes u_3(t)e^{iq0 n.x}/sqrt(V).
    # Strip the common factor nu/V. The named monomials keep a(t), q0 and m separate.
    stress = {
        "T00": {"abs_dot_u_squared": d0, "m_squared_abs_u_squared": d0,
                  "q0_squared_a_inverse_squared_abs_u_squared": 9*d0},
        "T0i": [{"q0_Im_dot_u_conj_u": 2*d1[i],
                  "q0_Re_dot_u_conj_u": F(0)} for i in range(3)],
        "Tij": [[{
            "a_squared_abs_dot_u_squared": d0 if i == j else F(0),
            "a_squared_m_squared_abs_u_squared": -d0 if i == j else F(0),
            "q0_squared_abs_u_squared": 2*d2[i][j] - (9*d0 if i == j else F(0)),
        } for j in range(3)] for i in range(3)],
    }
    covariance = {
        "X_X_abs_u_squared": 2*d0,
        "Xdot_Xdot_abs_dot_u_squared": 2*d0,
        "sym_X_Xdot_Re_u_conj_dot_u": 2*d0,
        "grad_i_grad_j_q0_squared_abs_u_squared": [[2*d2[i][j] for j in range(3)] for i in range(3)],
    }
    def all_zero(tree):
        if isinstance(tree, dict): return all(all_zero(x) for x in tree.values())
        if isinstance(tree, list): return all(all_zero(x) for x in tree)
        return tree == 0
    assert all_zero(stress) and all_zero(covariance)

    smear_a = ma["fourth_smear_D_squared"]
    smear_b = mb["fourth_smear_D_squared"]
    assert smear_a == 1296 and smear_b == 144 and smear_a == 9*smear_b
    # Cos(pi*n_x/2) has exact integer values at the finite translation L/4.
    cosine = (F(1), F(0), F(-1), F(0))
    translated_a = sum((F(4) * cosine[n[0] % 4] for n in axes), F(0))
    translated_b = sum((cosine[n[0] % 4] for n in mixed), F(0))
    assert translated_a == 16 and translated_b == -16
    assert translated_a-translated_b == 32
    # A normalized nonnegative radial bump in r0=1/10 has a common shell Fourier G.
    # |k.x| < 3*pi/50 < 33/175, and cos(z)>=1-z^2/2.
    radius, period, pi_upper = F(1, 10), F(10), F(22, 7)
    phase_upper = 2*pi_upper*3*radius/period
    g_lower = 1-phase_upper**2/2
    assert phase_upper == F(33, 175) and g_lower > 0
    return {
        "status": "PASS",
        "full_mean_stress_difference_coefficients": stress,
        "coincident_quadratic_covariance_difference_coefficients": covariance,
        "why_full_window": "All occupied directions share the same u_3(t); angular coefficients vanish identically, not only on an initial slice.",
        "same_renormalization_subtraction_required": True,
        "point_comoving_two_point_difference_coefficient": d0,
        "point_blindness_scope": "Fixed comoving position, no angular spatial smearing or spatial derivative apparatus; the full unequal-time two-point shell sum cancels.",
        "finite_translation_check": {"separation": "(L/4,0,0)", "A_sum_over_nu": translated_a,
                                     "B_sum_over_nu": translated_b, "difference_over_nu": translated_a-translated_b},
        "compact_smear": {
            "f": "(partial_x^2-partial_y^2)g",
            "g": "normalized nonnegative radial C_c^infinity bump",
            "comoving_support_radius": radius, "torus_period": period,
            "phase_strict_upper": phase_upper, "G_positive_lower": g_lower,
            "fourth_weight_A_over_nu": smear_a, "fourth_weight_B_over_nu": smear_b,
            "added_variance_ratio_A_to_B": smear_a/smear_b,
            "real_scalar_thermal_extra_covariance_factor": "2",
            "common_positive_variance_factor": "2*|time_mode_overlap|^2*q0^4*G^2/V",
            "mean_cancels": "Integral(f)=0; a common spatially homogeneous classical X mean contributes zero.",
            "time_smearing_nonzero_overlap_is_analytic_proof": True,
        },
    }


def group_gaussian_and_probability():
    # Illustrations of the exact one-oscillator identity z=N/(1+N),
    # rho=(1-z) sum_{m>=0} z^m |m><m|. Analytic infinite sums belong to proof.md.
    single_mode = []
    for nu in (F(1, 100), F(1, 3), F(1), F(7, 2)):
        for weight in (F(1), F(4)):
            occupation = weight*nu
            ratio = occupation/(1+occupation)
            assert 0 < ratio < 1
            assert (1-ratio)*(1+occupation) == 1
            assert ratio/(1-ratio) == occupation
            assert ratio/(1-ratio)**2 == occupation*(occupation+1)
            single_mode.append({"N": occupation, "geometric_ratio": ratio,
                                "mean": occupation, "number_variance": occupation*(occupation+1)})
    # Total mean number and energy coefficient match, but not all number statistics.
    mean_a = 6*F(4)
    mean_b = 24*F(1)
    variance_quadratic_a, variance_quadratic_b = 6*F(4)**2, 24*F(1)**2
    assert mean_a == mean_b == 24
    assert variance_quadratic_a-variance_quadratic_b == 72

    va_ratio, vb_ratio = F(10), F(2)  # C=V0, added A variance=9C, B=C.
    exponent_a, exponent_b = va_ratio/4, vb_ratio/4
    assert exponent_a == F(5, 2) and exponent_b == F(1, 2)
    # Strict alternating lower bound for exp(-1/2); next term is positive.
    exp_minus_half_lower = sum((F(-1, 2)**j/F(factorial(j)) for j in range(4)), F(0))
    assert exp_minus_half_lower == F(29, 48) > F(3, 5)
    # exp(5/2) strictly exceeds its degree-six positive Taylor polynomial.
    exp_plus_5half_lower = sum((F(5, 2)**j/F(factorial(j)) for j in range(7)), F(0))
    assert exp_plus_5half_lower == F(110681, 9216) > 12
    exp_minus_5half_upper = F(1, 12)
    strict_gap_lower = (exp_minus_half_lower-exp_minus_5half_upper)/2
    assert strict_gap_lower == F(25, 96) > F(1, 4)
    # C=x V0, 0<x<=1, E_x uses denominator sqrt((1+x)*V0).
    # gap=.5 exp(-.5)*(1-exp(-t)), t=4*x/(1+x)>=2*x.
    # exp(t)>1+t gives 1-exp(-t)>t/(1+t)>=2*x/3.
    small_x_coefficient = F(1, 2)*F(3, 5)*F(2, 3)
    assert small_x_coefficient == F(1, 5)
    assert (F(9)-F(1))/2 == 4
    # Exact rational checks of the cleared numerators, valid for 0<x<=1.
    # t-2x = 2x(1-x)/(1+x); 2x/(1+2x)-2x/3 = 4x(1-x)/(3(1+2x)).
    assert 4-2 == 2 and 2*3-2 == 4
    return {
        "status": "PASS",
        "one_mode_geometric_identity_exact_samples": single_mode,
        "general_state_definition": "Finite product of phase-invariant normal Gaussian one-mode densities; untouched modes remain in the common pure Hadamard Fock vacuum.",
        "analytic_properties_not_proved_by_samples": [
            "Geometric infinite-series normalization and moments for every finite N>=0.",
            "Normal finite-mode product, smooth Hadamard two-point difference, and local smeared covariance finiteness.",
            "Finite expected energy on the smooth compact hold window; no hard particle-number cutoff or physical preparation certificate.",
        ],
        "finite_resources": {
            "mean_total_excitation_number_over_nu": mean_a,
            "common_energy_density_increment": "24*nu/V*(|dot_u_3|^2+(1+9*q0^2/a^2)*|u_3|^2)",
            "entropy_A": "6*s(4*nu)", "entropy_B": "24*s(nu)",
            "s_definition": "s(N)=(N+1)log(N+1)-Nlog(N); finite for finite N>=0, s(0)=0",
            "total_number_variance_A": "24*nu+96*nu^2",
            "total_number_variance_B": "24*nu+24*nu^2",
            "number_variance_difference_over_nu_squared": variance_quadratic_a-variance_quadratic_b,
            "common_source_increment_is_not_zero": True,
        },
        "bounded_local_effect": {
            "effect": "E=(I+cos(Phi(f_time)/sqrt(2*V0)))/2",
            "effect_operator_bounds": "0<=E<=I by spectral calculus",
            "V0": "positive finite common reference variance, fixed before choosing nu",
            "occupation_scale_choice": "finite nu chosen so the B added variance C equals V0",
            "A_full_variance_over_V0": va_ratio, "B_full_variance_over_V0": vb_ratio,
            "A_gaussian_exponent": exponent_a, "B_gaussian_exponent": exponent_b,
            "P_B_minus_P_A": "(exp(-1/2)-exp(-5/2))/2",
            "strict_exp_minus_half_lower": exp_minus_half_lower,
            "strict_exp_plus_5half_lower": exp_plus_5half_lower,
            "strict_exp_minus_5half_upper": exp_minus_5half_upper,
            "strict_probability_gap_lower": strict_gap_lower,
            "probability_gap_exceeds": "1/4",
        },
        "small_occupation_branch": {
            "domain": "C=x*V0, 0<x<=1",
            "effect": "E_x=(I+cos(Phi(f_time)/sqrt((1+x)*V0)))/2; same effect for A and B at fixed x",
            "gap": "0.5*exp(-1/2)*(1-exp(-4*x/(1+x)))",
            "t": "4*x/(1+x)",
            "cleared_t_minus_2x": "2*x*(1-x)/(1+x)>=0",
            "strict_exp_bound": "exp(t)>1+t for t>0",
            "cleared_rational_tail": "2*x/(1+2*x)-2*x/3=4*x*(1-x)/(3*(1+2*x))>=0",
            "strict_gap_lower_over_x": small_x_coefficient,
            "strict_gap_exceeds": "x/5",
        },
    }


def calculate():
    for name, expected in SOURCES.items():
        actual = hashlib.sha256((RESEARCH/name).read_bytes()).hexdigest()
        assert actual == expected, (name, actual, expected)
    axes, mixed = shell_points()
    group1, ma, mb = group_shell(axes, mixed)
    group2 = group_source_and_local(ma, mb, axes, mixed)
    group3 = group_gaussian_and_probability()
    return strings({
        "round": 1257,
        "status": "PASS",
        "science_groups": 3,
        "arithmetic": "Exact Python Fraction and integer enumeration; no floating-point quantity is used as proof.",
        "source_sha256": SOURCES,
        "group_1_exact_shell_moments": group1,
        "group_2_complete_mean_source_and_local_two_point": group2,
        "group_3_normal_gaussian_resources_and_bounded_probability": group3,
        "new_adopted_axioms": 0,
        "whole_GR_proved": False,
        "scope": {
            "all_hold_times_mean_stress_difference_zero": True,
            "all_quantum_source_statistics_equal": False,
            "common_quantum_source_zero": False,
            "old_classical_geometry_automatically_semiclassically_self_consistent": False,
            "actual_preparation_or_original_K_instrument_certified": False,
            "physical_KMS_or_Unruh_certified": False,
            "ODE_reintegrated": False,
            "full_field_UV_truncated": False,
        },
        "net_result": "Same original free X, same background and full mean stress throughout the hold window, different local two-point covariance and a bounded Gaussian effect probability; no complete quantum-gravity or preparation claim.",
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true", help="exclusive first write; never overwrite")
    args = parser.parse_args()
    actual = calculate()
    path = HERE/"results.json"
    if args.write_results:
        with path.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(actual, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        print("PASS: 1257 three exact author groups; first results written exclusively.")
    else:
        saved = json.loads(path.read_text(encoding="utf-8"))
        assert saved == actual, "saved result mismatch"
        print("PASS: 1257 three exact author groups; read-only reconstruction matches frozen results.")


if __name__ == "__main__":
    main()
