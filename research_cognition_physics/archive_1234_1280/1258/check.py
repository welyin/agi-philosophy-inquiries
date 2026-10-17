"""1258: exact finite checks for an energy-bounded local Gaussian reference family.

Default reconstructs and compares results read-only. --write-results exclusively
creates the first results file. No ODE integration, mode truncation, quadrature,
Hadamard construction, or physical-probe implementation is performed here.
The all-mode estimates and norm-convergent Dyson argument are in proof.md.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from math import factorial
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parent.parent
SOURCES = {
    "archive_1234_/1247/proof.md": "f3a79f8f43aa235d2fa13cc87a232e2b64643fc1bade375bfee96c204ab4bf3d",
    "archive_1234_/1250/proof.md": "5e9f5e1495156f6bcd81c48a5046ac538ea72b2823a6cb7a9e3b7a59ff499d9f",
    "archive_1234_/1257/proof.md": "bae5e521022c6cea650e91f9e5c9cf9fe789936fac2d84460f455722dc1163f6",
    "archive_702_741/research_note_730.md": "e0ea20e66e44dfb0278efa85d1508435362cc21cd0baac3f4fbb092be1db4a83",
    "archive_990_1008/research_note_1004.md": "ae717955f940e15ff1861e56c4c02e067564d15767abe713bf4cc5db1548c5fd",
    "archive_1219_1233/1226/proof.md": "9d5558d399c0059d5f6d3ac8f0a077111f9702588edaadcfb0626458e1d42eb1",
}


def strings(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: strings(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [strings(item) for item in value]
    return value


def geometric_exp_upper(x):
    # For 0<x<1, exp(x)=sum x^n/n! < sum x^n=1/(1-x).
    assert 0 < x < 1
    return 1/(1-x)


def group_energy_modes_sources():
    # Strip a^3 and H; basis A=|u_dot|^2, M=|u|^2, Q=k^2*a^-2*|u|^2.
    e = [F(1), F(1), F(1)]
    derivative_over_H = [F(-3), F(3), F(1)]
    upper_margin = [3*x-y for x, y in zip(e, derivative_over_H)]
    lower_margin = [3*x+y for x, y in zip(e, derivative_over_H)]
    assert upper_margin == [6, 0, 2]
    assert lower_margin == [0, 6, 4]
    assert all(x >= 0 for x in upper_margin + lower_margin)

    center, left, right = F(17,8), F(2), F(9,4)
    max_delta = max(center-left, right-center)
    hmax = F(11,10)
    exponent = 3*hmax*max_delta
    assert exponent == F(33,80)
    exp_upper = geometric_exp_upper(exponent)
    assert exp_upper == F(80,47) < 2

    # A spatial pressure in a unit direction has coefficients (1,-1,2r-1),
    # with r = squared parallel wave-number / squared wave-number in [0,1].
    # rho-P = 2M+2(1-r)Q; rho+P = 2A+2rQ.
    for r in [F(0), F(1,3), F(1,2), F(1)]:
        pressure = [F(1), F(-1), 2*r-1]
        assert [x-y for x, y in zip(e, pressure)] == [0, 2, 2*(1-r)]
        assert [x+y for x, y in zip(e, pressure)] == [2, 0, 2*r]
    # Flux bound 2 sqrt(AQ)<=A+Q<=A+M+Q; null derivative bound
    # 2(sqrt(A)+sqrt(Q))^2<=4(A+Q)<=4 rho per occupation.
    # Exact square margins in variables x=sqrt(A), y=sqrt(Q).
    flux_margin_coefficients = [F(1), F(-2), F(1)]
    null_margin_coefficients = [F(2), F(-4), F(2)]
    assert null_margin_coefficients == [2*x for x in flux_margin_coefficients]

    return {
        "status": "PASS",
        "mode": "u_ddot+3H*u_dot+(1+k^2/a^2)*u=0; mass=1; V=1000",
        "physical_mode_energy": "e_k=a^3*(|u_dot|^2+(1+k^2/a^2)*|u|^2)",
        "energy_budget": "sum_k n_k*e_k(t_ref)<=E_star; n_k>=0 diagonal Gaussian, all modes retained",
        "t_ref": center,
        "energy_derivative_basis": ["|u_dot|^2", "|u|^2", "k^2*a^-2*|u|^2"],
        "e_dot_over_H_a3": derivative_over_H,
        "three_e_minus_derivative_margin": upper_margin,
        "three_e_plus_derivative_margin": lower_margin,
        "transport": {
            "statement": "|e_dot|<=3|H|e, hence e_k(t)<=B e_k(t_ref)",
            "window": [left, right], "H_abs_upper": hmax,
            "max_time_separation_from_reference": max_delta,
            "B_exponent_upper": exponent,
            "B_strict_geometric_upper": exp_upper,
            "B_coarse_strict_upper": 2,
        },
        "all_mode_covariance": {
            "z_k": "ell/sqrt(V) * exp(i k.x0) * h_hat(ell*k) * integral chi(s)u_k(tau+ell*s)ds",
            "variance_increment": "Delta V=2 sum_k n_k |z_k|^2",
            "Fourier_bound": "|h_hat(q)|<=min(C0,C2*|q|^2)",
            "mixed_Fourier_bound": "|h_hat(ell*k)|^2<=C0*C2*ell^2*|k|^2",
            "frequency_ratio": "|k|^2/(a^3*(1+|k|^2/a^2))<=1/a<=1/a_min",
            "D": "2*B*||chi||_1^2*C0*C2*E_star*ell^4/(V*a_min)",
            "cross_covariance": "|Delta W(F_i,F_j)|<=sqrt(Delta V_i Delta V_j)<=D for the common scale family",
            "zero_mode": "h_hat(0)=0; no zero-mode removal from the field",
            "finite_shell_upgrade": "fixed k shell gives O(E_star*ell^6), while all-mode uniform bound is O(E_star*ell^4)",
        },
        "complete_mean_excitation_source": {
            "energy_density": "0<=Delta rho(t)=sum_k n_k e_k(t)/(V*a(t)^3)<=B*E_star/(V*a_min^3)",
            "flux": "Euclidean norm of orthonormal-frame Delta T_0i <= Delta rho, including non-even occupation families",
            "pressure": "every spatial unit-direction pressure obeys |Delta T_vv|<=Delta rho; hence all orthonormal spatial components bounded by Delta rho",
            "pressure_margins": ["rho-P_v=2[M+(1-r)Q]", "rho+P_v=2[A+rQ], 0<=r<=1"],
            "flux_square_margin": flux_margin_coefficients,
            "null_square_margin": null_margin_coefficients,
            "null": "0<=Delta T_ab k^a k^b<=4*Delta rho*(k^hat0)^2 for null k",
            "renormalization": "state difference at a common subtraction prescription; not the reference vacuum source",
            "normality": "mass 1 and the Wronskian imply e_k>=omega_k>=1, so sum n_k<=E_star at t_ref",
            "Hadamard_scope": "finite energy alone does not imply Hadamard; retain Hadamard-qualified members if that qualification is required",
        },
    }


def group_bump_constants():
    # b(s)=C exp[-1/(1-s^2)] on (-1,1), normalized. On [-1/2,1/2]
    # integral unnormalized b >= exp(-4/3); exp(4/3)<(1/(1-2/3))^2=9.
    C_upper = geometric_exp_upper(F(2,3))**2
    assert C_upper == 9
    # v=1/(1-s^2). With |s|<=1, sup v^j exp(-v) <= (j/e)^j.
    # e>2 yields 4(4/e)^4+8(3/e)^3+2(2/e)^2 <93.
    polynomial_sup_upper = 4*F(4,2)**4 + 8*F(3,2)**3 + 2*F(2,2)**2
    assert polynomial_sup_upper == 93
    derivative_L1_upper = 2*polynomial_sup_upper*C_upper
    C0_upper = 2*derivative_L1_upper
    assert derivative_L1_upper == 1674
    assert C0_upper == 3348
    # Independently differentiate the bump logarithm with a rational test point.
    # log b'=-2s v^2 and log b''=-2v^2-8s^2 v^3;
    # b''/b=(log b')^2+log b''.
    derivative_samples = []
    for s in [F(0), F(1,4), F(1,2), F(3,4)]:
        v = 1/(1-s*s)
        logarithm_first = -2*s*v*v
        logarithm_second = -2*v*v-8*s*s*v**3
        value = logarithm_first**2 + logarithm_second
        formula = 4*s*s*v**4-8*s*s*v**3-2*v*v
        assert value == formula
        derivative_samples.append({"s": s, "v": v, "b_second_over_b": value})
    C2 = F(1)
    norm_chi = F(1)
    volume, amin = F(1000), F(1)
    D_coefficient = 2*2*norm_chi**2*C0_upper*C2/(volume*amin)
    assert D_coefficient == F(1674,125) < 14
    return {
        "status": "PASS",
        "b1": "C*exp(-1/(1-s^2)) for |s|<1, zero otherwise, integral b1=1",
        "normalization_strict_upper": C_upper,
        "normalization_argument": "integral[-1/2,1/2] exp(-1/(1-s^2))ds>=exp(-4/3), and exp(4/3)<9",
        "second_derivative": "b1''=C exp(-v)*(4s^2 v^4-8s^2 v^3-2v^2), v=(1-s^2)^-1",
        "second_derivative_rational_identity_samples": derivative_samples,
        "sup_abs_second_derivative_over_C_strict_upper": polynomial_sup_upper,
        "b1_second_derivative_L1_strict_upper": derivative_L1_upper,
        "spatial_base": "b(x)=b1(x1)b1(x2)b1(x3), normalized and nonnegative",
        "h": "partial_1^2 b - partial_2^2 b",
        "moments": "integral h=0 and integral x_i h=0 by compact-support integration by parts",
        "C0_L1_strict_upper": C0_upper,
        "C2_Fourier_bound": C2,
        "C2_qualification": "Fourier identity h_hat(q)=(q2^2-q1^2)b_hat(q), |b_hat|<=1; this is NOT a claim that the Taylor second-moment integral equals 1",
        "chi": "b1, with L1 norm 1",
        "D_strict_coefficient_over_E_star_ell4": D_coefficient,
        "D_coarse_strict_coefficient_over_E_star_ell4": 14,
        "no_numerical_quadrature": True,
    }, D_coefficient


def group_dyson_finite_budget(D_coefficient):
    # For r letters, Weyl expansion of products of cos has absolute weights
    # totaling 1. A positive covariance increment bounded by D per letter has
    # |delta variance|<=r^2 D. The common CCR phase is unchanged.
    # Exact coefficient checks supplement the analytic identity valid for all r:
    # sum_{n+m=r} (r^2/2)/(n!m!) = [L^r] L(1+2L) exp(2L).
    coefficient_checks = []
    for r in range(21):
        inverse_factorial_sum = sum((F(1, factorial(n)*factorial(r-n)) for n in range(r+1)), F(0))
        assert inverse_factorial_sum == F(2**r, factorial(r))
        left = F(r*r,2)*inverse_factorial_sum
        right = F(2**(r-1), factorial(r-1)) if r >= 1 else F(0)
        if r >= 2:
            right += F(2**(r-1), factorial(r-2))
        assert left == right
        assert r*r == r+r*(r-1)
        coefficient_checks.append({"r": r, "combined_coefficient": left})
    L = F(1,8)
    exp_upper = geometric_exp_upper(2*L)
    assert exp_upper == F(4,3)
    factor_upper = L*(1+2*L)*exp_upper
    assert factor_upper == F(5,24)
    tight_probability_coefficient = D_coefficient*factor_upper
    coarse_probability_coefficient = F(14)*factor_upper
    assert tight_probability_coefficient == F(279,100) < 3
    assert coarse_probability_coefficient == F(35,12) < 3

    ell, E_star = F(1,16), F(1)
    center = F(17,8)
    probe_radius = F(1,32)
    all_support = [center-probe_radius-ell, center+probe_radius+ell]
    assert all_support == [F(65,32), F(71,32)]
    assert F(2) < all_support[0] < all_support[1] < F(9,4)
    variance_upper = F(14)*E_star*ell**4
    probability_upper = F(3)*E_star*ell**4
    rho_upper = F(2)*E_star/F(1000)
    assert probability_upper == F(3,65536) < F(1,10000)
    assert rho_upper == F(1,500)
    # No matching-state normalization or variable sensor gain is used here.
    return {
        "status": "PASS",
        "bounded_test_generator": "lambda*c(tau)*M(tau) tensor cos(Phi(F_tau,ell)), ||M(tau)||<=1",
        "test_probe_scope": "norm-convergent bounded Dyson test; not an autonomous instantaneous implementation of time-smeared fields in the original mother action",
        "common_CCR": "Weyl product phase is state-independent and equal for both compared states",
        "single_word_difference": "at r field-cos letters, absolute expectation difference <=r^2 D/2",
        "Dyson_integrated_strength": "L=integral |lambda*c(tau)|dtau, independent of smear-internal chi",
        "all_order_identity": "sum_{n,m>=0} ((n+m)^2/2)*L^(n+m)/(n!m!) = L*(1+2L)*exp(2L)",
        "exact_identity_coefficient_checks": coefficient_checks,
        "all_order_proof": "binomial factorial identity and r^2=r+r(r-1); differentiate exp(2L)",
        "fixed_probe_envelope": "c(tau)=32*b1(32*(tau-17/8)), ||c||1=1",
        "lambda_and_L": L,
        "exp_two_L_strict_upper": exp_upper,
        "Dyson_factor_strict_upper": factor_upper,
        "probability_strict_coefficient_over_E_star_ell4": tight_probability_coefficient,
        "using_coarse_D_probability_coefficient": coarse_probability_coefficient,
        "reported_probability_strict_coefficient": 3,
        "finite_nonempty_contract": {
            "E_star": E_star, "ell": ell,
            "tau_center": center, "tau_radius": probe_radius,
            "all_field_time_support_envelope": all_support,
            "field_time_support_strictly_inside_hold": True,
            "spatial_support": "translated cube of side 2 ell on the period-10 torus; fixed patch chosen without wrap",
            "added_variance_strict_upper": variance_upper,
            "probe_probability_error_strict_upper": probability_upper,
            "requested_probability_budget": "1/10000",
            "excitation_energy_density_strict_upper": rho_upper,
            "orthonormal_flux_and_each_pressure_component": "bounded by the same Delta rho",
            "null_source_strict_upper_over_khat0_squared": 4*rho_upper,
            "nonvacuum_members_exist": "any finite-mode Gaussian occupation with 0<n_k<E_star/e_k(t_ref); no infinite energy required",
        },
        "finite_shell_interpretation": "For fixed shell and fixed gain, nonzero fixed probability separation under this shrinking menu requires growing energy at least of order ell^-6; a baseline-normalized effect additionally needs a baseline-variance lower bound.",
    }


def calculate():
    for name, expected in SOURCES.items():
        actual = hashlib.sha256((RESEARCH/name).read_bytes()).hexdigest()
        assert actual == expected, (name, actual, expected)
    group1 = group_energy_modes_sources()
    group2, coefficient = group_bump_constants()
    group3 = group_dyson_finite_budget(coefficient)
    return strings({
        "round": 1258, "status": "PASS", "science_groups": 3,
        "arithmetic": "Exact Python Fraction and integer identities only; finite algebra supports the analytic all-mode and all-Dyson-order arguments in proof.md.",
        "source_sha256": SOURCES,
        "group_1_uniform_mode_energy_covariance_and_full_excitation_source": group1,
        "group_2_explicit_bump_and_constants": group2,
        "group_3_bounded_Dyson_and_nonempty_finite_budget": group3,
        "new_adopted_axioms": 0,
        "scope": {
            "all_nonnegative_diagonal_occupation_modes_retained": True,
            "finite_energy_alone_implies_Hadamard": False,
            "all_TR_states_covered": False,
            "reference_baseline_or_KMS_selected": False,
            "reference_renormalized_source_included_in_excitation_bound": False,
            "original_classical_geometry_automatically_self_consistent": False,
            "actual_autonomous_probe_preparation_or_record_certified": False,
            "microscopic_continuity_or_minimal_length_inferred": False,
            "ODE_reintegrated": False,
            "QFT_UV_truncated": False,
        },
        "net_result": "A common finite energy bound uniformly controls local covariance, bounded finite-coupling test probabilities, and the complete mean excitation source in the declared Gaussian family, with an explicit nonempty finite menu; no unique physical reference or complete gravity claim.",
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true", help="exclusive first creation, never overwrite")
    args = parser.parse_args()
    result = calculate()
    path = HERE/"results.json"
    if args.write_results:
        with path.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        print("PASS: 1258 three exact author groups; first results created exclusively.")
    else:
        saved = json.loads(path.read_text(encoding="utf-8"))
        assert saved == result, "saved results mismatch"
        print("PASS: 1258 three exact author groups; read-only reconstruction matches saved results.")


if __name__ == "__main__":
    main()
