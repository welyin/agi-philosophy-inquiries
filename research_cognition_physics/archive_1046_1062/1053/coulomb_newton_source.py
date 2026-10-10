"""1053: actual Coulomb source dictionary and exact finite inequalities.

No Coulomb propagation is numerically simulated. Smooth Slater source moments
are quadrature diagnostics; strict task bounds use rational inequalities.
The default is read-only. --write-results creates a new result exclusively.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
RESULT = HERE / "results.json"
HISTORY = [
    "archive_894_922/research_note_896.md",
    "archive_894_922/research_note_899.md",
    "archive_956_989/research_note_986.md",
    "archive_1046_/research_note_1047.md",
    "archive_1046_/research_note_1050.md",
    "archive_1046_/1050/proof.md",
    "archive_1046_/1050/coulomb_natural_instrument_results.json",
    "archive_1046_/1050/research_round_1050_checks.json",
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def norm(a):
    return float(np.linalg.norm(a, 2))


def rational_certificate():
    s = F(1, 10**6)
    # S decreases with R. Never evaluate exp(-1000) as a floating zero.
    e20_lower = sum(F(20)**j / math.factorial(j) for j in range(50))
    e20_target = (1 + F(20) + F(20)**2 / 3) / s
    assert e20_lower > e20_target
    assert F(2) / (1 - s) < F(10, 7)**2
    assert (1 + s) / (1 - s)**3 < 4  # Lowdin center commutator < 2 S.
    assert F(55, 6)**2 > 84
    assert F(26, 45)**2 > F(1, 3)
    d_lower = F(221, 400)
    assert (1 - d_lower)**2 > F(1, 5)
    R, a, T, lb = F(1000), F(500), F(1000, 137), F(22, 7000)
    delta0 = F(4) * F(10, 7) / (3 * a) + 2 * s
    delta_natural = delta0 + T * (
        4 * F(55, 6) * F(26, 45) / (3 * a) + 3 / a**2)
    delta = delta_natural + T * 4 * lb * F(26, 45) / (3 * a)
    delta_simple = F(27, 250)
    eps = F(74239, 335650)
    old_eps = F(22, 7) / R * (F(20, 7) + T * F(37, 2) / 2)
    assert old_eps == eps
    assert delta == F(4156812479, 38839500000)
    assert delta < delta_simple
    bad = eps**2
    joint = d_lower - delta_simple - bad
    recover = (d_lower - delta_simple - 2 * bad) / 2
    assert joint > F(791, 2000)
    assert recover > F(1733, 10000)
    # Stronger exact rational bound is reported separately, not optimized.
    exact_joint = d_lower - delta - bad
    exact_recover = (d_lower - delta - 2 * bad) / 2
    vals = {
        "R": R, "a": a, "T": T,
        "overlap_strict_upper": s,
        "initial_source_error_strict_upper": delta0,
        "natural_source_error_strict_upper": delta_natural,
        "driven_source_error_strict_upper": delta,
        "driven_source_error_simple_upper": delta_simple,
        "d_strict_lower": d_lower,
        "old_instrument_error_strict_upper": eps,
        "bad_pointer_probability_strict_upper": bad,
        "pointer_total_variation_strict_upper": bad,
        "joint_source_half_test_strict_lower": joint,
        "record_recovery_half_test_strict_lower": recover,
        "joint_source_half_test_tighter_lower": exact_joint,
        "record_recovery_half_test_tighter_lower": exact_recover,
    }
    return {
        "exact": {key: str(value) for key, value in vals.items()},
        "decimal": {key: float(value) for key, value in vals.items()},
        "positive_Taylor_terms": 50,
        "e20_positive_Taylor_exceeds_required_threshold": True,
        "R1000_overlap_underflow_used": False,
        "all_inequalities_use_positive_rational_arithmetic": True,
    }


def orbital_moments(R, order):
    """Actual normalized right Slater orbital; a=R/2, angular axis along x."""
    t, wt = np.polynomial.laguerre.laggauss(order)
    z, wz = np.polynomial.legendre.leggauss(order)
    r = t[:, None] / 2
    a = R / 2
    weight = wt[:, None] * wz[None, :] * t[:, None]**2 / 4
    f = 1 / np.sqrt(r*r + a*a) - 1 / np.sqrt(
        r*r + R*R + 2*R*r*z[None, :] + a*a)
    return {"normalization": float(np.sum(weight)),
            "f": float(np.sum(weight*f)),
            "f2": float(np.sum(weight*f*f))}


def cross_scaled(R, order):
    """Prolate integral divided by exp(-R); never loses positive overlap."""
    u, wu = np.polynomial.laguerre.laggauss(order)
    v, wv = np.polynomial.legendre.leggauss(order)
    mu, nu = 1 + u[:, None]/R, v[None, :]
    a = R/2
    rR, rL = R*(mu-nu)/2, R*(mu+nu)/2
    f = 1/np.sqrt(rR*rR+a*a) - 1/np.sqrt(rL*rL+a*a)
    weight = R*R/4 * wu[:, None]*wv[None, :] * (mu*mu-nu*nu)
    return {"overlap_scaled": float(np.sum(weight)),
            "f_scaled": float(np.sum(weight*f)),
            "f2_scaled": float(np.sum(weight*f*f))}


def fermion_lift(a):
    creators = []
    for j in range(4):
        c = np.zeros((16, 16))
        for bits in range(16):
            if not ((bits >> j) & 1):
                c[bits | (1 << j), bits] = (-1)**((bits & ((1 << j)-1)).bit_count())
        creators.append(c)
    sector = [b for b in range(16) if b.bit_count() == 2]
    full = sum(a[i, j]*creators[i]@creators[j].T for i in range(4) for j in range(4))
    return full[np.ix_(sector, sector)], sector


def matrix_diagnostic(R=20., order=80):
    """A diagnostic of the analytic dictionary, not a new instrument parameter."""
    m, c = orbital_moments(R, order), cross_scaled(R, order)
    e = math.exp(-R)
    S = e*(1+R+R*R/3)
    G = np.array([[1., S], [S, 1.]])
    lam, vec = np.linalg.eigh(G)
    L = (vec*lam**(-.5))@vec.T
    raw_f = np.diag([-m["f"], m["f"]])
    raw_f2 = np.array([[m["f2"], e*c["f2_scaled"]],
                       [e*c["f2_scaled"], m["f2"]]])
    f1, f2 = L@raw_f@L, L@raw_f2@L
    q, sector = fermion_lift(np.diag([-.5, -.5, .5, .5]))
    fsum, _ = fermion_lift(np.kron(f1, np.eye(2)))
    tail, _ = fermion_lift(np.kron(f2-f1@f1, np.eye(2)))
    a, d = R/2, 1-1/math.sqrt(5)
    compressed = a/2*fsum
    second = (a/2)**2*(fsum@fsum+tail)
    gram = second-d*(compressed@q+q@compressed)+d*d*(q@q)
    eig = np.linalg.eigvalsh(gram)
    source = a*m["f"]/math.sqrt(1-S*S)
    err = norm(compressed-source*q)
    assert min(eig) > -2e-12
    assert err < 1e-12
    assert abs(c["f_scaled"]) < 1e-10
    assert abs(c["overlap_scaled"]-(1+R+R*R/3)) < 2e-10
    bound = 4*math.sqrt(2/(1-S))/(3*a)+d*S*math.sqrt(1+S)/(1-S)**1.5
    actual = math.sqrt(max(0., float(eig[-1])))
    assert actual <= bound
    return {
        "diagnostic_only_R": R,
        "quadrature_order": order,
        "overlap_floating_positive": S,
        "sector_binary_labels": sector,
        "Q_diagonal": np.diag(q).tolist(),
        "compressed_A_diagonal": np.diag(compressed).tolist(),
        "actual_compressed_source_coefficient": source,
        "compressed_A_equals_source_Q_residual": err,
        "full_source_error_gram_eigenvalues": eig.tolist(),
        "full_source_error_norm_diagnostic": actual,
        "analytic_initial_source_error_bound": bound,
        "cross_f_odd_symmetry_scaled_residual": abs(c["f_scaled"]),
        "off_code_variance_included": True,
    }


def actual_R1000_diagnostic():
    lower_order, higher_order = orbital_moments(1000., 48), orbital_moments(1000., 80)
    # S is strictly positive but tiny. Enclose the orthogonalization correction;
    # do not evaluate exp(-1000), set S=0, or call this quadrature certified.
    s_upper = 1e-6
    uncorrected = 500*higher_order["f"]
    correction = uncorrected*(1/math.sqrt(1-s_upper*s_upper)-1)
    convergence = max(abs(lower_order[k]-higher_order[k]) for k in higher_order)
    assert convergence < 2e-12
    assert abs(higher_order["normalization"]-1) < 3e-12
    return {
        "R": 1000, "a": 500,
        "quadrature_orders": [48, 80],
        "raw_right_orbital_f": higher_order["f"],
        "raw_right_orbital_f2": higher_order["f2"],
        "radial_angular_normalization": higher_order["normalization"],
        "coefficient_before_Lowdin": uncorrected,
        "positive_Lowdin_correction_upper_diagnostic": correction,
        "coefficient_formula": "500 * raw_f / sqrt(1-S(1000)^2), 0<S(1000)<1e-6",
        "quadrature_order_change_max": convergence,
        "quadrature_is_numerical_diagnostic_not_rigorous_error_bound": True,
        "actual_Coulomb_propagation_numerically_solved": False,
    }


def kernel_checks():
    # Algebraic derivative formula, with values only as consistency checks.
    a = 500.
    r = np.array([0., a/math.sqrt(2), a, 2*a])
    k = 1/np.sqrt(r*r+a*a)
    laplacian = -3*a*a/(r*r+a*a)**2.5
    density = 3*a*a/(4*math.pi*(r*r+a*a)**2.5)
    assert np.max(np.abs(laplacian+4*math.pi*density)) < 1e-20
    lip = r/(r*r+a*a)**1.5
    assert abs(lip[1]-2/(3*math.sqrt(3)*a*a)) < 1e-20
    # Exact radial CDF of w_a at infinity is 1: r^3/(r^2+a^2)^(3/2).
    return {"kernel": "(|r|^2+a^2)^(-1/2)",
            "positive_density": "3a^2/(4pi(|r|^2+a^2)^(5/2))",
            "radial_CDF": "r^3/(r^2+a^2)^(3/2)",
            "normalization_at_infinity": 1,
            "poisson_point_check_max_error": float(np.max(np.abs(laplacian+4*math.pi*density))),
            "fixed_nuclei_symmetric_source_difference": 0,
            "kernel_values": k.tolist()}


def run():
    return {
        "round": 1053,
        "new_conditional_physical_source_calibration_groups": 1,
        "new_cognitive_axioms": 0,
        "rational_certificate": rational_certificate(),
        "kernel": kernel_checks(),
        "actual_R1000_orbital_integral": actual_R1000_diagnostic(),
        "finite_CAR_source_diagnostic": matrix_diagnostic(),
        "scope": {
            "same_actual_1050_H_and_pulse_and_preparation": True,
            "all_six_dimensional_inputs_and_passive_references": True,
            "source": "electronic rest-mass leading-G Newton constraint dual coefficient",
            "full_actual_U_source_dictionary": "|| A U_F V - U_F V d Q || <= delta",
            "positive_overlap_not_replaced_by_underflow_zero": True,
            "joint_unormalized_pointer_and_source": True,
            "record_decoder_obstruction": "one use of this binary pointer, same CPTP decoder; not arbitrary histories or repeated sampling",
            "Newton_constraint_adopted_not_generated": True,
            "total_apparatus_source_matched": False,
            "physical_finite_G_backreaction_remainder_certified": False,
            "gravity_detector_constructed": False,
            "whole_M4_or_roadmap_complete": False,
        },
        "historical_sha256": {rel: digest(ROOT/rel) for rel in HISTORY},
    }


def compare(actual, expected, path="root"):
    if isinstance(actual, dict):
        assert actual.keys() == expected.keys(), path
        for key in actual:
            compare(actual[key], expected[key], path+"."+key)
    elif isinstance(actual, list):
        assert len(actual) == len(expected), path
        for i, value in enumerate(actual):
            compare(value, expected[i], f"{path}[{i}]")
    elif isinstance(actual, float):
        assert math.isclose(actual, expected, rel_tol=3e-11, abs_tol=2e-12), (path, actual, expected)
    else:
        assert actual == expected, (path, actual, expected)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    output = run()
    if args.write_results:
        with RESULT.open("x", encoding="utf-8") as f:
            json.dump(output, f, ensure_ascii=False, indent=2)
            f.write("\n")
    else:
        compare(output, json.loads(RESULT.read_text(encoding="utf-8")))
    print(json.dumps({"round": 1053, "all_passed": True,
                      "mode": "exclusive_write" if args.write_results else "read_only",
                      "source_delta_upper": output["rational_certificate"]["decimal"]["driven_source_error_strict_upper"]}))


if __name__ == "__main__":
    main()
