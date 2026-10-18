"""Independent audit of round 1291; no author implementation is imported.

Exact checks: mass-shell Ward vertices, Wick weights, Beta constants, fixed
compact-smearing budget, and the frozen round-1290 vacuum bound.
Finite oscillator matrices are diagnostics only and explicitly retain a cutoff.
Default execution compares without writing. --write creates results exclusively.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents
            if (p / "research_cognition_physics").is_dir()
            and (p / "scripts").is_dir())
STAGE = HERE.parent


def rat(x):
    x = Q(x)
    return f"{x.numerator}/{x.denominator}"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def shell_vertices(m, p, ep, q, eq):
    """Build the tensor by upper/lower metric conversion, using exact rationals."""
    sign = (-1, 1, 1, 1)
    upper_p = tuple(map(Q, (ep, *p)))
    upper_q = tuple(map(Q, (eq, *q)))
    lower_p = [sign[j] * upper_p[j] for j in range(4)]
    lower_q = [sign[j] * upper_q[j] for j in range(4)]
    assert sum(upper_p[j] * lower_p[j] for j in range(4)) == -Q(m) ** 2
    assert sum(upper_q[j] * lower_q[j] for j in range(4)) == -Q(m) ** 2
    inner = sum(upper_p[j] * lower_q[j] for j in range(4))
    pair, scatter = [], []
    for i in range(4):
        tp, ts = [], []
        for j in range(4):
            sym = (lower_p[i] * lower_q[j] + lower_q[i] * lower_p[j]) / 2
            metric = sign[i] if i == j else 0
            tp.append(-sym + metric * (inner - Q(m) ** 2) / 2)
            ts.append(sym - metric * (inner + Q(m) ** 2) / 2)
        pair.append(tp)
        scatter.append(ts)
    pair_ward = [sum((upper_p[i] + upper_q[i]) * pair[i][j]
                     for i in range(4)) for j in range(4)]
    scatter_ward = [sum((upper_p[i] - upper_q[i]) * scatter[i][j]
                        for i in range(4)) for j in range(4)]
    assert pair_ward == [0] * 4
    assert scatter_ward == [0] * 4
    assert max(abs(x) for row in pair for x in row) <= 3 * Q(ep) * Q(eq)
    assert max(abs(x) for row in scatter for x in row) <= 3 * Q(ep) * Q(eq)
    return pair, scatter, pair_ward, scatter_ward


def exact_algebra():
    cases = [
        (3, (0, 4, 0), 5, (0, 0, 4), 5),
        (12, (0, 5, 0), 13, (0, 0, 9), 15),
        (1, (0, Q(3, 4), 0), Q(5, 4),
         (Q(4, 3), 0, 0), Q(5, 3)),
        (5, (12, 0, 0), 13, (-12, 0, 0), 13),
    ]
    vertices = []
    for m, p, ep, q, eq in cases:
        pair, sc, wp, ws = shell_vertices(m, p, ep, q, eq)
        vertices.append({
            "mass": rat(m),
            "p": [rat(x) for x in p], "q": [rat(x) for x in q],
            "delta_energy": rat(Q(ep) - Q(eq)),
            "s00": rat(sc[0][0]),
            "pair_Ward": [rat(x) for x in wp],
            "scattering_Ward": [rat(x) for x in ws],
        })
    assert vertices[0]["delta_energy"] == "0/1"
    assert Q(vertices[0]["s00"]) > 0

    # Gamma half-integers carry a common sqrt(pi). This computes the two
    # radial integrals independently through their Euler-Beta values.
    def gamma_half_coeff(n):
        return Q(math.factorial(2 * n), 4 ** n * math.factorial(n))
    radial0_pi_coeff = gamma_half_coeff(1) * gamma_half_coeff(2) / (2 * 6)
    radial1 = Q(math.factorial(1) * math.factorial(1), 2 * math.factorial(3))
    assert radial0_pi_coeff == Q(1, 32)
    assert radial1 == Q(1, 12)
    J0_inverse_pi_coeff = radial0_pi_coeff / 2
    J1_inverse_pi2_coeff = radial1 / 2
    assert J0_inverse_pi_coeff == Q(1, 64)
    assert J1_inverse_pi2_coeff == Q(1, 24)

    # P_j is built recursively by integration by parts, not copied as P3/P4.
    polynomials = [[Q(1)]]
    for j in range(1, 5):
        prev = polynomials[-1]
        now = [j * x for x in prev] + [Q(1)]
        polynomials.append(now)
    assert polynomials[3] == list(map(Q, [6, 6, 3, 1]))
    assert polynomials[4] == list(map(Q, [24, 24, 12, 4, 1]))
    for j, coeffs in enumerate(polynomials):
        derivative = [(k + 1) * coeffs[k + 1] for k in range(j)]
        # beta=1: dP_j/dM - P_j = -M^j.
        difference = [(derivative[k] if k < j else 0) - coeffs[k]
                      for k in range(j + 1)]
        assert difference == [0] * j + [-1]

    wick_rows = []
    for n in (Q(1, 7), Q(2, 3), Q(3), Q(11, 2)):
        # Thermal geometric factorial moments.
        N = n
        N2 = 2 * n * n + n
        adag2a2 = N2 - N
        a2adag2 = N2 + 3 * N + 2
        assert adag2a2 == 2 * n * n
        assert a2adag2 == 2 * (n + 1) ** 2
        assert N2 - N * N == n * (n + 1)
        wick_rows.append({
            "n": rat(n), "creation_annihilation_sum": rat(a2adag2 + adag2a2),
            "number_variance": rat(N2 - N * N),
        })
    return {
        "exact_vertex_cases": vertices,
        "zero_energy_transfer_nonzero_vertex": True,
        "J0_coefficient_of_inverse_pi_ell3": rat(J0_inverse_pi_coeff),
        "J1_coefficient_of_inverse_pi2_ell4": rat(J1_inverse_pi2_coeff),
        "beta1_P3_coefficients_ascending": [rat(x) for x in polynomials[3]],
        "beta1_P4_coefficients_ascending": [rat(x) for x in polynomials[4]],
        "exact_single_mode_Wick_moments": wick_rows,
    }


def compact_budget():
    # Full box width 1/8: ||u'||TV=16. Time has six boxes, each spatial
    # direction has four; the common C-infinity mollifier has radius 1/8.
    width = Q(1, 8)
    one_derivative_TV = 2 / width
    derivative_bounds = {}
    for r in (0, 2, 4, 6):
        # Exact total variation of the r-fold finite-difference measure.
        weights = [Q((-1) ** j * math.comb(r, j)) / width ** r
                   for j in range(r + 1)]
        bound = sum(abs(x) for x in weights)
        assert bound == one_derivative_TV ** r
        derivative_bounds[r] = bound
    B6 = derivative_bounds[6]
    C1 = 1 + 2 * 3 * derivative_bounds[2]
    C1 += 3 * derivative_bounds[4] + 6 * derivative_bounds[2] ** 2
    assert C1 == 591361 == 769 ** 2
    assert 6 * width / 2 + Q(1, 8) == Q(1, 2)
    assert 4 * width / 2 + Q(1, 8) == Q(3, 8)

    M = 2 ** 20
    # e>2 and M>256: exp(-M)<2^-256. Thus d>1/2 and n_M<2*x.
    x = Q(1, 2 ** 256)
    P3 = sum(Q(math.comb(3, j) * math.factorial(j) * M ** (3 - j))
             for j in range(4))
    P4 = sum(Q(math.comb(4, j) * math.factorial(j) * M ** (4 - j))
             for j in range(5))
    vacuum = Q(9, 1120 * 3 ** 4 * 4) * B6 ** 2 / (2 * M) ** 4
    old = json.loads((STAGE / "1290/results.json").read_text(encoding="utf-8"))
    frozen_vacuum = Q(old["exact"]["compact_example"]["variance_upper_using_pi_gt_3"])
    assert vacuum == frozen_vacuum == Q(1, 40320 * 2 ** 36)

    # All are outward upper bounds using pi>3 and d>1/2.
    mean = x * P3 / 9
    scattering = 2 * C1 ** 2 * x * (P4 / 192 + P3 / 216)
    pair_multiplier = 1 + 4 * x + 8 * x * x
    full_variance = pair_multiplier * vacuum + scattering
    # Match the author's simpler published <2 vacuum guard as well.
    coarse_full_variance = 2 * vacuum + scattering
    assert pair_multiplier < 2
    assert mean < Q(1, 10 ** 59)
    assert scattering < Q(1, 10 ** 43)
    assert full_variance < coarse_full_variance < Q(1, 10 ** 15)
    assert full_variance + mean ** 2 < Q(1, 10 ** 15)

    return {
        "test": "f00=g(t) v(x1)v(x2)v(x3); all other components zero",
        "time_boxes": 6, "boxes_each_spatial_axis": 4,
        "box_full_width": rat(width), "mollifier_radius": "1/8",
        "time_support": ["-1/2", "1/2"],
        "each_spatial_support": ["-3/8", "3/8"],
        "beta": "1/1", "ell": "1/1", "mass": M,
        "B0": "1/1", "B6_upper": rat(B6), "C1_upper": rat(C1),
        "exp_minus_mass_upper": rat(x),
        "beta1_P3": rat(P3), "beta1_P4": rat(P4),
        "vacuum_from_frozen_1290_recomputed": rat(vacuum),
        "thermal_mean_upper": rat(mean),
        "scattering_variance_upper": rat(scattering),
        "pair_multiplier_upper": rat(pair_multiplier),
        "full_variance_upper": rat(full_variance),
        "coarse_full_variance_upper": rat(coarse_full_variance),
        "uncentered_second_moment_upper": rat(full_variance + mean ** 2),
        "targets": {"absolute_mean": "1/10^59", "scattering": "1/10^43",
                    "variance_and_second_moment": "1/10^15"},
        "strict_coarse_variance_margin": rat(Q(1, 10 ** 15) - coarse_full_variance),
        "fixed_UV_mass_domain_still_required": True,
    }


def oscillator_diagnostic():
    rows = []
    for q, alpha, b in ((0.5, 1 / 3 + 2j / 5, -2 / 7),
                       (1 / 3, -1 / 5 + 1j / 2, 3 / 4)):
        n_infinite = q / (1 - q)
        expected_var = (2 * abs(alpha) ** 2
                        * ((n_infinite + 1) ** 2 + n_infinite ** 2)
                        + b ** 2 * n_infinite * (1 + n_infinite))
        for dim in (32, 64):
            a = np.zeros((dim, dim), dtype=complex)
            for j in range(1, dim):
                a[j - 1, j] = math.sqrt(j)
            adag = a.conj().T
            number = adag @ a
            observable = alpha * (a @ a) + alpha.conjugate() * (adag @ adag) + b * number
            weights = (1 - q) * q ** np.arange(dim) / (1 - q ** dim)
            mean = float(np.sum(weights * np.diag(observable)).real)
            second = float(np.sum(weights * np.diag(observable @ observable)).real)
            var = second - mean * mean
            ccr_expect = float(np.sum(weights * np.diag(a @ adag - adag @ a)).real)
            assert np.max(np.abs(observable - observable.conj().T)) < 1e-13
            assert var > 0
            if dim == 64:
                assert abs(var - expected_var) < 1e-11
            rows.append({
                "geometric_q": q, "dimension": dim,
                "alpha_real": alpha.real, "alpha_imag": alpha.imag, "b": b,
                "finite_normalized_mean": mean, "finite_variance": var,
                "infinite_Wick_variance": expected_var,
                "absolute_variance_difference": abs(var - expected_var),
                "infinite_geometric_tail_probability": q ** dim,
                "thermal_CCR_expectation": ccr_expect,
            })
    return {
        "scope": "Finite oscillator matrix diagnostic only; not a field truncation or proof of an operator-norm tail.",
        "boundary_CCR_not_exact": True,
        "rows": rows,
    }


def calculate():
    files = [
        ("round1291_admission", STAGE / "_shared/round1291_admission.md"),
        ("author_proof", HERE / "proof.md"),
        ("author_check", HERE / "check.py"),
        ("author_results", HERE / "results.json"),
        ("independent_check", Path(__file__).resolve()),
        ("round1290_proof", STAGE / "1290/proof.md"),
        ("round1290_results", STAGE / "1290/results.json"),
    ]
    return {
        "schema": "round1291-independent-v1",
        "status": "PASS_WITH_SCOPE",
        "independent_check_groups": 1,
        "source_bindings": {label: {"path": str(path.relative_to(ROOT)).replace("\\", "/"),
                                   "sha256": sha(path)} for label, path in files},
        "exact_algebra": exact_algebra(),
        "compact_budget": compact_budget(),
        "oscillator_diagnostic": oscillator_diagnostic(),
        "scope": {
            "fixed_beta_clock_spacetime_tests": True,
            "complete_minimal_KG_stress_vacuum_Wick": True,
            "thermal_mean_retained": True,
            "gapless_number_preserving_channel_retained": True,
            "same_beta_product_reference_only_at_zero_portal_fixed_background": True,
            "infinite_volume_reference_is_only_locally_normal": True,
            "finite_total_thermal_reservoir_preparation_claimed": False,
            "nonzero_portal_or_curved_geometry_certificate": False,
            "actual_instrument_constructed": False,
        },
    }


def match(now, saved, key="$"):
    if type(now) is not type(saved):
        raise AssertionError(f"type mismatch: {key}")
    if isinstance(now, dict):
        if now.keys() != saved.keys():
            raise AssertionError(f"key mismatch: {key}")
        for k in now:
            match(now[k], saved[k], f"{key}.{k}")
    elif isinstance(now, list):
        if len(now) != len(saved):
            raise AssertionError(f"length mismatch: {key}")
        for j, (x, y) in enumerate(zip(now, saved)):
            match(x, y, f"{key}[{j}]")
    elif isinstance(now, float):
        if not math.isclose(now, saved, rel_tol=1e-10, abs_tol=1e-12):
            raise AssertionError(f"float mismatch: {key}")
    elif now != saved:
        raise AssertionError(f"value mismatch: {key}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = calculate()
    dest = HERE / "independent_results.json"
    if args.write:
        with dest.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, sort_keys=True,
                      allow_nan=False)
            stream.write("\n")
        receipt = "CREATED_EXCLUSIVELY"
    else:
        match(result, json.loads(dest.read_text(encoding="utf-8")))
        receipt = "READ_ONLY_MATCH"
    print(json.dumps({"status": result["status"], "receipt": receipt,
                      "result_sha256": sha(dest)}, ensure_ascii=False))


if __name__ == "__main__":
    main()

