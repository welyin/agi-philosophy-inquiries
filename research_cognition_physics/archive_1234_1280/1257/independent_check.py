"""Independent finite-shell audit for round 1257.

The integer shell is enumerated from a cube; no author module is imported.
The default compares an existing saved result without modifying any file.
Only --write performs exclusive first creation of independent_results.json.
Analytic field regularity, local smear existence and actual permissions are
outside this finite exact-arithmetic audit and belong to the written proof.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def vector_sum(weights, degree):
    return {
        index: sum(
            weight * math.prod(n[i] for i in index)
            for n, weight in weights.items()
        )
        for index in itertools.product(range(3), repeat=degree)
    }


def reduced_tensor(tensor):
    return {",".join(map(str, index)): str(value)
            for index, value in tensor.items() if value != 0}


def calculate():
    shell = tuple(
        n for n in itertools.product(range(-3, 4), repeat=3)
        if sum(x*x for x in n) == 9
    )
    axes = tuple(n for n in shell if sum(x != 0 for x in n) == 1)
    remaining = tuple(n for n in shell if n not in axes)
    assert len(shell) == 30 and len(axes) == 6 and len(remaining) == 24
    assert all(sorted(abs(x) for x in n) == [1, 2, 2] for n in remaining)

    # Occupations are normalized by one common positive parameter nu.
    a = {n: Q(4 if n in axes else 0) for n in shell}
    b = {n: Q(0 if n in axes else 1) for n in shell}
    assert all(a[n] == a[tuple(-x for x in n)] for n in shell)
    assert all(b[n] == b[tuple(-x for x in n)] for n in shell)
    assert all(value >= 0 for value in (*a.values(), *b.values()))
    ma = {d: vector_sum(a, d) for d in range(3)}
    mb = {d: vector_sum(b, d) for d in range(3)}
    assert ma == mb
    assert ma[0][()] == 24
    assert all(x == 0 for x in ma[1].values())
    assert all(ma[2][i, j] == (72 if i == j else 0)
               for i in range(3) for j in range(3))

    # T_00, all T_0i, and all T_ij are built from the common scalar
    # mode/time coefficients times these zeroth/first/second moments.
    # Check the entire ordered second tensor, not just pressure/trace.
    stress_basis_difference = {
        "occupation_count": ma[0][()] - mb[0][()],
        "energy_spatial_trace": sum(ma[2][i, i]-mb[2][i, i]
                                   for i in range(3)),
        **{"flux_"+str(i): ma[1][(i,)]-mb[1][(i,)] for i in range(3)},
        **{"gradient_"+str(i)+str(j): ma[2][i, j]-mb[2][i, j]
           for i in range(3) for j in range(3)},
    }
    assert all(x == 0 for x in stress_basis_difference.values())
    fourth_a, fourth_b = vector_sum(a, 4), vector_sum(b, 4)
    fourth_delta = {idx: fourth_a[idx]-fourth_b[idx] for idx in fourth_a}
    assert fourth_delta[0, 0, 0, 0] == 384
    assert fourth_delta[0, 0, 1, 1] == -192

    even_pairs = {min(n, tuple(-x for x in n)) for n in shell}
    assert len(even_pairs) == 15
    assert sum(n in axes for n in even_pairs) == 3

    # The local spatial operator is (d_x^2-d_y^2) applied to a radial bump.
    # On this fixed shell the common Fourier bump value and time overlap
    # factor out, leaving the following integer angular weight.
    angular_weight = {n: Q((n[0]*n[0]-n[1]*n[1])**2) for n in shell}
    wa = sum(a[n]*angular_weight[n] for n in shell)
    wb = sum(b[n]*angular_weight[n] for n in shell)
    assert (wa, wb) == (1296, 144) and wa/wb == 9
    assert wa - wb == 1152

    # Let the unexcited reference variance be V0>0, choose nu so B adds V0.
    # The exact centered Gaussian variances are then 10 V0 and 2 V0.
    # With lambda^2 V0=1/2, cos expectations are exp(-5/2), exp(-1/2).
    variance_a, variance_b = Q(1)+wa/wb, Q(2)
    lambda_squared_v0 = Q(1, 2)
    exponent_a = lambda_squared_v0*variance_a/2
    exponent_b = lambda_squared_v0*variance_b/2
    assert (exponent_a, exponent_b) == (Q(5, 2), Q(1, 2))

    # Strict rational bounds use different Taylor arguments on the two sides.
    # Alternating odd partial sum is a strict lower bound at x=1/2.
    exp_minus_half_lower = sum((-Q(1, 2))**j / math.factorial(j)
                              for j in range(4))
    assert exp_minus_half_lower == Q(29, 48)
    # Positive partial sum is below exp(5/2), and already exceeds 12.
    exp_five_halves_lower = sum(Q(5, 2)**j / math.factorial(j)
                               for j in range(7))
    assert exp_five_halves_lower > 12
    exp_minus_five_halves_upper = Q(1, 12)
    gap_lower = (exp_minus_half_lower-exp_minus_five_halves_upper)/2
    assert gap_lower == Q(25, 96) > Q(1, 4)
    minimax_lower = gap_lower/2
    assert minimax_lower == Q(25, 192) > Q(1, 8)

    # A same-point comoving or angularly equal menu has weight one:
    # it cannot detect this occupation difference. The actual anisotropic
    # local smear above is essential, not supplied by equal stress alone.
    constant_weight_difference = sum(a.values())-sum(b.values())
    assert constant_weight_difference == 0

    # A radius 1/10 bump on L=10 has |k|r=3*pi/50 < 1/5.
    # This proves the real Fourier bump value is positive for nonnegative
    # radial bump data; it does not compute the temporal mode function.
    kr_upper = Q(3, 50)*Q(22, 7)
    cos_guard = 1-kr_upper*kr_upper/2
    assert kr_upper < Q(1, 5) and cos_guard > Q(49, 50)

    return {
        "status": "PASS",
        "independent_groups": 2,
        "independent_check_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "group_complete_source_moments": {
            "enumeration_cube": [-3, 3],
            "shell_norm_squared": 9,
            "shell_count": len(shell),
            "axis_count": len(axes),
            "remaining_count": len(remaining),
            "opposite_pairs": len(even_pairs),
            "axis_opposite_pairs": 3,
            "remaining_opposite_pairs": 12,
            "occupations_per_nu": {"axis_in_A": "4", "remaining_in_B": "1"},
            "common_occupation_count_per_nu": str(ma[0][()]),
            "common_first_moments": {str(i): str(ma[1][(i,)]) for i in range(3)},
            "common_complete_second_tensor": {
                str(i)+str(j): str(ma[2][i, j])
                for i in range(3) for j in range(3)
            },
            "stress_basis_difference": {k: str(v) for k, v in stress_basis_difference.items()},
            "common_nonzero_spatial_trace_per_nu": "216",
            "nonzero_fourth_difference": reduced_tensor(fourth_delta),
            "scope": "All mean quadratic stress components coincide on a common isotropic background with identical shell mode functions; the common quantum source is not zero or erased.",
        },
        "group_local_bounded_gaussian_effect": {
            "spatial_operator": "d_x^2-d_y^2 on one common radial compact bump",
            "angular_weights_per_nu": {"A": str(wa), "B": str(wb)},
            "ratio_A_to_B_added_variance": str(wa/wb),
            "total_variance_in_V0": {"A": str(variance_a), "B": str(variance_b)},
            "lambda_squared_V0": str(lambda_squared_v0),
            "negative_exponent_magnitudes": {"A": str(exponent_a), "B": str(exponent_b)},
            "exp_minus_half_strict_lower": str(exp_minus_half_lower),
            "exp_five_halves_partial_sum_strict_lower": str(exp_five_halves_lower),
            "exp_minus_five_halves_strict_upper": str(exp_minus_five_halves_upper),
            "probability_gap_strict_lower": str(gap_lower),
            "gap_margin_above_quarter": str(gap_lower-Q(1, 4)),
            "mean_source_only_minimax_error_strict_lower": str(minimax_lower),
            "same_point_comoving_difference": str(constant_weight_difference),
            "local_bump_radius": "1/10",
            "kr_upper_from_pi_lt_22_over_7": str(kr_upper),
            "bump_Fourier_positive_cos_guard": str(cos_guard),
            "scope": "Exact centered Gaussian bounded-effect probabilities, conditional on the stated common finite spacetime smear and preparation/effect permissions; not a finite-coupling detector construction.",
        },
        "not_certified": [
            "Hadamard existence and propagation by finite arithmetic",
            "temporal mode integration or its nonzero-overlap constants",
            "automatic internal preparation or actual implementation of the terminal effect",
            "two states satisfying an unproved KMS or TR reference role",
            "unchanged stress fluctuations or complete quantum geometry",
            "original classical FLRW solving the equation after silently adding the common quantum source",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true",
                        help="create the independent result once, refusing overwrite")
    args = parser.parse_args()
    actual = calculate()
    path = HERE / "independent_results.json"
    if args.write:
        with path.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(actual, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        action = "FIRST_WRITE"
    else:
        expected = json.loads(path.read_text(encoding="utf-8"))
        assert expected == actual, "independent saved result mismatch"
        action = "READ_ONLY_MATCH"
    print(json.dumps({"status": "PASS", "action": action,
                      "independent_groups": actual["independent_groups"],
                      "gap_strict_lower": "25/96",
                      "shell_counts": [6, 24]}))


if __name__ == "__main__":
    main()
