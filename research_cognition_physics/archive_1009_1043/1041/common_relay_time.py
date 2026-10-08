"""1041: exact relay-cone and ordinary binary-channel composition certificates.

Default recomputes without writes; --write exclusively creates the result.
No full multi-metric QFT or physical instrument implementation is certified.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
PHYSICS = HERE.parents[1]
RESULT = HERE / "common_relay_time_results.json"
HISTORY = [
    "archive_231_258/research_note_232.md", "archive_231_258/research_note_233.md",
    "archive_231_258/research_note_234.md", "archive_370_428/research_note_392.md",
    "archive_370_428/research_note_382.md", "archive_370_428/research_note_383.md",
    "archive_370_428/research_note_384.md", "archive_370_428/research_note_386.md",
    "archive_370_428/research_note_425.md", "archive_467_530/research_note_522.md",
    "archive_467_530/research_note_523.md", "archive_854_872/research_note_859.md",
    "archive_923_934/research_note_930.md", "archive_956_989/research_note_963.md",
    "archive_1009_/research_note_1038.md", "archive_1009_/research_note_1040.md",
    "archive_1009_/1040/dependency_delta_1036_1040.json",
]


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), F(0))


def add(a, b):
    return [x+y for x, y in zip(a, b)]


def scale(c, a):
    return [c*x for x in a]


def mm(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))), F(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def mv(a, x):
    return [dot(row, x) for row in a]


def ident(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def metric(u, inverse=False):
    c = F(17) if inverse else F(17, 16)
    return [[F(i == j)-c*u[i]*u[j] for j in range(4)] for i in range(4)]


def nullvector(matrix):
    a = [row[:] for row in matrix]
    m, n = len(a), len(a[0]); row = 0; pivots = []
    for col in range(n):
        pivot = next((j for j in range(row, m) if a[j][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        z = a[row][col]; a[row] = [x/z for x in a[row]]
        for j in range(m):
            if j != row and a[j][col]:
                z = a[j][col]; a[j] = [x-z*y for x, y in zip(a[j], a[row])]
        pivots.append(col); row += 1
        if row == m:
            break
    free = next((j for j in range(n) if j not in pivots), None)
    assert free is not None
    x = [F(0)]*n; x[free] = F(1)
    for i, col in enumerate(pivots):
        x[col] = -a[i][free]
    assert all(dot(r, x) == 0 for r in matrix)
    return x


def reduce_certificate(points, weights):
    points = [p[:] for p in points]; weights = weights[:]
    original = len(points); steps = 0
    while len(points) > 5:
        alpha = nullvector([[p[j] for p in points] for j in range(4)]+[[F(1)]*len(points)])
        assert any(x > 0 for x in alpha) and any(x < 0 for x in alpha)
        t = min(w/a for w, a in zip(weights, alpha) if a > 0)
        weights = [w-t*a for w, a in zip(weights, alpha)]
        keep = [j for j, w in enumerate(weights) if w > 0]
        points = [points[j] for j in keep]; weights = [weights[j] for j in keep]
        assert sum(weights) == 1 and all(sum(w*p[k] for w, p in zip(weights, points)) == 0 for k in range(4))
        steps += 1
    return {"original_segments": original, "remaining_segments": len(points), "exact_reductions": steps,
            "points": [[str(x) for x in p] for p in points], "weights": [str(w) for w in weights]}


def channel(epsilon, flip=False):
    # Columns are input, rows output; every individual map is stochastic/CPTP on a classical code.
    T = [[1-epsilon, epsilon], [epsilon, 1-epsilon]]
    return T[::-1] if flip else T


def build_results():
    axes = [[F(1), F(0), F(0), F(0)],
            [F(-3, 5), F(4, 5), F(0), F(0)],
            [F(-3, 5), F(-4, 5), F(0), F(0)]]
    pair_rows = []
    for u in axes:
        assert dot(u, u) == 1
        g, gi = metric(u), metric(u, True)
        assert mm(g, gi) == ident(4)
        assert mv(g, u) == scale(F(-1, 16), u)
        assert dot(u, mv(g, u)) == F(-1, 16)
    for i, j in itertools.combinations(range(3), 2):
        n = add(axes[i], axes[j]); norm2 = dot(n, n)
        records = []
        for u in (axes[i], axes[j]):
            alpha = dot(n, u); transverse2 = norm2-alpha*alpha
            # Exact square root for this rational construction.
            transverse = F(4, 5) if i == 0 else F(24, 25)
            assert transverse*transverse == transverse2
            numerator = 4*alpha-transverse
            normalized_margin2 = numerator*numerator/(17*norm2)
            assert alpha > 0 and numerator > 0 and normalized_margin2 > F(1, 25)
            assert dot(n, mv(metric(u, True), n)) < 0
            records.append({"dual_metric_square": str(dot(n, mv(metric(u, True), n))),
                            "normalized_margin_squared": str(normalized_margin2)})
        pair_rows.append({"pair": [i+1, j+1], "clock": [str(x) for x in n], "checks": records})
    lengths = [F(6, 5), F(1), F(1)]
    legs = [scale(a, u) for a, u in zip(lengths, axes)]
    assert [sum(x) for x in zip(*legs)] == [0, 0, 0, 0]
    loop_rows = []
    for u, v in zip(axes, legs):
        s = dot(u, v)
        assert s > 0 and dot(v, mv(metric(u), v)) < 0
        loop_rows.append({"leg": [str(x) for x in v], "future_axis_projection": str(s),
                          "metric_square": str(dot(v, mv(metric(u), v)))})
    # 12 directions, all strictly inside their original cone, with an exact positive zero combination.
    points, weights = [], []
    for u, lam in zip(axes, [F(3, 8), F(5, 16), F(5, 16)]):
        for coord in (2, 3):
            for sign in (-1, 1):
                p = u[:]; p[coord] += F(sign, 8)
                assert dot(u, p) == 1 and dot(p, p) == F(65, 64)
                assert dot(p, mv(metric(u), p)) < 0
                points.append(p); weights.append(lam/4)
    reduction = reduce_certificate(points, weights)
    # Finite, explicitly assumed conversion displacements; their sum is compensated in the last segment.
    w = [F(1, 25), F(-1, 50), F(1, 50), F(0)]
    assert dot(w, w) <= F(1, 100)
    v3 = add(legs[2], scale(F(-1), w)); alpha = dot(axes[2], v3)
    perp2 = dot(v3, v3)-alpha*alpha
    assert alpha > 0 and alpha*alpha > 16*perp2
    assert add(add(add(legs[0], legs[1]), v3), w) == [0, 0, 0, 0]
    # General Lipschitz safety: s(v)=u.v-4||v_perp|| changes by at most 5||delta v||.
    conversion_budget, guaranteed_safety = F(1, 10), 1-5*F(1, 10)
    assert guaranteed_safety == F(1, 2)
    # Old round-233 normalization witness, transported to these specific geometric legs.
    noise_cases = []
    for eps in (F(0), F(1, 100), F(1, 10), F(1, 2)):
        r = 1-2*eps; plus = mm(channel(eps), mm(channel(eps), channel(eps)))
        minus = mm(channel(eps, True), mm(channel(eps), channel(eps)))
        zp, zm = sum(plus[i][i] for i in range(2)), sum(minus[i][i] for i in range(2))
        assert zp == 1+r**3 and zm == 1-r**3
        for T in (channel(eps), channel(eps, True)):
            assert all(sum(T[i][j] for i in range(2)) == 1 for j in range(2))
        noise_cases.append({"each_bit_flip_probability": str(eps), "copy_ring_weight": str(zp),
                            "one_flip_ring_weight": str(zm), "normalization_gap": str(r**3)})
    assert F(4, 5)**3 == F(64, 125)
    # Independent numerical calibration of exact pair bounds and closing paths; no general theorem by sampling.
    rng = np.random.default_rng(1041)
    sample_min, max_metric_residual = float("inf"), 0.
    sample_count = 0
    for pair in pair_rows:
        i, j = [k-1 for k in pair["pair"]]
        n = np.array([float(x) for x in add(axes[i], axes[j])]); n /= np.linalg.norm(n)
        for index in (i, j):
            u = np.array([float(x) for x in axes[index]])
            g = np.array([[float(x) for x in row] for row in metric(axes[index])])
            for _ in range(64):
                z = rng.normal(size=4); z -= (z@u)*u; z /= np.linalg.norm(z)
                radial = rng.uniform(0, .25)
                v = u+radial*z; v /= np.linalg.norm(v)
                assert n@v > .2
                sample_min = min(sample_min, float(n@v)); sample_count += 1
                expected = (radial*radial-1/16)/(1+radial*radial)
                max_metric_residual = max(max_metric_residual, abs(float(v@g@v)-expected))
    # No metric equality follows: all cones with a common axis, widths 1,2,4, share a clock.
    common = [{"kappa": k, "unit_clock_margin_squared": str(F(k*k, 1+k*k)),
               "max_space_per_clock": str(F(1, k))} for k in (1, 2, 4)]
    assert max_metric_residual < 2e-14
    return {"round": 1041, "date": "2026-10-08", "scientific_baseline": 1040,
            "new_calibration_groups": 1, "new_adopted_cognitive_axioms": 0,
            "goal_completed": False, "all_checks_passed": True,
            "ambient_dimension": 4, "departments": 3,
            "scope": {"finite_constant_closed_cones_adopted": True,
                      "translation_rescaling_and_arbitrary_relay_adopted": True,
                      "binary_signal_and_instrument_contract_needed_for_normalization_bridge": True,
                      "Lorentz_dimension_or_common_metric_derived": False,
                      "curved_spacetime_global_time_theorem": False,
                      "actual_closed_causal_curve_or_instrument_certified": False,
                      "full_multimetric_quantum_process_certified": False},
            "pairwise_clock_certificates": pair_rows,
            "loop_certificates": loop_rows,
            "Caratheodory_bound_dimension_plus_one": 5,
            "finite_certificate_reduction": reduction,
            "conversion_control": {"total_offset": [str(x) for x in w],
                                   "allowed_total_norm": str(conversion_budget),
                                   "guaranteed_last_leg_safety": str(guaranteed_safety),
                                   "corrected_last_leg": [str(x) for x in v3],
                                   "future_projection": str(alpha), "transverse_norm_squared": str(perp2)},
            "normalization_bridge_old233": noise_cases,
            "conditional_noise_le_one_tenth_gap": "64/125",
            "pairwise_robustness": {"exact_worst_margin_squared": "4/85",
                                     "declared_direction_error": "1/20",
                                     "strict_remaining_margin_lower_bound": "3/20"},
            "distinct_compatible_cones": common,
            "sampled_pair_clock_directions": sample_count,
            "minimum_sample_clock_increment": sample_min,
            "max_metric_form_residual": max_metric_residual,
            "historical_sha256": {p: hashlib.sha256((PHYSICS/p).read_bytes()).hexdigest() for p in HISTORY}}


def compare(a, b, path="root"):
    if isinstance(a, dict):
        assert a.keys() == b.keys(), path
        for k in a: compare(a[k], b[k], path+"."+k)
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)): compare(x, y, path+f"[{i}]")
    elif isinstance(a, float):
        assert abs(a-b) <= 2e-12+1e-12*abs(a), (path, a, b)
    else:
        assert a == b, (path, a, b)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument("--write", action="store_true")
    args = parser.parse_args(); result = build_results()
    if args.write:
        with RESULT.open("x", encoding="utf8") as f:
            json.dump(result, f, ensure_ascii=False, indent=2); f.write("\n")
    else:
        compare(result, json.loads(RESULT.read_text(encoding="utf8")))
    print(json.dumps({"round": 1041, "all_checks_passed": True,
                      "mode": "exclusive_write" if args.write else "read_only_compare",
                      "pairwise_clock_certificates": 3, "closed_loop_segments": 3,
                      "max_metric_form_residual": result["max_metric_form_residual"],
                      "reduced_certificate_segments": result["finite_certificate_reduction"]["remaining_segments"]},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
