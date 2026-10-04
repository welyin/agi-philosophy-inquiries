"""Round 416: a paid finite steering alphabet can isotropize large-scale cost.

The additive position sector, its dimension, rotations, real-amplitude body-axis
translations, and their costs are inputs. This is not a generation of space.
Only this round's new JSON can be saved, using exclusive creation.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
import json
from pathlib import Path
import platform
import unittest

import numpy as np


TARGET = Path(__file__).with_name("finite_steering_budget_audit_results.json")
TAU = F(1)
M = 8
RHO = F(3, 4)
DELTA = F(18, 25)
INV = (1, 0, 3, 2)
I3 = tuple(tuple(F(i == j) for j in range(3)) for i in range(3))
RX = ((F(1), F(0), F(0)), (F(0), F(3, 5), F(-4, 5)),
      (F(0), F(4, 5), F(3, 5)))
RZ = ((F(3, 5), F(-4, 5), F(0)), (F(4, 5), F(3, 5), F(0)),
      (F(0), F(0), F(1)))


def transpose(a):
    return tuple(zip(*a))


ROT = (RX, transpose(RX), RZ, transpose(RZ))


def matmul(a, b):
    return tuple(tuple(sum(x * y for x, y in zip(row, col)) for col in zip(*b))
                 for row in a)


def exact_solve(a, b):
    rows = [list(map(F, row)) + [F(v)] for row, v in zip(a, b)]
    for j in range(len(rows)):
        pivot = next((i for i in range(j, len(rows)) if rows[i][j]), None)
        if pivot is None:
            raise ValueError("singular")
        rows[j], rows[pivot] = rows[pivot], rows[j]
        divisor = rows[j][j]
        rows[j] = [x / divisor for x in rows[j]]
        for i in range(len(rows)):
            if i != j:
                scale = rows[i][j]
                rows[i] = [x - scale * y for x, y in zip(rows[i], rows[j])]
    return tuple(row[-1] for row in rows)


def squared_norm(x):
    return sum(F(v) ** 2 for v in x)


def inverse_word(word):
    return tuple(INV[i] for i in reversed(word))


def rotation_word(word):
    q = I3
    for symbol in word:
        q = matmul(q, ROT[symbol])
    return q


def power_word(k, plus):
    return (plus if k >= 0 else plus + 1,) * abs(k)


def angle_power(k):
    c, s = F(1), F(0)
    step_s = F(4 if k >= 0 else -4, 5)
    for _ in range(abs(k)):
        c, s = c * F(3, 5) - s * step_s, s * F(3, 5) + c * step_s
    return c, s


@lru_cache(None)
def dictionary():
    found = {}
    for k, ell in product(range(-M, M + 1), repeat=2):
        ck, sk = angle_power(k)
        cl, sl = angle_power(ell)
        direction = (cl, ck * sl, sk * sl)
        word = power_word(k, 0) + power_word(ell, 2)
        if direction not in found or len(word) < len(found[direction]):
            found[direction] = word
    return tuple((v, found[v]) for v in sorted(found))


def atan_interval(x, terms=20):
    # Alternating-series rational enclosure, 0 < x < 1.
    x = F(x)
    partial = sum((-1) ** j * x ** (2 * j + 1) / (2 * j + 1)
                  for j in range(terms))
    other = partial + (-1) ** terms * x ** (2 * terms + 1) / (2 * terms + 1)
    return min(partial, other), max(partial, other)


def scale_interval(k, interval):
    values = (k * interval[0], k * interval[1])
    return min(values), max(values)


def angular_net_certificate():
    # Machin: pi/4 = 4 atan(1/5) - atan(1/239).
    a, b = atan_interval(F(1, 5)), atan_interval(F(1, 239))
    pi = (16 * a[0] - 4 * b[1], 16 * a[1] - 4 * b[0])
    c = atan_interval(F(1, 7))
    alpha = (pi[0] / 4 + c[0], pi[1] / 4 + c[1])
    # alpha = atan(4/3) = pi/4 + atan(1/7).
    angular = []
    for k in range(-M, M + 1):
        lo, hi = scale_interval(k, alpha)
        turn = int(np.floor(float((lo + hi) / (2 * (pi[0] + pi[1])))))
        tlo, thi = scale_interval(2 * turn, pi)
        angular.append((lo - thi, hi - tlo))
    angular.sort(key=lambda x: x[0])
    in_circle = all(0 <= lo <= hi < 2 * pi[0] for lo, hi in angular)
    ordered = all(angular[j][1] < angular[j + 1][0] for j in range(len(angular) - 1))
    gaps = [angular[j + 1][1] - angular[j][0] for j in range(len(angular) - 1)]
    gaps.append(angular[0][1] + 2 * pi[1] - angular[-1][0])
    # Last included cosine term is negative, hence this is a lower bound.
    factorial = 1
    cosine_lower = F(1)
    for j in range(1, 10):
        factorial *= (2 * j - 1) * (2 * j)
        cosine_lower += (-1) ** j * DELTA ** (2 * j) / factorial
    return dict(ordered=ordered, in_circle=in_circle, gap_upper=max(gaps), cosine_lower=cosine_lower,
                delta=DELTA, alpha_interval=alpha)


def decompose(target):
    """Find a route, not the exact optimal cost C(target).

    Float candidates are checked by a second, exact rational solve. The general
    convex-hull theorem is proved in the note and is not inferred from this search.
    """
    target = tuple(map(F, target))
    if not squared_norm(target):
        return ()
    entries = dictionary()
    vectors = np.array([v for v, _ in entries], dtype=float)
    x = np.array(target, dtype=float)
    dots = vectors @ (x / np.linalg.norm(x))
    for index in np.argsort(-dots)[:5]:
        v, word = entries[index]
        pivot = next(j for j in range(3) if v[j])
        coefficient = target[pivot] / v[pivot]
        if coefficient >= 0 and all(coefficient * v[j] == target[j] for j in range(3)):
            return ((int(index), coefficient),)
    nearest = np.argsort(-dots)[:48]
    best = None
    for selected in combinations(nearest, 3):
        matrix = vectors[list(selected)].T
        if abs(np.linalg.det(matrix)) < 1e-11:
            continue
        coefficients = np.linalg.solve(matrix, x)
        if min(coefficients) < -1e-12:
            continue
        if best is not None and sum(coefficients) >= best[0]:
            continue
        exact_matrix = tuple(tuple(entries[int(index)][0][j] for index in selected)
                             for j in range(3))
        exact_coefficients = exact_solve(exact_matrix, target)
        if min(exact_coefficients) < 0:
            continue
        best = (sum(exact_coefficients), tuple((int(index), value)
                                             for index, value in zip(selected, exact_coefficients)
                                             if value))
    if best is None:
        raise RuntimeError("The finite illustrative route search found no admissible cone")
    return best[1]


def route(plan, scale=1):
    steps = []
    for index, amount in plan:
        word = dictionary()[index][1]
        steps.extend(("r", symbol) for symbol in word)
        steps.append(("t", F(scale) * amount))
        steps.extend(("r", symbol) for symbol in inverse_word(word))
    return tuple(steps)


def inverse_route(steps):
    return tuple((kind, INV[value] if kind == "r" else -value)
                 for kind, value in reversed(steps))


def route_cost(steps):
    return sum(TAU if kind == "r" else abs(value) for kind, value in steps)


def se_matrix(step, exact=False):
    kind, value = step
    result = [[F(i == j) for j in range(4)] for i in range(4)]
    if kind == "r":
        for i in range(3):
            for j in range(3):
                result[i][j] = ROT[value][i][j]
    else:
        result[0][3] = F(value)
    return tuple(map(tuple, result)) if exact else np.array(result, dtype=float)


def execute(steps, exact=False):
    if exact:
        state = tuple(tuple(F(i == j) for j in range(4)) for i in range(4))
        for step in steps:
            state = matmul(state, se_matrix(step, True))
        return state
    state = np.eye(4)
    for step in steps:
        state = state @ se_matrix(step)
    return state


def short_returning_rotation_words():
    return tuple(word for length in range(3) for word in product(range(4), repeat=length)
                 if rotation_word(word) == I3)


def transverse_route(epsilon):
    # Rz(alpha), body translation 5epsilon/4, inverse Rz, body -3epsilon/4.
    return (("r", 2), ("t", F(5, 4) * epsilon), ("r", 3), ("t", F(-3, 4) * epsilon))


def as_text(value):
    value = F(value)
    return f"{value.numerator}/{value.denominator}"


@lru_cache(None)
def report():
    entries = dictionary()
    certificate = angular_net_certificate()
    targets = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1),
               (2, -1, 3), (-3, 2, 1), (-2, -3, -1), (3, 1, -2), (0, 0, 0))
    plans = [decompose(p) for p in targets]
    decomposition_ok = True
    translation_bound_ok = True
    max_endpoint = max_posture = 0.
    endpoint_table = []
    for target, plan in zip(targets, plans):
        reconstructed = tuple(sum(amount * entries[index][0][j] for index, amount in plan)
                              for j in range(3))
        decomposition_ok &= reconstructed == tuple(map(F, target))
        amount = sum(a for _, a in plan)
        translation_bound_ok &= amount ** 2 * RHO ** 2 <= squared_norm(target)
        steps = route(plan)
        state = execute(steps)
        max_endpoint = max(max_endpoint, float(np.linalg.norm(state[:3, 3] - target)))
        max_posture = max(max_posture, float(np.linalg.norm(state[:3, :3] - np.eye(3))))
        endpoint_table.append(dict(target=list(target),directions=len(plan),
                                   translation_cost=float(amount),
                                   rotation_calls=sum(kind == "r" for kind, _ in steps),
                                   total_cost=float(route_cost(steps)),
                                   normalized_translation_ratio=(float(amount) / np.linalg.norm(target)
                                                                 if squared_norm(target) else 1.)))
    exact_state = execute(route(plans[3]), True)
    exact_route_ok = exact_state == tuple(tuple(F(i == j) if j < 3 or i == 3
                                              else F(targets[3][i]) for j in range(4))
                                         for i in range(4))
    all_reverse_cost = True
    reverse_endpoint = concat_endpoint = 0.
    for target, plan in zip(targets, plans):
        steps = route(plan)
        back = inverse_route(steps)
        all_reverse_cost &= route_cost(steps) == route_cost(back)
        reverse_endpoint = max(reverse_endpoint, float(np.linalg.norm(execute(back)[:3, 3] + target)))
        concat_endpoint = max(concat_endpoint, float(np.linalg.norm(execute(steps + back) - np.eye(4))))
    pair_steps = route(plans[3]) + route(plans[4])
    expected_pair = np.array(targets[3]) + np.array(targets[4])
    pair_error = float(np.linalg.norm(execute(pair_steps)[:3, 3] - expected_pair))

    scale_rows = []
    for scale in (1, 16, 256, 4096):
        max_surplus = 0.
        max_normalized_error = 0.
        all_budget_bounds = True
        for target, plan in zip(targets, plans):
            steps = route(plan, scale)
            cost = route_cost(steps)
            translation = F(scale) * sum(a for _, a in plan)
            turn_cost = cost - translation
            all_budget_bounds &= turn_cost <= 2 * 4 * TAU * (2 * M)
            all_budget_bounds &= translation ** 2 * RHO ** 2 <= scale ** 2 * squared_norm(target)
            max_surplus = max(max_surplus, float(cost / scale) - np.linalg.norm(target))
            state = execute(steps)
            max_normalized_error = max(max_normalized_error,
                                       float(np.linalg.norm(state[:3, 3] / scale - target)))
        scale_rows.append(dict(scale=scale,maximum_route_upper_surplus=max_surplus,
                               theorem_upper_for_this_window=(1 / float(RHO) - 1)
                               * max(np.linalg.norm(p) for p in targets) + 128 / scale,
                               maximum_normalized_endpoint_error=max_normalized_error,
                               exact_budget_components_within_bounds=bool(all_budget_bounds)))

    epsilons = (F(1, 2), F(1, 8), F(1, 32), F(1, 128))
    transverse = []
    for epsilon in epsilons:
        steps = transverse_route(epsilon)
        state = execute(steps, True)
        transverse.append(dict(epsilon=as_text(epsilon),cost=as_text(route_cost(steps)),
                               proven_optimum=as_text(2 * TAU + 2 * epsilon),
                               exact_endpoint_and_posture=(state ==
                                  ((F(1),F(0),F(0),F(0)),(F(0),F(1),F(0),epsilon),
                                   (F(0),F(0),F(1),F(0)),(F(0),F(0),F(0),F(1))))))
    separation_points = [TAU / (4 * i) for i in range(1, 13)]
    minimum_pair_cost = min(2 * TAU + 2 * abs(a - b)
                            for a, b in combinations(separation_points, 2))

    return dict(round=416,scientific_baseline=414,
        scope="Paid finite steering with dense axis orbit has Euclidean large-scale position budget; the input geometry and microscopic cost topology remain distinct",
        model=dict(dimension=3,rotation_cosine="3/5",rotation_sine="4/5",rotation_cost=as_text(TAU),
                   body_axis="e1",translation_cost="abs(t)",final_posture="identity"),
        alphabet=dict(exact_orthogonality=all(matmul(r, transpose(r)) == I3 for r in ROT),
                      exact_inverse_pairs=all(matmul(ROT[i], ROT[INV[i]]) == I3 for i in range(4)),
                      rotations_noncommute=matmul(RX, RZ) != matmul(RZ, RX),
                      irrationality_and_density_proved_analytically=True),
        finite_net=dict(maximum_each_axis_exponent=M,distinct_directions=len(entries),maximum_word_length=2*M,
                        exact_direction_norms=all(squared_norm(v) == 1 for v, _ in entries),
                        exact_words_match_directions=all(tuple(row[0] for row in rotation_word(w)) == v
                                                        for v, w in entries),
                        interval_order_verified=certificate["ordered"],
                        exact_circle_representatives=certificate["in_circle"],
                        certified_gap_upper=float(certificate["gap_upper"]),
                        chosen_delta=as_text(DELTA),gap_strictly_below_delta=certificate["gap_upper"] < DELTA,
                        cosine_lower=float(certificate["cosine_lower"]),
                        certified_inner_ball_radius=as_text(RHO),
                        cosine_strictly_above_inner_radius=certificate["cosine_lower"] > RHO),
        exact_decomposition=dict(all_targets_exact=bool(decomposition_ok),
                                 all_translation_bounds=bool(translation_bound_ok),
                                 maximum_directions=max(map(len, plans))),
        actual_routes=dict(samples=endpoint_table,maximum_endpoint_error=max_endpoint,
                           maximum_posture_error=max_posture,one_full_rational_se3_route_exact=exact_route_ok,
                           example_target=[1,1,1],
                           example_chunks=[dict(outward_rotation_symbols=list(entries[index][1]),
                                                body_translation=as_text(amount),
                                                return_rotation_symbols=list(inverse_word(entries[index][1])))
                                           for index, amount in plans[3]],
                           route_cost_is_upper_bound_not_computed_optimum=True),
        composition=dict(inverse_costs_match=all_reverse_cost,inverse_endpoint_error=reverse_endpoint,
                         concatenation_identity_error=concat_endpoint,two_endpoint_sum_error=pair_error),
        rescaled_route_bounds=scale_rows,
        microscopic_boundary=dict(short_returning_rotation_words=[list(w) for w in short_returning_rotation_words()],
                                  transverse_samples=transverse,
                                  origin_cost=0,transverse_limit_cost=2.,euclidean_continuity=False,
                                  separated_sample_count=len(separation_points),
                                  maximum_sample_origin_cost=as_text(max(2*TAU+2*x for x in separation_points)),
                                  minimum_sample_pair_cost=as_text(minimum_pair_cost),
                                  infinite_separated_sequence_proved_analytically=True,
                                  cost_metric_is_proper=False),
        input_position_dimension_is_derived=False,euclidean_rotation_action_is_input=True,
        arbitrary_rotation_is_a_free_primitive=False,steering_cost_omitted=False,
        euclidean_sublevel_compactness_proved=True,uniform_limit_is_on_declared_euclidean_windows=True,
        universal_rate_from_density_alone=False,all_round385_conditions_satisfied=False,
        proper_pointed_gh_convergence_claimed=False,cognition_to_space_completed=False)


class Audit(unittest.TestCase):
    def test_01_exact_finite_rotation_alphabet(self):
        r = report()["alphabet"]
        self.assertTrue(r["exact_orthogonality"] and r["exact_inverse_pairs"] and r["rotations_noncommute"])

    def test_02_global_finite_net_has_rational_interval_certificate(self):
        r = report()["finite_net"]
        self.assertTrue(r["exact_direction_norms"] and r["exact_words_match_directions"])
        self.assertTrue(r["interval_order_verified"] and r["exact_circle_representatives"]
                        and r["gap_strictly_below_delta"])
        self.assertTrue(r["cosine_strictly_above_inner_radius"])

    def test_03_exact_endpoint_decompositions_and_translation_bound(self):
        r = report()["exact_decomposition"]
        self.assertTrue(r["all_targets_exact"] and r["all_translation_bounds"])
        self.assertLessEqual(r["maximum_directions"], 4)

    def test_04_actual_se3_execution_and_returned_posture(self):
        r = report()["actual_routes"]
        self.assertLess(r["maximum_endpoint_error"], 1e-11)
        self.assertLess(r["maximum_posture_error"], 1e-12)
        self.assertTrue(r["one_full_rational_se3_route_exact"])

    def test_05_route_reversal_and_composition(self):
        r = report()["composition"]
        self.assertTrue(r["inverse_costs_match"])
        for key in ("inverse_endpoint_error", "concatenation_identity_error", "two_endpoint_sum_error"):
            self.assertLess(r[key], 1e-10)

    def test_06_large_scale_budget_brackets_and_real_endpoints(self):
        for r in report()["rescaled_route_bounds"]:
            self.assertTrue(r["exact_budget_components_within_bounds"])
            self.assertLessEqual(r["maximum_route_upper_surplus"], r["theorem_upper_for_this_window"])
            self.assertLess(r["maximum_normalized_endpoint_error"], 1e-11)

    def test_07_positive_transverse_toll_and_exact_short_routes(self):
        r = report()["microscopic_boundary"]
        self.assertEqual(len(r["short_returning_rotation_words"]), 5)
        for item in r["transverse_samples"]:
            self.assertTrue(item["exact_endpoint_and_posture"])
            self.assertEqual(item["cost"], item["proven_optimum"])

    def test_08_same_bounded_cost_ball_has_uniformly_separated_points(self):
        r = report()["microscopic_boundary"]
        self.assertLessEqual(F(r["maximum_sample_origin_cost"]), F(5, 2))
        self.assertGreater(F(r["minimum_sample_pair_cost"]), 2 * TAU)
        self.assertFalse(r["euclidean_continuity"] or r["cost_metric_is_proper"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    run = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not run.wasSuccessful():
        raise SystemExit(1)
    result = dict(report(),checks=dict(run=run.testsRun,failures=len(run.failures),errors=len(run.errors)),
                  runtime=dict(python=platform.python_version(),numpy=np.__version__))
    if args.write_results:
        with TARGET.open("x", encoding="utf-8", newline="\n") as out:
            out.write(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))
