"""Round 414: fixed finite translation costs and their polyhedral limit.

The input is a specified lattice action and a fixed symmetric alphabet, not a
derivation of spatial dimension or an implementation of the cognitive axioms.
Only this module's new result may be saved, with exclusive creation.
"""
import argparse
from collections import deque
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations, product
import json
from pathlib import Path
import platform
import unittest

import numpy as np


TARGET = Path(__file__).with_name("word_budget_limit_audit_results.json")
STEPS = ((2, 0, 0), (-2, 0, 0), (3, 0, 0), (-3, 0, 0),
         (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
DUAL = tuple((F(a, 3), F(b), F(c)) for a, b, c in product((-1, 1), repeat=3))


def bfs(steps, radius):
    zero = (0,) * len(steps[0])
    distance = {zero: 0}
    queue = deque([zero])
    while queue:
        p = queue.popleft()
        if distance[p] == radius:
            continue
        for s in steps:
            q = tuple(x + y for x, y in zip(p, s))
            if q not in distance:
                distance[q] = distance[p] + 1
                queue.append(q)
    return distance


def axis_cost(x):
    n = abs(int(x))
    return 0 if n == 0 else 2 if n == 1 else (n + 2) // 3


def word_cost(p):
    x, y, z = p
    return axis_cost(x) + abs(y) + abs(z)


def limit_norm(p):
    x, y, z = map(F, p)
    return abs(x) / 3 + abs(y) + abs(z)


def dual_norm(p):
    return max(sum(a * F(x) for a, x in zip(row, p)) for row in DUAL)


def solve_exact(matrix, rhs):
    a = [list(map(F, row)) + [F(b)] for row, b in zip(matrix, rhs)]
    n = len(a)
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            raise ValueError("singular")
        a[j], a[pivot] = a[pivot], a[j]
        div = a[j][j]
        a[j] = [v / div for v in a[j]]
        for i in range(n):
            if i != j:
                mul = a[i][j]
                a[i] = [v - mul * w for v, w in zip(a[i], a[j])]
    return tuple(a[i][-1] for i in range(n))


def polytope_vertices():
    vertices = set()
    for rows in combinations(DUAL, 3):
        try:
            p = solve_exact(rows, (1, 1, 1))
        except ValueError:
            continue
        if dual_norm(p) <= 1:
            vertices.add(p)
    return vertices


def signed_permutations():
    out = []
    for perm in permutations(range(3)):
        for signs in product((-1, 1), repeat=3):
            matrix = np.zeros((3, 3), dtype=int)
            for i in range(3):
                matrix[i, perm[i]] = signs[i]
            out.append(matrix)
    return out


def rational_action(matrix, p):
    # A = D matrix D^{-1}; D = diag(3,1,1). These are limit isometries,
    # not all isometries of the original lattice or the original alphabet.
    diagonal = (F(3), F(1), F(1))
    return tuple(diagonal[i] * sum(F(int(matrix[i, j])) * F(p[j]) / diagonal[j]
                                  for j in range(3)) for i in range(3))


def nearest_integer(value):
    # Exact round-to-nearest, with ties rounded upwards.
    shifted = F(value) + F(1, 2)
    return shifted.numerator // shifted.denominator


def exact_text(value):
    value = F(value)
    return f"{value.numerator}/{value.denominator}"


@lru_cache(None)
def report():
    one = bfs(((2,), (-2,), (3,), (-3,)), 20)
    one_errors = sum(axis_cost(p[0]) != n for p, n in one.items())
    ball = bfs(STEPS, 12)
    three_errors = sum(word_cost(p) != n for p, n in ball.items())
    gaps = [F(n) - limit_norm(p) for p, n in ball.items()]
    vertices = polytope_vertices()
    expected_vertices = {(F(a), F(b), F(c)) for a, b, c in
                         ((3,0,0),(-3,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))}

    # An independent half-space enumeration identifies conv(S) exactly.
    hull_correct = vertices == expected_vertices
    generators_contained = all(dual_norm(p) <= 1 for p in STEPS)
    vertices_are_generators = all(p in STEPS for p in vertices)
    norm_errors = sum(limit_norm(p) != dual_norm(p) for p in ball)

    matrices = signed_permutations()
    matrix_keys = {tuple(m.flat) for m in matrices}
    closure = all(tuple((a @ b).flat) in matrix_keys for a in matrices for b in matrices)
    polytope_preserved = all({rational_action(m, p) for p in vertices} == vertices for m in matrices)
    generic_point = (F(3, 7), F(2, 7), F(4, 7))
    orbit = {rational_action(m, generic_point) for m in matrices}

    # Map actual scaled lattice positions to fixed continuum test points.
    # The universal rounding bound is proved in the note, not inferred here.
    values = (F(-4, 5), F(-3, 10), F(1, 5), F(7, 10))
    targets = list(product(values, repeat=3))
    scaled = []
    for scale in (4, 8, 16, 32):
        grid = [tuple(nearest_integer(scale * x) for x in p) for p in targets]
        max_placement = max(limit_norm(tuple(F(x, scale) - y for x, y in zip(q, p)))
                            for p, q in zip(targets, grid))
        max_distortion = F(0)
        for i, p in enumerate(targets):
            for j, q in enumerate(targets):
                lattice_difference = tuple(a - b for a, b in zip(grid[i], grid[j]))
                continuum_difference = tuple(a - b for a, b in zip(p, q))
                discrepancy = abs(F(word_cost(lattice_difference), scale) - limit_norm(continuum_difference))
                max_distortion = max(max_distortion, discrepancy)
        scaled.append(dict(scale=scale,placement_error=exact_text(max_placement),
                           placement_bound=exact_text(F(7, 6 * scale)),
                           pair_distance_error=exact_text(max_distortion),
                           pair_distance_bound=exact_text(F(4, scale))))

    # A genuine qubit rotation need not be a symmetry of this task budget.
    theta = np.pi / 4
    rotation = np.array([[np.cos(theta), -np.sin(theta), 0.],
                         [np.sin(theta), np.cos(theta), 0.], [0., 0., 1.]])
    direction = np.array([1., 0., 0.])
    rotated = rotation @ direction
    sigma = (np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),
             np.array([[1,0],[0,-1]],complex))
    u = np.diag(np.exp(np.array([-1j, 1j]) * theta / 2))
    e = (np.eye(2) + sigma[0]) / 2
    target_e = (np.eye(2) + sum(v*s for v,s in zip(rotated,sigma))) / 2
    covariance = float(np.linalg.norm(u @ e @ u.conj().T - target_e))

    return dict(round=414,scientific_baseline=412,
        scope="A fixed finite lattice translation alphabet yields a proper word budget with a polyhedral continuum limit; this budget cannot supply the minimal continuous budget-preserving redirection required by rounds 385 and 412",
        model=dict(group="Z^3",steps=[list(s) for s in STEPS],step_cost=1,
                   limit_norm="abs(x)/3+abs(y)+abs(z)",sharp_uniform_gap="5/3"),
        independent_bfs=dict(one_dimensional_radius=20,one_dimensional_points=len(one),one_dimensional_mismatches=one_errors,
                             three_dimensional_radius=12,three_dimensional_points=len(ball),three_dimensional_mismatches=three_errors),
        polytope=dict(exact_vertex_enumeration_correct=hull_correct,generators_contained=generators_contained,
                      vertices_are_generators=vertices_are_generators,dual_norm_mismatches=norm_errors,vertex_count=len(vertices)),
        bounded_difference=dict(minimum=exact_text(min(gaps)),maximum=exact_text(max(gaps)),
                                maximizing_point=[1,0,0],analytic_bound_all_integer_points=True),
        scaled_lattice_comparison=scaled,
        limit_symmetry=dict(linear_isometry_count=len(matrix_keys),multiplication_closed=closure,
                            exact_vertex_set_preserved=polytope_preserved,generic_orbit_size=len(orbit),
                            original_discrete_alphabet_isometry_count_claimed=False),
        direction_rotation=dict(angle=theta,unit_direction_budget=1.,rotated_direction_budget=float(np.sum(abs(rotated))),
                                budget_ratio=float(np.sqrt(2)),qubit_covariance_error=covariance),
        fixed_finite_translation_alphabet_is_model_input=True,endpoint_group_identity_is_model_input=True,
        polynomial_growth_is_cognitive_axiom=False,schreier_quotient_automatically_a_group=False,
        arbitrary_finite_control_alphabet_ruled_out=False,random_walk_limit_identified_with_word_budget=False,
        finite_symmetry_count_selects_spatial_dimension=False,three_dimensional_topology_ruled_out=False,
        cognition_to_gr_refuted=False,new_cognitive_axiom_adopted=False)


class Audit(unittest.TestCase):
    def test_01_one_dimensional_word_cost_independent_bfs(self):
        self.assertEqual(report()["independent_bfs"]["one_dimensional_mismatches"],0)
        self.assertEqual(axis_cost(1),2)

    def test_02_complete_three_dimensional_word_ball(self):
        self.assertEqual(report()["independent_bfs"]["three_dimensional_mismatches"],0)
        self.assertGreater(report()["independent_bfs"]["three_dimensional_points"],5000)

    def test_03_exact_convex_hull_and_dual_norm(self):
        r=report()["polytope"]
        self.assertTrue(r["exact_vertex_enumeration_correct"] and r["generators_contained"] and r["vertices_are_generators"])
        self.assertEqual(r["dual_norm_mismatches"],0)
        self.assertEqual(r["vertex_count"],6)

    def test_04_sharp_bounded_difference(self):
        r=report()["bounded_difference"]
        self.assertEqual(F(r["minimum"]),0)
        self.assertEqual(F(r["maximum"]),F(5,3))

    def test_05_uniform_scaled_lattice_comparison(self):
        for r in report()["scaled_lattice_comparison"]:
            self.assertLessEqual(F(r["placement_error"]),F(r["placement_bound"]))
            self.assertLessEqual(F(r["pair_distance_error"]),F(r["pair_distance_bound"]))

    def test_06_all_limit_linear_isometries_form_a_finite_group(self):
        r=report()["limit_symmetry"]
        self.assertEqual(r["linear_isometry_count"],48)
        self.assertEqual(r["generic_orbit_size"],48)
        self.assertTrue(r["multiplication_closed"] and r["exact_vertex_set_preserved"])

    def test_07_qubit_rotation_does_not_preserve_this_budget(self):
        r=report()["direction_rotation"]
        self.assertLess(r["qubit_covariance_error"],1e-14)
        self.assertAlmostEqual(r["rotated_direction_budget"],np.sqrt(2),places=14)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    run=unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not run.wasSuccessful():
        raise SystemExit(1)
    result=dict(report(),checks=dict(run=run.testsRun,failures=len(run.failures),errors=len(run.errors)),
                runtime=dict(python=platform.python_version(),numpy=np.__version__))
    if args.write_results:
        with TARGET.open("x",encoding="utf-8",newline="\n") as out:
            out.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(result,ensure_ascii=False,indent=2))
