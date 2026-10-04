"""Round 183: self-duality, homogeneity and continuous pure transitivity differ.

A regular pentagon is self-dual but not homogeneous. A four-dimensional ball
has all three single-system properties, without being a complex qubit. No
locally tomographic homogeneous composite for that ball is asserted here.
"""

import argparse
import json
from pathlib import Path
import unittest

import numpy as np
from cone_duality_audit import ray_fixing_constraints


def pentagon_rays():
    angles = 2 * np.pi * np.arange(5) / 5
    radius = 1 / np.sqrt(np.cos(np.pi / 5))
    return np.column_stack((np.ones(5), radius * np.cos(angles), radius * np.sin(angles)))


def lorentz_to_state(vector):
    vector = np.asarray(vector, dtype=float)
    norm_squared = vector @ vector
    if norm_squared >= 1:
        raise ValueError("An interior target is required.")
    gamma = 1 / np.sqrt(1 - norm_squared)
    boost = np.block([[np.array([[gamma]]), gamma * vector[None, :]],
                      [gamma * vector[:, None], np.eye(len(vector))
                       + gamma ** 2 / (gamma + 1) * np.outer(vector, vector)]])
    return boost / gamma


def spin_product(a, b):
    return np.r_[a[0] * b[0] + a[1:] @ b[1:], a[0] * b[1:] + b[0] * a[1:]]


class JordanScopeAuditTests(unittest.TestCase):
    def test_pentagon_dual_facets_are_its_own_rays(self):
        rays = pentagon_rays()
        gram = rays @ rays.T
        self.assertGreaterEqual(gram.min(), -1e-15)
        for i in range(5):
            zeros = {(i + 2) % 5, (i + 3) % 5}
            self.assertEqual(set(np.flatnonzero(np.abs(gram[i]) < 1e-14)), zeros)
            self.assertEqual(np.linalg.matrix_rank(rays[list(zeros)]), 2)

    def test_pentagon_connected_automorphisms_only_scale(self):
        rays = pentagon_rays()
        equations = ray_fixing_constraints(rays)
        self.assertEqual(equations.shape[1] - np.linalg.matrix_rank(equations, tol=1e-10), 1)
        expansion = np.linalg.solve(rays[:3].T, rays[3])
        self.assertTrue(np.all(np.abs(expansion) > .1))

    def test_lorentz_duality_and_explicit_separator_outside_cone(self):
        rng = np.random.default_rng(183)
        for dimension in (2, 3, 4, 5):
            for _ in range(12):
                x, y = rng.normal(size=(2, dimension))
                a, b = np.r_[np.linalg.norm(x) + .1, x], np.r_[np.linalg.norm(y), y]
                self.assertGreaterEqual(a @ b, -1e-15)
                outside = np.r_[np.linalg.norm(x) - .2, x]
                separator = np.r_[1., -x / np.linalg.norm(x)]
                self.assertAlmostEqual(outside @ separator, -.2)

    def test_lorentz_boosts_connect_interior_states_and_have_positive_inverses(self):
        r = np.array([.2, -.3, .1, .4])
        s = np.array([-.1, .2, .3, -.2])
        tr, ts = lorentz_to_state(r), lorentz_to_state(s)
        metric = np.diag([1, -1, -1, -1, -1])
        for vector, transform in ((r, tr), (s, ts)):
            np.testing.assert_allclose(transform[:, 0], np.r_[1, vector], atol=2e-16)
            np.testing.assert_allclose(transform.T @ metric @ transform,
                                       (1 - vector @ vector) * metric, atol=5e-16)
            np.testing.assert_allclose(np.linalg.inv(transform),
                                       lorentz_to_state(-vector) / (1 - vector @ vector), atol=5e-16)
        np.testing.assert_allclose(ts @ np.linalg.inv(tr) @ np.r_[1, r], np.r_[1, s], atol=3e-16)

    def test_four_ball_rotations_continuously_connect_pure_states(self):
        for angle in np.linspace(0, np.pi, 9):
            rotation = np.eye(4)
            rotation[:2, :2] = [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
            np.testing.assert_allclose(rotation.T @ rotation, np.eye(4), atol=3e-16)
            state = rotation @ np.array([1., 0., 0., 0.])
            self.assertAlmostEqual(np.linalg.norm(state), 1)

    def test_spin_product_satisfies_jordan_identity_with_integer_witnesses(self):
        rng = np.random.default_rng(283)
        for _ in range(30):
            a, b = rng.integers(-3, 4, size=(2, 5))
            aa = spin_product(a, a)
            np.testing.assert_array_equal(spin_product(aa, spin_product(a, b)),
                                           spin_product(a, spin_product(aa, b)))
            self.assertGreaterEqual(aa[0] ** 2 - aa[1:] @ aa[1:], 0)

    def test_capacity_two_and_conditional_rank_four_dimension_obstruction(self):
        axis = np.array([1., 0., 0., 0.])
        states = [np.r_[1, axis], np.r_[1, -axis]]
        effects = [state / 2 for state in states]
        np.testing.assert_array_equal(np.array(effects) @ np.array(states).T, np.eye(2))
        np.testing.assert_array_equal(sum(effects), np.r_[1, np.zeros(4)])
        self.assertNotIn(5 ** 2, (10, 16, 28))
        # This last mismatch applies only after additionally requiring a
        # rank-four simple Jordan composite; it does not define a composite.


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(JordanScopeAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = {
        "round": 183,
        "pentagon_self_dual": True,
        "pentagon_homogeneous": False,
        "four_ball_self_dual": True,
        "four_ball_homogeneous": True,
        "four_ball_continuous_pure_transitivity": True,
        "four_ball_normalized_dimension": 4,
        "four_ball_linear_dimension": 5,
        "four_ball_perfect_distinguishability_capacity": 2,
        "locally_tomographic_joint_dimension_required": 25,
        "rank_four_simple_jordan_dimensions": [10, 16, 28],
        "four_ball_full_compositional_countertheory_constructed": False,
        "homogeneity_and_transitivity_are_derived_from_cognition": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("jordan_scope_audit_results.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
