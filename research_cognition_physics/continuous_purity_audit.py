"""Round 186: a finite classical pure-state path has a sharp 1/2 barrier.

The barrier concerns complete finite-simplex states in total variation distance.
It is not a test against every hidden-variable model or infinite phase space.
"""

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import unittest

import numpy as np


def classical_impurity(p):
    return 1 - max(p)


def real_rotation(theta):
    return np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])


class ContinuousPurityAuditTests(unittest.TestCase):
    def test_distance_to_nearest_classical_pure_state_has_exact_formula(self):
        for p in ([F(1, 2), F(1, 3), F(1, 6)], [F(3, 4), F(1, 4)], [F(0), F(0), F(1)]):
            distances = [sum(abs(value - int(i == j)) for i, value in enumerate(p)) / 2
                         for j in range(len(p))]
            self.assertEqual(min(distances), classical_impurity(p))

    def test_direct_and_extra_label_routes_attain_but_do_not_avoid_half_barrier(self):
        for dimension in (2, 3, 8):
            for first, second in zip(range(dimension - 1), range(1, dimension)):
                impurities = []
                for step in range(21):
                    t = F(step, 20)
                    p = [F(0)] * dimension
                    p[first], p[second] = 1 - t, t
                    impurities.append(classical_impurity(p))
                self.assertEqual(max(impurities), F(1, 2))

    def test_continuous_stochastic_interpolation_is_not_reversible(self):
        swap = np.array([[0., 1.], [1., 0.]])
        for t in (.1, .25, .75, .9):
            transition = (1 - t) * np.eye(2) + t * swap
            self.assertLess(np.linalg.inv(transition).min(), 0)
        self.assertEqual(np.linalg.matrix_rank((np.eye(2) + swap) / 2), 1)

    def test_real_quantum_path_is_pure_and_reversible_throughout(self):
        p0 = np.diag([1., 0.])
        for theta in np.linspace(0, np.pi / 2, 21):
            rotation = real_rotation(theta)
            rho = rotation @ p0 @ rotation.T
            np.testing.assert_allclose(rho @ rho, rho, atol=3e-16)
            np.testing.assert_allclose(rotation.T @ rho @ rotation, p0, atol=5e-16)
        np.testing.assert_allclose(real_rotation(np.pi / 2) @ p0 @ real_rotation(np.pi / 2).T,
                                   np.diag([0., 1.]), atol=1e-16)

    def test_one_reading_can_hide_difference_between_pure_rotation_and_classical_mixing(self):
        plus_effect = np.array([[.5, .5], [.5, .5]])
        for theta in (.2, np.pi / 4, 1.):
            rotation = real_rotation(theta)
            pure = rotation @ np.diag([1., 0.]) @ rotation.T
            mixed = np.diag(np.diag(pure))
            np.testing.assert_array_equal(np.diag(pure), np.diag(mixed))
            self.assertAlmostEqual(np.trace(plus_effect @ mixed), .5)
            self.assertAlmostEqual(np.trace(plus_effect @ pure), (1 + np.sin(2 * theta)) / 2)

    def test_sampling_bound_needs_an_explicit_speed_bound(self):
        # For p(t)=(1-t,t), total-variation speed L is 1. An odd interval
        # count misses t=1/2 and attains the lower bound 1/2-L*Delta/2.
        for intervals in (3, 5, 9):
            samples = [classical_impurity([1 - F(k, intervals), F(k, intervals)])
                       for k in range(intervals + 1)]
            self.assertEqual(max(samples), F(1, 2) - F(1, 2 * intervals))
        # A transition compressed entirely between observed times is invisible
        # to the pure endpoints; continuity alone supplies no speed bound.
        samples = [[F(1), F(0)], [F(0), F(1)]]
        self.assertEqual(max(classical_impurity(p) for p in samples), 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(ContinuousPurityAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = {
        "round": 186,
        "finite_classical_path_worst_distance_to_pure_states_lower_bound": "1/2",
        "bound_is_attainable": True,
        "more_finite_classical_labels_avoid_bound": False,
        "real_quantum_rotation_stays_pure": True,
        "continuous_probabilities_imply_continuous_reversible_pure_motion": False,
        "sampled_impurity_bound_with_TV_speed_L_and_spacing_Delta": "max(0,1/2-L*Delta/2)",
        "finite_sampled_controls_without_speed_bound_certify_continuous_purity": False,
        "scope": "Complete finite classical simplex, not arbitrary hidden-state or infinite-dimensional models",
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("continuous_purity_audit_results.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
