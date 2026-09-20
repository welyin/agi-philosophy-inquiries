"""Sharp longitudinal retention bound for the declared general classical class."""

import argparse
import json
import math
import platform
import unittest
from pathlib import Path

import numpy as np

from outcome_updates import WINDOW_ATTENUATION as ALPHA, partial_branch_matrix, validate_unit_interval
from randomized_retention import nonselective_matrix
from retention_tradeoff import fresh_branch_matrix, quadrature
from operational_capacity import PREPARATION_HALF_WIDTH as HALF_WIDTH


def optimal_longitudinal_gain(contrast):
    validate_unit_interval(contrast)
    return 1 - (1 - ALPHA) * contrast


def fixed_probe_error_floor(contrast, probe_contrast=0.5):
    validate_unit_interval(contrast)
    validate_unit_interval(probe_contrast)
    return probe_contrast * ALPHA * (1 - ALPHA) * contrast / 2


def necessary_bound_slacks(matrix, outcome, contrast):
    """Two necessary inequalities whose sum bounds the cosine coefficient.

    One uses a whole legal input window, the other pointwise kernel positivity.
    Neither is by itself a full feasibility criterion.
    """
    sign = 2 * outcome - 1
    a, b, _ = matrix[1]
    preparation_slack = ALPHA * (1 + contrast * ALPHA) / 2 - (sign * a + ALPHA * b)
    pointwise_slack = (1 - contrast) / 2 - (b - sign * a)
    return np.array([preparation_slack, pointwise_slack])


def disk_positive_candidate(outcome, contrast):
    """A summary-cone-positive matrix that has no compatible pointwise kernel."""
    sign = 2 * outcome - 1
    return np.array([[1, sign * contrast, 0],
                     [sign * ALPHA**2 * contrast, 1, 0],
                     [0, 0, math.sqrt(1 - ALPHA**2 * contrast**2)]]) / 2


def direct_baseline_window_summary(center, outcome, contrast):
    nodes, weights = quadrature(32)
    x = center + HALF_WIDTH * nodes
    sign = 2 * outcome - 1
    idle_mass = (1 - contrast) / 2
    active_mass = contrast * (1 + sign * np.cos(x)) / 2
    # Actual kernel: idle leaves x; active emits mu_0 or mu_pi.
    first = idle_mass + active_mass
    cosine = idle_mass * np.cos(x) + active_mass * sign * ALPHA
    sine = idle_mass * np.sin(x)
    return np.array([weights @ first, weights @ cosine, weights @ sine])


def comparison_rows():
    rows = []
    for eta in (0, 0.005, 0.1, 0.5, 1):
        baseline = sum(partial_branch_matrix(0, mark, eta) for mark in (0, 1))
        randomized = nonselective_matrix(eta)
        candidate = disk_positive_candidate(1, eta)
        rows.append({"contrast": eta,
                     "general_longitudinal_upper_attained": optimal_longitudinal_gain(eta),
                     "baseline_longitudinal_gain": float(baseline[1, 1]),
                     "baseline_tangent_gain": float(baseline[2, 2]),
                     "round_9_longitudinal_gain": float(randomized[1, 1]),
                     "round_9_tangent_gain": float(randomized[2, 2]),
                     "fixed_half_contrast_probe_worst_input_error_floor": fixed_probe_error_floor(eta),
                     "disk_candidate_pointwise_slack": float(necessary_bound_slacks(candidate, 1, eta)[1])})
    return rows


class LongitudinalBoundTests(unittest.TestCase):
    def test_baseline_attains_both_necessary_inequalities(self):
        for eta in (0, 0.005, 0.5, 1):
            for mark in (0, 1):
                matrix = partial_branch_matrix(0, mark, eta)
                np.testing.assert_allclose(necessary_bound_slacks(matrix, mark, eta), 0, atol=1e-14)
                self.assertAlmostEqual(matrix[1, 1], optimal_longitudinal_gain(eta) / 2)

    def test_actual_classical_kernel_integration_matches_attaining_matrix(self):
        for center in (-1.2, 0, 0.9, math.pi):
            initial = np.array([1, ALPHA * math.cos(center), ALPHA * math.sin(center)])
            for eta in (0, 0.1, 0.5, 1):
                for mark in (0, 1):
                    direct = direct_baseline_window_summary(center, mark, eta)
                    np.testing.assert_allclose(direct, partial_branch_matrix(0, mark, eta) @ initial, atol=1e-14)

    def test_slack_identity_is_independent_of_drift_and_transverse_coefficient(self):
        rng = np.random.default_rng(1301)
        for _ in range(100):
            eta = float(rng.random())
            matrix = rng.normal(size=(3, 3))
            for mark in (0, 1):
                summed = necessary_bound_slacks(matrix, mark, eta).sum()
                self.assertAlmostEqual(summed, (1 + ALPHA) * (optimal_longitudinal_gain(eta) / 2 - matrix[1, 1]))

    def test_previously_admissible_families_respect_the_bound(self):
        for eta in (0, 0.01, 0.5, 0.99, 1):
            for mark in (0, 1):
                for longitudinal in (0, 0.3, 1):
                    slacks = necessary_bound_slacks(fresh_branch_matrix(0, mark, eta, longitudinal), mark, eta)
                    self.assertGreaterEqual(float(slacks.min()), -1e-14)

    def test_disk_positive_candidate_really_preserves_summary_cone(self):
        angles = np.linspace(0, 2 * math.pi, 257)
        states = np.column_stack((np.ones_like(angles), ALPHA * np.cos(angles), ALPHA * np.sin(angles)))
        for eta in (0, 0.1, 0.5, 1):
            for mark in (0, 1):
                after = states @ disk_positive_candidate(mark, eta).T
                self.assertLessEqual(float(np.max(np.linalg.norm(after[:, 1:], axis=1) - ALPHA * after[:, 0])), 1e-14)

    def test_same_candidate_violates_actual_kernel_bound_at_a_specific_point(self):
        for eta in (0.005, 0.5, 1):
            for mark in (0, 1):
                sign = 2 * mark - 1
                after = disk_positive_candidate(mark, eta) @ [1, -sign, 0]
                violation = abs(after[1]) - after[0]
                self.assertAlmostEqual(violation, eta * (1 - ALPHA**2) / 2)
                self.assertGreater(violation, 0)

    def test_unknown_drift_cannot_cancel_both_opposite_input_errors(self):
        eta = 0.5
        for drift in np.linspace(-0.5, 0.5, 51):
            for gain in (0, 0.5, optimal_longitudinal_gain(eta)):
                errors = [abs(drift + sign * ALPHA * (gain - 1)) / 4 for sign in (-1, 1)]
                self.assertGreaterEqual(max(errors), fixed_probe_error_floor(eta) - 1e-14)

    def test_best_longitudinal_mechanism_does_not_optimize_transverse_task(self):
        row = next(row for row in comparison_rows() if row["contrast"] == 0.5)
        self.assertGreater(row["baseline_longitudinal_gain"], row["round_9_longitudinal_gain"])
        self.assertLess(row["baseline_tangent_gain"], row["round_9_tangent_gain"])

    def test_positive_linear_floor_excludes_uniform_quadratic_disturbance(self):
        for eta in (0.1, 0.01, 0.001, 0.0001):
            self.assertAlmostEqual(fixed_probe_error_floor(eta) / eta, ALPHA * (1 - ALPHA) / 4)
        self.assertGreater(fixed_probe_error_floor(0.0001) / 0.0001**2, 25)

    def test_zero_contrast_allows_identity(self):
        np.testing.assert_array_equal(sum(partial_branch_matrix(0, mark, 0) for mark in (0, 1)), np.eye(3))
        self.assertEqual(fixed_probe_error_floor(0), 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(LongitudinalBoundTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": tests.testsRun, "failures": len(tests.failures), "errors": len(tests.errors)},
               "half_width": HALF_WIDTH, "alpha": ALPHA,
               "sharp_longitudinal_loss_slope": 1 - ALPHA,
               "comparison_rows": comparison_rows(),
               "scope": "Analytic optimum of signed unconditional longitudinal gain, over periodic classical branch kernels with exact three-moment closure and full original-window preparation preservation. Not an optimum of conditional tangent retention. Numerical checks support the written proof; they are not an exhaustive kernel search."}
    if args.write_results:
        Path(__file__).with_name("longitudinal_bound_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
