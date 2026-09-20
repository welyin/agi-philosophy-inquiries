"""Exact irrational rotation cannot be a finite time-homogeneous Markov model.

The theorem is analytic (note 19). Cyclic and finite-prefix constructions make
its exactness, horizon and finite-state assumptions explicit.
"""

import argparse
import cmath
import json
import math
import platform
import unittest
from pathlib import Path

import numpy as np

from bounded_responses import translation_matrix
from outcome_updates import WINDOW_ATTENUATION as ALPHA


def rotation_probabilities(angle, repeats, contrast=0.5):
    return (1 + contrast * ALPHA * np.cos(angle * np.arange(repeats + 1))) / 2


def cycle_matrix(size, advance=1):
    if isinstance(size, bool) or not isinstance(size, int) or size < 1:
        raise ValueError("A positive finite state count is required.")
    matrix = np.zeros((size, size))
    matrix[(np.arange(size) + advance) % size, np.arange(size)] = 1
    return matrix


def simulate_repeated_action(matrix, response, initial, repeats):
    state = np.asarray(initial, dtype=float)
    values = []
    for index in range(repeats + 1):
        values.append(float(response @ state))
        if index < repeats:
            state = matrix @ state
    return np.asarray(values)


def finite_prefix_model(angle, horizon, contrast=0.5):
    """A counter encoded as horizon+1 states; exactly fits just one finite prefix."""
    size = horizon + 1
    matrix = np.zeros((size, size))
    matrix[np.minimum(np.arange(size) + 1, horizon), np.arange(size)] = 1
    response = rotation_probabilities(angle, horizon, contrast)
    initial = np.eye(size)[:, 0]
    return matrix, response, initial


def nearest_peripheral_root(angle, state_budget):
    target = cmath.exp(1j * angle)
    candidates = []
    for denominator in range(1, state_budget + 1):
        numerator = round(angle / (2 * math.pi) * denominator)
        candidate = cmath.exp(2j * math.pi * numerator / denominator)
        candidates.append((abs(candidate - target), denominator, numerator))
    error, denominator, numerator = min(candidates)
    return {"state_budget": state_budget, "nearest_root_order_candidate": denominator,
            "numerator": numerator, "complex_distance": error}


def recurrence_residual(values, coefficients):
    degree = len(coefficients) - 1
    return np.array([np.dot(coefficients, values[start:start + degree + 1])
                     for start in range(len(values) - degree)])


class FiniteRotationObstructionTests(unittest.TestCase):
    def test_rotation_has_three_dimensional_linear_predictor(self):
        angle = 1.0
        matrix = translation_matrix(angle, (1.0,))
        values = simulate_repeated_action(matrix, np.array([0.5, 0.25, 0]), [1, ALPHA, 0], 100)
        np.testing.assert_allclose(values, rotation_probabilities(angle, 100), atol=4e-15)
        self.assertLess(matrix.min(), 0)  # This predictor is not a probability transition.

    def test_required_rotation_eigenvalues_are_present(self):
        matrix = translation_matrix(1, (1.0,))
        eigenvalues = np.linalg.eigvals(matrix)
        for value in (1, cmath.exp(1j), cmath.exp(-1j)):
            self.assertLess(float(np.min(abs(eigenvalues - value))), 1e-14)

    def test_analytic_cubic_annihilates_entire_sampled_rotation_sequence(self):
        angle = 1.0
        coefficient = 1 + 2 * math.cos(angle)
        polynomial = [-1, coefficient, -coefficient, 1]
        residual = recurrence_residual(rotation_probabilities(angle, 100), polynomial)
        self.assertLess(float(np.max(abs(residual))), 2e-15)
        roots = np.roots(polynomial[::-1])
        self.assertLess(float(np.min(abs(roots - cmath.exp(1j * angle)))), 1e-14)

    def test_cayley_hamilton_recurrence_for_an_independent_markov_model(self):
        rng = np.random.default_rng(1901)
        matrix = rng.random((5, 5))
        matrix /= matrix.sum(axis=0)
        values = simulate_repeated_action(matrix, rng.random(5), np.ones(5) / 5, 30)
        residual = recurrence_residual(values, np.poly(matrix)[::-1])
        self.assertLess(float(np.max(abs(residual))), 1e-14)

    def test_rational_rotation_has_an_exact_finite_cyclic_realization(self):
        for size, advance in ((3, 1), (5, 2), (7, 3), (16, 1)):
            matrix = cycle_matrix(size, advance)
            response = (1 + 0.5 * ALPHA * np.cos(2 * math.pi * np.arange(size) / size)) / 2
            actual = simulate_repeated_action(matrix, response, np.eye(size)[:, 0], 100)
            expected = rotation_probabilities(2 * math.pi * advance / size, 100)
            np.testing.assert_allclose(actual, expected, atol=3e-14)

    def test_cyclic_peripheral_spectrum_has_the_expected_roots(self):
        for size in (3, 5, 7):
            eigenvalues = np.linalg.eigvals(cycle_matrix(size))
            np.testing.assert_allclose(eigenvalues**size, 1, atol=1e-14)

    def test_finite_prefix_is_exact_then_fails_beyond_its_horizon(self):
        matrix, response, initial = finite_prefix_model(1, 24)
        actual = simulate_repeated_action(matrix, response, initial, 100)
        expected = rotation_probabilities(1, 100)
        np.testing.assert_allclose(actual[:25], expected[:25], atol=0)
        self.assertGreater(float(np.max(abs(actual[25:] - expected[25:]))), 0.3)
        np.testing.assert_array_equal(matrix.sum(axis=0), 1)

    def test_finite_root_approximations_improve_but_do_not_certify_impossibility(self):
        rows = [nearest_peripheral_root(1, budget) for budget in (4, 8, 16, 32, 64)]
        errors = [row["complex_distance"] for row in rows]
        self.assertTrue(all(a >= b for a, b in zip(errors, errors[1:])))
        self.assertGreater(min(errors), 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(FiniteRotationObstructionTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    matrix, response, initial = finite_prefix_model(1, 24)
    actual = simulate_repeated_action(matrix, response, initial, 100)
    expected = rotation_probabilities(1, 100)
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": tests.testsRun, "failures": len(tests.failures), "errors": len(tests.errors)},
               "exact_symbolic_angle_radians": 1,
               "finite_root_diagnostics": [nearest_peripheral_root(1, m) for m in (4, 8, 16, 32, 64)],
               "finite_prefix_control": {"states": 25, "exact_horizon": 24,
                                         "max_error_through_horizon": float(np.max(abs(actual[:25] - expected[:25]))),
                                         "max_error_to_repeat_100": float(np.max(abs(actual - expected)))},
               "scope": "No finite homogeneous Markov transition and fixed response can reproduce an undamped irrational rotation for every repeat count. The irrational angle is a mathematical exact input, not a claim of experimentally infinite precision. Finite prefixes and rational rotations admit finite classical realizations. No exclusion of continuous classical latent states is implied."}
    if args.write_results:
        Path(__file__).with_name("finite_rotation_obstruction_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
