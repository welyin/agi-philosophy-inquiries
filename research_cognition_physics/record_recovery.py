"""Algebraic reconstruction, affine recovery obstruction, and disk identity."""

import argparse
import json
import math
import platform
import unittest
from pathlib import Path

import numpy as np

from marked_event_limit import ALPHA, KAPPA, active_matrix


def normalized_output(matrix, state):
    result = matrix @ state
    if result[0] <= 0:
        raise ValueError("Normalization needs positive branch weight.")
    return result / result[0]


def inverse_witness(outcome):
    sign = 2 * outcome - 1
    legal_output = np.array([1, -sign * ALPHA, 0])
    inferred = normalized_output(np.linalg.inv(active_matrix(outcome)), legal_output)
    return {"outcome": outcome, "legal_output": legal_output.tolist(),
            "normalized_inverse": inferred.tolist(),
            "excess_over_allowed_radius": float(np.linalg.norm(inferred[1:]) - ALPHA)}


def recovery_curvature(outcome):
    """Required corrected branch is non-affine even on three collinear inputs."""
    matrix = active_matrix(outcome)
    center = np.array([1, 0, 0])
    plus, minus = np.array([1, ALPHA, 0]), np.array([1, -ALPHA, 0])
    target = lambda state: (matrix @ state)[0] * state
    return target(plus) + target(minus) - 2 * target(center)


def boundary_fit(effect):
    """Fit a linear map to a purported nondisturbing branch on the disk boundary.

    Nonconstant effects generate a second harmonic that no 3x3 map can fit.
    The analytic argument in note 12, not this finite fit, proves the theorem.
    """
    angles = 2 * math.pi * np.arange(128) / 128
    states = np.column_stack((np.ones(128), ALPHA * np.cos(angles), ALPHA * np.sin(angles)))
    required = (states @ effect)[:, None] * states
    coefficients = np.linalg.lstsq(states, required, rcond=None)[0]
    residual = required - states @ coefficients
    harmonic_cos = 2 / len(states) * np.cos(2 * angles) @ required[:, 1:]
    harmonic_sin = 2 / len(states) * np.sin(2 * angles) @ required[:, 1:]
    return coefficients.T, residual, harmonic_cos, harmonic_sin


class RecordRecoveryTests(unittest.TestCase):
    def test_inverting_actual_image_recovers_input_summary(self):
        for outcome in (0, 1):
            matrix = active_matrix(outcome)
            for angle in np.linspace(0, 2 * math.pi, 65):
                initial = np.array([1, ALPHA * math.cos(angle), ALPHA * math.sin(angle)])
                after = normalized_output(matrix, initial)
                recovered = normalized_output(np.linalg.inv(matrix), after)
                np.testing.assert_allclose(recovered, initial, atol=1e-14)

    def test_same_inverse_is_not_positive_on_the_full_preparation_cone(self):
        for outcome in (0, 1):
            row = inverse_witness(outcome)
            self.assertAlmostEqual(np.linalg.norm(row["normalized_inverse"][1:]), 1)
            self.assertGreater(row["excess_over_allowed_radius"], 0.01)

    def test_required_conditional_restoration_violates_mixture_affinity(self):
        for outcome in (0, 1):
            np.testing.assert_allclose(recovery_curvature(outcome),
                                       [0, (2 * outcome - 1) * KAPPA * ALPHA**2, 0], atol=1e-14)

    def test_constant_effect_allows_only_the_expected_boundary_map(self):
        matrix, residual, cosine, sine = boundary_fit(np.array([0.4, 0, 0]))
        np.testing.assert_allclose(matrix, 0.4 * np.eye(3), atol=1e-14)
        self.assertLess(float(np.max(abs(residual))), 1e-14)

    def test_nonconstant_effect_has_the_analytic_second_harmonic(self):
        for effect in (np.array([0.5, KAPPA / 2, 0]), np.array([0.5, 0.1, -0.2])):
            _, residual, cosine, sine = boundary_fit(effect)
            b, c = effect[1:]
            np.testing.assert_allclose(cosine, ALPHA**2 / 2 * np.array([b, -c]), atol=1e-14)
            np.testing.assert_allclose(sine, ALPHA**2 / 2 * np.array([c, b]), atol=1e-14)
            self.assertGreater(float(np.max(abs(residual))), 0.01)

    def test_input_correlated_backup_is_a_different_operation_scope(self):
        # With a source label already supplied, re-preparing that known source
        # is possible. That does not define a source-independent affine inverse.
        first, second = np.array([1, ALPHA, 0]), np.array([1, -ALPHA, 0])
        prior = (first + second) / 2
        matrix = active_matrix(1)
        correct_posterior = ((matrix @ first)[0] * first + (matrix @ second)[0] * second)
        correct_posterior /= correct_posterior[0]
        self.assertGreater(float(np.linalg.norm(correct_posterior - prior)), 0.1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RecordRecoveryTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": tests.testsRun, "failures": len(tests.failures), "errors": len(tests.errors)},
               "inverse_witnesses": [inverse_witness(outcome) for outcome in (0, 1)],
               "conditional_recovery_affinity_defects": [recovery_curvature(outcome).tolist() for outcome in (0, 1)],
               "scope": "Inverting a known conditional summary is inference. Deterministic correction uses affine positive channels on the same disk, without an input-correlated backup. Boundary extremality proves even average identity recovery would force constant record effects. Tests illustrate but do not replace the proof."}
    if args.write_results:
        Path(__file__).with_name("record_recovery_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
