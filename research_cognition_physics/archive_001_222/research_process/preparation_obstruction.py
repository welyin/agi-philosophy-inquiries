"""A classical kernel passes first-moment checks but violates preparation closure."""

import argparse
import json
import math
import platform
import unittest
from pathlib import Path

import numpy as np

from bounded_responses import harmonic_features, probability_effect, translation_matrix
from operational_capacity import PREPARATION_HALF_WIDTH
from outcome_updates import WINDOW_ATTENUATION, validate_unit_interval


def candidate_matrix(setting, outcome, contrast):
    validate_unit_interval(contrast)
    probability_effect(setting, outcome)
    local = np.array([[0.5, (2 * outcome - 1) * contrast / 2, 0],
                      [0, 0, 0], [0, 0, math.sqrt(1 - contrast**2) / 2]])
    rotation = translation_matrix(setting, (1.0,))
    return rotation @ local @ rotation.T


def candidate_atoms(scores, setting, outcome, contrast):
    """Actual atomic kernel; deliberately NOT an admissible preparation update."""
    validate_unit_interval(contrast)
    probability_effect(setting, outcome)
    relative = np.asarray(scores, dtype=float) - setting
    if not np.all(np.isfinite(relative)):
        raise ValueError("Finite scores are required.")
    angle = np.arctan2(math.sqrt(1 - contrast**2) * np.sin(relative),
                       np.cos(relative) + (2 * outcome - 1) * contrast)
    locations = setting + np.stack((angle, math.pi - angle), axis=-1)
    probability = (1 + (2 * outcome - 1) * contrast * np.cos(relative)) / 2
    weights = np.broadcast_to(probability[..., None] / 2, locations.shape)
    return locations, weights


def output_harmonics(center=0.0, contrast=0.5, order=40):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    locations, masses = candidate_atoms(center + PREPARATION_HALF_WIDTH * nodes, 0, 1, contrast)
    return np.einsum("i,ij,ijk->k", weights / 2, masses,
                     harmonic_features(locations, (1.0, 2.0)))


def obstruction_data(contrast=0.5):
    validate_unit_interval(contrast)
    half_width = PREPARATION_HALF_WIDTH
    alpha = WINDOW_ATTENUATION
    q2 = 1 - contrast**2
    probability = (1 + contrast * alpha) / 2
    output_half_width = 2 * math.atan(math.sqrt((1 - contrast) / (1 + contrast))
                                    * math.tan(half_width / 2))
    harmonics = output_harmonics(contrast=contrast)
    beta = math.sin(2 * half_width) / (2 * half_width)
    lower_bound = 1 - q2 * half_width**2 / (3 * probability * (1 + contrast * math.cos(half_width)))
    return {
        "contrast": contrast,
        "input_half_width": half_width,
        "output_component_half_width": output_half_width,
        "output_support_components": [[-output_half_width, output_half_width],
                                      [math.pi - output_half_width, math.pi + output_half_width]],
        "branch_probability": probability,
        "conditional_first_moment_summary": (harmonics[:3] / harmonics[0]).tolist(),
        "conditional_second_cosine": float(harmonics[3] / harmonics[0]),
        "analytic_second_cosine_lower_bound": lower_bound,
        "allowed_preparation_second_harmonic_bound": beta,
        "second_harmonic_violation": float(harmonics[3] / harmonics[0] - beta),
        "formal_tangent_gain": math.sqrt(q2),
        "admissible": False if contrast == 0.5 else None,
    }


class PreparationObstructionTests(unittest.TestCase):
    def test_candidate_uses_fixed_effect(self):
        for outcome in (0, 1):
            effect = probability_effect(0.7, outcome)
            effect[1:] *= 0.5
            np.testing.assert_allclose(candidate_matrix(0.7, outcome, 0.5)[0], effect, atol=1e-14)

    def test_candidate_is_positive_on_all_sampled_summary_radii(self):
        for radius in (0, WINDOW_ATTENUATION / 2, WINDOW_ATTENUATION):
            states = harmonic_features(np.linspace(-math.pi, math.pi, 501), (1.0,))
            states[:, 1:] *= radius
            for outcome in (0, 1):
                outputs = states @ candidate_matrix(0, outcome, 0.5).T
                self.assertGreater(float(outputs[:, 0].min()), 0)
                self.assertLessEqual(float(np.max(np.linalg.norm(outputs[:, 1:], axis=1)
                                                  - WINDOW_ATTENUATION * outputs[:, 0])), 1e-14)

    def test_actual_kernel_has_valid_classical_mass(self):
        scores = np.linspace(-math.pi, math.pi, 401)
        total = np.zeros_like(scores)
        for outcome in (0, 1):
            _, weights = candidate_atoms(scores, 0, outcome, 0.5)
            self.assertGreaterEqual(float(weights.min()), 0)
            total += weights.sum(axis=1)
        np.testing.assert_allclose(total, 1)

    def test_actual_atomic_kernel_matches_first_moment_matrix(self):
        scores = np.linspace(-math.pi, math.pi, 307)
        for setting in (0, 0.7):
            for outcome in (0, 1):
                locations, weights = candidate_atoms(scores, setting, outcome, 0.5)
                moments = np.einsum("ij,ijk->ik", weights, harmonic_features(locations, (1.0,)))
                np.testing.assert_allclose(moments, harmonic_features(scores, (1.0,))
                                           @ candidate_matrix(setting, outcome, 0.5).T, atol=1e-14)

    def test_actual_kernel_meets_universal_pointwise_moment_bound(self):
        scores = np.linspace(-math.pi, math.pi, 301)
        outputs = harmonic_features(scores, (1.0,)) @ candidate_matrix(0, 1, 0.5).T
        self.assertLessEqual(float(np.max(np.linalg.norm(outputs[:, 1:], axis=1) - outputs[:, 0])), 1e-14)

    def test_formal_tangent_bound_is_saturated(self):
        contrast = 0.5
        angle = math.acos(-contrast)
        output = candidate_matrix(0, 1, contrast) @ harmonic_features(angle, (1.0,))
        self.assertAlmostEqual(float(output[2]), float(output[0]))
        self.assertAlmostEqual(candidate_matrix(0, 1, contrast)[2, 2] * 2, math.sqrt(1 - contrast**2))

    def test_input_window_has_legal_but_uninformative_output_first_moments(self):
        data = obstruction_data()
        np.testing.assert_allclose(data["conditional_first_moment_summary"], [1, 0, 0], atol=1e-14)

    def test_complete_output_support_cannot_contain_an_original_width_window(self):
        data = obstruction_data()
        self.assertGreater(data["output_component_half_width"], 0)
        self.assertLess(data["output_component_half_width"], PREPARATION_HALF_WIDTH)
        first, second = data["output_support_components"]
        self.assertGreater(second[0] - first[1], 2 * PREPARATION_HALF_WIDTH)

    def test_second_harmonic_proves_failure_beyond_first_moments(self):
        data = obstruction_data()
        self.assertGreater(data["analytic_second_cosine_lower_bound"], data["allowed_preparation_second_harmonic_bound"])
        self.assertGreater(data["conditional_second_cosine"], data["analytic_second_cosine_lower_bound"])
        self.assertGreater(data["second_harmonic_violation"], 0.02)

    def test_independent_quadrature_orders_agree(self):
        np.testing.assert_allclose(output_harmonics(order=16), output_harmonics(order=48), atol=1e-14)

    def test_endpoint_no_read_symmetrization_preserves_original_window_width(self):
        data = obstruction_data(0)
        self.assertAlmostEqual(data["output_component_half_width"], PREPARATION_HALF_WIDTH)
        self.assertAlmostEqual(data["conditional_second_cosine"], data["allowed_preparation_second_harmonic_bound"])

    def test_invalid_contrast_is_rejected(self):
        for contrast in (-1, 1.1, math.nan):
            with self.assertRaises(ValueError):
                candidate_matrix(0, 1, contrast)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(PreparationObstructionTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    results = {
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "automated_checks": {"run": tests.testsRun, "failures": len(tests.failures), "errors": len(tests.errors)},
        "fixed_contrast_witness": obstruction_data(),
        "scope": "This particular atomic kernel fails preparation closure despite positive first-moment maps and a valid pointwise classical kernel. It does not exclude other kernels realizing the same first-moment matrix. The general optimum remains unresolved.",
    }
    if args.write_results:
        Path(__file__).with_name("preparation_obstruction_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
