"""Sharp preparation-noncontextual visibility for a disk, with one final read.

Only preparation noncontextuality is imposed. The continuous construction has
no claimed implementation of the fresh selective update instruments.
"""

import argparse
import json
import math
import platform
import unittest
from pathlib import Path

import numpy as np

from outcome_updates import WINDOW_ATTENUATION as ALPHA


THRESHOLD = 2 / math.pi


def unit(angle):
    return np.array([math.cos(angle), math.sin(angle)])


def preparation_density(state, angles):
    state = np.asarray(state, dtype=float)
    if state.shape != (2,) or np.linalg.norm(state) > 1 + 1e-13:
        raise ValueError("Use a normalized disk state of radius at most one.")
    angles = np.asarray(angles)
    return (1 + state[0] * np.cos(angles) + state[1] * np.sin(angles)) / (2 * math.pi)


def response_scale(visibility):
    if not 0 <= visibility <= THRESHOLD + 1e-14:
        raise ValueError("This noncontextual construction requires visibility <= 2/pi.")
    return min(1.0, math.pi * visibility / 2)


def continuous_probability(state, setting, visibility):
    """Independent quadrature on the two response semicircles."""
    scale = response_scale(visibility)
    nodes, weights = np.polynomial.legendre.leggauss(24)
    positive_angles = setting + nodes * math.pi / 2
    negative_angles = setting + math.pi + nodes * math.pi / 2
    positive_mass = float(weights @ preparation_density(state, positive_angles) * math.pi / 2)
    negative_mass = float(weights @ preparation_density(state, negative_angles) * math.pi / 2)
    return (1 + scale) * positive_mass / 2 + (1 - scale) * negative_mass / 2


def four_label_encoding(state):
    return (1 + np.array([[1, 0], [0, 1], [-1, 0], [0, -1]]) @ np.asarray(state)) / 4


def four_label_response(setting, visibility):
    if not 0 <= visibility <= 0.5:
        raise ValueError("This four-label construction is sufficient only up to visibility 1/2.")
    directions = np.array([[1, 0], [0, 1], [-1, 0], [0, -1]])
    return (1 + 2 * visibility * directions @ unit(setting)) / 2


def support_function(weights, vectors, settings):
    directions = np.stack((np.cos(settings), np.sin(settings)), axis=-1)
    return np.abs(directions @ np.asarray(vectors).T) @ np.asarray(weights)


def diagnostics():
    points = np.linspace(0, 2 * math.pi, 32768, endpoint=False)
    finite_rows = []
    for labels in (4, 8, 16, 32):
        angles = 2 * math.pi * np.arange(labels) / labels
        vectors = np.stack((np.cos(angles), np.sin(angles)), axis=1)
        support = support_function(np.full(labels, 1 / labels), vectors, points)
        finite_rows.append({"labels": labels, "necessary_visibility_ceiling_for_this_encoding": float(support.min()),
                            "support_max": float(support.max()), "angular_mean": float(support.mean())})
    return {"sharp_single_read_visibility_threshold": THRESHOLD,
            "original_radius": ALPHA, "original_contrast_threshold_for_all_disk_preparations": THRESHOLD / ALPHA,
            "four_preparation_witness_visibility_threshold": 1 / math.sqrt(2),
            "four_label_sufficient_visibility": 0.5,
            "example_missed_by_four_preparation_witness": {"contrast": 0.68, "visibility": 0.68 * ALPHA,
                                                         "four_preparation_success": (1 + 0.68 * ALPHA / math.sqrt(2)) / 2},
            "finite_encoding_support_diagnostics_not_general_optima": finite_rows}


class AffinePreparationThresholdTests(unittest.TestCase):
    def test_continuous_preparations_are_positive_and_normalized(self):
        points = np.linspace(0, 2 * math.pi, 4096, endpoint=False)
        for state in ([0, 0], [0.3, -0.4], [1, 0], unit(0.7)):
            values = preparation_density(state, points)
            self.assertGreaterEqual(values.min(), -1e-15)
            self.assertAlmostEqual(values.mean() * 2 * math.pi, 1)

    def test_actual_mixtures_have_unique_affine_density(self):
        points = np.linspace(-1, 7, 101)
        first, second = unit(0.3), unit(2.1)
        mixture = 0.3 * preparation_density(first, points) + 0.7 * preparation_density(second, points)
        np.testing.assert_allclose(mixture, preparation_density(0.3 * first + 0.7 * second, points), atol=1e-15)

    def test_opposite_diameter_recipes_have_identical_hidden_density(self):
        points = np.linspace(-2, 8, 105)
        for angle in (0, 0.2, 1.3):
            direction = unit(angle)
            np.testing.assert_allclose((preparation_density(direction, points) + preparation_density(-direction, points)) / 2,
                                       preparation_density([0, 0], points), atol=1e-15)

    def test_quadrature_reproduces_all_direction_single_read_probabilities(self):
        for visibility in (0, 0.2, 0.5, THRESHOLD):
            for state in ([0.2, -0.7], unit(0.7), unit(-1)):
                for setting in (0, 0.5, 1.9, -2):
                    expected = (1 + visibility * np.dot(state, unit(setting))) / 2
                    self.assertAlmostEqual(continuous_probability(state, setting, visibility), expected, places=13)

    def test_response_above_threshold_is_rejected(self):
        with self.assertRaises(ValueError):
            response_scale(THRESHOLD + 0.001)

    def test_rotation_is_exact_before_the_final_read(self):
        points = np.linspace(0, 2 * math.pi, 93)
        angle, offset = 0.2, 0.7
        np.testing.assert_allclose(preparation_density(unit(angle + offset), points),
                                   preparation_density(unit(angle), points - offset), atol=1e-15)

    def test_four_label_model_is_positive_affine_and_exact_below_half(self):
        for state in ([0.3, 0.5], unit(0.2), unit(1.2)):
            encoding = four_label_encoding(state)
            self.assertGreaterEqual(encoding.min(), 0)
            self.assertAlmostEqual(encoding.sum(), 1)
            for setting in (0, 0.5, 1.5):
                response = four_label_response(setting, 0.5)
                self.assertGreaterEqual(response.min(), 0)
                self.assertLessEqual(response.max(), 1)
                self.assertAlmostEqual(response @ encoding, (1 + 0.5 * unit(setting) @ state) / 2)

    def test_general_finite_encoding_angular_support_average(self):
        rng = np.random.default_rng(2202)
        angles = rng.uniform(0, 2 * math.pi, 7)
        radii = rng.uniform(0.1, 1, 7)
        vectors = np.stack((np.cos(angles), np.sin(angles)), axis=1) * radii[:, None]
        vectors = np.concatenate((vectors, -vectors))
        weights = np.tile(rng.dirichlet(np.ones(7)) / 2, 2)
        settings = np.linspace(0, 2 * math.pi, 32768, endpoint=False)
        values = support_function(weights, vectors, settings)
        self.assertAlmostEqual(values.mean(), THRESHOLD * np.sum(weights * np.linalg.norm(vectors, axis=1)), places=8)
        self.assertLess(values.mean(), THRESHOLD)

    def test_symmetric_finite_encodings_do_not_reach_continuous_threshold(self):
        for labels in (4, 8, 16):
            angles = 2 * math.pi * np.arange(labels) / labels
            vectors = np.stack((np.cos(angles), np.sin(angles)), axis=1)
            support = support_function(np.full(labels, 1 / labels), vectors, np.array([0.0, math.pi / labels]))
            self.assertLess(support.min(), THRESHOLD)
            self.assertGreater(support.max(), support.min())

    def test_single_four_preparation_witness_is_not_complete(self):
        visibility = 0.68 * ALPHA
        self.assertGreater(visibility, THRESHOLD)
        self.assertLess((1 + visibility / math.sqrt(2)) / 2, 0.75)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(AffinePreparationThresholdTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": checks.testsRun, "failures": len(checks.failures), "errors": len(checks.errors)},
               **diagnostics(),
               "scope": "Sharp 2/pi visibility threshold for preparation-noncontextual affine encoding of the entire disk with arbitrary bounded, possibly measurement-contextual, single-read responses. Continuous circle model attains it; no finite-label affine preparation encoding attains the threshold. No fresh selective-update kernel is constructed here."}
    if args.write_results:
        Path(__file__).with_name("affine_preparation_threshold_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
