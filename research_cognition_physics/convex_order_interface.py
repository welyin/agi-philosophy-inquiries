"""Convex-order diagnostics for the fixed affine circle preparation encoding.

The existence equivalence is proved in note 29 using Strassen's theorem. Tests
of selected convex functions provide necessary conditions, not sufficiency.
"""

import argparse
import json
import math
import platform
import unittest
from pathlib import Path

import numpy as np

from atomic_reset_kernel import construction_threshold
from compensated_angle_kernel import preconvolution_density, warped_angle
from outcome_updates import WINDOW_ATTENUATION as ALPHA


def angular_vectors(points=32768):
    angles = (np.arange(points) + 0.5) * 2 * math.pi / points
    return angles, np.stack((np.cos(angles), np.sin(angles)), axis=1)


def convex_expectations(contrast, function, points=32768):
    angles, vectors = angular_vectors(points)
    source = float(np.mean(preconvolution_density(angles, contrast) * function(ALPHA * vectors)))
    target = float(np.mean(function(vectors)))
    return source, target


def absolute_factor(contrast):
    if not 0 <= contrast <= 1:
        raise ValueError("Use contrast in [0,1].")
    return math.sqrt(1 - contrast**2) + contrast * math.asin(contrast)


def absolute_threshold():
    lower, upper = 0.0, 0.5
    for _ in range(56):
        middle = (lower + upper) / 2
        if ALPHA * absolute_factor(middle) <= 1:
            lower = middle
        else:
            upper = middle
    return lower, upper


def quadratic_x_expectation(contrast):
    q = math.sqrt(1 - contrast**2)
    t = contrast / (1 + q)
    return ALPHA**2 * (1 + t**2 * (1 + 2 * q)) / 2


def diagnostics():
    rows = []
    for eta in (0.083, 0.11, 0.14, 0.145, 0.15, 0.17):
        rows.append({"contrast": eta,
                     "absolute_projection_source": 2 * ALPHA * absolute_factor(eta) / math.pi,
                     "absolute_projection_target": 2 / math.pi,
                     "quadratic_x_source": quadratic_x_expectation(eta), "quadratic_x_target": 0.5})
    return rows


class ConvexOrderInterfaceTests(unittest.TestCase):
    def test_two_marginals_have_unit_mass_and_zero_mean(self):
        angles, vectors = angular_vectors()
        for eta in (0, 0.083, 0.15):
            weights = preconvolution_density(angles, eta)
            self.assertAlmostEqual(weights.mean(), 1, places=13)
            np.testing.assert_allclose(np.mean(weights[:, None] * ALPHA * vectors, axis=0), 0, atol=1e-13)

    def test_output_label_pushforward_matches_inner_circle_law(self):
        angles, _ = angular_vectors()
        for eta in (0.03, 0.15):
            for power in (1, 2, 4):
                direct = 0.0
                for mark in (0, 1):
                    weight = (1 + (2 * mark - 1) * eta * np.cos(angles)) / 2
                    direct += np.mean(weight * (ALPHA * np.cos(warped_angle(angles, mark, eta)))**power)
                source, _ = convex_expectations(eta, lambda vectors: vectors[:, 0]**power)
                self.assertAlmostEqual(direct, source, places=12)

    def test_absolute_projection_formula_matches_independent_angular_integral(self):
        for eta in (0, 0.083, 0.145):
            source, target = convex_expectations(eta, lambda vectors: np.abs(vectors[:, 0]))
            self.assertAlmostEqual(source, 2 * ALPHA * absolute_factor(eta) / math.pi, places=8)
            self.assertAlmostEqual(target, 2 / math.pi, places=8)

    def test_quadratic_formula_matches_independent_integral(self):
        for eta in (0, 0.1, 0.17):
            source, target = convex_expectations(eta, lambda vectors: vectors[:, 0]**2)
            self.assertAlmostEqual(source, quadratic_x_expectation(eta), places=13)
            self.assertAlmostEqual(target, 0.5, places=13)

    def test_absolute_witness_excludes_a_kernel_even_when_quadratic_passes(self):
        self.assertGreater(ALPHA * absolute_factor(0.15), 1)
        self.assertLess(quadratic_x_expectation(0.15), 0.5)

    def test_known_constructive_range_passes_selected_convex_tests(self):
        functions = [lambda vectors: np.abs(vectors[:, 0]), lambda vectors: vectors[:, 0]**2,
                     lambda vectors: np.maximum(vectors[:, 0] - 0.3, 0), lambda vectors: np.linalg.norm(vectors, axis=1)]
        for function in functions:
            source, target = convex_expectations(construction_threshold(), function)
            self.assertLessEqual(source, target + 1e-8)

    def test_threshold_bracket_and_strict_monotonicity(self):
        lower, upper = absolute_threshold()
        self.assertGreater(lower, 0.144)
        self.assertLess(upper, 0.145)
        self.assertGreater(absolute_factor(0.145), absolute_factor(0.144))
        self.assertAlmostEqual(ALPHA * absolute_factor((lower + upper) / 2), 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ConvexOrderInterfaceTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": checks.testsRun, "failures": len(checks.failures), "errors": len(checks.errors)},
               "absolute_projection_threshold_floating_bracket": absolute_threshold(), "necessary_test_diagnostics": diagnostics(),
               "scope": "For the fixed affine circle encoding, full positive kernels are equivalent to convex order between alpha*n_psi with density R_eta and the uniform unit circle. Selected test functions are only necessary checks. Extension to arbitrary encodings is not asserted in this round."}
    if args.write_results:
        Path(__file__).with_name("convex_order_interface_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
