"""A universal preparation-noncontextual bound with all exact silent rotations.

Rotation invariance is derived for the affine preparation scores, not assumed
for the complete hidden variables. The proof covers arbitrary score radii.
"""

import argparse
import json
import math
import platform
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from convex_order_interface import absolute_factor, absolute_threshold
from outcome_updates import WINDOW_ATTENUATION as ALPHA
from stopped_readout_witness import alpha_interval


def radial_max_average(radius, contrast):
    if not 0 <= radius <= 1 or not 0 <= contrast <= 1:
        raise ValueError("Use radius and contrast in [0,1].")
    if radius <= contrast:
        return contrast
    return 2 / math.pi * (math.sqrt(radius**2 - contrast**2) + contrast * math.asin(contrast / radius))


def radial_jensen_gap(radius, contrast):
    return ALPHA * radial_max_average(radius, contrast) - 2 * radius / math.pi


def backward_score(score, outcome, contrast):
    score = np.asarray(score, dtype=float)
    sign = 2 * outcome - 1
    q = math.sqrt(1 - contrast**2)
    denominator = 1 + sign * contrast * score[..., 0]
    value = ALPHA * np.stack((score[..., 0] + sign * contrast, q * score[..., 1]), axis=-1) / denominator[..., None]
    return denominator / 2, value


def rational_exclusion_margin(contrast=Fraction(29, 200)):
    alpha_low, _ = alpha_interval()
    # G(eta)=1+eta^2/2+eta^4/24+eta^6/80+..., all remaining terms positive.
    return alpha_low * (1 + contrast**2 / 2 + contrast**4 / 24) - 1


def radial_diagnostics():
    return [{"contrast": eta, "radii_and_gaps": [{"radius": radius, "jensen_gap": radial_jensen_gap(radius, eta)}
                                                 for radius in (0, 0.1, 0.3, 0.7, 0.99, 1)]}
            for eta in (0.11, 0.14, 0.145, 0.15)]


class RotationNoncontextualBoundTests(unittest.TestCase):
    def test_backward_score_matches_affine_branch_coefficients(self):
        score = np.array([0.3, -0.4])
        state = np.array([0.2, 0.6])
        eta = 0.145
        for mark in (0, 1):
            weight, backward = backward_score(score, mark, eta)
            sign = 2 * mark - 1
            output_weight = (1 + sign * ALPHA * eta * state[0]) / 2
            output_vector = np.array([ALPHA * state[0] + sign * eta, ALPHA * math.sqrt(1 - eta**2) * state[1]]) / 2
            self.assertAlmostEqual(weight * (1 + backward @ state), output_weight + score @ output_vector)

    def test_sum_of_branch_absolute_scores_is_the_maximum_identity(self):
        angles = np.linspace(-3, 3, 47)
        for radius in (0, 0.3, 1):
            scores = radius * np.stack((np.cos(angles), np.sin(angles)), axis=1)
            accumulated = np.zeros(len(angles))
            for mark in (0, 1):
                weight, backward = backward_score(scores, mark, 0.145)
                accumulated += weight * np.abs(backward[:, 0])
            np.testing.assert_allclose(accumulated, ALPHA * np.maximum(np.abs(scores[:, 0]), 0.145), atol=1e-14)

    def test_radial_angular_integral_formula(self):
        angles = (np.arange(65536) + 0.5) * 2 * math.pi / 65536
        for radius in (0, 0.1, 0.3, 1):
            actual = np.mean(np.maximum(np.abs(radius * np.cos(angles)), 0.145))
            self.assertAlmostEqual(actual, radial_max_average(radius, 0.145), places=8)

    def test_each_radius_obeys_the_same_lower_ratio_bound(self):
        for eta in (0.03, 0.11, 0.145):
            for radius in np.linspace(0, 1, 101):
                self.assertGreaterEqual(radial_max_average(radius, eta) + 1e-14,
                                        2 * radius / math.pi * absolute_factor(eta))

    def test_rational_margin_proves_exclusion_at_point_145(self):
        self.assertGreater(rational_exclusion_margin(), 0)
        self.assertGreater(ALPHA * absolute_factor(0.145), 1)

    def test_all_radial_components_violate_beyond_threshold(self):
        for radius in np.linspace(0, 1, 1001):
            self.assertGreater(radial_jensen_gap(radius, 0.145), 0)

    def test_endpoint_equality_requires_unit_score_radius(self):
        lower, upper = absolute_threshold()
        eta = (lower + upper) / 2
        self.assertAlmostEqual(radial_jensen_gap(1, eta), 0, places=13)
        for radius in (0, 0.1, 0.8, 0.99, 0.9999):
            self.assertGreater(radial_jensen_gap(radius, eta), 0)

    def test_exact_grid_rotation_has_zero_conditional_score_variance(self):
        size, advance = 16, 3
        angles = np.arange(size) * 2 * math.pi / size
        scores = np.stack((np.cos(angles), np.sin(angles)), axis=1)
        angle = advance * 2 * math.pi / size
        rotation = np.array([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]])
        for incoming in range(size):
            outgoing = (incoming + advance) % size
            np.testing.assert_allclose(scores[incoming], rotation.T @ scores[outgoing], atol=1e-14)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RotationNoncontextualBoundTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": checks.testsRun, "failures": len(checks.failures), "errors": len(checks.errors)},
               "universal_necessary_threshold_floating_bracket": absolute_threshold(),
               "rationally_excluded_contrast": 0.145, "rational_margin": str(rational_exclusion_margin()),
               "rational_margin_decimal": float(rational_exclusion_margin()), "radial_diagnostics": radial_diagnostics(),
               "scope": "For full preparation-noncontextual models respecting mixing, conditional preparations and every exact silent rotation, affine scores have a rotation-invariant base law. Jensen then forces alpha*(sqrt(1-eta^2)+eta*asin(eta))<=1 for arbitrary hidden spaces and score radii. This necessary bound is not claimed attainable. At equality unit score radius and a sign-preserving conditional coupling are necessary."}
    if args.write_results:
        Path(__file__).with_name("rotation_noncontextual_bound_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
