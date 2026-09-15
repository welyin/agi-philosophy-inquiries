"""A positive continuous kernel for sufficiently weak fresh-read instruments.

Poisson convolution supplies the first-harmonic retention. Bounded corrections
match the complete selective affine update while retaining a positive margin.
Numerical angle grids verify integrals; they are not exact finite-label models.
"""

import argparse
import json
import math
import platform
import unittest
from pathlib import Path

import numpy as np

from bounded_responses import translation_matrix
from outcome_updates import WINDOW_ATTENUATION as ALPHA
from retention_tradeoff import fresh_branch_matrix


def poisson_kernel(offset, alpha=ALPHA):
    if not 0 < alpha < 1:
        raise ValueError("Use retention strictly between zero and one.")
    # Equivalent stable denominator: (1-alpha)^2 + 4*alpha*sin(offset/2)^2.
    return (1 - alpha**2) / ((1 - alpha)**2 + 4 * alpha * np.sin(np.asarray(offset) / 2)**2)


def sufficient_contrast(alpha=ALPHA):
    minimum = (1 - alpha) / (1 + alpha)
    linear = 1 + 2 * alpha
    # Stable positive root of 2*alpha*eta^2 + (1+2*alpha)*eta = minimum.
    return 2 * minimum / (linear + math.sqrt(linear**2 + 8 * alpha * minimum))


def certified_density_margin(contrast, alpha=ALPHA):
    if not 0 <= contrast <= 1:
        raise ValueError("Contrast must lie in [0,1].")
    shrink = contrast**2 / (1 + math.sqrt(1 - contrast**2))
    return ((1 - alpha) / (1 + alpha) - (1 + 2 * alpha) * contrast - 2 * alpha * shrink) / (4 * math.pi)


def kernel_density(output_angle, input_angle, outcome, contrast, setting=0.0):
    """Density with respect to output angle Lebesgue measure."""
    if outcome not in (0, 1) or not 0 <= contrast <= 1:
        raise ValueError("Use a binary outcome and contrast in [0,1].")
    outgoing = np.asarray(output_angle) - setting
    incoming = np.asarray(input_angle) - setting
    sign = 2 * outcome - 1
    tangent_loss = -contrast**2 / (1 + math.sqrt(1 - contrast**2))
    return (poisson_kernel(outgoing - incoming)
            + sign * contrast * (np.cos(outgoing) + 2 * ALPHA * np.cos(incoming))
            + 2 * ALPHA * tangent_loss * np.sin(outgoing) * np.sin(incoming)) / (4 * math.pi)


def summary_density(summary, angles):
    weight, u, v = np.asarray(summary, dtype=float)
    return (weight + (u * np.cos(angles) + v * np.sin(angles)) / ALPHA) / (2 * math.pi)


def integrate_output_density(summary, output_angles, outcome, contrast, setting=0.0, nodes=4096):
    incoming = np.arange(nodes) * 2 * math.pi / nodes
    kernel = kernel_density(np.asarray(output_angles)[:, None], incoming[None, :], outcome, contrast, setting)
    return kernel @ summary_density(summary, incoming) * (2 * math.pi / nodes)


def diagnostics():
    threshold = sufficient_contrast()
    contrasts = (0, 0.001, threshold, 0.01, 0.5)
    return {"original_radius": ALPHA, "poisson_kernel_minimum": (1 - ALPHA) / (1 + ALPHA),
            "sufficient_maximum_contrast_for_full_sequential_model": threshold,
            "selected_contrast": 0.001,
            "uniform_kernel_density_lower_bounds": [{"contrast": eta, "density_lower_bound": certified_density_margin(eta)}
                                                     for eta in contrasts],
            "single_read_threshold_from_round22": 2 / (math.pi * ALPHA),
            "excluded_contrast_from_five_read_witness": 0.5,
            "interpretation": "A negative lower bound is only failure of this sufficient certificate; it is not proof that this kernel or every other kernel is negative."}


class NoncontextualUpdateKernelTests(unittest.TestCase):
    def test_poisson_kernel_normalization_and_first_harmonic(self):
        angles = np.arange(4096) * 2 * math.pi / 4096
        values = poisson_kernel(angles)
        self.assertAlmostEqual(values.mean(), 1, places=12)
        self.assertAlmostEqual(np.mean(values * np.cos(angles)), ALPHA, places=12)
        self.assertAlmostEqual(np.mean(values * np.sin(angles)), 0, places=12)

    def test_sufficient_contrast_has_a_positive_uniform_margin(self):
        threshold = sufficient_contrast()
        self.assertGreater(threshold, 0.001)
        self.assertGreater(certified_density_margin(0.001), 0)
        self.assertGreaterEqual(certified_density_margin(threshold), 0)

    def test_pointwise_kernel_matches_uniform_lower_bound(self):
        grid = np.arange(128) * 2 * math.pi / 128
        for eta in (0, 0.001, sufficient_contrast()):
            for outcome in (0, 1):
                values = kernel_density(grid[:, None], grid[None, :], outcome, eta, 0.3)
                self.assertGreaterEqual(values.min() + 1e-14, certified_density_margin(eta))

    def test_complete_instrument_preserves_mass_at_each_input(self):
        grid = np.arange(4096) * 2 * math.pi / 4096
        incoming = np.array([0.0, 0.3, 1.1, 2.9, -0.6])
        for eta in (0, 0.001, sufficient_contrast()):
            integrals = sum(kernel_density(grid[:, None], incoming[None, :], outcome, eta, 0.4).mean(axis=0) * 2 * math.pi
                            for outcome in (0, 1))
            np.testing.assert_allclose(integrals, 1, atol=3e-12)

    def test_outcome_mass_is_bounded_and_has_correct_hidden_response(self):
        grid = np.arange(4096) * 2 * math.pi / 4096
        incoming = np.array([0.0, 0.2, 1.5, 3.1])
        eta, setting = 0.001, 0.4
        for mark in (0, 1):
            integrals = kernel_density(grid[:, None], incoming[None, :], mark, eta, setting).mean(axis=0) * 2 * math.pi
            expected = (1 + (2 * mark - 1) * 2 * ALPHA * eta * np.cos(incoming - setting)) / 2
            np.testing.assert_allclose(integrals, expected, atol=3e-12)
            self.assertGreater(integrals.min(), 0)
            self.assertLess(integrals.max(), 1)

    def test_selective_kernel_reproduces_entire_output_density(self):
        outgoing = np.linspace(-1, 7, 23)
        for summary in ([1, 0, 0], [1, ALPHA, 0], [1, 0, ALPHA], [0.4, 0.1, -0.2]):
            for eta in (0, 0.001, sufficient_contrast()):
                for setting in (0, 0.7):
                    for mark in (0, 1):
                        actual = integrate_output_density(summary, outgoing, mark, eta, setting)
                        expected_summary = fresh_branch_matrix(setting, mark, eta) @ summary
                        np.testing.assert_allclose(actual, summary_density(expected_summary, outgoing), atol=4e-12)

    def test_rotation_preserves_unique_affine_preparation_encoding(self):
        summary = np.array([1, 0.3, -0.4])
        angles = np.linspace(0, 7, 37)
        angle = 0.7
        np.testing.assert_allclose(summary_density(translation_matrix(angle, (1.0,)) @ summary, angles),
                                   summary_density(summary, angles - angle), atol=1e-15)

    def test_equivalent_mixture_recipes_remain_identical_after_update(self):
        angles = np.linspace(0, 2 * math.pi, 19)
        first, second = np.array([1, ALPHA, 0]), np.array([1, 0, ALPHA])
        mixed = (integrate_output_density(first, angles, 1, 0.001)
                 + integrate_output_density(second, angles, 1, 0.001)) / 2
        np.testing.assert_allclose(mixed, integrate_output_density((first + second) / 2, angles, 1, 0.001), atol=1e-14)

    def test_zero_strength_read_is_not_silently_replaced_by_identity(self):
        summary = np.array([1, 0.4, 0.3])
        output = sum(fresh_branch_matrix(0, mark, 0) @ summary for mark in (0, 1))
        np.testing.assert_allclose(output, [1, 0.4 * ALPHA, 0.3 * ALPHA], atol=1e-14)
        self.assertGreater(np.linalg.norm(output - summary), 0)

    def test_certificate_failure_does_not_get_reported_as_general_impossibility(self):
        self.assertLess(certified_density_margin(0.01), 0)
        self.assertLess(sufficient_contrast(), 0.01)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(NoncontextualUpdateKernelTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": checks.testsRun, "failures": len(checks.failures), "errors": len(checks.errors)},
               **diagnostics(),
               "scope": "Continuous preparation-noncontextual model exactly intertwines every fresh branch for contrasts up to the sufficient bound, with arbitrary rotations and arbitrary finite adaptive protocols. The full conditional hidden density, not only first moments, equals the unique affine preparation encoding. No claim of optimal contrast or finite-label realization. The original fresh eta=0 read retains its alpha disturbance."}
    if args.write_results:
        Path(__file__).with_name("noncontextual_update_kernel_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
