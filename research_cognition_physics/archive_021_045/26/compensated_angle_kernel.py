"""A normalized, positive angle-warp kernel for a larger noncontextual range.

The correction has zero zeroth and first Fourier moments. Its quadratic onset
avoids the linear positivity cost of the previous additive kernel.
"""

import argparse
import json
import math
import platform
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from noncontextual_update_kernel import poisson_kernel, summary_density
from outcome_updates import WINDOW_ATTENUATION as ALPHA
from retention_tradeoff import fresh_branch_matrix


def parameters(contrast):
    if not 0 <= contrast < 1:
        raise ValueError("Use contrast in [0,1).")
    q = math.sqrt(1 - contrast**2)
    t = contrast / (1 + q)
    return q, t, (ALPHA * t)**2


def warped_angle(outgoing, outcome, contrast):
    q, _, _ = parameters(contrast)
    sign = 2 * outcome - 1
    return np.arctan2(q * np.sin(outgoing), np.cos(outgoing) + sign * contrast)


def preconvolution_density(angles, contrast):
    q, _, _ = parameters(contrast)
    cosine = np.cos(angles)
    return q**3 * ((1 - contrast * cosine)**-2 + (1 + contrast * cosine)**-2) / 2


def compensation(angles, contrast):
    """g - 1, with angles relative to the read direction."""
    q, _, r = parameters(contrast)
    z = r * np.exp(2j * np.asarray(angles))
    return 2 * np.real(z / (1 - z) + 2 * q * z / (1 - z)**2)


def maximum_compensation(contrast):
    q, _, r = parameters(contrast)
    return 2 * (r / (1 - r) + 2 * q * r / (1 - r)**2)


def positivity_bracket(contrast):
    return (1 - ALPHA) / (1 + ALPHA) - maximum_compensation(contrast)


def certificate_zero_bracket():
    lower, upper = 0.05, 0.07
    for _ in range(48):
        middle = (lower + upper) / 2
        if positivity_bracket(middle) >= 0:
            lower = middle
        else:
            upper = middle
    return lower, upper


def rational_certificate(ceiling=Fraction(59, 1000)):
    upper_alpha = Fraction(95, 96) + Fraction(1, 30720)
    lower_q = 1 - ceiling**2 / 2 - ceiling**4 / 2
    upper_r = (upper_alpha * ceiling / (1 + lower_q))**2
    upper_compensation = 2 * upper_r / (1 - upper_r) + 4 * upper_r / (1 - upper_r)**2
    lower_poisson = (1 - upper_alpha) / (1 + upper_alpha)
    return lower_poisson - upper_compensation


def kernel_density(output_angle, input_angle, outcome, contrast, setting=0.0):
    if outcome not in (0, 1):
        raise ValueError("Use a binary outcome.")
    outgoing, incoming = np.asarray(output_angle) - setting, np.asarray(input_angle) - setting
    weight = (1 + (2 * outcome - 1) * contrast * np.cos(outgoing)) / 2
    center = warped_angle(outgoing, outcome, contrast)
    return weight * (poisson_kernel(incoming - center) - compensation(incoming, contrast)) / (2 * math.pi)


def integrate_density(summary, outgoing, outcome, contrast, setting=0.0, nodes=4096):
    incoming = np.arange(nodes) * 2 * math.pi / nodes
    kernel = kernel_density(np.asarray(outgoing)[:, None], incoming[None, :], outcome, contrast, setting)
    return kernel @ summary_density(summary, incoming) * 2 * math.pi / nodes


class CompensatedAngleKernelTests(unittest.TestCase):
    def test_rational_certificate_is_strictly_positive(self):
        self.assertGreater(rational_certificate(), 0)
        self.assertGreater(positivity_bracket(0.059), 0)

    def test_compensation_has_no_mass_or_first_harmonic(self):
        grid = np.arange(4096) * 2 * math.pi / 4096
        for eta in (0, 0.03, 0.059):
            values = compensation(grid, eta)
            self.assertAlmostEqual(values.mean(), 0, places=14)
            self.assertAlmostEqual(np.mean(values * np.cos(grid)), 0, places=14)
            self.assertAlmostEqual(np.mean(values * np.sin(grid)), 0, places=14)

    def test_angular_mixture_matches_independent_poisson_convolution(self):
        grid = np.arange(4096) * 2 * math.pi / 4096
        output_points = np.linspace(-0.5, 6.5, 19)
        eta = 0.059
        convolution = poisson_kernel(output_points[:, None] - grid[None, :]) @ preconvolution_density(grid, eta) / len(grid)
        np.testing.assert_allclose(convolution, 1 + compensation(output_points, eta), atol=3e-12)

    def test_angle_warp_pushforward_has_the_claimed_density(self):
        angles = np.linspace(-2.9, 2.9, 101)
        eta = 0.059
        q, _, _ = parameters(eta)
        total = np.zeros_like(angles)
        for sign in (-1, 1):
            inverse_cosine = (np.cos(angles) - sign * eta) / (1 - sign * eta * np.cos(angles))
            source_weight = (1 + sign * eta * inverse_cosine) / 2
            jacobian_inverse = (1 + sign * eta * inverse_cosine) / q
            total += source_weight * jacobian_inverse
        np.testing.assert_allclose(total, preconvolution_density(angles, eta), atol=1e-14)

    def test_exact_compensation_maximum(self):
        grid = np.linspace(-math.pi, math.pi, 1001)
        for eta in (0.001, 0.03, 0.059):
            self.assertLessEqual(compensation(grid, eta).max(), maximum_compensation(eta) + 1e-15)
            self.assertAlmostEqual(compensation(0, eta), maximum_compensation(eta))

    def test_uniform_positivity_certificate_for_both_outcomes(self):
        grid = np.arange(128) * 2 * math.pi / 128
        for eta in (0, 0.03, 0.059):
            lower_density = (1 - eta) * positivity_bracket(eta) / (4 * math.pi)
            for mark in (0, 1):
                values = kernel_density(grid[:, None], grid[None, :], mark, eta, 0.3)
                self.assertGreaterEqual(values.min() + 1e-14, lower_density)

    def test_instrument_preserves_total_mass_for_each_hidden_input(self):
        grid = np.arange(4096) * 2 * math.pi / 4096
        incoming = np.array([-0.7, 0.1, 0.4, 1.3, 2.2, 3.4])
        for eta in (0, 0.03, 0.059):
            mass = sum(kernel_density(grid[:, None], incoming[None, :], mark, eta, 0.4).mean(axis=0) * 2 * math.pi
                       for mark in (0, 1))
            np.testing.assert_allclose(mass, 1, atol=3e-12)

    def test_entire_selective_output_density_matches_original_fresh_update(self):
        outgoing = np.linspace(-1, 7, 17)
        for summary in ([1, 0, 0], [1, ALPHA, 0], [1, 0, ALPHA], [0.4, 0.1, -0.2]):
            for eta in (0.001, 0.03, 0.059):
                for setting in (0, 0.7):
                    for mark in (0, 1):
                        actual = integrate_density(summary, outgoing, mark, eta, setting)
                        expected = summary_density(fresh_branch_matrix(setting, mark, eta) @ summary, outgoing)
                        np.testing.assert_allclose(actual, expected, atol=3e-12)

    def test_actual_preparation_mixtures_share_identical_hidden_output(self):
        points = np.linspace(0, 6, 17)
        first, second = np.array([1, ALPHA, 0]), np.array([1, -ALPHA, 0])
        actual = (integrate_density(first, points, 1, 0.059) + integrate_density(second, points, 1, 0.059)) / 2
        np.testing.assert_allclose(actual, integrate_density([1, 0, 0], points, 1, 0.059), atol=1e-14)

    def test_zero_contrast_recovers_the_previous_poisson_update(self):
        outgoing, incoming = np.meshgrid(np.linspace(0, 6, 12), np.linspace(0, 6, 13))
        np.testing.assert_allclose(kernel_density(outgoing, incoming, 0, 0), poisson_kernel(incoming - outgoing) / (4 * math.pi), atol=1e-11)

    def test_compensation_starts_quadratically(self):
        for eta in (1e-4, 2e-4):
            self.assertAlmostEqual(maximum_compensation(eta) / eta**2, 1.5 * ALPHA**2, places=6)

    def test_numerical_certificate_root_does_not_claim_global_optimality(self):
        lower, upper = certificate_zero_bracket()
        self.assertGreater(lower, 0.059)
        self.assertLess(upper, 0.06)
        self.assertGreaterEqual(positivity_bracket(lower), -1e-16)
        self.assertLessEqual(positivity_bracket(upper), 1e-16)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CompensatedAngleKernelTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": checks.testsRun, "failures": len(checks.failures), "errors": len(checks.errors)},
               "rationally_certified_all_strength_ceiling": 0.059,
               "uniform_bracket_rational_lower_bound": float(rational_certificate()),
               "exact_formula_margin_at_0_059": positivity_bracket(0.059),
               "sufficient_certificate_zero_floating_bracket": certificate_zero_bracket(),
               "quadratic_compensation_coefficient": 1.5 * ALPHA**2,
               "scope": "All original fresh branches with contrast at most 0.059, arbitrary directions and finite adaptive protocols have a positive continuous preparation-noncontextual implementation. Exact full-density identity, not only summary closure. The certificate root is not a general impossibility boundary, and integration grids are not finite hidden-label realizations."}
    if args.write_results:
        Path(__file__).with_name("compensated_angle_kernel_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
