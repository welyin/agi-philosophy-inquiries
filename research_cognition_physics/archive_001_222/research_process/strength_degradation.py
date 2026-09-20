"""A positive rotation mixture transfers canonical feasibility to lower strengths."""

import argparse
import json
import math
import platform
import unittest
from pathlib import Path

import numpy as np

from compensated_angle_kernel import preconvolution_density


def parameters(contrast):
    if not 0 <= contrast < 1:
        raise ValueError("Use contrast in [0,1).")
    q = math.sqrt(1 - contrast**2)
    return q, contrast / (1 + q)


def certified_ceiling_condition(upper):
    # If T^4+4T^2<=1, all lower strengths have the positive kernel below.
    _, t = parameters(upper)
    return 1 - 4 * t**2 - t**4


def coefficients(lower, upper):
    if not 0 <= lower < upper < 1:
        raise ValueError("For densities use lower<upper; equality uses an identity atom.")
    q, t = parameters(lower)
    big_q, big_t = parameters(upper)
    return q / big_q, (t / big_t)**2, big_q


def poisson(angles, radius):
    return (1 - radius**2) / (1 - 2 * radius * np.cos(angles) + radius**2)


def mixing_density(angles, lower, upper, nodes=256):
    """Density of a rotation angle; exact definition uses the integral over u."""
    ratio, radius, big_q = coefficients(lower, upper)
    abscissas, weights = np.polynomial.legendre.leggauss(nodes)
    u = (abscissas + 1) / 2
    smaller = radius * u**(2 * big_q)
    angle = 2 * np.asarray(angles)
    integrated = np.sum(poisson(angle[..., None], smaller) * (weights / 2), axis=-1)
    return (ratio * poisson(angle, radius) - (ratio - 1) * integrated) / (2 * math.pi)


def positive_relative_bound(lower, upper):
    _, t = parameters(lower)
    _, big_t = parameters(upper)
    return 1 - 4 * t**2 / ((1 + t**2) * (1 - big_t**2))


def fourier_multiplier(k, lower, upper):
    if lower == upper:
        return 1.0
    if upper == 0:
        raise ValueError("Only the equal zero case is allowed.")
    q, t = parameters(lower)
    big_q, big_t = parameters(upper)
    return (t / big_t)**(2 * k) * (1 + 2 * k * q) / (1 + 2 * k * big_q)


class StrengthDegradationTests(unittest.TestCase):
    def test_current_research_range_satisfies_global_positive_condition(self):
        self.assertGreater(certified_ceiling_condition(0.145), 0.97)
        threshold = math.sqrt((math.sqrt(5) - 1) / 2)
        self.assertAlmostEqual(certified_ceiling_condition(threshold), 0, places=13)

    def test_poisson_comparison_bound(self):
        grid = np.linspace(-math.pi, math.pi, 301)
        for r in (0.1, 0.8, 0.99):
            for s in (0, r / 2, r):
                ratio = poisson(grid, s) / poisson(grid, r)
                self.assertLessEqual(float(ratio.max()), (1 + r) / (1 - r) + 1e-10)

    def test_density_is_positive_and_normalized(self):
        grid = (np.arange(4096) + 0.5) * 2 * math.pi / 4096
        for lower in (0, 0.05, 0.12, 0.143):
            value = mixing_density(grid, lower, 0.144)
            self.assertGreater(float(value.min()), 0)
            self.assertAlmostEqual(float(value.mean() * 2 * math.pi), 1, places=10)
            _, radius, _ = coefficients(lower, 0.144)
            bound = positive_relative_bound(lower, 0.144) * poisson(2 * grid, radius) / (2 * math.pi)
            self.assertGreaterEqual(float(np.min(value - bound)), -1e-11)

    def test_fourier_multipliers_match_independent_integrals(self):
        grid = (np.arange(4096) + 0.5) * 2 * math.pi / 4096
        density = mixing_density(grid, 0.12, 0.144)
        for k in range(1, 7):
            actual = np.mean(density * np.cos(2 * k * grid)) * 2 * math.pi
            self.assertAlmostEqual(float(actual), fourier_multiplier(k, 0.12, 0.144), places=11)
        self.assertAlmostEqual(float(np.mean(density * np.cos(grid)) * 2 * math.pi), 0, places=12)

    def test_convolution_recovers_entire_lower_strength_density(self):
        grid = (np.arange(8192) + 0.5) * 2 * math.pi / 8192
        for lower in (0, 0.05, 0.12):
            kernel = mixing_density(grid, lower, 0.144)
            for output in (0, 0.3, 0.9, 1.5):
                actual = np.mean(kernel * preconvolution_density(output - grid, 0.144)) * 2 * math.pi
                self.assertAlmostEqual(float(actual), float(preconvolution_density(output, lower)), places=10)

    def test_degradation_multiplies_transitively_in_fourier_space(self):
        for k in (1, 2, 5, 11):
            self.assertAlmostEqual(fourier_multiplier(k, 0.05, 0.12) * fourier_multiplier(k, 0.12, 0.144),
                                   fourier_multiplier(k, 0.05, 0.144), places=14)

    def test_endpoints_have_explicit_identity_and_uniform_cases(self):
        for k in (1, 2, 9):
            self.assertEqual(fourier_multiplier(k, 0.144, 0.144), 1)
            self.assertEqual(fourier_multiplier(k, 0, 0.144), 0)
        grid = np.linspace(0, 6, 31)
        np.testing.assert_allclose(mixing_density(grid, 0, 0.144), 1 / (2 * math.pi), atol=1e-14)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(StrengthDegradationTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": checks.testsRun, "failures": len(checks.failures), "errors": len(checks.errors)},
               "sufficient_upper_strength_limit": math.sqrt((math.sqrt(5) - 1) / 2),
               "ceiling_condition_margin_at_0_145": certified_ceiling_condition(0.145),
               "relative_positive_bounds_to_0_144": [{"lower": low, "bound": positive_relative_bound(low, 0.144)}
                                                   for low in (0, 0.05, 0.12, 0.143)],
               "scope": "For the canonical circle encoding and upper strength E<=sqrt((sqrt(5)-1)/2), R_eta is a positive rotation mixture of R_E for every eta<=E. Rotating any endpoint martingale coupling therefore gives exact canonical instruments for the full lower-strength family. This does not assume convex order between two different distributions on the same-radius circle, and does not assert arbitrary encodings can be degraded."}
    if args.write_results:
        Path(__file__).with_name("strength_degradation_results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
