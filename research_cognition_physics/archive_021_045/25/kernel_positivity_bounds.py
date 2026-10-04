"""Sharper global bounds for the specific Poisson-plus-linear kernel of round 24."""

import argparse
import json
import math
import platform
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from noncontextual_update_kernel import kernel_density, poisson_kernel, sufficient_contrast
from outcome_updates import WINDOW_ATTENUATION as ALPHA


def linear_threshold(alpha=ALPHA):
    return (1 - alpha**2) / (1 + 2 * alpha**2)**1.5


def improved_sufficient_contrast(alpha=ALPHA):
    linear = linear_threshold(alpha)
    correction = 2 * alpha / (1 + 2 * alpha)
    return 2 * linear / (1 + math.sqrt(1 + 4 * correction * linear))


def density_lower_bound(contrast, alpha=ALPHA):
    if not 0 <= contrast <= 1:
        raise ValueError("Use contrast in [0,1].")
    loss = contrast**2 / (1 + math.sqrt(1 - contrast**2))
    remainder = linear_threshold(alpha) - contrast - 2 * alpha / (1 + 2 * alpha) * loss
    # For a negative remainder the largest B is needed, reversing the shortcut.
    factor = abs(2 * alpha - 1) if remainder >= 0 else 1 + 2 * alpha
    return factor * remainder / (4 * math.pi)


def witness_angles(alpha=ALPHA):
    difference = math.acos(-alpha / 2)
    first = (1 + 2 * alpha) * math.cos(difference / 2)
    second = (2 * alpha - 1) * math.sin(difference / 2)
    midpoint = math.atan2(second, first) + math.pi
    return midpoint + difference / 2, midpoint - difference / 2


def witness_zero_bracket():
    outgoing, incoming = witness_angles()
    lower, upper = 0.004, 0.005
    for _ in range(48):
        middle = (lower + upper) / 2
        if kernel_density(outgoing, incoming, 1, middle) >= 0:
            lower = middle
        else:
            upper = middle
    return lower, upper


def rational_negative_certificate():
    lower_alpha = Fraction(95, 96)
    upper_alpha = lower_alpha + Fraction(1, 30720)
    eta = Fraction(41, 10000)
    upper_poisson = (1 - lower_alpha**2) / (1 + lower_alpha + lower_alpha**2)
    # Angles output=3*pi/2, input=5*pi/6. sqrt(3)>433/250 and 1-q<=eta^2.
    return upper_poisson - Fraction(433, 250) * lower_alpha * eta + upper_alpha * eta**2


class KernelPositivityBoundsTests(unittest.TestCase):
    def test_rational_negative_certificate_at_simple_fixed_angles(self):
        self.assertLess(rational_negative_certificate(), 0)
        self.assertLess(kernel_density(3 * math.pi / 2, 5 * math.pi / 6, 1, 0.0041), 0)

    def test_analytic_linear_ratio_minimum(self):
        cosines = np.linspace(-1, 1, 10001)
        ratios = (1 - ALPHA**2) / ((1 + ALPHA**2 - 2 * ALPHA * cosines)
                                                   * np.sqrt(1 + 4 * ALPHA**2 + 4 * ALPHA * cosines))
        self.assertGreaterEqual(ratios.min() + 1e-15, linear_threshold())
        c = -ALPHA / 2
        self.assertAlmostEqual((1 - ALPHA**2) / ((1 + ALPHA**2 - 2 * ALPHA * c)
                                                   * math.sqrt(1 + 4 * ALPHA**2 + 4 * ALPHA * c)), linear_threshold())

    def test_refined_certificate_improves_the_previous_one(self):
        self.assertGreater(improved_sufficient_contrast(), 2 * sufficient_contrast())
        self.assertGreater(density_lower_bound(improved_sufficient_contrast()), 0)

    def test_geometric_linear_amplitude_and_sine_product_bounds(self):
        rng = np.random.default_rng(2501)
        outgoing, incoming = rng.uniform(-math.pi, math.pi, (2, 2000))
        cosine = np.cos(outgoing - incoming)
        amplitude = np.sqrt(1 + 4 * ALPHA**2 + 4 * ALPHA * cosine)
        self.assertTrue(np.all(np.abs(np.cos(outgoing) + 2 * ALPHA * np.cos(incoming)) <= amplitude + 1e-14))
        self.assertTrue(np.all(np.sin(outgoing) * np.sin(incoming) <= (1 + cosine) / 2 + 1e-14))

    def test_global_density_bound_on_independent_angles(self):
        grid = np.linspace(0, 2 * math.pi, 151)
        for contrast in (0.001, 0.004, improved_sufficient_contrast()):
            values = kernel_density(grid[:, None], grid[None, :], 1, contrast)
            self.assertGreaterEqual(values.min() + 1e-14, density_lower_bound(contrast))

    def test_explicit_pair_really_makes_the_old_kernel_negative(self):
        outgoing, incoming = witness_angles()
        self.assertLess(kernel_density(outgoing, incoming, 1, 0.0041), -1e-6)
        self.assertGreater(kernel_density(outgoing, incoming, 1, 0.004), 0)

    def test_displayed_zero_bracket_is_only_a_point_witness(self):
        outgoing, incoming = witness_angles()
        lower, upper = witness_zero_bracket()
        self.assertGreaterEqual(kernel_density(outgoing, incoming, 1, lower), -1e-17)
        self.assertLessEqual(kernel_density(outgoing, incoming, 1, upper), 1e-17)
        self.assertGreater(lower, improved_sufficient_contrast())
        self.assertLess(upper, 0.0041)

    def test_point_witness_uses_the_exact_linear_worst_case(self):
        outgoing, incoming = witness_angles()
        self.assertAlmostEqual(math.cos(outgoing - incoming), -ALPHA / 2)
        amplitude = math.sqrt(1 + 2 * ALPHA**2)
        self.assertAlmostEqual(math.cos(outgoing) + 2 * ALPHA * math.cos(incoming), -amplitude)
        self.assertAlmostEqual(poisson_kernel(outgoing - incoming), amplitude * linear_threshold())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(KernelPositivityBoundsTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    outgoing, incoming = witness_angles()
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": checks.testsRun, "failures": len(checks.failures), "errors": len(checks.errors)},
               "old_sufficient_contrast": sufficient_contrast(), "new_sufficient_contrast": improved_sufficient_contrast(),
               "linear_part_exact_threshold": linear_threshold(),
               "point_witness": {"outgoing": outgoing, "incoming": incoming, "contrast": 0.0041,
                                 "negative_density": float(kernel_density(outgoing, incoming, 1, 0.0041)),
                                 "floating_zero_bracket": witness_zero_bracket()},
               "simple_angle_rational_upper_bound_on_4pi_density": float(rational_negative_certificate()),
               "scope": "The sufficient certificate is global. The negative point excludes this particular round-24 kernel, not all preparation-noncontextual implementations. The point zero is not asserted to be the global kernel threshold."}
    if args.write_results:
        Path(__file__).with_name("kernel_positivity_bounds_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
