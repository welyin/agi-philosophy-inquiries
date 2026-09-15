"""Grid-specific obstruction, reserve bounds, and checks for the limiting argument."""

import argparse
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I, SCALE, pi_interval, sin_interval, cos_interval
from convex_order_interface import absolute_factor, absolute_threshold
from outcome_updates import WINDOW_ATTENUATION as ALPHA
from polygon_martingale_certificate import source_weights, quantize_source


def factor_interval(contrast, terms=32):
    """G(eta)=1+integral_0^eta asin(t) dt, with an explicit positive tail."""
    eta = I.exact(contrast)
    if eta.lo < 0 or eta.hi >= SCALE:
        raise ValueError("Use 0 <= contrast < 1.")
    result = I.exact(1)
    for n in range(1, terms + 1):
        coefficient = Fraction(math.comb(2 * n - 2, n - 1), 4**(n - 1) * (2 * n - 1) * (2 * n))
        result += I.exact(coefficient) * eta**(2 * n)
    # All omitted coefficients are positive and at most one.
    tail = eta**(2 * terms + 2) / (1 - eta**2)
    return I(result.lo, result.hi + tail.hi)


def grid_certificate(size, contrast):
    if size < 4 or size % 4:
        raise ValueError("This axis-aligned identity requires N divisible by four.")
    alpha = 4 * sin_interval(I.rational(1, 4))
    cosine = cos_interval(pi_interval() / size)
    gap = cosine - alpha * factor_interval(contrast)
    reserve = gap / cosine
    return {"size": size, "contrast": str(contrast),
            "projection_budget_interval": gap.floats(),
            "uniform_reserve_upper_bound_interval": reserve.floats(),
            "grid_excluded": gap.hi < 0,
            "budget_integer_interval": [str(gap.lo), str(gap.hi)],
            "interval_denominator": str(SCALE)}


def grid_ceiling(size):
    """Floating diagnostic root, not an independently certified feasibility claim."""
    lower, upper = 0.0, sum(absolute_threshold()) / 2
    for _ in range(60):
        middle = (lower + upper) / 2
        if ALPHA * absolute_factor(middle) <= math.cos(math.pi / size):
            lower = middle
        else:
            upper = middle
    return (lower + upper) / 2


def quantization_error_bounds(size, radius):
    h = math.pi / size
    if radius * math.cos(h) < ALPHA:
        raise ValueError("Polygon does not contain the source circle.")
    return {"center_mass_upper_bound": 1 - ALPHA / radius,
            "source_rms_displacement_upper_bound": math.sqrt(radius**2 - ALPHA**2),
            "target_rms_displacement": math.sqrt(1 - np.sinc(1 / size)**2)}


class PolygonLimitTests(unittest.TestCase):
    def test_grid_target_absolute_moment_has_exact_cosine_loss(self):
        for size in (32, 64, 128, 256, 512, 1024):
            angle = np.arange(size) * 2 * math.pi / size
            value = np.mean(np.abs(np.sinc(1 / size) * np.cos(angle)))
            self.assertAlmostEqual(value, 2 / math.pi * math.cos(math.pi / size), places=14)

    def test_source_quantization_preserves_absolute_read_axis_pointwise(self):
        size, radius = 512, 0.989635
        angles = np.linspace(-4, 8, 3001)
        left, lower, upper, center = quantize_source(angles, size, radius)
        vertices = radius * np.cos(np.arange(size) * 2 * math.pi / size)
        value = lower * np.abs(vertices[left]) + upper * np.abs(vertices[(left + 1) % size])
        np.testing.assert_allclose(value, ALPHA * np.abs(np.cos(angles)), atol=2e-15)
        self.assertGreaterEqual(float(center.min()), 0)

    def test_source_absolute_moment_matches_analytic_projection(self):
        for size, eta, radius in ((256, .144, .9897), (512, .1445, .989635), (1024, .1447, .9896206)):
            weights = source_weights(size, eta, radius)
            vertices = radius * np.cos(np.arange(size) * 2 * math.pi / size)
            self.assertAlmostEqual(float(weights[:-1] @ np.abs(vertices)),
                                   2 * ALPHA / math.pi * absolute_factor(eta), places=13)

    def test_series_interval_agrees_with_independent_elementary_function(self):
        for eta in (Fraction(0), Fraction(1, 10), Fraction(1447, 10000), Fraction(1, 2)):
            lo, hi = factor_interval(eta).floats()
            self.assertLess(abs((lo + hi) / 2 - absolute_factor(float(eta))), 5e-16)

    def test_512_grid_is_strictly_excluded_at_new_1024_strength(self):
        certificate = grid_certificate(512, Fraction(1447, 10000))
        self.assertTrue(certificate["grid_excluded"])
        self.assertLess(int(certificate["budget_integer_interval"][1]), 0)
        self.assertFalse(grid_certificate(1024, Fraction(1447, 10000))["grid_excluded"])

    def test_computed_grid_ceiling_has_inverse_square_distance(self):
        endpoint = sum(absolute_threshold()) / 2
        expected = math.pi**2 / (2 * ALPHA * math.asin(endpoint))
        for size in (512, 1024, 2048):
            coefficient = (endpoint - grid_ceiling(size)) * size**2
            self.assertLess(abs(coefficient / expected - 1), 0.001)

    def test_center_and_quantization_errors_vanish_in_controlled_sequence(self):
        previous = None
        for size in (128, 256, 512, 1024):
            radius = ALPHA * (1 + (math.pi / size)**2)
            bounds = quantization_error_bounds(size, radius)
            actual_center = source_weights(size, .144, radius)[-1]
            self.assertGreater(actual_center, 0)
            self.assertLessEqual(actual_center, bounds["center_mass_upper_bound"])
            if previous is not None:
                self.assertLess(bounds["center_mass_upper_bound"], previous["center_mass_upper_bound"] * .251)
                self.assertLess(bounds["source_rms_displacement_upper_bound"], previous["source_rms_displacement_upper_bound"] * .501)
            previous = bounds

    def test_martingale_variance_identity_for_barycentric_quantization(self):
        size, radius = 512, .989635
        angles = np.linspace(0, 2 * math.pi, 3001)
        left, lower, upper, center = quantize_source(angles, size, radius)
        points = radius * np.stack((np.cos(np.arange(size) * 2 * math.pi / size),
                                    np.sin(np.arange(size) * 2 * math.pi / size)), axis=1)
        y = ALPHA * np.stack((np.cos(angles), np.sin(angles)), axis=1)
        displacement = (lower * np.sum((points[left] - y)**2, axis=1)
                        + upper * np.sum((points[(left + 1) % size] - y)**2, axis=1)
                        + center * ALPHA**2)
        np.testing.assert_allclose(displacement, radius**2 * (1 - center) - ALPHA**2, atol=8e-16)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(PolygonLimitTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    endpoint = sum(absolute_threshold()) / 2
    report = {"grid_certificates": [grid_certificate(n, Fraction(e)) for n, e in
                                    ((256, "0.144"), (512, "0.1445"), (512, "0.1447"), (1024, "0.1447"))],
              "diagnostic_grid_ceilings": {str(n): grid_ceiling(n) for n in (128, 256, 512, 1024, 2048)},
              "asymptotic_endpoint_gap_times_size_squared": math.pi**2 / (2 * ALPHA * math.asin(endpoint)),
              "error_bounds_1024": quantization_error_bounds(1024, .9896206),
              "automated_checks": {"run": checks.testsRun, "failures": len(checks.failures), "errors": len(checks.errors)},
              "scope": "Grid exclusions apply only to this finite intermediate geometry. Compactness and martingale closure are proved in note 37, not by these numerical tests. No sequence reaching eta_B and no endpoint coupling have yet been established."}
    if args.write_results:
        Path(__file__).with_name("polygon_limit_analysis_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
