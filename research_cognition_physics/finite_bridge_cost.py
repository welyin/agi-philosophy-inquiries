"""Certified quantitative cost of a finite nonnegative bridge near the endpoint."""

import argparse
import json
import math
import unittest
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I, SCALE, pi_interval, sin_interval, cos_interval
from polygon_limit_analysis import factor_interval
from compensated_angle_kernel import preconvolution_density
from convex_order_interface import absolute_factor, absolute_threshold
from outcome_updates import WINDOW_ATTENUATION as ALPHA


SOURCE_DENSITY_LOWER = Fraction(623, 2000)
TARGET_DENSITY_UPPER = Fraction(8, 25)
BAND_WIDTH = Fraction(1, 10)
SEGMENT_GROWTH = Fraction(7, 5)


@lru_cache(maxsize=1)
def constants_certificate():
    pi = pi_interval()
    alpha = 4 * sin_interval(I.rational(1, 4))
    endpoint_upper = Fraction(29, 200)
    q = (1 - I.exact(endpoint_upper)**2).sqrt()
    margins = {
        "endpoint_below_0.145": alpha * factor_interval(endpoint_upper) - 1,
        "alpha_below_one": 1 - alpha,
        "pi_above_three": pi - 3,
        "asin_endpoint_below_0.15": sin_interval(I.rational(3, 20)) - endpoint_upper,
        "band_inside_source_support": alpha - BAND_WIDTH,
        "source_density_above_0.3115": q**3 / (pi * alpha) - SOURCE_DENSITY_LOWER,
        "target_density_below_0.32": TARGET_DENSITY_UPPER - 1 / (pi * (1 - I.exact(BAND_WIDTH)**2).sqrt()),
        "twenty_one_labels_cannot_contain_source_circle": alpha - cos_interval(pi / 21),
    }
    ratio = (TARGET_DENSITY_UPPER - SOURCE_DENSITY_LOWER) / SOURCE_DENSITY_LOWER
    assert ratio == Fraction(17, 623) and ratio < Fraction(1, 36)
    return {"verified": all(value.lo > 0 for value in margins.values()),
            "strict_margin_lower_bounds_exact": {name: str(Fraction(value.lo, SCALE)) for name, value in margins.items()},
            "strict_margin_diagnostics": {name: value.lo / SCALE for name, value in margins.items()},
            "source_density_lower": str(SOURCE_DENSITY_LOWER),
            "target_density_upper": str(TARGET_DENSITY_UPPER),
            "band_width": str(BAND_WIDTH), "segment_growth_upper": str(SEGMENT_GROWTH),
            "density_ratio": str(ratio), "geometric_label_minimum": 22}


def gap_lower_bound(labels):
    """Necessary eta_B-eta bound for a canonical bridge with <= labels factors."""
    if type(labels) is not int or labels < 1:
        raise ValueError("Use a positive integer factor count.")
    segments = labels // 2 + 1
    # Exact rational form of 303 / (250000 * ((7/5)^segments - 1)^2).
    return Fraction(303 * 5**(2 * segments), 250000 * (7**segments - 5**segments)**2)


def necessary_label_count(gap_upper):
    """If eta_B-eta <= gap_upper, return a necessary, not sufficient, label count."""
    gap_upper = Fraction(gap_upper)
    if gap_upper <= 0:
        raise ValueError("A finite positive gap upper bound is required.")
    labels = 22
    while gap_lower_bound(labels) > gap_upper:
        labels += 1
    return labels


def source_projection_density(x, contrast):
    """Density of Y_x for |x|<alpha; floating values are only diagnostics."""
    x = np.asarray(x)
    if np.any(np.abs(x) >= ALPHA):
        raise ValueError("Evaluate strictly inside the source interval.")
    return preconvolution_density(np.arccos(x / ALPHA), contrast) / (math.pi * np.sqrt(ALPHA**2 - x**2))


def source_call(threshold, contrast):
    """Independent angular quadrature for E[(Y_x-threshold)_+], threshold>=0."""
    if threshold < 0:
        raise ValueError("Use the nonnegative half-axis.")
    if threshold >= ALPHA:
        return 0.0
    nodes, weights = np.polynomial.legendre.leggauss(64)
    end = math.acos(threshold / ALPHA)
    angles = (nodes + 1) * end / 2
    return float(np.sum(weights * (ALPHA * np.cos(angles) - threshold)
                        * preconvolution_density(angles, contrast)) * end / (2 * math.pi))


def target_call(threshold):
    if threshold < 0:
        raise ValueError("Use the nonnegative half-axis.")
    if threshold >= 1:
        return 0.0
    return (math.sqrt(1 - threshold**2) - threshold * math.acos(threshold)) / math.pi


class FiniteBridgeCostTests(unittest.TestCase):
    def test_all_uniform_constants_have_integer_certificates(self):
        self.assertTrue(constants_certificate()["verified"])

    def test_projected_density_matches_independent_angular_band_integrals(self):
        nodes, weights = np.polynomial.legendre.leggauss(48)
        width = float(BAND_WIDTH)
        x = (nodes + 1) * width / 2
        end = math.asin(width / ALPHA)
        angle = math.pi / 2 - (nodes + 1) * end / 2
        for eta in (0, .1, .14473104667, .145):
            density_mass = float(weights @ source_projection_density(x, eta) * width / 2)
            angle_mass = float(weights @ preconvolution_density(angle, eta) * end / (2 * math.pi))
            self.assertAlmostEqual(density_mass, angle_mass, places=14)
            self.assertGreater(float(source_projection_density(x, eta).min()), float(SOURCE_DENSITY_LOWER))

    def test_hinge_integrals_match_projection_moment_and_uniform_arc(self):
        for eta in (0, .1, .14473104667):
            self.assertAlmostEqual(source_call(0, eta), ALPHA * absolute_factor(eta) / math.pi, places=14)
        nodes, weights = np.polynomial.legendre.leggauss(64)
        for threshold in (0, .05, .5, .99):
            end = math.acos(threshold)
            angles = (nodes + 1) * end / 2
            integrated = float(weights @ (np.cos(angles) - threshold) * end / (2 * math.pi))
            self.assertAlmostEqual(target_call(threshold), integrated, places=14)

    def test_saved_feasible_arc_bridges_obey_projected_call_sandwich(self):
        for size, eta in ((32, 0), (2048, .14473104667)):
            points = np.sinc(1 / size) * np.cos(np.arange(size) * 2 * math.pi / size)
            for threshold in (0, .001, .003, .01, .05, .1, .5, .9, .99):
                middle = float(np.maximum(points - threshold, 0).mean())
                self.assertLessEqual(source_call(threshold, eta), middle + 2e-15)
                self.assertLessEqual(middle, target_call(threshold) + 2e-15)

    def test_empty_segments_of_known_bridges_satisfy_growth_bound(self):
        for size, eta in ((32, 0), (2048, .14473104667)):
            points = np.sinc(1 / size) * np.cos(np.arange(size) * 2 * math.pi / size)
            knots = sorted(set(float(x) for x in points if 1e-12 < x < float(BAND_WIDTH)))
            endpoints = [0] + knots + [float(BAND_WIDTH)]
            deficit = (1 - ALPHA * absolute_factor(eta)) / math.pi
            additive = math.sqrt(8 * deficit / float(2 * SOURCE_DENSITY_LOWER - TARGET_DENSITY_UPPER))
            for left, right in zip(endpoints, endpoints[1:]):
                self.assertLessEqual(right, float(SEGMENT_GROWTH) * left + additive + 1e-13)

    def test_x_projection_all_hinges_pass_at_endpoint_but_origin_fails_above(self):
        endpoint = sum(absolute_threshold()) / 2
        for eta in (0, .1, endpoint):
            for threshold in np.linspace(0, 1, 41):
                self.assertLessEqual(source_call(threshold, eta), target_call(threshold) + 3e-15)
        self.assertGreater(source_call(0, .145), target_call(0))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(FiniteBridgeCostTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    requested_gaps = [Fraction(1, 10**exponent) for exponent in (6, 9, 12, 18, 24, 60)]
    data = {"constants_certificate": constants_certificate(),
            "necessary_counts": [{"endpoint_gap_at_most": str(gap),
                                   "necessary_labels_at_least": necessary_label_count(gap)} for gap in requested_gaps],
            "fixed_count_gap_bounds": [{"labels": count, "endpoint_gap_lower_bound_exact": str(gap_lower_bound(count)),
                                        "endpoint_gap_lower_bound_float": float(gap_lower_bound(count))}
                                       for count in (22, 32, 64, 128)],
            "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
            "scope": "Necessary bounds for finite nonnegative product decompositions of the canonical circle coupling. Neither sufficient label counts nor a general nonendpoint encoding bound; the logarithmic lower bound is not claimed sharp."}
    if args.write_results:
        Path(__file__).with_name("finite_bridge_cost_results.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(data, indent=2))


if __name__ == "__main__":
    main()
