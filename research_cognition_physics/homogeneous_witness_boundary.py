"""Exact exhaustion of degree-one convex witnesses, plus endpoint reduction checks."""

import argparse
import json
import math
import platform
import unittest
from pathlib import Path

import numpy as np

from compensated_angle_kernel import preconvolution_density, warped_angle
from convex_order_interface import absolute_factor, absolute_threshold
from outcome_updates import WINDOW_ATTENUATION as ALPHA
from refined_cyclic_transport import SHIFTS, reverse_masses, resonance


def projection_factor(contrast, direction):
    if not 0 <= contrast < 1:
        raise ValueError("Use contrast in [0,1).")
    q = math.sqrt(1 - contrast**2)
    b = contrast * np.abs(np.cos(np.asarray(direction)))
    return q + b * np.arctan(b / q)


def projection_gap(contrast, direction):
    """Target minus source expectation of |n_direction dot vector|."""
    return 2 / math.pi * (1 - ALPHA * projection_factor(contrast, direction))


def circle_expectations(contrast, function, points=65536):
    angles = (np.arange(points) + 0.5) * 2 * math.pi / points
    vectors = np.stack((np.cos(angles), np.sin(angles)), axis=-1)
    source = np.mean(preconvolution_density(angles, contrast) * function(ALPHA * vectors))
    target = np.mean(function(vectors))
    return float(source), float(target)


def halfcircle_summary(contrast, points=65536):
    angles = -math.pi / 2 + (np.arange(points) + 0.5) * math.pi / points
    weights = preconvolution_density(angles, contrast)
    source_mean = ALPHA * np.array([np.mean(weights * np.cos(angles)), np.mean(weights * np.sin(angles))])
    target_mean = np.array([np.mean(np.cos(angles)), np.mean(np.sin(angles))])
    return {"source_mass": float(np.mean(weights)), "target_mass": 1.0,
            "source_mean": source_mean.tolist(), "target_mean": target_mean.tolist()}


def crossing_budget(contrast):
    """E[|X_x| 1_{X_x Y_x < 0}] for any admissible canonical coupling."""
    return (1 - ALPHA * absolute_factor(contrast)) / math.pi


def shifted_projection_scan(contrast, points=8192):
    """Finite exploratory scan, explicitly not an all-convex sufficiency test."""
    angles = (np.arange(points) + 0.5) * 2 * math.pi / points
    weights = preconvolution_density(angles, contrast)
    offsets = np.linspace(0.0, 0.98, 50)
    smallest, location = math.inf, None
    for direction in np.linspace(0, math.pi / 2, 33):
        values = np.abs(np.cos(angles - direction))
        gaps = np.mean(np.maximum(values[:, None] - offsets, 0)
                       - weights[:, None] * np.maximum(ALPHA * values[:, None] - offsets, 0), axis=0)
        index = int(np.argmin(gaps))
        if gaps[index] < smallest:
            smallest = float(gaps[index])
            location = {"direction": float(direction), "offset": float(offsets[index])}
    # Both laws centrally symmetric: a one-sided hinge gap is half this gap.
    return {"contrast": contrast, "directions": 33, "offsets": 50, "quadrature_points": points,
            "smallest_sampled_gap": smallest, "location": location,
            "status": "finite scan only; near-zero quadrature values have no certified sign"}


class HomogeneousWitnessBoundaryTests(unittest.TestCase):
    def test_directional_formula_matches_independent_circle_integrals(self):
        for eta in (0, 0.1182, 0.144, 0.145):
            for theta in (0, 0.3, 0.7, math.pi / 2):
                direction = np.array([math.cos(theta), math.sin(theta)])
                source, target = circle_expectations(eta, lambda v: np.abs(v @ direction))
                self.assertAlmostEqual(source, 2 * ALPHA * float(projection_factor(eta, theta)) / math.pi, places=8)
                self.assertAlmostEqual(target, 2 / math.pi, places=8)

    def test_directional_formula_matches_original_output_label_pushforward(self):
        angles = (np.arange(65536) + 0.5) * 2 * math.pi / 65536
        for theta in (0.2, 0.8):
            value = 0.0
            for outcome in (0, 1):
                f = (1 + (2 * outcome - 1) * 0.144 * np.cos(angles)) / 2
                value += np.mean(f * ALPHA * np.abs(np.cos(warped_angle(angles, outcome, 0.144) - theta)))
            self.assertAlmostEqual(float(value), 2 * ALPHA * float(projection_factor(0.144, theta)) / math.pi, places=8)

    def test_read_axis_is_the_global_maximum_and_recovers_the_old_bound(self):
        for eta in (0.02, 0.1182, 0.144):
            directions = np.linspace(-math.pi, math.pi, 301)
            self.assertLessEqual(float(np.max(projection_factor(eta, directions))), absolute_factor(eta) + 1e-14)
            self.assertAlmostEqual(float(projection_factor(eta, 0)), absolute_factor(eta), places=14)
            self.assertAlmostEqual(float(projection_factor(eta, math.pi / 2)), math.sqrt(1 - eta**2), places=14)

    def test_positive_sums_of_absolute_projections_cannot_improve_bound(self):
        directions = np.array([0.0, 0.2, 0.7, 1.4, 2.0])
        coefficients = np.array([0.3, 0.2, 0.7, 0.1, 0.8])
        vectors = np.stack((np.cos(directions), np.sin(directions)), axis=-1)
        source, target = circle_expectations(0.144, lambda v: np.abs(v @ vectors.T) @ coefficients)
        exact_gap = float(projection_gap(0.144, directions) @ coefficients)
        self.assertAlmostEqual(target - source, exact_gap, places=8)
        self.assertGreater(exact_gap, 0)

    def test_symmetrizing_a_non_even_support_function_keeps_expectations(self):
        slopes = np.array([[1, 0], [0.2, 1.3], [-0.7, -0.4]])
        phi = lambda v: np.max(v @ slopes.T, axis=1)
        even_phi = lambda v: (phi(v) + phi(-v)) / 2
        np.testing.assert_allclose(circle_expectations(0.144, phi), circle_expectations(0.144, even_phi), atol=1e-13)

    def test_endpoint_equality_is_confined_to_read_axis_projections(self):
        eta = sum(absolute_threshold()) / 2
        self.assertAlmostEqual(float(projection_gap(eta, 0)), 0, places=13)
        for theta in (0.01, 0.1, 0.7, math.pi / 2):
            self.assertGreater(float(projection_gap(eta, theta)), 0)

    def test_endpoint_halfcircle_means_match_but_below_endpoint_do_not(self):
        eta = sum(absolute_threshold()) / 2
        summary = halfcircle_summary(eta)
        self.assertAlmostEqual(summary["source_mass"], 1, places=13)
        np.testing.assert_allclose(summary["source_mean"], summary["target_mean"], atol=2e-10)
        below = halfcircle_summary(0.14)
        self.assertLess(below["source_mean"][0], below["target_mean"][0])

    def test_homogeneous_order_does_not_imply_full_convex_order_control(self):
        # Unrelated control pair: mu=1/4(delta_-2+delta_2)+1/2 delta_0,
        # nu=1/2(delta_-1+delta_1), on the horizontal line.
        mu = np.array([[-2, 0], [0, 0], [2, 0]])
        nu = np.array([[-1, 0], [1, 0]])
        weights = np.array([0.25, 0.5, 0.25])
        for theta in (0, 0.4, 1.2):
            direction = np.array([math.cos(theta), math.sin(theta)])
            self.assertAlmostEqual(float(weights @ np.abs(mu @ direction)), float(np.mean(np.abs(nu @ direction))))
        witness = lambda v: np.maximum(np.abs(v[:, 0]) - 1, 0)
        self.assertEqual(float(weights @ witness(mu)), 0.5)
        self.assertEqual(float(np.mean(witness(nu))), 0)

    def test_constructed_kernel_obeys_the_exact_weighted_sign_crossing_budget(self):
        angles = (np.arange(65536) + 0.5) * 2 * math.pi / 65536
        for eta in (0, 0.05, 0.1182):
            a, u, c, _ = reverse_masses(angles, eta)
            crossing = np.zeros_like(angles)
            for shift, mass in zip(SHIFTS, (a, u, u, c)):
                cosine = np.cos(angles + shift)
                crossing += mass * np.abs(cosine) * (cosine * np.cos(angles) < 0)
            actual = float(np.mean(crossing)) + resonance(0, eta)[1] / math.pi
            self.assertAlmostEqual(actual, crossing_budget(eta), places=8)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(HomogeneousWitnessBoundaryTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    eta = sum(absolute_threshold()) / 2
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": checks.testsRun, "failures": len(checks.failures), "errors": len(checks.errors)},
               "homogeneous_test_class_threshold_approximation": eta,
               "endpoint_direction_gaps": [{"direction": theta, "gap": float(projection_gap(eta, theta))}
                                            for theta in (0, 0.1, 0.5, 1.0, math.pi / 2)],
               "endpoint_halfcircle_summary": halfcircle_summary(eta),
               "weighted_sign_crossing_budget": [{"contrast": e, "budget": crossing_budget(e)}
                                                  for e in (0.1182, 0.14, eta)],
               "exploratory_shifted_projection_scans": [shifted_projection_scan(e) for e in (0.1182, 0.14, eta)],
               "scope": "For the fixed circle encoding, every finite convex positively degree-one homogeneous function passes exactly when alpha*G(eta)<=1. The proof uses planar even support-function representation and an exact directional integral. This is not full convex order. At the upper endpoint, arbitrary preparation-noncontextual models must project to the canonical angle encoding, and existence for that strength reduces to convex order of the positive halfcircles. Finite shifted-hinge scans are not sufficiency certificates."}
    if args.write_results:
        Path(__file__).with_name("homogeneous_witness_boundary_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
