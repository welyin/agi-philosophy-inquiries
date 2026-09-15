"""Checks supporting the endpoint obstruction to finite positive factorizations."""

import argparse
import json
import math
import unittest
from pathlib import Path

import numpy as np

from outcome_updates import WINDOW_ATTENUATION as ALPHA
from compensated_angle_kernel import preconvolution_density
from convex_order_interface import absolute_threshold
from polygon_martingale_certificate import quantize_source, source_weights


def positive_axis_band_mass(contrast, width):
    """Numerical integral for 0 < Y_x < width, with Y=alpha*n_psi."""
    if not 0 < width <= ALPHA:
        raise ValueError("Use 0 < width <= alpha.")
    nodes, weights = np.polynomial.legendre.leggauss(32)
    # Shift from the pole, avoiding cancellation in acos near width=0.
    angle = math.asin(width / ALPHA)
    offsets = (nodes + 1) * angle / 2
    return float(np.sum(weights * preconvolution_density(math.pi / 2 - offsets, contrast)) * angle / (2 * math.pi))


def positive_axis_band_lower_bound(contrast, width):
    """Analytic bound q^3 * width / (pi*alpha); evaluation here is diagnostic float."""
    return (1 - contrast**2)**1.5 * width / (math.pi * ALPHA)


def zero_strength_finite_bridge(size=32, radius=.995):
    """Explicit finite bridge control at eta=0, not at the critical endpoint."""
    beta = np.sinc(1 / size)
    if radius * math.cos(math.pi / size) < ALPHA or radius >= beta:
        raise ValueError("The chosen geometry does not permit this simple control.")
    weights = source_weights(size, 0, radius)
    conditional = np.full((size + 1, size), (1 - radius / beta) / size)
    conditional[np.arange(size), np.arange(size)] += radius / beta
    conditional[-1] = 1 / size
    return weights[:, None] * conditional


class FiniteBridgeObstructionTests(unittest.TestCase):
    def test_positive_band_lower_bound_matches_independent_integrals(self):
        endpoint = sum(absolute_threshold()) / 2
        for eta in (0, .14473, endpoint):
            for width in (1e-5, .001, .01, .1):
                observed = positive_axis_band_mass(eta, width)
                lower = positive_axis_band_lower_bound(eta, width)
                self.assertGreaterEqual(observed, lower * (1 - 1e-12))
                self.assertGreater(observed, 0)

    def test_four_arc_control_saturates_projection_without_full_source_support(self):
        # Each target quadrant is uniform. Its barycenter is a finite source atom.
        half_width = math.pi / 4
        centers = (np.arange(4) + .5) * math.pi / 2
        means = np.sinc(1 / 4) * np.stack((np.cos(centers), np.sin(centers)), axis=1)
        self.assertAlmostEqual(float(np.mean(np.abs(means[:, 0]))), 2 / math.pi, places=14)
        self.assertGreater(float(np.min(np.abs(means[:, 0]))), .6)
        nodes, weights = np.polynomial.legendre.leggauss(32)
        target = centers[:, None] + half_width * nodes
        integrated = np.stack((np.cos(target) @ weights, np.sin(target) @ weights), axis=1) / 2
        np.testing.assert_allclose(means, integrated, atol=2e-15)

    def test_finite_bridge_remains_possible_strictly_below_endpoint(self):
        size, radius = 32, .995
        matrix = zero_strength_finite_bridge(size, radius)
        weights = source_weights(size, 0, radius)
        angles = np.arange(size) * 2 * math.pi / size
        targets = np.sinc(1 / size) * np.stack((np.cos(angles), np.sin(angles)), axis=1)
        self.assertGreater(float(matrix.min()), 0)
        np.testing.assert_allclose(matrix.sum(axis=0), 1 / size, atol=2e-16)
        psi = np.linspace(-1, 8, 311)
        left, lower, upper, center = quantize_source(psi, size, radius)
        conditional = matrix / weights[:, None]
        probabilities = (lower[:, None] * conditional[left] + upper[:, None] * conditional[(left + 1) % size]
                         + center[:, None] * conditional[-1])
        expected = ALPHA * np.stack((np.cos(psi), np.sin(psi)), axis=1)
        np.testing.assert_allclose(probabilities @ targets, expected, atol=2e-15)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(FiniteBridgeObstructionTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    endpoint = sum(absolute_threshold()) / 2
    data = {"endpoint_band_diagnostics": [{"width": width,
                                            "source_mass_quadrature": positive_axis_band_mass(endpoint, width),
                                            "analytic_lower_bound_float": positive_axis_band_lower_bound(endpoint, width)}
                                           for width in (1e-5, .001, .01, .1)],
            "zero_strength_control_bridge_states": 32,
            "automated_checks": {"run": checks.testsRun, "errors": 0, "failures": 0},
            "scope": "Note 41 proves that an endpoint coupling, if it exists, has no finite nonnegative product decomposition, and fixed factor count cannot approach the endpoint. This is not an exclusion of infinite classical kernels or a computed convergence rate. Numerical band values are diagnostics, not certified root evaluations."}
    if args.write_results:
        Path(__file__).with_name("finite_bridge_obstruction_results.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(data, indent=2))


if __name__ == "__main__":
    main()
