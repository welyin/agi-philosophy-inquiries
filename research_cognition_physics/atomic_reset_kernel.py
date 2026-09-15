"""An exact angle-retention/reset instrument with a sharp construction threshold.

The kernel is a measure with a deterministic atom and a continuous reset part.
It implements the same full affine preparation densities as earlier rounds.
"""

import argparse
import json
import math
import platform
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from compensated_angle_kernel import preconvolution_density, warped_angle
from noncontextual_update_kernel import summary_density
from outcome_updates import WINDOW_ATTENUATION as ALPHA
from retention_tradeoff import fresh_branch_matrix


def construction_threshold(alpha=ALPHA):
    if not 0 < alpha < 1:
        raise ValueError("Use preparation attenuation strictly between zero and one.")
    square = alpha**2
    return math.sqrt(2 * (1 - square) / (1 + 2 * square + math.sqrt(1 + 8 * square)))


def atom_weight(incoming, outcome, contrast, setting=0.0):
    q = math.sqrt(1 - contrast**2)
    sign = 2 * outcome - 1
    return ALPHA * q**3 / (2 * (1 - sign * contrast * np.cos(np.asarray(incoming) - setting))**2)


def atom_target(incoming, outcome, contrast, setting=0.0):
    relative = np.asarray(incoming) - setting
    q = math.sqrt(1 - contrast**2)
    sign = 2 * outcome - 1
    return np.arctan2(q * np.sin(relative), np.cos(relative) - sign * contrast) + setting


def reset_weight(incoming, contrast, setting=0.0):
    return 1 - ALPHA * preconvolution_density(np.asarray(incoming) - setting, contrast)


def minimum_reset_weight(contrast):
    return 1 - ALPHA * (1 + contrast**2) / math.sqrt(1 - contrast**2)


def reset_output_density(outgoing, outcome, contrast, setting=0.0):
    return (1 + (2 * outcome - 1) * contrast * np.cos(np.asarray(outgoing) - setting)) / (4 * math.pi)


def rational_083_certificate():
    upper_alpha = Fraction(95, 96) + Fraction(1, 30720)
    ceiling = Fraction(83, 1000)
    return 1 - ceiling**2 - (upper_alpha * (1 + ceiling**2))**2


def pushed_output_density(summary, outgoing, outcome, contrast, setting=0.0, nodes=4096):
    """Independent atom change of variables plus reset integral."""
    outgoing = np.asarray(outgoing)
    relative = outgoing - setting
    sign = 2 * outcome - 1
    original_angle = warped_angle(relative, outcome, contrast) + setting
    jacobian = math.sqrt(1 - contrast**2) / (1 + sign * contrast * np.cos(relative))
    atom_density = atom_weight(original_angle, outcome, contrast, setting) * summary_density(summary, original_angle) * jacobian
    incoming = np.arange(nodes) * 2 * math.pi / nodes
    reset_mass = np.mean(reset_weight(incoming, contrast, setting) * summary_density(summary, incoming)) * 2 * math.pi
    return atom_density + reset_mass * reset_output_density(outgoing, outcome, contrast, setting)


def weak_output_moment(summary, outcome, contrast, setting, frequency, sine=False, nodes=4096):
    grid = np.arange(nodes) * 2 * math.pi / nodes
    function = np.sin if sine else np.cos
    reset_integral = np.mean(reset_output_density(grid, outcome, contrast, setting) * function(frequency * grid)) * 2 * math.pi
    response = (atom_weight(grid, outcome, contrast, setting) * function(frequency * atom_target(grid, outcome, contrast, setting))
                + reset_weight(grid, contrast, setting) * reset_integral)
    return float(np.mean(summary_density(summary, grid) * response) * 2 * math.pi)


class AtomicResetKernelTests(unittest.TestCase):
    def test_exact_threshold_and_rational_safe_ceiling(self):
        self.assertGreater(construction_threshold(), 0.083)
        self.assertLess(construction_threshold(), 0.084)
        self.assertGreater(rational_083_certificate(), 0)
        self.assertAlmostEqual(minimum_reset_weight(construction_threshold()), 0, places=14)

    def test_all_mode_weights_are_positive_up_to_threshold(self):
        grid = np.linspace(-math.pi, math.pi, 301)
        for eta in (0, 0.03, 0.083, construction_threshold()):
            self.assertGreaterEqual(reset_weight(grid, eta).min(), -1e-14)
            for mark in (0, 1):
                self.assertGreaterEqual(atom_weight(grid, mark, eta).min(), 0)

    def test_atom_and_reset_masses_normalize_pointwise(self):
        grid = np.linspace(0, 7, 77)
        for eta in (0, 0.03, 0.083):
            total = reset_weight(grid, eta, 0.4) + sum(atom_weight(grid, mark, eta, 0.4) for mark in (0, 1))
            np.testing.assert_allclose(total, 1, atol=1e-14)

    def test_inverse_angle_map_has_the_required_forward_direction(self):
        incoming = np.linspace(-3, 3, 71)
        for mark in (0, 1):
            target = atom_target(incoming, mark, 0.083)
            recovered = warped_angle(target, mark, 0.083)
            np.testing.assert_allclose(np.exp(1j * recovered), np.exp(1j * incoming), atol=1e-14)

    def test_reset_sampler_joint_density_normalizes(self):
        grid = np.arange(4096) * 2 * math.pi / 4096
        for eta in (0, 0.083):
            masses = [np.mean(reset_output_density(grid, mark, eta, 0.2)) * 2 * math.pi for mark in (0, 1)]
            np.testing.assert_allclose(masses, [0.5, 0.5], atol=1e-14)

    def test_entire_output_density_matches_original_fresh_instrument(self):
        outgoing = np.linspace(-2, 8, 37)
        for summary in ([1, 0, 0], [1, ALPHA, 0], [1, 0, ALPHA], [0.4, 0.1, -0.2]):
            for eta in (0, 0.03, 0.083, construction_threshold()):
                for setting in (0, 0.7):
                    for mark in (0, 1):
                        expected = summary_density(fresh_branch_matrix(setting, mark, eta) @ summary, outgoing)
                        np.testing.assert_allclose(pushed_output_density(summary, outgoing, mark, eta, setting), expected, atol=1e-14)

    def test_direct_measure_integrals_have_correct_higher_moments(self):
        summary = np.array([1, 0.3, -0.4])
        for setting in (0, 0.7):
            for mark in (0, 1):
                expected_summary = fresh_branch_matrix(setting, mark, 0.083) @ summary
                self.assertAlmostEqual(weak_output_moment(summary, mark, 0.083, setting, 0), expected_summary[0])
                self.assertAlmostEqual(weak_output_moment(summary, mark, 0.083, setting, 1), expected_summary[1] / (2 * ALPHA))
                self.assertAlmostEqual(weak_output_moment(summary, mark, 0.083, setting, 1, True), expected_summary[2] / (2 * ALPHA))
                for frequency in (2, 3, 6):
                    self.assertAlmostEqual(weak_output_moment(summary, mark, 0.083, setting, frequency), 0, places=13)

    def test_preparation_mixtures_still_have_identical_full_hidden_outputs(self):
        outgoing = np.linspace(0, 6, 53)
        first, second = np.array([1, ALPHA, 0]), np.array([1, 0, ALPHA])
        mixed = (pushed_output_density(first, outgoing, 1, 0.083) + pushed_output_density(second, outgoing, 1, 0.083)) / 2
        np.testing.assert_allclose(mixed, pushed_output_density((first + second) / 2, outgoing, 1, 0.083), atol=1e-14)

    def test_above_threshold_the_reset_part_is_a_negative_measure(self):
        self.assertLess(reset_weight(0, 0.09), 0)
        self.assertGreater(sum(atom_weight(0, mark, 0.09) for mark in (0, 1)), 1)

    def test_zero_strength_read_retains_original_disturbance(self):
        grid = np.array([0.2, 0.7, 2.3])
        np.testing.assert_allclose(reset_weight(grid, 0), 1 - ALPHA)
        for mark in (0, 1):
            np.testing.assert_allclose(atom_weight(grid, mark, 0), ALPHA / 2)
            np.testing.assert_allclose(atom_target(grid, mark, 0), grid)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(AtomicResetKernelTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": checks.testsRun, "failures": len(checks.failures), "errors": len(checks.errors)},
               "exact_threshold_for_this_kernel_formula": construction_threshold(),
               "rational_safe_ceiling": 0.083, "rational_squared_margin": float(rational_083_certificate()),
               "minimum_reset_weight_at_0_083": minimum_reset_weight(0.083),
               "minimum_reset_weight_at_0_09": minimum_reset_weight(0.09),
               "scope": "A continuous hidden angle with deterministic outcome-dependent atoms and a nonnegative reset density reproduces every full selective affine preparation density up to the explicit threshold. This threshold is exact for this formula, sufficient for general implementations, and not a universal optimum. All finite adaptive fresh-read/rotation protocols are included; no original preparation width or zero-strength disturbance is changed."}
    if args.write_results:
        Path(__file__).with_name("atomic_reset_kernel_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
