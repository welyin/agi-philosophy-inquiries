"""Identify noise and signal from disjoint record blocks, not hidden states."""

import argparse
import json
import math
import platform
import unittest
from pathlib import Path

import numpy as np

from shrinking_preparations import branch_matrix, parameters


def decay_integral(rate, duration):
    if rate < 0 or duration < 0 or not math.isfinite(rate + duration):
        raise ValueError("Finite nonnegative rate and duration are required.")
    return -math.expm1(-rate * duration) / rate if rate else duration


def block_cross_moment(length=1.0, gap=1.0, noise=1.0, signal=1.0):
    """Raw second moment of two disjoint record increments of equal length."""
    if signal <= 0 or not math.isfinite(signal):
        raise ValueError("Signal must be finite and positive.")
    rate = noise / 6
    integral = decay_integral(rate, length)
    decay_integral(rate, gap)  # Validate gap.
    return signal**2 * math.exp(-rate * gap) * integral**2


def finite_block_cross_moment(count, block_steps, gap_steps, noise=1.0, signal=1.0):
    if block_steps < 1 or gap_steps < 0:
        raise ValueError("Positive block size and nonnegative gap are required.")
    _, alpha, eta = parameters(count, noise, signal)
    geometric_sum = (-math.expm1(block_steps * math.log(alpha)) / (1 - alpha)
                     if alpha < 1 else float(block_steps))
    return eta**2 / count * alpha**(gap_steps + 1) * geometric_sum**2


def infer_parameters(adjacent_moment, separated_moment, length=1.0, gap=1.0):
    """Population-moment identification, not error-free estimation from one run."""
    if not (0 < separated_moment <= adjacent_moment and length > 0 and gap > 0):
        raise ValueError("Positive moments with nonincreasing correlation and positive lengths are required.")
    rate = -math.log(separated_moment / adjacent_moment) / gap
    signal_squared = adjacent_moment / decay_integral(rate, length)**2
    return {"noise": 6 * rate, "signal": math.sqrt(signal_squared),
            "time_unit_invariant_noise_to_signal_squared": 6 * rate / signal_squared}


def block_tree(count, block_steps, gap_steps, initial, noise=1.0, signal=1.0):
    total_steps = 2 * block_steps + gap_steps
    if total_steps > 16:
        raise ValueError("Full record enumeration is limited to 16 steps.")
    states = np.asarray(initial, dtype=float)[None, :]
    first = second = np.zeros(1)
    for index in range(total_steps):
        states = np.concatenate([states @ branch_matrix(mark, count, noise, signal).T for mark in (0, 1)])
        first_increment = 1 if index < block_steps else 0
        second_increment = 1 if index >= block_steps + gap_steps else 0
        first = np.concatenate((first - first_increment, first + first_increment))
        second = np.concatenate((second - second_increment, second + second_increment))
    return states[:, 0], first / math.sqrt(count), second / math.sqrt(count)


def population_rows():
    rows = []
    for noise in (0, 1, 3, 6):
        adjacent = block_cross_moment(gap=0, noise=noise)
        separated = block_cross_moment(gap=1, noise=noise)
        rows.append({"noise": noise, "signal": 1,
                     "adjacent_block_raw_moment": adjacent,
                     "one_unit_gap_raw_moment": separated,
                     "inferred_from_population_moments": infer_parameters(adjacent, separated)})
    return rows


class RecordIdentifiabilityTests(unittest.TestCase):
    def test_nonoverlapping_block_formula_matches_full_tree(self):
        for initial in ([1, 0, 0], [1, 0.4, 0.3], [1, -0.7, 0.1]):
            for block_steps, gap_steps in ((2, 0), (3, 2), (4, 3)):
                probabilities, first, second = block_tree(16, block_steps, gap_steps, initial)
                self.assertAlmostEqual(probabilities.sum(), 1)
                actual = probabilities @ (first * second)
                self.assertAlmostEqual(actual, finite_block_cross_moment(16, block_steps, gap_steps))

    def test_continuous_block_formula_matches_independent_double_quadrature(self):
        nodes, weights = np.polynomial.legendre.leggauss(32)
        length, gap, noise, signal = 0.7, 0.4, 3, 1.3
        first = length * (nodes + 1) / 2
        second = length + gap + length * (nodes + 1) / 2
        kernel = signal**2 * np.exp(-noise / 6 * (second[None, :] - first[:, None]))
        integrated = length**2 / 4 * weights @ kernel @ weights
        self.assertAlmostEqual(integrated, block_cross_moment(length, gap, noise, signal))

    def test_discrete_formula_converges_to_continuous_record_statistic(self):
        target = block_cross_moment()
        errors = [abs(finite_block_cross_moment(n, n, n) - target) for n in (20, 100, 1000)]
        self.assertTrue(errors[0] > errors[1] > errors[2])
        self.assertLess(errors[-1], 0.0001)

    def test_two_lag_moments_identify_both_parameters(self):
        for noise in (0, 1, 3, 6):
            for signal in (0.5, 1, 2):
                adjacent = block_cross_moment(gap=0, noise=noise, signal=signal)
                separated = block_cross_moment(gap=1, noise=noise, signal=signal)
                inferred = infer_parameters(adjacent, separated)
                self.assertAlmostEqual(inferred["noise"], noise)
                self.assertAlmostEqual(inferred["signal"], signal)

    def test_single_lag_cannot_identify_both_parameters(self):
        # Distinct noise values can match one block moment by changing signal.
        target = block_cross_moment(noise=1)
        other_signal = math.sqrt(target / block_cross_moment(noise=6))
        self.assertAlmostEqual(block_cross_moment(noise=6, signal=other_signal), target)
        self.assertNotAlmostEqual(block_cross_moment(gap=2, noise=6, signal=other_signal),
                                  block_cross_moment(gap=2, noise=1))

    def test_time_relabeling_preserves_a_dimensionless_parameter(self):
        noise, signal, clock_factor = 3, 1.4, 7
        transformed_noise, transformed_signal = noise / clock_factor, signal / math.sqrt(clock_factor)
        self.assertAlmostEqual(noise / signal**2, transformed_noise / transformed_signal**2)
        original = block_cross_moment(0.8, 0.6, noise, signal)
        transformed = block_cross_moment(0.8 * clock_factor, 0.6 * clock_factor,
                                         transformed_noise, transformed_signal)
        self.assertAlmostEqual(transformed, clock_factor * original)

    def test_fixed_axis_records_cannot_identify_transverse_initial_state(self):
        first, _, _ = block_tree(16, 4, 2, [1, 0.3, 0.6])
        second, _, _ = block_tree(16, 4, 2, [1, 0.3, -0.6])
        np.testing.assert_array_equal(first, second)

    def test_unphysical_population_moments_are_rejected(self):
        for moments in ((0, 0), (1, -0.1), (1, 1.1)):
            with self.assertRaises(ValueError):
                infer_parameters(*moments)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RecordIdentifiabilityTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    finite = []
    for count in (20, 100, 1000, 10000):
        adjacent = finite_block_cross_moment(count, count, 0)
        separated = finite_block_cross_moment(count, count, count)
        finite.append({"count": count, "adjacent_raw_moment": adjacent,
                       "separated_raw_moment": separated,
                       "continuous_model_inference_from_discrete_moments": infer_parameters(adjacent, separated)})
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": tests.testsRun, "failures": len(tests.failures), "errors": len(tests.errors)},
               "population_parameter_comparison": population_rows(), "finite_step_bias": finite,
               "scope": "Fixed direction, normalized record noise and declared operation-time unit. Two population raw block moments identify nu and positive sigma. This is structural identification, not a finite-sample accuracy claim. A single fixed-axis record law is independent of initial transverse summary v."}
    if args.write_results:
        Path(__file__).with_name("record_identifiability_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
