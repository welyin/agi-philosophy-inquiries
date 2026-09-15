"""Exact fixed-axis record equivalence to a binary hidden Markov model."""

import argparse
import cmath
import json
import math
import platform
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from shrinking_preparations import branch_matrix, parameters


def hidden_branch_matrix(outcome, count, noise=1.0, signal=1.0):
    """Emit from the current bit, then independently flip it for the next step."""
    _, alpha, eta = parameters(count, noise, signal)
    flip = (1 - alpha) / 2
    transition = np.array([[1 - flip, flip], [flip, 1 - flip]])
    sign = 2 * outcome - 1
    emission = np.diag([(1 + sign * eta) / 2, (1 - sign * eta) / 2])
    return transition @ emission


def word_probability(marks, initial, count, noise=1.0, signal=1.0, hidden=False):
    state = (np.array([(1 + initial[1]) / 2, (1 - initial[1]) / 2]) if hidden
             else np.asarray(initial, dtype=float).copy())
    matrix = hidden_branch_matrix if hidden else branch_matrix
    for mark in marks:
        state = matrix(mark, count, noise, signal) @ state
    return float(state.sum() if hidden else state[0])


def finite_characteristic(frequency, count, duration=1.0, noise=1.0, signal=1.0, initial_u=0.0):
    steps = round(count * duration)
    if not math.isclose(steps / count, duration):
        raise ValueError("Duration must lie on the discrete time grid.")
    transfer = sum(cmath.exp(1j * frequency * (2 * mark - 1) / math.sqrt(count))
                   * hidden_branch_matrix(mark, count, noise, signal) for mark in (0, 1))
    initial = np.array([(1 + initial_u) / 2, (1 - initial_u) / 2])
    return complex(np.sum(np.linalg.matrix_power(transfer, steps) @ initial))


def telegraph_characteristic(frequency, duration=1.0, noise=1.0, signal=1.0, initial_u=0.0):
    """Characteristic function of sigma*integral(X dt)+independent Brownian."""
    if duration < 0 or noise < 0 or signal <= 0 or abs(initial_u) > 1:
        raise ValueError("Use nonnegative duration/noise, positive signal, and |u0| <= 1.")
    rate = noise / 12  # Symmetric bit-flip rate; mean decay is twice this rate.
    delta = cmath.sqrt(rate**2 - signal**2 * frequency**2)
    ratio = cmath.sinh(delta * duration) / delta if abs(delta) > 1e-12 else duration + delta**2 * duration**3 / 6
    return cmath.exp(-(rate + frequency**2 / 2) * duration) * (
        cmath.cosh(delta * duration) + (rate + 1j * frequency * signal * initial_u) * ratio)


def direct_latent_word_probability(marks, initial_u, count, noise=1.0, signal=1.0):
    """Independent enumeration of bit paths, including hidden final flips."""
    _, alpha, eta = parameters(count, noise, signal)
    flip = (1 - alpha) / 2
    total = 0.0
    for bits in product((-1, 1), repeat=len(marks) + 1):
        probability = (1 + bits[0] * initial_u) / 2
        for index, mark in enumerate(marks):
            probability *= (1 + (2 * mark - 1) * eta * bits[index]) / 2
            probability *= (1 - flip) if bits[index + 1] == bits[index] else flip
        total += probability
    return total


def characteristic_rows():
    rows = []
    for noise in (0, 1, 6):
        for frequency in (0.5, 1.0, 2.0):
            target = telegraph_characteristic(frequency, noise=noise, initial_u=0.4)
            errors = [{"count": count,
                       "absolute_error": abs(finite_characteristic(frequency, count, noise=noise, initial_u=0.4) - target)}
                      for count in (32, 128, 512, 2048)]
            rows.append({"noise": noise, "frequency": frequency,
                         "limiting_characteristic_real_imag": [target.real, target.imag],
                         "finite_record_characteristic_errors": errors})
    return rows


class TwoStateFilterTests(unittest.TestCase):
    def test_hidden_instrument_is_a_positive_normalized_classical_model(self):
        for noise in (0, 1, 6):
            matrices = [hidden_branch_matrix(mark, 16, noise) for mark in (0, 1)]
            self.assertGreaterEqual(min(matrix.min() for matrix in matrices), 0)
            np.testing.assert_allclose(sum(matrices).sum(axis=0), [1, 1], atol=1e-14)

    def test_similarity_transformation_equals_entire_relevant_block(self):
        transform = np.array([[1, 1], [1, -1]])
        for noise in (0, 1, 6):
            for mark in (0, 1):
                hidden = transform @ hidden_branch_matrix(mark, 32, noise) @ (transform / 2)
                np.testing.assert_allclose(hidden, branch_matrix(mark, 32, noise)[:2, :2], atol=1e-14)

    def test_every_short_word_matches_for_multiple_initial_preparations(self):
        for initial in ([1, 0, 0], [1, 0.4, 0.3], [1, -0.7, 0.2]):
            for marks in product((0, 1), repeat=7):
                self.assertAlmostEqual(word_probability(marks, initial, 32),
                                       word_probability(marks, initial, 32, hidden=True))

    def test_independent_latent_path_sum_matches_matrix_recursion(self):
        for marks in product((0, 1), repeat=4):
            direct = direct_latent_word_probability(marks, 0.4, 32)
            self.assertAlmostEqual(direct, word_probability(marks, [1, 0.4, 0.3], 32))

    def test_conditional_bit_mean_matches_conditional_summary(self):
        state = np.array([1, 0.4, 0.3])
        probabilities = np.array([0.7, 0.3])
        for mark in (1, 0, 1, 1, 0, 0, 1):
            state = branch_matrix(mark, 64) @ state
            state /= state[0]
            probabilities = hidden_branch_matrix(mark, 64) @ probabilities
            probabilities /= probabilities.sum()
            self.assertAlmostEqual(state[1], probabilities[0] - probabilities[1])

    def test_hidden_flip_rate_converges_with_factor_two_accounted_for(self):
        for noise in (1, 3, 6):
            _, alpha, _ = parameters(10000, noise)
            self.assertAlmostEqual(10000 * (1 - alpha) / 2, noise / 12, places=4)

    def test_finite_characteristic_matches_full_record_enumeration(self):
        count, frequency, initial = 8, 1.3, [1, 0.4, 0.3]
        enumerated = sum(word_probability(marks, initial, count) *
                         cmath.exp(1j * frequency * sum(2 * mark - 1 for mark in marks) / math.sqrt(count))
                         for marks in product((0, 1), repeat=count))
        self.assertAlmostEqual(abs(enumerated - finite_characteristic(frequency, count, initial_u=0.4)), 0)

    def test_limiting_characteristic_matches_discrete_model(self):
        for row in characteristic_rows():
            errors = [item["absolute_error"] for item in row["finite_record_characteristic_errors"]]
            self.assertTrue(all(first > second for first, second in zip(errors, errors[1:])))
            self.assertLess(errors[-1], 0.001)

    def test_zero_flip_limit_is_the_expected_two_gaussian_mixture(self):
        for frequency in (0, 0.5, 1, 2):
            expected = math.exp(-frequency**2 / 2) * (math.cos(frequency) + 0.4j * math.sin(frequency))
            self.assertAlmostEqual(abs(telegraph_characteristic(frequency, noise=0, initial_u=0.4) - expected), 0)

    def test_characteristic_is_normalized_and_has_conjugate_symmetry(self):
        for noise in (0, 1, 6):
            self.assertAlmostEqual(abs(telegraph_characteristic(0, noise=noise) - 1), 0)
            for frequency in (noise / 12, 0.3, 1.7):
                positive = telegraph_characteristic(frequency, noise=noise, initial_u=0.4)
                negative = telegraph_characteristic(-frequency, noise=noise, initial_u=0.4)
                self.assertLessEqual(abs(positive), 1 + 1e-14)
                self.assertAlmostEqual(abs(negative - positive.conjugate()), 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(TwoStateFilterTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": tests.testsRun, "failures": len(tests.failures), "errors": len(tests.errors)},
               "continuous_hidden_flip_rate": "nu/12", "characteristic_comparisons": characteristic_rows(),
               "scope": "Exact equivalence of all finite fixed-axis binary record probabilities and u predictions. The auxiliary bit is an alternative realization, not an identified original latent variable. It does not reproduce arbitrary rotated probes or v-dependent protocols. Continuous equivalence uses the two-state Wonham filter, with no quantum postulates."}
    if args.write_results:
        Path(__file__).with_name("two_state_filter_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
