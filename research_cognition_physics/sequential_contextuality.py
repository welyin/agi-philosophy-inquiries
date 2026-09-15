"""Sequential fresh reads can obstruct a valid single-read affine model.

Count recursion is exact for a declared finite read block. It neither assumes
independent repeated observations nor postselects a favorable output branch.
"""

import argparse
import json
import math
import platform
import unittest
from itertools import combinations, product
from fractions import Fraction
from pathlib import Path

import numpy as np

from affine_preparation_threshold import THRESHOLD, four_label_encoding, four_label_response
from outcome_updates import WINDOW_ATTENUATION as ALPHA
from preparation_contextuality import SIGNS, square_preparations
from retention_tradeoff import fresh_branch_matrix


def count_transfers(reads, contrast, setting=0.0):
    if not isinstance(reads, int) or reads < 1 or not 0 <= contrast <= 1:
        raise ValueError("Use a positive integer read count and contrast in [0,1].")
    minus, plus = [fresh_branch_matrix(setting, outcome, contrast) for outcome in (0, 1)]
    transfers = np.eye(3)[None, :, :]
    for _ in range(reads):
        updated = np.zeros((len(transfers) + 1, 3, 3))
        updated[:-1] += minus @ transfers
        updated[1:] += plus @ transfers
        transfers = updated
    return transfers


def majority_bias_row(reads, contrast, setting=0.0):
    """Ties are guessed by an independent fair coin, so their bias is zero."""
    decisions = np.sign(2 * np.arange(reads + 1) - reads)
    return decisions @ count_transfers(reads, contrast, setting)[:, 0, :]


def block_success(reads, contrast):
    return float((1 + ALPHA * majority_bias_row(reads, contrast)[1] / math.sqrt(2)) / 2)


def three_read_gain(contrast, alpha=ALPHA):
    return contrast * (1 + alpha + alpha**2 - alpha * contrast**2) / 2


def odd_product_gain(indices, contrast, alpha=ALPHA):
    """Indices start at zero: bit flips occur after each emission."""
    if not indices or len(indices) % 2 == 0:
        raise ValueError("This formula is for nonempty odd products.")
    exponent = indices[0] + sum(indices[j + 1] - indices[j] for j in range(1, len(indices), 2))
    return contrast**len(indices) * alpha**exponent


def five_read_gain(contrast, alpha=ALPHA):
    singles = sum(odd_product_gain((index,), contrast, alpha) for index in range(5))
    triples = sum(odd_product_gain(indices, contrast, alpha) for indices in combinations(range(5), 3))
    quintuple = odd_product_gain(tuple(range(5)), contrast, alpha)
    return (3 * singles - triples + 3 * quintuple) / 8


def enumerate_word_transfers(reads, contrast, setting=0):
    matrices = [fresh_branch_matrix(setting, mark, contrast) for mark in (0, 1)]
    words = []
    for marks in product((0, 1), repeat=reads):
        transfer = np.eye(3)
        for mark in marks:
            transfer = matrices[mark] @ transfer
        words.append((marks, transfer))
    return words


def diagnostics():
    rows = []
    for contrast in (0.1, 0.2, 0.3, 0.5):
        for reads in (1, 3, 5, 11, 31, 101, 301):
            gain = float(majority_bias_row(reads, contrast)[1])
            score = (1 + ALPHA * gain / math.sqrt(2)) / 2
            rows.append({"contrast": contrast, "reads": reads, "majority_effective_contrast": gain,
                         "effective_visibility": ALPHA * gain, "four_preparation_success": score,
                         "violates_four_preparation_bound": bool(score > 0.75),
                         "violates_all_disk_visibility_bound": bool(ALPHA * gain > THRESHOLD)})
    three = three_read_gain(0.5)
    five = five_read_gain(0.5)
    return {"preparation_radius": ALPHA, "selected_contrast": 0.5,
            "single_read_visibility": 0.5 * ALPHA, "four_label_single_read_model_is_available": bool(0.5 * ALPHA <= 0.5),
            "three_read_effective_contrast": three, "three_read_effective_visibility": ALPHA * three,
            "five_read_effective_contrast": five, "five_read_four_preparation_success": block_success(5, 0.5),
            "five_read_hidden_parity_tv_lower_bound": max(0.0, 4 * block_success(5, 0.5) - 3),
            "declared_majority_protocol_scan_not_global_optimization": rows}


class SequentialContextualityTests(unittest.TestCase):
    def test_five_read_violation_has_a_rational_margin_certificate(self):
        lower_alpha = Fraction(49, 50)
        lower_gain = (Fraction(3, 2) * sum(lower_alpha**index for index in range(5))
                      - Fraction(5, 4) + Fraction(3, 32) * lower_alpha**2) / 8
        self.assertGreater(Fraction(95, 96), lower_alpha)
        self.assertGreater(lower_gain, Fraction(3, 4))
        self.assertGreater((lower_alpha * Fraction(3, 4))**2, Fraction(1, 2))

    def test_count_recursion_matches_complete_word_enumeration(self):
        reads, contrast = 5, 0.5
        enumerated = np.zeros((reads + 1, 3, 3))
        for marks, transfer in enumerate_word_transfers(reads, contrast, 0.4):
            enumerated[sum(marks)] += transfer
        np.testing.assert_allclose(count_transfers(reads, contrast, 0.4), enumerated, atol=1e-14)

    def test_count_distributions_are_positive_and_normalized(self):
        for contrast in (0, 0.2, 0.5, 1):
            effects = count_transfers(31, contrast)[:, 0, :]
            for state in square_preparations():
                probabilities = effects @ state
                self.assertGreaterEqual(probabilities.min(), -1e-15)
                self.assertAlmostEqual(probabilities.sum(), 1)

    def test_total_transfer_is_nonselective_channel_power(self):
        contrast, reads = 0.5, 11
        channel = sum(fresh_branch_matrix(0, mark, contrast) for mark in (0, 1))
        np.testing.assert_allclose(count_transfers(reads, contrast).sum(axis=0), np.linalg.matrix_power(channel, reads), atol=1e-14)

    def test_majority_bias_has_no_offset_or_transverse_component(self):
        for reads in (1, 3, 5, 11, 30):
            row = majority_bias_row(reads, 0.5)
            self.assertAlmostEqual(row[0], 0)
            self.assertAlmostEqual(row[2], 0)
            self.assertLessEqual(abs(row[1]), 1)

    def test_three_read_gain_matches_independent_odd_moment_formula(self):
        for contrast in (0, 0.2, 0.5, 1):
            self.assertAlmostEqual(majority_bias_row(3, contrast)[1], three_read_gain(contrast))

    def test_five_read_gain_matches_independent_odd_moment_formula(self):
        for contrast in (0, 0.2, 0.5, 1):
            self.assertAlmostEqual(majority_bias_row(5, contrast)[1], five_read_gain(contrast))

    def test_majority_five_boolean_fourier_identity(self):
        for signs in product((-1, 1), repeat=5):
            polynomial = (3 * sum(signs) - sum(math.prod(signs[i] for i in indices) for indices in combinations(range(5), 3))
                          + 3 * math.prod(signs)) / 8
            self.assertEqual(polynomial, np.sign(sum(signs)))

    def test_four_label_single_read_model_exists_at_selected_contrast(self):
        for state in square_preparations():
            encoding = four_label_encoding(state[1:] / ALPHA)
            for axis in (0, math.pi / 2):
                probability = four_label_response(axis, ALPHA * 0.5) @ encoding
                self.assertAlmostEqual(probability, fresh_branch_matrix(axis, 1, 0.5)[0] @ state)

    def test_three_reads_cross_general_disk_threshold(self):
        self.assertLess(0.5 * ALPHA, THRESHOLD)
        self.assertGreater(three_read_gain(0.5) * ALPHA, THRESHOLD)

    def test_five_reads_violate_finite_four_preparation_bound(self):
        self.assertLess(block_success(1, 0.5), 0.75)
        self.assertGreater(block_success(5, 0.5), 0.75)

    def test_rotated_block_yields_same_success_on_both_requested_bits(self):
        success = []
        for state, signs in zip(square_preparations(), SIGNS):
            for axis in (0, 1):
                row = majority_bias_row(5, 0.5, axis * math.pi / 2)
                success.append((1 + signs[axis] * row @ state) / 2)
        np.testing.assert_allclose(success, block_success(5, 0.5), atol=1e-14)

    def test_parity_mixtures_agree_for_each_complete_word(self):
        parity = np.prod(SIGNS, axis=1)
        states = square_preparations()
        for axis in (0, math.pi / 2):
            for _, transfer in enumerate_word_transfers(5, 0.5, axis):
                probabilities = states @ transfer[0]
                self.assertAlmostEqual(probabilities[parity == 1].mean(), probabilities[parity == -1].mean())

    def test_more_reads_need_not_improve_a_majority_summary(self):
        self.assertGreater(block_success(31, 0.5), block_success(301, 0.5))

    def test_complete_records_preserve_earlier_information_by_marginalization(self):
        prefix_effects = {marks: transfer[0] for marks, transfer in enumerate_word_transfers(3, 0.5)}
        marginalized = {marks: np.zeros(3) for marks in prefix_effects}
        for marks, transfer in enumerate_word_transfers(5, 0.5):
            marginalized[marks[:3]] += transfer[0]
        for marks in prefix_effects:
            np.testing.assert_allclose(marginalized[marks], prefix_effects[marks], atol=1e-14)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SequentialContextualityTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": checks.testsRun, "failures": len(checks.failures), "errors": len(checks.errors)},
               **diagnostics(),
               "scope": "At contrast 0.5 a four-label preparation-noncontextual model reproduces all single reads, but the original fresh update gives three-read majority visibility > 2/pi and a five-read four-preparation score > 3/4. This excludes any preparation-noncontextual extension to those sequential protocols, not ordinary contextual classical implementations. Scans are for majority postprocessing only and are not optimal read strategies."}
    if args.write_results:
        Path(__file__).with_name("sequential_contextuality_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
