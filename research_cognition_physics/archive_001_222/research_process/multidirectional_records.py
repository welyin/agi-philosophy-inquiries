"""Two read axes defeat a common binary encoding; predictive dimension is 3."""

import argparse
import json
import math
import platform
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from bounded_responses import translation_matrix
from outcome_updates import WINDOW_ATTENUATION as ALPHA
from retention_tradeoff import fresh_branch_matrix
from shrinking_preparations import branch_matrix, parameters, record_mean_and_tv_floor


def rotated_weak_branch(setting, outcome, count, noise=1.0, signal=1.0):
    rotation = translation_matrix(setting, (1.0,))
    return rotation @ branch_matrix(outcome, count, noise, signal) @ rotation.T


def record_distribution(initial, settings, contrast=0.5):
    """Predeclared directions on a single system, with every result retained."""
    probabilities = []
    for marks in product((0, 1), repeat=len(settings)):
        state = np.asarray(initial, dtype=float)
        for setting, mark in zip(settings, marks):
            state = fresh_branch_matrix(setting, mark, contrast) @ state
        probabilities.append(state[0])
    return np.asarray(probabilities)


def preparation_protocol_table(contrast=0.5):
    preparations = np.array([[1, ALPHA, 0], [1, -ALPHA, 0], [1, 0, ALPHA]]).T
    effects = np.stack(([1, 0, 0], fresh_branch_matrix(0, 1, contrast)[0],
                        fresh_branch_matrix(math.pi / 2, 1, contrast)[0]))
    return effects @ preparations


def two_axis_block_mean_table(count, noise=1.0, signal=1.0):
    """Rows are normalization and mean records from separate fixed-axis runs."""
    mean, _ = record_mean_and_tv_floor(count, noise, signal)
    return np.array([[1, 1, 1], [mean, -mean, 0], [0, 0, mean]])


def witnesses():
    first, second = [1, 0, ALPHA], [1, 0, -ALPHA]
    same_axis = [record_distribution(state, (0, 0)) for state in (first, second)]
    crossed = [record_distribution(state, (0, math.pi / 2)) for state in (first, second)]
    table = preparation_protocol_table()
    return {"contrast": 0.5, "same_axis_record_tv": float(abs(same_axis[0] - same_axis[1]).sum() / 2),
            "crossed_axis_record_tv": float(abs(crossed[0] - crossed[1]).sum() / 2),
            "crossed_axis_distributions": [row.tolist() for row in crossed],
            "preparation_protocol_table": table.tolist(), "table_rank": int(np.linalg.matrix_rank(table)),
            "table_determinant": float(np.linalg.det(table)),
            "weak_block_rank_diagnostics": [
                {"count": count, "determinant": float(np.linalg.det(two_axis_block_mean_table(count)))}
                for count in (100, 1000, 10000)]}


class MultidirectionalRecordTests(unittest.TestCase):
    def test_same_axis_cannot_see_opposite_transverse_preparations(self):
        for depth in (1, 3, 6):
            np.testing.assert_allclose(record_distribution([1, 0, ALPHA], (0,) * depth),
                                       record_distribution([1, 0, -ALPHA], (0,) * depth), atol=1e-14)

    def test_second_axis_gives_the_analytic_complete_record_gap(self):
        for eta in (0.1, 0.5, 0.9):
            first = record_distribution([1, 0, ALPHA], (0, math.pi / 2), eta)
            second = record_distribution([1, 0, -ALPHA], (0, math.pi / 2), eta)
            self.assertAlmostEqual(first.sum(), 1)
            self.assertAlmostEqual(second.sum(), 1)
            self.assertAlmostEqual(abs(first - second).sum() / 2, eta * ALPHA**2 * math.sqrt(1 - eta**2))

    def test_probability_table_has_analytic_rank_three_minor(self):
        for eta in (0.001, 0.5, 1):
            table = preparation_protocol_table(eta)
            self.assertEqual(np.linalg.matrix_rank(table), 3)
            self.assertAlmostEqual(np.linalg.det(table), -eta**2 * ALPHA**2 / 2)

    def test_zero_contrast_has_no_single_read_rank_witness(self):
        self.assertEqual(np.linalg.matrix_rank(preparation_protocol_table(0)), 1)

    def test_rotated_weak_update_observes_the_intended_coordinate(self):
        for angle in (0, 0.7, math.pi / 2):
            for mark in (0, 1):
                _, _, eta = parameters(128)
                expected = np.array([1, (2 * mark - 1) * eta * math.cos(angle),
                                     (2 * mark - 1) * eta * math.sin(angle)]) / 2
                np.testing.assert_allclose(rotated_weak_branch(angle, mark, 128)[0], expected, atol=1e-14)

    def test_alternating_and_adaptive_protocols_normalize(self):
        for adaptive in (False, True):
            total = 0.0
            for marks in product((0, 1), repeat=6):
                state = np.array([1, 0.2, 0.5])
                for step, mark in enumerate(marks):
                    setting = math.pi / 2 * (marks[step - 1] if adaptive and step else step % 2)
                    state = fresh_branch_matrix(setting, mark, 0.5) @ state
                total += state[0]
            self.assertAlmostEqual(total, 1)

    def test_finite_duration_weak_blocks_keep_a_nonzero_rank_minor(self):
        for count in (100, 1000, 10000):
            table = two_axis_block_mean_table(count)
            self.assertEqual(np.linalg.matrix_rank(table), 3)
            self.assertLess(np.linalg.det(table), -1.6)

    def test_block_mean_rows_match_exact_matrix_recursion(self):
        count = 32
        _, alpha, eta = parameters(count)
        for setting, spatial_coordinate in ((0, 0), (math.pi / 2, 1)):
            nonselective = sum(rotated_weak_branch(setting, mark, count) for mark in (0, 1))
            for column, initial in enumerate(([1, alpha, 0], [1, -alpha, 0], [1, 0, alpha])):
                state = np.asarray(initial, dtype=float)
                mean = 0.0
                effect = rotated_weak_branch(setting, 1, count)[0] - rotated_weak_branch(setting, 0, count)[0]
                for _ in range(count):
                    mean += effect @ state / math.sqrt(count)
                    state = nonselective @ state
                self.assertAlmostEqual(mean, two_axis_block_mean_table(count)[spatial_coordinate + 1, column])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(MultidirectionalRecordTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": tests.testsRun, "failures": len(tests.failures), "errors": len(tests.errors)},
               **witnesses(),
               "scope": "One preparation encoding must serve all subsequently chosen protocols. Rank 3 excludes any common two-hidden-state Markov realization and gives the minimal homogeneous linear predictive dimension. It does not prove that three classical hidden states suffice or that all classical models fail."}
    if args.write_results:
        Path(__file__).with_name("multidirectional_records_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
