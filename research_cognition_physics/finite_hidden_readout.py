"""Finite classical read instruments and a controlled approximation of rotations.

Exactness concerns the allowed predictive records, not the original full latent
distribution. The grid size is a constructive upper bound, not a general optimum.
"""

import argparse
import cmath
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
from shrinking_preparations import parameters


def minimum_grid_size(radius):
    if not 0 <= radius <= 1:
        raise ValueError("Use a normalized preparation radius.")
    if radius == 1:
        return None
    size = max(3, math.ceil(math.pi / math.acos(radius)))
    while size > 3 and math.cos(math.pi / (size - 1)) >= radius:
        size -= 1
    return size


def embedding(size):
    if isinstance(size, bool) or not isinstance(size, int) or size < 3:
        raise ValueError("The regular polygon needs at least three labels.")
    angles = 2 * math.pi * np.arange(size) / size
    return np.stack((np.ones(size), np.cos(angles), np.sin(angles)))


def encode_summary(summary, size):
    """Positive polygon decomposition: uniform center plus two adjacent vertices."""
    summary = np.asarray(summary, dtype=float)
    if summary.shape != (3,) or not np.all(np.isfinite(summary)) or summary[0] < 0:
        raise ValueError("Use a finite subnormalized positive summary.")
    embedding(size)
    weight = summary[0]
    if weight == 0:
        if np.linalg.norm(summary[1:]) > 0:
            raise ValueError("Zero mass has zero spatial moments.")
        return np.zeros(size)
    radius = np.linalg.norm(summary[1:]) / weight
    if radius == 0:
        return np.full(size, weight / size)
    angle = math.atan2(summary[2], summary[1]) % (2 * math.pi)
    delta = 2 * math.pi / size
    left = math.floor(angle / delta) % size
    offset = angle - left * delta
    first = radius * math.sin(delta - offset) / math.sin(delta)
    second = radius * math.sin(offset) / math.sin(delta)
    center = 1 - first - second
    if min(first, second, center) < -2e-13:
        raise ValueError("The summary lies outside this polygon.")
    probabilities = np.full(size, max(center, 0) / size)
    probabilities[left] += max(first, 0)
    probabilities[(left + 1) % size] += max(second, 0)
    return weight * probabilities / probabilities.sum()


def fresh_matrix(setting, outcome, contrast, radius=ALPHA):
    matrix = fresh_branch_matrix(setting, outcome, contrast)
    matrix[1:] *= radius / ALPHA
    return matrix


def read_kernel(size, setting, outcome, contrast, radius=ALPHA):
    if radius > math.cos(math.pi / size) + 1e-14:
        raise ValueError("This grid cannot contain every allowed preparation.")
    outputs = fresh_matrix(setting, outcome, contrast, radius) @ embedding(size)
    return np.column_stack([encode_summary(outputs[:, index], size) for index in range(size)])


def grid_rotation(size, angle):
    embedding(size)
    if not math.isfinite(angle):
        raise ValueError("The rotation angle must be finite.")
    delta = 2 * math.pi / size
    position = (angle % (2 * math.pi)) / delta
    advance = math.floor(position)
    fraction = position - advance
    result = np.zeros((size, size))
    columns = np.arange(size)
    result[(columns + advance) % size, columns] += 1 - fraction
    result[(columns + advance + 1) % size, columns] += fraction
    return result


def rotation_multiplier(size, angle):
    delta = 2 * math.pi / size
    position = (angle % (2 * math.pi)) / delta
    advance = math.floor(position)
    fraction = position - advance
    value = (1 - fraction) * cmath.exp(1j * advance * delta) + fraction * cmath.exp(1j * (advance + 1) * delta)
    return value, fraction


def record_tv_bound(size, rotations):
    return min(1.0, rotations * math.pi**2 / (4 * size**2))


def sufficient_grid_size(rotations, tolerance, radius=ALPHA):
    if rotations < 0 or not 0 < tolerance < 1:
        raise ValueError("Use a nonnegative rotation budget and tolerance in (0,1).")
    minimum = minimum_grid_size(radius)
    if minimum is None:
        return None
    return max(minimum, math.ceil(math.pi / 2 * math.sqrt(rotations / tolerance)))


def protocol_distribution(initial, operations, size, hidden):
    states = (encode_summary(initial, size) if hidden else np.asarray(initial, dtype=float))[None, :]
    for operation in operations:
        if operation[0] == "rotate":
            matrix = grid_rotation(size, operation[1]) if hidden else translation_matrix(operation[1], (1.0,))
            states = states @ matrix.T
        elif operation[0] == "read":
            _, setting, contrast = operation
            matrices = [read_kernel(size, setting, mark, contrast) if hidden else fresh_matrix(setting, mark, contrast)
                        for mark in (0, 1)]
            states = np.concatenate([states @ matrix.T for matrix in matrices])
        else:
            raise ValueError("Unknown operation.")
    return states.sum(axis=1) if hidden else states[:, 0]


def adaptive_distribution(initial, size, hidden):
    probabilities = []
    for marks in product((0, 1), repeat=4):
        state = encode_summary(initial, size) if hidden else np.asarray(initial, dtype=float)
        for step, mark in enumerate(marks):
            angle = 0.4 + 0.2 * sum(marks[:step])
            setting = math.pi / 2 * (sum(marks[:step]) % 2)
            rotation = grid_rotation(size, angle) if hidden else translation_matrix(angle, (1.0,))
            reading = read_kernel(size, setting, mark, 0.5) if hidden else fresh_matrix(setting, mark, 0.5)
            state = reading @ rotation @ state
        probabilities.append(state.sum() if hidden else state[0])
    return np.asarray(probabilities)


def diagnostics():
    rows = []
    for size in (22, 32, 64, 128, 256):
        multiplier, fraction = rotation_multiplier(size, 1.0)
        delta = 2 * math.pi / size
        rows.append({"labels": size, "grid_inradius": math.cos(math.pi / size),
                     "rotation_one_radian_single_moment_error": abs(multiplier - cmath.exp(1j)),
                     "single_moment_error_bound": fraction * (1 - fraction) * delta**2 / 2,
                     "one_roundtrip_moment_gain": abs(multiplier)**2,
                     "hundred_roundtrips_moment_gain": abs(multiplier)**200,
                     "all_record_tv_bound_for_100_rotations": record_tv_bound(size, 100)})
    scaling = []
    for count in (100, 1000, 10000, 100000):
        _, alpha, _ = parameters(count)
        scaling.append({"weak_experiment_steps": count, "noise": 1, "preparation_radius": alpha,
                        "sufficient_labels_for_exact_reads": minimum_grid_size(alpha),
                        "leading_sqrt_count_formula": math.pi * math.sqrt(3 * count)})
    return {"rotation_diagnostics": rows, "shrinking_preparation_read_resources": scaling}


class FiniteHiddenReadoutTests(unittest.TestCase):
    def test_twenty_two_labels_contain_the_original_preparation_disk(self):
        self.assertEqual(minimum_grid_size(ALPHA), 22)
        self.assertGreaterEqual(math.cos(math.pi / 22), ALPHA)
        self.assertLess(math.cos(math.pi / 21), ALPHA)
        self.assertIsNone(minimum_grid_size(1))

    def test_positive_encoding_recovers_boundary_and_interior_summaries(self):
        for radius in (0, 0.5 * ALPHA, ALPHA):
            for angle in np.linspace(0, 2 * math.pi, 65):
                summary = np.array([1, radius * math.cos(angle), radius * math.sin(angle)])
                encoded = encode_summary(summary, 22)
                self.assertGreaterEqual(encoded.min(), 0)
                np.testing.assert_allclose(embedding(22) @ encoded, summary, atol=2e-14)

    def test_insufficient_grid_has_a_specific_preparation_witness(self):
        angle = math.pi / 21
        with self.assertRaises(ValueError):
            encode_summary([1, ALPHA * math.cos(angle), ALPHA * math.sin(angle)], 21)

    def test_all_read_instruments_are_positive_and_normalized(self):
        for angle, eta in product((0, 0.3, math.pi / 2), (0, 0.5, 1)):
            matrices = [read_kernel(22, angle, mark, eta) for mark in (0, 1)]
            self.assertGreaterEqual(min(matrix.min() for matrix in matrices), 0)
            np.testing.assert_allclose(sum(matrices).sum(axis=0), 1, atol=1e-14)

    def test_exact_intertwining_including_sharp_and_rotated_reads(self):
        for angle, eta, mark in product((0, 0.3, math.pi / 2), (0, 0.5, 1), (0, 1)):
            np.testing.assert_allclose(embedding(22) @ read_kernel(22, angle, mark, eta),
                                       fresh_matrix(angle, mark, eta) @ embedding(22), atol=2e-14)

    def test_complete_multidirectional_records_are_exact_without_silent_rotations(self):
        operations = [("read", angle, eta) for angle, eta in ((0, 0.5), (0.7, 0.2), (math.pi / 2, 1), (-0.4, 0.7))]
        for initial in ([1, ALPHA, 0], [1, 0, ALPHA], [1, 0.3, -0.2]):
            np.testing.assert_allclose(protocol_distribution(initial, operations, 22, True),
                                       protocol_distribution(initial, operations, 22, False), atol=2e-14)

    def test_different_hidden_encodings_of_same_moments_have_identical_read_records(self):
        first = np.ones(22) / 22
        second = np.zeros(22)
        second[0] = second[11] = 0.5
        np.testing.assert_allclose(embedding(22) @ first, embedding(22) @ second, atol=1e-14)
        for setting, outcome in ((0.3, 1), (1.2, 0), (-0.2, 1)):
            matrix = read_kernel(22, setting, outcome, 0.5)
            first, second = matrix @ first, matrix @ second
            self.assertAlmostEqual(first.sum(), second.sum())

    def test_grid_rotation_has_the_claimed_harmonic_multiplier(self):
        size, angle = 32, 0.73
        complex_row = embedding(size)[1] + 1j * embedding(size)[2]
        multiplier, _ = rotation_multiplier(size, angle)
        np.testing.assert_allclose(complex_row @ grid_rotation(size, angle), multiplier * complex_row, atol=2e-14)

    def test_actual_preparation_mixtures_need_not_use_a_unique_hidden_encoding(self):
        first, second = np.array([1, ALPHA, 0]), np.array([1, 0, ALPHA])
        mixture = (encode_summary(first, 22) + encode_summary(second, 22)) / 2
        representative = encode_summary((first + second) / 2, 22)
        self.assertGreater(np.linalg.norm(mixture - representative), 0.1)
        np.testing.assert_allclose(embedding(22) @ mixture, embedding(22) @ representative, atol=1e-14)
        for setting, outcome in ((0.2, 1), (1.0, 0), (-0.5, 1)):
            kernel = read_kernel(22, setting, outcome, 0.5)
            mixture, representative = kernel @ mixture, kernel @ representative
            self.assertAlmostEqual(mixture.sum(), representative.sum())

    def test_uniform_interpolation_error_and_irreversible_roundtrip(self):
        for size, angle in product((22, 32, 64), (0.1, 0.7, 1.0, -0.3)):
            value, fraction = rotation_multiplier(size, angle)
            delta = 2 * math.pi / size
            self.assertLessEqual(abs(value - cmath.exp(1j * angle)), fraction * (1 - fraction) * delta**2 / 2 + 1e-14)
            self.assertAlmostEqual(abs(value)**2, 1 - 4 * fraction * (1 - fraction) * math.sin(delta / 2)**2)
            self.assertLess(abs(value), 1)

    def test_grid_aligned_rotations_are_exact_permutations(self):
        size = 32
        rotation = grid_rotation(size, 6 * math.pi / size)
        reverse = grid_rotation(size, -6 * math.pi / size)
        np.testing.assert_allclose(reverse @ rotation, np.eye(size), atol=2e-14)

    def test_complete_record_error_respects_bound_with_interleaved_rotations(self):
        operations = [("rotate", 0.7), ("read", 0, 0.5), ("rotate", -0.3),
                      ("read", 0.8, 0.9), ("rotate", 1.0), ("read", math.pi / 2, 0.4)]
        for size in (22, 32, 64):
            exact = protocol_distribution([1, 0.3, 0.8], operations, size, False)
            approximate = protocol_distribution([1, 0.3, 0.8], operations, size, True)
            self.assertAlmostEqual(approximate.sum(), 1)
            self.assertLessEqual(abs(exact - approximate).sum() / 2, record_tv_bound(size, 3) + 1e-14)

    def test_adaptive_record_error_respects_the_same_budget(self):
        exact = adaptive_distribution([1, 0.3, 0.8], 22, False)
        approximate = adaptive_distribution([1, 0.3, 0.8], 22, True)
        self.assertAlmostEqual(approximate.sum(), 1)
        self.assertLessEqual(abs(exact - approximate).sum() / 2, record_tv_bound(22, 4) + 1e-14)

    def test_constructive_resource_budget_meets_requested_record_accuracy(self):
        size = sufficient_grid_size(100, 0.01)
        self.assertEqual(size, 158)
        self.assertLessEqual(record_tv_bound(size, 100), 0.01)
        self.assertGreater(record_tv_bound(size - 1, 100), 0.01)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(FiniteHiddenReadoutTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": tests.testsRun, "failures": len(tests.failures), "errors": len(tests.errors)},
               "original_preparation_radius": ALPHA, "sufficient_labels_for_exact_fresh_reads": minimum_grid_size(ALPHA),
               "sufficient_labels_for_100_rotations_tv_0_01": sufficient_grid_size(100, 0.01),
               **diagnostics(),
               "scope": "All directions and strengths of the lambda=1 fresh-read family are exact at the operational three-moment level. Actual preparation mixtures mix hidden distributions linearly, but operationally equivalent recipes need not have identical hidden distributions. Silent rotations are approximated with complete-record TV <= L*pi^2/(4*M^2), even with adaptive settings, for at most L rotations. This is a constructive grid bound, not the optimal hidden-label count over all classical models or finite-precision bits."}
    if args.write_results:
        Path(__file__).with_name("finite_hidden_readout_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
