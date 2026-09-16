"""Round 127: deterministic approximate intake with exact host preservation.

The negative orientation branch uses the known physical approximate transpose
T_d(X)=(tr(X) I + X.T)/(d+1). All real Kraus branches are retained. Ideal
control is assumed; this is not a compilation into the old noisy primitives.
"""

import argparse
import json
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from common_orientation_structure import encode_state
from encoded_composition_audit import independent_encoding
from flagged_subject_joining import joining_kraus, random_state
from quantum_interface_audit import PAULI_X, PAULI_Y, PAULI_Z


def trace_distance(first, second):
    return float(np.abs(np.linalg.eigvalsh(first - second)).sum() / 2)


def transpose_kraus(dimension):
    if not isinstance(dimension, int) or dimension < 2:
        raise ValueError("A target dimension of at least two is required.")
    result = []
    for i in range(dimension):
        diagonal = np.zeros((dimension, dimension))
        diagonal[i, i] = np.sqrt(2 / (dimension + 1))
        result.append(diagonal)
        for j in range(i + 1, dimension):
            off_diagonal = np.zeros((dimension, dimension))
            off_diagonal[i, j] = off_diagonal[j, i] = 1 / np.sqrt(dimension + 1)
            result.append(off_diagonal)
    return tuple(result)


def approximate_transpose(matrix):
    dimension = len(matrix)
    return (np.trace(matrix) * np.eye(dimension) + matrix.T) / (dimension + 1)


def deterministic_joining_kraus(host_dimension=2, newcomer_dimension=2):
    """First Kraus label h=0; all subsequent labels h=1. No rejected result."""
    initial = joining_kraus((host_dimension, newcomer_dimension))
    result = [initial[(0, 0)]]
    for k in transpose_kraus(newcomer_dimension):
        repair = np.kron(np.eye(2 * host_dimension), k)
        result.append(repair @ initial[(0, 1)])
    return tuple(result)


def apply_kraus(kraus, density):
    return sum(k @ density @ k.T for k in kraus)


def join_output(host, newcomer):
    return apply_kraus(deterministic_joining_kraus(len(host), len(newcomer)),
                       independent_encoding((host, newcomer)))


def logical_newcomer(newcomer):
    dimension = len(newcomer)
    return ((dimension + 2) * newcomer + np.trace(newcomer) * np.eye(dimension)) / (2 * (dimension + 1))


def error_bound(dimension):
    if not isinstance(dimension, int) or dimension < 2:
        raise ValueError("A target dimension of at least two is required.")
    return Fraction(dimension - 1, 2 * (dimension + 1))


def encoded_host_marginal(output, host_dimension, newcomer_dimension):
    combined = 2 * host_dimension
    return np.einsum("abcb->ac", output.reshape(combined, newcomer_dimension, combined, newcomer_dimension))


class ApproximateSubjectJoiningTests(unittest.TestCase):
    def test_real_transpose_kraus_form_a_channel_and_match_known_map(self):
        rng = np.random.default_rng(127)
        for d in (2, 3, 4):
            ks = transpose_kraus(d)
            self.assertEqual(len(ks), d * (d + 1) // 2)
            np.testing.assert_allclose(sum(k.T @ k for k in ks), np.eye(d), atol=5e-16)
            # Non-Hermitian matrix units check the entire linear map.
            for i in range(d):
                for j in range(d):
                    matrix = np.zeros((d, d), dtype=complex)
                    matrix[i, j] = 1j
                    np.testing.assert_allclose(apply_kraus(ks, matrix), approximate_transpose(matrix), atol=3e-16)
            density = random_state(rng, d)
            self.assertGreaterEqual(np.linalg.eigvalsh(apply_kraus(ks, density)).min(), 0)

    def test_full_join_is_trace_preserving_with_no_selected_success_branch(self):
        for da, db in ((2, 2), (3, 2), (2, 3)):
            ks = deterministic_joining_kraus(da, db)
            np.testing.assert_allclose(sum(k.T @ k for k in ks), np.eye(4 * da * db), atol=1e-15)
            self.assertTrue(all(np.isrealobj(k) for k in ks))

    def test_unknown_mixed_input_formula_and_exact_host_preservation(self):
        rng = np.random.default_rng(227)
        for da, db in ((2, 2), (3, 2), (2, 3), (3, 3)):
            for _ in range(3):
                host, new = random_state(rng, da), random_state(rng, db)
                actual = join_output(host, new)
                np.testing.assert_allclose(actual, encode_state(np.kron(host, logical_newcomer(new))), atol=3e-16)
                np.testing.assert_allclose(encoded_host_marginal(actual, da, db), encode_state(host), atol=3e-16)

    def test_uniform_error_bound_is_attained_by_every_tested_pure_newcomer(self):
        rng = np.random.default_rng(327)
        for d in (2, 3, 4):
            host = random_state(rng, 2)
            for _ in range(5):
                vector = rng.normal(size=d) + 1j * rng.normal(size=d)
                vector /= np.linalg.norm(vector)
                new = np.outer(vector, vector.conj())
                error = trace_distance(join_output(host, new), encode_state(np.kron(host, new)))
                self.assertAlmostEqual(error, float(error_bound(d)), places=14)
            mixed = random_state(rng, d)
            self.assertLessEqual(trace_distance(join_output(host, mixed), encode_state(np.kron(host, mixed))),
                                 float(error_bound(d)) + 1e-15)

    def test_kept_orientation_flag_has_zero_and_one_third_conditional_errors(self):
        host = (np.eye(2) + .6 * PAULI_Y) / 2
        new = (np.eye(2) + PAULI_X) / 2
        target = encode_state(np.kron(host, new))
        source = independent_encoding((host, new))
        ks = deterministic_joining_kraus()
        branches = (ks[0] @ source @ ks[0].T, apply_kraus(ks[1:], source))
        for branch in branches:
            self.assertAlmostEqual(np.trace(branch), .5)
        self.assertAlmostEqual(trace_distance(2 * branches[0], target), 0.)
        self.assertAlmostEqual(trace_distance(2 * branches[1], target), 1 / 3)
        self.assertAlmostEqual(sum(trace_distance(branch, target / 2) for branch in branches), 1 / 6)

    def test_qubit_old_distinguishability_contracts_by_two_thirds_on_every_axis(self):
        host = np.array([[.6, .1j], [-.1j, .4]])
        for pauli in (PAULI_X, PAULI_Y, PAULI_Z):
            first, second = [(np.eye(2) + sign * pauli) / 2 for sign in (-1, 1)]
            self.assertAlmostEqual(trace_distance(join_output(host, first), join_output(host, second)), 2 / 3)

    def test_retained_environment_is_an_isometric_extension_of_the_actual_channel(self):
        ks = deterministic_joining_kraus()
        isometry = np.vstack(ks)  # Four environment labels, then eight active dimensions.
        np.testing.assert_allclose(isometry.T @ isometry, np.eye(16), atol=1e-15)
        rng = np.random.default_rng(427)
        raw = rng.normal(size=(16, 16))
        original = raw @ raw.T
        original /= np.trace(original)
        whole = isometry @ original @ isometry.T
        np.testing.assert_allclose(isometry.T @ whole @ isometry, original, atol=3e-16)
        active = np.einsum("abad->bd", whole.reshape(4, 8, 4, 8))
        np.testing.assert_allclose(active, apply_kraus(ks, original), atol=1e-16)

    def test_invalid_dimensions_and_exact_resource_numbers(self):
        for d in (0, 1, -2, 2.5):
            with self.assertRaises(ValueError):
                transpose_kraus(d)
            with self.assertRaises(ValueError):
                error_bound(d)
        self.assertEqual(error_bound(2), Fraction(1, 6))
        self.assertEqual(error_bound(3), Fraction(1, 4))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(ApproximateSubjectJoiningTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 127,
        "input": "One independent standard real encoding per unknown subject",
        "all_outcomes_accepted": True,
        "logical_output": "rho_A tensor (((d_B+2) rho_B + I)/(2(d_B+1)))",
        "host_state_preserved_exactly": True,
        "qubit_newcomer_worst_trace_distance_exact": "1/6",
        "qubit_newcomer_distinguishability_factor_exact": "2/3",
        "qubit_negative_orientation_conditional_error_exact": "1/3",
        "orientation_flagged_error_against_same_uniform_flag_exact": "1/6",
        "all_fine_environment_records_claimed_to_obey_same_error_bound": False,
        "whole_environment_retained_in_isometric_dilation": True,
        "qubit_kraus_count": 4,
        "qubit_dilation_environment_dimension": 4,
        "qubit_input_dimension": 16,
        "qubit_active_output_dimension": 8,
        "qubit_additional_pure_rebits_sufficient_for_abstract_orthogonal_extension": 1,
        "pure_ancilla_initialization_and_flag_copying_not_free": True,
        "compiled_into_old_noisy_primitives": False,
        "global_optimality_claimed_in_this_round": False,
        "rows": [{"newcomer_dimension": d, "worst_error_exact": str(error_bound(d))}
                 for d in (2, 3, 4, 8)],
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("approximate_subject_joining_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
