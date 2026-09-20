"""Round 214: Jordan-Lie compatibility is stronger than sequential execution.

With declared symmetric observable product and Lie derivation bracket, the
complexified ordered product is associative exactly when its two real
associator contributions cancel. Matrix examples and real/classical controls
audit scope; no cognitive origin of this compatibility is assumed.
"""

import argparse
import json
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from composite_question_symmetry import QUESTIONS
from real_interface_record_bound import IDENTITY, PAULI_X, PAULI_Y, PAULI_Z


def jordan(first, second):
    return (first @ second + second @ first) / 2


def bracket(first, second, scale=1):
    return scale * (first @ second - second @ first) / (2j)


def ordered_product(first, second, scale=1):
    return jordan(first, second) + 1j * bracket(first, second, scale)


def jordan_associator(first, second, third):
    return jordan(jordan(first, second), third) - jordan(first, jordan(second, third))


def ordered_associator(first, second, third, scale=1):
    return (ordered_product(ordered_product(first, second, scale), third, scale)
            - ordered_product(first, ordered_product(second, third, scale), scale))


def compatibility_residual(first, second, third, scale=1):
    return jordan_associator(first, second, third) + bracket(second, bracket(first, third, scale), scale)


def left_jordan_commutator(first, second, third):
    return jordan(first, jordan(second, third)) - jordan(second, jordan(first, third))


class ObservableActionCompatibilityTests(unittest.TestCase):
    def test_square_polarization_recovers_the_declared_jordan_product(self):
        first = PAULI_X + 2 * PAULI_Z
        second = PAULI_Y - PAULI_Z
        polarized = ((first + second) @ (first + second) - first @ first - second @ second) / 2
        np.testing.assert_array_equal(polarized, jordan(first, second))

    def test_lie_bracket_is_skew_jacobi_and_a_jordan_derivation(self):
        basis = (IDENTITY, PAULI_X, PAULI_Y, PAULI_Z)
        for first, second, third in product(basis, repeat=3):
            np.testing.assert_array_equal(bracket(first, second), -bracket(second, first))
            jacobi = (bracket(first, bracket(second, third)) + bracket(second, bracket(third, first))
                      + bracket(third, bracket(first, second)))
            np.testing.assert_array_equal(jacobi, np.zeros((2, 2)))
            np.testing.assert_array_equal(bracket(first, jordan(second, third)),
                                          jordan(bracket(first, second), third) + jordan(second, bracket(first, third)))

    def test_ordered_associator_equals_the_real_compatibility_residual_exactly(self):
        first = PAULI_X + PAULI_Y
        second = 2 * PAULI_X - PAULI_Z
        third = PAULI_Z + PAULI_Y
        for scale in (-2, -1, -0.5, 0, 0.5, 1, 2):
            np.testing.assert_array_equal(ordered_associator(first, second, third, scale),
                                          compatibility_residual(first, second, third, scale))
            np.testing.assert_array_equal(ordered_associator(first, second, third, scale),
                                          (1 - scale ** 2) * jordan_associator(first, second, third))

    def test_associativity_selects_two_orientations_not_one_sign(self):
        for scale in (-2, -1, 0, 1, 2):
            residual = ordered_associator(PAULI_X, PAULI_X, PAULI_Z, scale)
            np.testing.assert_array_equal(residual, (1 - scale ** 2) * PAULI_Z)
        np.testing.assert_array_equal(ordered_product(PAULI_X, PAULI_Z, 1), PAULI_X @ PAULI_Z)
        np.testing.assert_array_equal(ordered_product(PAULI_X, PAULI_Z, -1), PAULI_Z @ PAULI_X)

    def test_commutator_compatibility_holds_on_the_entire_two_qubit_basis(self):
        basis = (np.eye(4),) + QUESTIONS
        for first, second, third in product(basis, repeat=3):
            dynamic = (bracket(first, bracket(second, third)) - bracket(second, bracket(first, third)))
            np.testing.assert_array_equal(dynamic, -left_jordan_commutator(first, second, third))
            np.testing.assert_array_equal(compatibility_residual(first, second, third), np.zeros((4, 4)))

    def test_general_integer_hermitian_examples_do_not_depend_on_pauli_dimension(self):
        random = np.random.default_rng(214)
        for dimension in (2, 3, 4):
            matrices = []
            for _ in range(3):
                raw = random.integers(-2, 3, size=(dimension, dimension)) + 1j * random.integers(-2, 3, size=(dimension, dimension))
                matrices.append(raw + raw.conj().T)
            np.testing.assert_array_equal(compatibility_residual(*matrices), np.zeros((dimension, dimension)))
            np.testing.assert_array_equal(ordered_product(matrices[0], matrices[1]).conj().T,
                                          ordered_product(matrices[1], matrices[0]))

    def test_real_symmetric_observables_do_not_close_under_the_complex_bracket(self):
        np.testing.assert_array_equal(jordan(PAULI_X, PAULI_Z).imag, np.zeros((2, 2)))
        outside = bracket(PAULI_X, PAULI_Z)
        np.testing.assert_array_equal(outside, -PAULI_Y)
        self.assertGreater(np.linalg.norm(outside.imag), 1)
        np.testing.assert_array_equal(jordan_associator(PAULI_X, PAULI_X, PAULI_Z), PAULI_Z)

    def test_classical_diagonal_algebra_with_zero_bracket_is_associative(self):
        first = np.diag([1, 2, -1])
        second = np.diag([2, 0, 3])
        third = np.diag([-1, 1, 2])
        np.testing.assert_array_equal(bracket(first, second, 0), np.zeros((3, 3)))
        np.testing.assert_array_equal(ordered_associator(first, second, third, 0), np.zeros((3, 3)))
        np.testing.assert_array_equal(ordered_product(first, second, 0), first @ second)

    def test_associative_real_channel_execution_does_not_imply_jordan_associativity(self):
        gates = ((IDENTITY + PAULI_X) / 2, (IDENTITY + PAULI_Z) / 2, PAULI_X)
        density = np.diag([1., 0.])
        left_product = (gates[2] @ gates[1]) @ gates[0]
        right_product = gates[2] @ (gates[1] @ gates[0])
        np.testing.assert_array_equal(left_product @ density @ left_product.T,
                                      right_product @ density @ right_product.T)
        np.testing.assert_array_equal(jordan_associator(PAULI_X, PAULI_X, PAULI_Z), PAULI_Z)

    def test_tensor_jordan_and_lie_products_need_cross_terms(self):
        basis = (PAULI_X, PAULI_Y, PAULI_Z)
        for first, second, third, fourth in product(basis, repeat=4):
            left, right = np.kron(first, second), np.kron(third, fourth)
            expected_jordan = (np.kron(jordan(first, third), jordan(second, fourth))
                               - np.kron(bracket(first, third), bracket(second, fourth)))
            expected_bracket = (np.kron(bracket(first, third), jordan(second, fourth))
                                + np.kron(jordan(first, third), bracket(second, fourth)))
            np.testing.assert_array_equal(jordan(left, right), expected_jordan)
            np.testing.assert_array_equal(bracket(left, right), expected_bracket)
        actual = jordan(np.kron(PAULI_X, PAULI_X), np.kron(PAULI_Z, PAULI_Z))
        np.testing.assert_array_equal(actual, -np.kron(PAULI_Y, PAULI_Y))
        np.testing.assert_array_equal(actual.imag, np.zeros((4, 4)))
        self.assertGreater(np.linalg.norm(actual), 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(ObservableActionCompatibilityTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 214,
        "antecedent_rounds": [80, 85, 183, 189, 204, 210, 211, 213],
        "scope": "Jordan product plus Lie derivation bracket; algebraic compatibility audit, not a physical reconstruction",
        "jordan_product": "a circle b = (ab+ba)/2 in matrix controls",
        "bracket": "B(a,b)=(ab-ba)/(2i) in complex matrix controls",
        "ordered_product": "a star b = a circle b + i B(a,b), complex bilinear extension",
        "ordered_associator": "(a circle b) circle c - a circle (b circle c) + B(b,B(a,c))",
        "associativity_iff_compatibility_given_lie_derivation_axioms": True,
        "equivalent_derivation_identity": "[D_a,D_b]=-[L_a,L_b]",
        "scaled_matrix_bracket_associator": "(1-kappa^2) times Jordan associator",
        "noncommutative_matrix_compatible_scales": [-1, 1],
        "absolute_time_scale_or_orientation_selected": False,
        "classical_diagonal_zero_bracket_passes": True,
        "real_symmetric_complex_bracket_is_closed": False,
        "all_possible_real_dynamical_correspondences_excluded_here": False,
        "real_channel_composition_alone_implies_observable_compatibility": False,
        "tensor_cross_term_example": "(X tensor X) circle (Z tensor Z) = -Y tensor Y, a real matrix",
        "extra_observable_permissions_or_resource_costs_derived": False,
        "C_star_norm_or_full_quantum_theory_reconstructed": False,
        "compatibility_derived_from_cognition": False,
        "quantum_incompleteness_proved": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("observable_action_compatibility_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()