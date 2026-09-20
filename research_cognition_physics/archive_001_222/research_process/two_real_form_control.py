"""Round 208: two active real forms on one declared complex carrier.

The generators have an exact commutator certificate for su(d). Complex
matrices, active relative frame control and closure permissions are inputs,
not consequences of cognition. No finite-cost universality is claimed.
"""

import argparse
import json
import unittest
from pathlib import Path

import numpy as np

from common_orientation_structure import encode_state, orientation
from complex_control_from_reference import real_lift


def matrix_unit(dimension, row, column):
    matrix = np.zeros((dimension, dimension), dtype=complex)
    matrix[row, column] = 1
    return matrix


def real_generator(dimension, row, column):
    return matrix_unit(dimension, row, column) - matrix_unit(dimension, column, row)


def imaginary_generator(dimension, row, column):
    return 1j * (matrix_unit(dimension, row, column) + matrix_unit(dimension, column, row))


def commutator(left, right):
    return left @ right - right @ left


def phase_frame(dimension, angle=None):
    frame = np.eye(dimension, dtype=complex)
    frame[0, 0] = 1j if angle is None else np.exp(1j * angle)
    return frame


def flow(generator, parameter):
    eigenvalues, eigenvectors = np.linalg.eigh(1j * generator)
    return (eigenvectors * np.exp(-1j * parameter * eigenvalues)) @ eigenvectors.conj().T


def frame_shadow(density, frame):
    coordinates = frame.conj().T @ density @ frame
    return frame @ coordinates.real @ frame.conj().T


def generated_basis(dimension):
    real_part = [real_generator(dimension, row, column)
                 for row in range(dimension) for column in range(row + 1, dimension)]
    frame = phase_frame(dimension)
    star = {column: frame @ real_generator(dimension, 0, column) @ frame.conj().T
            for column in range(1, dimension)}
    imaginary_part = list(star.values())
    imaginary_part.extend(-commutator(real_generator(dimension, 0, row), star[column])
                          for row in range(1, dimension)
                          for column in range(row + 1, dimension))
    diagonal_part = [commutator(real_generator(dimension, 0, column), star[column]) / 2
                     for column in range(1, dimension)]
    return real_part + imaginary_part + diagonal_part


def gram_matrix(basis):
    return np.array([[np.vdot(left, right).real for right in basis] for left in basis])


class TwoRealFormControlTests(unittest.TestCase):
    def test_phase_star_generators_and_commutators_are_exact(self):
        for dimension in range(2, 7):
            frame = phase_frame(dimension)
            for column in range(1, dimension):
                original = real_generator(dimension, 0, column)
                shifted = frame @ original @ frame.conj().T
                np.testing.assert_array_equal(shifted, imaginary_generator(dimension, 0, column))
                expected = 2j * (matrix_unit(dimension, 0, 0)
                                 - matrix_unit(dimension, column, column))
                np.testing.assert_array_equal(commutator(original, shifted), expected)
            for row in range(1, dimension):
                for column in range(row + 1, dimension):
                    actual = -commutator(real_generator(dimension, 0, row),
                                         imaginary_generator(dimension, 0, column))
                    np.testing.assert_array_equal(actual, imaginary_generator(dimension, row, column))

    def test_constructive_basis_has_exact_nonsingular_gram_certificate(self):
        for dimension in range(2, 7):
            basis = generated_basis(dimension)
            self.assertEqual(len(basis), dimension * dimension - 1)
            off_diagonal_count = dimension * (dimension - 1)
            expected = 2 * np.eye(len(basis))
            expected[off_diagonal_count:, off_diagonal_count:] = (
                np.eye(dimension - 1) + np.ones((dimension - 1, dimension - 1)))
            np.testing.assert_array_equal(gram_matrix(basis), expected)
            for generator in basis:
                np.testing.assert_array_equal(generator.conj().T, -generator)
                self.assertEqual(np.trace(generator), 0)

    def test_generic_nonzero_phase_isolates_the_same_imaginary_generator(self):
        for angle in (0.03, 0.4, 1.2, 2.8, -0.7):
            frame = phase_frame(4, angle)
            for column in range(1, 4):
                original = real_generator(4, 0, column)
                shifted = frame @ original @ frame.conj().T
                isolated = (shifted - np.cos(angle) * original) / np.sin(angle)
                np.testing.assert_allclose(isolated, imaginary_generator(4, 0, column), atol=1e-14)

    def test_real_frame_endpoints_add_no_new_generators(self):
        for sign in (-1, 1):
            frame = np.diag([sign, 1, 1, 1])
            for row in range(4):
                for column in range(row + 1, 4):
                    original = real_generator(4, row, column)
                    factor = sign if row == 0 else 1
                    np.testing.assert_array_equal(frame @ original @ frame.T, factor * original)

    def test_each_family_preserves_its_own_real_state_cone(self):
        raw = np.array([[1., 2., 0.], [0., 1., -1.], [2., 0., 1.]])
        density = raw @ raw.T / np.trace(raw @ raw.T)
        rotation = flow(real_generator(3, 0, 1), 0.37)
        for frame in (np.eye(3), phase_frame(3)):
            encoded = frame @ density @ frame.conj().T
            gate = frame @ rotation @ frame.conj().T
            output = gate @ encoded @ gate.conj().T
            np.testing.assert_allclose(frame_shadow(output, frame), output, atol=1e-15)

    def test_active_second_family_can_leave_the_first_real_cone(self):
        initial = np.diag([1., 0.])
        gate = flow(imaginary_generator(2, 0, 1), np.pi / 4)
        output = gate @ initial @ gate.conj().T
        expected = np.array([[0.5, -0.5j], [0.5j, 0.5]])
        np.testing.assert_allclose(output, expected, atol=5e-16)
        self.assertGreater(np.linalg.norm(output - output.real), 0.7)
        np.testing.assert_allclose(frame_shadow(output, phase_frame(2)), output, atol=1e-15)

    def test_passive_relabeling_preserves_control_dimension_and_records(self):
        dimension = 3
        frame = phase_frame(dimension)
        basis = [real_generator(dimension, row, column)
                 for row in range(dimension) for column in range(row + 1, dimension)]
        changed = [frame @ generator @ frame.conj().T for generator in basis]
        np.testing.assert_array_equal(gram_matrix(changed), gram_matrix(basis))
        density = np.diag([1., 0., 0.])
        effect = np.diag([0., 1., 0.])
        gate = flow(basis[0], 0.31)
        original_probability = np.trace(effect @ gate @ density @ gate.conj().T)
        changed_gate = frame @ gate @ frame.conj().T
        changed_probability = np.trace((frame @ effect @ frame.conj().T) @ changed_gate
                                       @ (frame @ density @ frame.conj().T) @ changed_gate.conj().T)
        self.assertAlmostEqual(original_probability.real, changed_probability.real)

    def test_real_lift_reproduces_an_alternating_control_record(self):
        dimension = 3
        density = np.diag([1., 0., 0.])
        first = flow(real_generator(dimension, 0, 1), 0.47)
        second = flow(imaginary_generator(dimension, 0, 2), 0.39)
        gate = second @ first
        lifted = real_lift(gate)
        np.testing.assert_allclose(lifted.T @ lifted, np.eye(2 * dimension), atol=1e-15)
        np.testing.assert_array_equal(lifted @ orientation(dimension), orientation(dimension) @ lifted)
        real_output = lifted @ encode_state(density) @ lifted.T
        complex_output = gate @ density @ gate.conj().T
        np.testing.assert_allclose(real_output, encode_state(complex_output), atol=1e-15)
        for index in range(dimension):
            effect = matrix_unit(dimension, index, index)
            self.assertAlmostEqual(np.trace(real_output @ real_lift(effect)).real,
                                   np.trace(complex_output @ effect).real)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(TwoRealFormControlTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 208,
        "antecedent_rounds": [50, 54, 57, 80, 85, 87, 88, 171, 204, 207],
        "scope": "Two active orthogonal-control frames on the same declared finite complex carrier",
        "inputs": ["complex positive matrices and trace probabilities", "all real SO(d) flows",
                   "second frame D_phi=diag(exp(i phi),1,...,1)",
                   "both families executable on the same carrier with signed continuous parameters",
                   "finite-product closure for arbitrarily accurate control"],
        "lie_algebra_if_sin_phi_nonzero": "su(d)",
        "lie_algebra_if_sin_phi_zero": "so(d)",
        "exact_certificate_phase": "pi/2",
        "basis_gram_determinant": "d * 2**(d*(d-1))",
        "dimension_rows": [{"d": dimension, "one_real_form": dimension * (dimension - 1) // 2,
                            "two_active_forms": len(generated_basis(dimension))}
                           for dimension in range(2, 7)],
        "passive_coordinate_change_adds_capability": False,
        "first_real_cone_stable_under_all_cross_frame_controls": False,
        "same_process_has_common_J_real_representation": True,
        "finite_gate_cost_or_uniform_small_angle_budget_proved": False,
        "all_general_instruments_or_composition_rules_derived": False,
        "relative_frame_necessity_derived_from_cognition": False,
        "quantum_necessity_or_incompleteness_proved": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("two_real_form_control_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()