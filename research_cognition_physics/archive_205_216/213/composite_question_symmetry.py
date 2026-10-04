"""Round 213: full SO(n) question symmetry is not a composite quantum axiom.

Two declared complex qubits have a fifteen-dimensional traceless question
space. Their commutator action is linear, faithful and conservative, but
covariant under unitary-induced rotations rather than all SO(15). An exact
positive-cone obstruction prevents treating the full composite as a ball.
"""

import argparse
import json
import math
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from common_orientation_structure import encode_state
from complex_control_from_reference import real_lift
from question_action_covariance import quarter_turn
from real_interface_record_bound import IDENTITY, PAULI_X, PAULI_Y, PAULI_Z


PAULI = {"I": IDENTITY, "X": PAULI_X, "Y": PAULI_Y, "Z": PAULI_Z}
LABELS = tuple(first + second for first, second in product(PAULI, repeat=2)
               if first + second != "II")
QUESTIONS = tuple(np.kron(PAULI[label[0]], PAULI[label[1]]) for label in LABELS)


def coefficients(matrix):
    return np.array([np.trace(question @ matrix).real / 4 for question in QUESTIONS])


def observable(coordinates):
    if len(coordinates) != 15:
        raise ValueError("The two-qubit question vector has fifteen components.")
    return sum(coordinate * question for coordinate, question in zip(coordinates, QUESTIONS))


def density_from_answers(coordinates):
    return (np.eye(4) + observable(coordinates)) / 4


def basis_vector(label):
    vector = np.zeros(15)
    vector[LABELS.index(label)] = 1
    return vector


def action_generator(coordinates):
    hamiltonian = observable(coordinates) / 2
    return np.column_stack([coefficients(-1j * (hamiltonian @ question - question @ hamiltonian))
                            for question in QUESTIONS])


def pauli_unitary(label, angle):
    question = QUESTIONS[LABELS.index(label)]
    return math.cos(angle / 2) * np.eye(4) - 1j * math.sin(angle / 2) * question


def induced_rotation(unitary):
    return np.column_stack([coefficients(unitary @ question @ unitary.conj().T) for question in QUESTIONS])


def coordinate_rotation(first_label, second_label, angle):
    rotation = np.eye(15)
    first = LABELS.index(first_label)
    second = LABELS.index(second_label)
    rotation[first, first] = rotation[second, second] = math.cos(angle)
    rotation[first, second] = -math.sin(angle)
    rotation[second, first] = math.sin(angle)
    return rotation


class CompositeQuestionSymmetryTests(unittest.TestCase):
    def test_pauli_coordinates_are_an_exact_orthogonal_complete_traceless_basis(self):
        self.assertEqual(len(QUESTIONS), 15)
        gram = np.array([[np.trace(first @ second).real for second in QUESTIONS] for first in QUESTIONS])
        np.testing.assert_array_equal(gram, 4 * np.eye(15))
        for question in QUESTIONS:
            np.testing.assert_array_equal(question.conj().T, question)
            np.testing.assert_array_equal(question @ question, np.eye(4))
            self.assertEqual(np.trace(question), 0)

    def test_each_action_generator_is_skew_and_exactly_conserves_its_question(self):
        for axis in np.eye(15):
            generator = action_generator(axis)
            np.testing.assert_array_equal(generator.T, -generator)
            np.testing.assert_array_equal(generator @ axis, np.zeros(15))
        first = np.arange(15, dtype=int) - 7
        second = np.arange(15, dtype=int) % 3 - 1
        np.testing.assert_array_equal(action_generator(2 * first - 3 * second),
                                      2 * action_generator(first) - 3 * action_generator(second))
        np.testing.assert_array_equal(action_generator(first) @ first, np.zeros(15))

    def test_generator_map_is_injective_by_an_exact_gram_certificate(self):
        generators = [action_generator(axis) for axis in np.eye(15)]
        gram = np.array([[np.sum(first * second) for second in generators] for first in generators])
        np.testing.assert_array_equal(gram, 8 * np.eye(15))

    def test_action_commutators_close_and_satisfy_the_induced_jacobi_relation(self):
        choices = (basis_vector("XI"), basis_vector("ZI"), basis_vector("XZ"), basis_vector("IY"))
        for first in choices:
            for second in choices:
                left = action_generator(first)
                right = action_generator(second)
                np.testing.assert_array_equal(left @ right - right @ left, action_generator(left @ second))

    def test_actual_unitary_rotations_preserve_the_cone_and_covariant_action(self):
        random = np.random.default_rng(213)
        for _ in range(8):
            raw = random.normal(size=(4, 4)) + 1j * random.normal(size=(4, 4))
            density = raw @ raw.conj().T
            density /= np.trace(density)
            coordinates = np.array([np.trace(density @ question).real for question in QUESTIONS])
            question = random.integers(-2, 3, size=15)
            unitary = pauli_unitary("XY", 0.37) @ pauli_unitary("ZI", -0.51)
            rotation = induced_rotation(unitary)
            np.testing.assert_allclose(rotation.T @ rotation, np.eye(15), atol=8e-16)
            self.assertAlmostEqual(np.linalg.det(rotation), 1)
            output = unitary @ density @ unitary.conj().T
            np.testing.assert_allclose(density_from_answers(rotation @ coordinates), output, atol=3e-16)
            self.assertGreater(np.linalg.eigvalsh(output).min(), 0)
            np.testing.assert_allclose(action_generator(rotation @ question),
                                       rotation @ action_generator(question) @ rotation.T, atol=2e-15)

    def test_so15_covariance_fails_on_a_commuting_to_noncommuting_pair(self):
        first = basis_vector("XI")
        second = basis_vector("IX")
        third = basis_vector("ZI")
        rotation = quarter_turn(15, LABELS.index("IX"), LABELS.index("ZI"))
        np.testing.assert_array_equal(rotation @ first, first)
        np.testing.assert_array_equal(rotation @ second, third)
        np.testing.assert_array_equal(action_generator(first) @ second, np.zeros(15))
        self.assertEqual(np.linalg.norm(action_generator(first) @ third), 1)
        defect = action_generator(rotation @ first) - rotation @ action_generator(first) @ rotation.T
        self.assertEqual(np.linalg.norm(defect @ (rotation @ second)), 1)

    def test_norm_preserving_so15_rotation_sends_a_valid_state_outside_the_cone(self):
        coordinates = basis_vector("XI")
        initial = density_from_answers(coordinates)
        rotation = coordinate_rotation("XI", "IX", math.pi / 4)
        changed = density_from_answers(rotation @ coordinates)
        expected = np.array([(1 - np.sqrt(2)) / 4, 1 / 4, 1 / 4, (1 + np.sqrt(2)) / 4])
        np.testing.assert_allclose(rotation.T @ rotation, np.eye(15), atol=2e-16)
        self.assertAlmostEqual(np.linalg.det(rotation), 1)
        np.testing.assert_array_equal(np.linalg.eigvalsh(initial), [0., 0., 0.5, 0.5])
        np.testing.assert_allclose(np.linalg.eigvalsh(changed), expected, atol=2e-16)
        self.assertAlmostEqual(np.trace(changed @ changed), np.trace(initial @ initial))
        self.assertLess(np.linalg.eigvalsh(changed)[0], -0.1)
        negative_x = np.array([1., -1.]) / np.sqrt(2)
        vector = np.kron(negative_x, negative_x)
        effect = np.outer(vector, vector)
        self.assertAlmostEqual(np.trace(effect @ changed).real, (1 - np.sqrt(2)) / 4)

    def test_fourth_spectral_moment_detects_what_coordinate_norm_misses(self):
        first = observable(basis_vector("XI"))
        second = observable((basis_vector("XI") + basis_vector("IX")) / np.sqrt(2))
        self.assertEqual(np.trace(first @ first), 4)
        self.assertAlmostEqual(np.trace(second @ second).real, 4)
        self.assertEqual(np.trace(np.linalg.matrix_power(first, 4)), 4)
        self.assertAlmostEqual(np.trace(np.linalg.matrix_power(second, 4)).real, 8)

    def test_local_real_disk_and_single_qubit_subalgebras_remain_available(self):
        labels = ("XI", "YI", "ZI")
        indices = [LABELS.index(label) for label in labels]
        for label in labels:
            generator = action_generator(basis_vector(label))
            external = [index for index in range(15) if index not in indices]
            np.testing.assert_array_equal(generator[np.ix_(external, indices)], np.zeros((12, 3)))
        radius = np.sqrt(0.2 ** 2 + 0.7 ** 2)
        real_state = density_from_answers(0.2 * basis_vector("XI") + 0.7 * basis_vector("ZI"))
        np.testing.assert_allclose(real_state.imag, 0, atol=0)
        np.testing.assert_allclose(np.linalg.eigvalsh(real_state),
                                   np.repeat([(1 - radius) / 4, (1 + radius) / 4], 2), atol=1e-16)

    def test_global_real_J_lift_reproduces_the_composite_action_and_records(self):
        vector = np.array([1., 0., 0., 1.]) / np.sqrt(2)
        density = np.outer(vector, vector)
        unitary = pauli_unitary("YX", 0.43) @ pauli_unitary("IZ", 0.27)
        lifted = real_lift(unitary)
        output = lifted @ encode_state(density) @ lifted.T
        expected = unitary @ density @ unitary.conj().T
        np.testing.assert_allclose(lifted.T @ lifted, np.eye(8), atol=5e-16)
        np.testing.assert_allclose(output, encode_state(expected), atol=2e-16)
        for question in QUESTIONS:
            effect = (np.eye(4) + question) / 2
            self.assertAlmostEqual(np.trace(real_lift(effect) @ output).real, np.trace(effect @ expected).real)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(CompositeQuestionSymmetryTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 213,
        "antecedent_rounds": [57, 80, 85, 171, 183, 189, 210, 211, 212],
        "scope": "Scope audit inside the declared standard complex two-qubit theory, not a cognitive derivation",
        "complex_carrier_dimension": 4,
        "traceless_question_dimension": 15,
        "question_labels": list(LABELS),
        "generator_definition": "D(x)y = coefficients(-i[Q(x),Q(y)]/2)",
        "linear_skew_and_own_question_preserving": True,
        "injective_generator_gram_exact": "8 I_15",
        "commutator_closure": "[D(x),D(y)]=D(D(x)y)",
        "covariance_group_used": "unitary-induced adjoint SU(4) rotations (global phase invisible)",
        "full_SO15_covariance": False,
        "commuting_pair_witness": {"fixed": "XI", "before": "IX", "after": "ZI", "defect_on_rotated_input_norm": "1"},
        "cone_witness": {
            "input": "(I+XI)/4",
            "rotated": "(I+(XI+IX)/sqrt(2))/4",
            "same_coordinate_norm": True,
            "same_purity_trace_square": True,
            "minimum_eigenvalue_exact": "(1-sqrt(2))/4",
            "minimum_eigenvalue_diagnostic": (1 - np.sqrt(2)) / 4,
            "negative_effect": "projector onto |-X> tensor |-X>",
            "question_fourth_moments": [4, 8]
        },
        "full_composite_state_space_is_a_ball": False,
        "local_real_disk_and_complex_qubit_interfaces_retained": True,
        "common_J_real_carrier_dimension": 8,
        "round211_all_SO_n_condition_can_be_imposed_on_every_quantum_composite": False,
        "round211_conditional_theorem_refuted": False,
        "weaker_group_covariance_uniquely_selects_quantum_theory": False,
        "cognitive_necessity_or_quantum_incompleteness_proved": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("composite_question_symmetry_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()