"""Round 212: an internal reference restores joint SO(4) covariance.

A unit real reference t specifies D_t(x) by a four-dimensional volume form.
The rule is linear in x, skew, answer-preserving and jointly covariant in
(x,t). Real Clifford controls implement it. Keeping t in the full state is
not equivalent to a reference-independent map on a four-ball.
"""

import argparse
import json
import math
import unittest
from itertools import combinations, permutations
from pathlib import Path

import numpy as np

from question_action_covariance import cross_generator, quarter_turn
from question_algebra_orientation import question_family, question_sum, state_from_answers


def volume_tensor():
    volume = np.zeros((4, 4, 4, 4), dtype=np.int64)
    for order in permutations(range(4)):
        inversions = sum(order[first] > order[second] for first in range(4)
                         for second in range(first + 1, 4))
        volume[order] = (-1) ** inversions
    return volume


VOLUME = volume_tensor()


def reference_generator(question, reference, scale=1):
    question = np.asarray(question)
    reference = np.asarray(reference)
    if question.shape != (4,) or reference.shape != (4,):
        raise ValueError("Question and reference must each have four coordinates.")
    return -scale * np.einsum("ijkl,k,l->ij", VOLUME, question, reference)


def wedge_generator(question, reference):
    return np.outer(question, reference) - np.outer(reference, question)


def generator_frequency(generator):
    return math.sqrt(float(np.sum(generator * generator)) / 2)


def vector_flow(generator, parameter):
    frequency = generator_frequency(generator)
    first = parameter * np.sinc(parameter * frequency / np.pi)
    second = parameter * parameter * np.sinc(parameter * frequency / (2 * np.pi)) ** 2 / 2
    return np.eye(4) + first * generator + second * (generator @ generator)


def real_spin_generator(generator):
    questions = question_family(4)
    return sum(generator[row, column] * questions[row] @ questions[column] / 2
               for row, column in combinations(range(4), 2))


def real_spin_flow(generator, parameter):
    frequency = generator_frequency(generator)
    spin = real_spin_generator(generator)
    return (math.cos(parameter * frequency / 2) * np.eye(len(spin))
            + parameter * np.sinc(parameter * frequency / (2 * np.pi)) * spin)


def block_diagonal(first, second):
    block_size = len(first)
    result = np.zeros((2 * block_size, 2 * block_size))
    result[:block_size, :block_size] = first
    result[block_size:, block_size:] = second
    return result


def correlation_witness(correlation_sign):
    if correlation_sign not in (-1, 1):
        raise ValueError("Correlation sign must be plus or minus one.")
    questions = question_family(4)
    first_axis = np.eye(4)[0]
    first = state_from_answers(questions, correlation_sign * first_axis) / 2
    second = state_from_answers(questions, -correlation_sign * first_axis) / 2
    return block_diagonal(first, second)


def target_marginal(joint):
    return joint[:8, :8] + joint[8:, 8:]


class InternalReferenceCovarianceTests(unittest.TestCase):
    def test_real_volume_definition_reduces_to_the_old_three_plane_exactly(self):
        reference = np.array([0, 0, 0, 1], dtype=np.int64)
        for question in (np.array([1, 2, 3, 4]), np.array([-2, 0, 1, -1])):
            expected = np.zeros((4, 4), dtype=np.int64)
            expected[:3, :3] = cross_generator(question[:3])
            np.testing.assert_array_equal(reference_generator(question, reference), expected)

    def test_linearity_skewness_and_both_preserved_vectors_are_exact(self):
        question = np.array([2, -1, 3, 1], dtype=np.int64)
        reference = np.array([1, 2, 0, -1], dtype=np.int64)
        other = np.array([-1, 1, 2, 0], dtype=np.int64)
        generator = reference_generator(question, reference)
        np.testing.assert_array_equal(generator.T, -generator)
        np.testing.assert_array_equal(generator @ question, np.zeros(4, dtype=int))
        np.testing.assert_array_equal(generator @ reference, np.zeros(4, dtype=int))
        np.testing.assert_array_equal(reference_generator(3 * question - 2 * other, reference),
                                      3 * generator - 2 * reference_generator(other, reference))
        np.testing.assert_array_equal(reference_generator(question, -reference), -generator)

    def test_joint_covariance_holds_for_all_integer_plane_quarter_turns(self):
        question = np.array([1, 2, -1, 3], dtype=np.int64)
        reference = np.array([-1, 1, 2, 0], dtype=np.int64)
        for row, column in combinations(range(4), 2):
            rotation = quarter_turn(4, row, column)
            np.testing.assert_array_equal(reference_generator(rotation @ question, rotation @ reference),
                                          rotation @ reference_generator(question, reference) @ rotation.T)

    def test_joint_so4_covariance_and_reflection_sign_beyond_integer_examples(self):
        random = np.random.default_rng(212)
        for _ in range(12):
            rotation, _ = np.linalg.qr(random.normal(size=(4, 4)))
            if np.linalg.det(rotation) < 0:
                rotation[:, 0] *= -1
            question = random.normal(size=4)
            reference = random.normal(size=4)
            reference /= np.linalg.norm(reference)
            np.testing.assert_allclose(reference_generator(rotation @ question, rotation @ reference),
                                       rotation @ reference_generator(question, reference) @ rotation.T, atol=3e-15)
        reflection = np.diag([-1, 1, 1, 1])
        np.testing.assert_allclose(reference_generator(reflection @ question, reflection @ reference),
                                   -reflection @ reference_generator(question, reference) @ reflection.T, atol=1e-15)

    def test_unrotated_reference_does_not_satisfy_the_round211_hypothesis(self):
        question, reference = np.eye(4)[0], np.eye(4)[3]
        rotation = quarter_turn(4, 0, 3)
        changed_question_only = reference_generator(rotation @ question, reference)
        changed_rule = rotation @ reference_generator(question, reference) @ rotation.T
        np.testing.assert_array_equal(changed_question_only, np.zeros((4, 4)))
        self.assertAlmostEqual(np.linalg.norm(changed_question_only - changed_rule), np.sqrt(2))

    def test_frequency_cubic_identity_and_aligned_question_degeneracy(self):
        for question, reference in ((np.array([1, 2, 3, 4]), np.array([1, 0, -1, 2])),
                                    (np.array([0, 0, 0, 2]), np.array([0, 0, 0, 1]))):
            generator = reference_generator(question, reference)
            squared = int(question @ question) * int(reference @ reference) - int(question @ reference) ** 2
            np.testing.assert_array_equal(generator @ generator @ generator, -squared * generator)
            self.assertEqual(int(np.sum(generator * generator)), 2 * squared)
        np.testing.assert_array_equal(reference_generator(reference, reference), np.zeros((4, 4)))

    def test_wedge_term_is_covariant_but_cannot_preserve_its_question(self):
        question = np.array([1., 2., 0., 3.])
        reference = np.eye(4)[3]
        wedge = wedge_generator(question, reference)
        expected = question * np.dot(reference, question) - reference * np.dot(question, question)
        np.testing.assert_array_equal(wedge @ question, expected)
        self.assertGreater(np.linalg.norm(expected), 1)
        for row, column in combinations(range(4), 2):
            rotation = quarter_turn(4, row, column)
            np.testing.assert_array_equal(wedge_generator(rotation @ question, rotation @ reference),
                                          rotation @ wedge @ rotation.T)

    def test_vector_and_real_spin_flows_match_all_four_answers_and_are_reversible(self):
        questions = question_family(4)
        random = np.random.default_rng(1212)
        for _ in range(12):
            question = random.normal(size=4)
            reference = random.normal(size=4)
            reference /= np.linalg.norm(reference)
            coordinates = random.normal(size=4)
            coordinates /= 1 + np.linalg.norm(coordinates)
            parameter = 0.71
            generator = reference_generator(question, reference)
            rotation = vector_flow(generator, parameter)
            spin = real_spin_generator(generator)
            gate = real_spin_flow(generator, parameter)
            frequency = generator_frequency(generator)
            np.testing.assert_allclose(spin @ spin, -(frequency ** 2) * np.eye(8) / 4, atol=2e-15)
            np.testing.assert_allclose(gate.T @ gate, np.eye(8), atol=2e-15)
            np.testing.assert_allclose(rotation.T @ rotation, np.eye(4), atol=2e-15)
            np.testing.assert_allclose(vector_flow(generator, -parameter) @ rotation, np.eye(4), atol=2e-15)
            density = state_from_answers(questions, coordinates)
            output = gate @ density @ gate.T
            np.testing.assert_allclose(output, state_from_answers(questions, rotation @ coordinates), atol=3e-16)
            observable = question_sum(questions, question)
            np.testing.assert_allclose(gate.T @ observable @ gate, observable, atol=3e-15)

    def test_same_reference_and_target_marginals_have_opposite_later_answers(self):
        positive, negative = correlation_witness(1), correlation_witness(-1)
        np.testing.assert_array_equal(target_marginal(positive), target_marginal(negative))
        self.assertEqual(np.trace(positive[:8, :8]), 0.5)
        self.assertEqual(np.trace(negative[:8, :8]), 0.5)
        self.assertGreaterEqual(np.linalg.eigvalsh(positive).min(), 0)
        self.assertGreaterEqual(np.linalg.eigvalsh(negative).min(), 0)
        question, reference = np.eye(4)[2], np.eye(4)[3]
        first = real_spin_flow(reference_generator(question, reference), math.pi / 2)
        second = real_spin_flow(reference_generator(question, -reference), math.pi / 2)
        gate = block_diagonal(first, second)
        effect = (np.eye(8) + question_family(4)[1]) / 2
        positive_probability = np.trace(effect @ target_marginal(gate @ positive @ gate.T))
        negative_probability = np.trace(effect @ target_marginal(gate @ negative @ gate.T))
        self.assertAlmostEqual(positive_probability, 1)
        self.assertAlmostEqual(negative_probability, 0)
        for joint in (positive, negative):
            output = gate @ joint @ gate.T
            np.testing.assert_allclose(gate.T @ output @ gate, joint, atol=2e-16)

    def test_complete_flag_control_is_positive_linear_and_keeps_external_correlations(self):
        question, reference = np.eye(4)[2], np.eye(4)[3]
        first = real_spin_flow(reference_generator(question, reference), 0.43)
        second = real_spin_flow(reference_generator(question, -reference), 0.43)
        gate = block_diagonal(first, second)
        expanded = np.kron(gate, np.eye(2))
        random = np.random.default_rng(2212)
        raw = random.normal(size=(32, 32))
        joint = raw @ raw.T
        joint /= np.trace(joint)
        output = expanded @ joint @ expanded.T
        self.assertGreater(np.linalg.eigvalsh(output).min(), 0)
        np.testing.assert_allclose(expanded.T @ output @ expanded, joint, atol=2e-16)
        initial_external = np.trace(joint.reshape(16, 2, 16, 2), axis1=0, axis2=2)
        final_external = np.trace(output.reshape(16, 2, 16, 2), axis1=0, axis2=2)
        np.testing.assert_allclose(final_external, initial_external, atol=3e-16)
        other = np.eye(32) / 32
        mixed = 0.3 * joint + 0.7 * other
        np.testing.assert_allclose(expanded @ mixed @ expanded.T,
                                   0.3 * output + 0.7 * expanded @ other @ expanded.T, atol=2e-16)

    def test_arbitrary_coherent_reference_marginal_need_not_be_preserved(self):
        question, reference = np.eye(4)[2], np.eye(4)[3]
        first = real_spin_flow(reference_generator(question, reference), math.pi / 2)
        second = real_spin_flow(reference_generator(question, -reference), math.pi / 2)
        gate = block_diagonal(first, second)
        coherent_reference = np.ones((2, 2)) / 2
        initial = np.kron(coherent_reference, np.eye(8) / 8)
        output = gate @ initial @ gate.T
        reduced = np.trace(output.reshape(2, 8, 2, 8), axis1=1, axis2=3)
        np.testing.assert_allclose(reduced, np.eye(2) / 2, atol=2e-16)
        self.assertAlmostEqual(np.linalg.norm(reduced - coherent_reference, ord=2), 0.5)
        np.testing.assert_allclose(gate.T @ output @ gate, initial, atol=2e-16)

    def test_zero_mean_reference_does_not_imply_zero_finite_action(self):
        question, reference = np.eye(4)[2], np.eye(4)[3]
        generator = reference_generator(question, reference)
        np.testing.assert_array_equal((generator + reference_generator(question, -reference)) / 2,
                                      np.zeros((4, 4)))
        parameter = 0.7
        averaged = (vector_flow(generator, parameter) + vector_flow(-generator, parameter)) / 2
        np.testing.assert_allclose(averaged, np.diag([math.cos(parameter), math.cos(parameter), 1, 1]), atol=2e-16)
        self.assertGreater(np.linalg.norm(averaged - np.eye(4)), 0.3)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(InternalReferenceCovarianceTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 212,
        "antecedent_rounds": [1, 80, 85, 106, 175, 176, 209, 210, 211],
        "scope": "SO(4)-covariant real question controls with an internal unit-vector reference",
        "generator": "D_t(x)_ij = -kappa * sum_kl epsilon_ijkl x_k t_l",
        "classification": "with unit t as the only extra directional datum, linear x and SO(4) covariance allow wedge plus dual wedge; own-question preservation leaves only dual wedge",
        "joint_covariance": "D_(Rt)(Rx)=R D_t(x) R^T",
        "reference_independent_SO4_covariance_holds": False,
        "D_x_zero_and_D_t_zero": True,
        "nonzero_for_every_nonzero_question_at_fixed_reference": False,
        "squared_frequency": "kappa^2*(norm(x)^2*norm(t)^2-dot(x,t)^2)",
        "real_spin_generator": "one_half * sum_(i<j) D_ij Q_i Q_j",
        "real_target_dimension": 8,
        "two_reference_value_control_dimension": 16,
        "continuous_all_reference_model_is_finite_memory_hardware": False,
        "equal_marginal_correlation_witness": {
            "reference_values": ["+e4", "-e4"],
            "reference_probabilities": ["1/2", "1/2"],
            "target_marginal": "I_8/8 in both preparations",
            "fixed_control_question": "e3",
            "parameter": "pi/2",
            "final_e2_positive_probabilities": ["1", "0"],
            "minimum_same_marginal_predictor_worst_probability_error": "1/2"
        },
        "classical_reference_distribution_preserved_but_joint_correlations_may_change": True,
        "arbitrary_coherent_reference_marginal_preserved": False,
        "zero_mean_reference_can_have_nontrivial_average_finite_action": True,
        "reference_internalization_implies_reference_independence": False,
        "full_finite_control_extension_is_real_orthogonal": True,
        "closed_autonomous_energy_and_clock_implementation_proved": False,
        "cognitive_necessity_or_quantum_incompleteness_proved": False,
        "round211_theorem_refuted": False,
        "new_composite_axioms_derived": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("internal_reference_covariance_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()