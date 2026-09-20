"""Round 210: reconstruct an intrinsic J from a declared sharp-question family.

All construction matrices are real. Three pairwise anticommuting symmetric
involutions generate a central orthogonal complex structure; two and four
questions give different centers. Matrix probability, full-background
sequential neutrality and permission completeness are hypotheses, not
derived cognitive laws.
"""

import argparse
import json
import math
import unittest
from itertools import permutations
from pathlib import Path

import numpy as np

from ancillary_operation_equivalence import J
from real_interface_record_bound import shared_memory_state


REAL_X = np.array([[0, 1], [1, 0]], dtype=np.int64)
REAL_Z = np.diag([1, -1]).astype(np.int64)
REAL_J = J.astype(np.int64)


def tensor_factors(*factors):
    output = np.ones((1, 1), dtype=np.int64)
    for factor in factors:
        output = np.kron(output, factor)
    return output


def question_family(count, multiplicity=1):
    if not isinstance(multiplicity, int) or multiplicity < 1:
        raise ValueError("Multiplicity must be a positive integer.")
    identity = np.eye(2, dtype=np.int64)
    if count == 2:
        questions = (REAL_Z, REAL_X)
    elif count == 3:
        questions = (np.kron(REAL_Z, identity), np.kron(REAL_X, identity),
                     np.kron(REAL_J, REAL_J))
    elif count == 4:
        questions = (tensor_factors(REAL_Z, identity, identity),
                     tensor_factors(REAL_X, REAL_Z, identity),
                     tensor_factors(REAL_X, REAL_X, REAL_Z),
                     tensor_factors(REAL_X, REAL_X, REAL_X))
    else:
        raise ValueError("This certificate implements two, three or four questions.")
    copies = np.eye(multiplicity, dtype=np.int64)
    return tuple(np.kron(question, copies) for question in questions)


def ordered_product(matrices):
    output = np.eye(len(matrices[0]), dtype=np.int64)
    for matrix in matrices:
        output = output @ matrix
    return output


def word_basis(questions):
    identity = np.eye(len(questions[0]), dtype=np.int64)
    words = []
    for mask in range(1 << len(questions)):
        word = identity.copy()
        for index, question in enumerate(questions):
            if mask & (1 << index):
                word = word @ question
        words.append(word)
    return tuple(words)


def central_word_indices(questions):
    return [index for index, word in enumerate(word_basis(questions))
            if all(np.array_equal(word @ question, question @ word) for question in questions)]


def symmetric_word_indices(questions):
    return [index for index, word in enumerate(word_basis(questions))
            if np.array_equal(word, word.T)]


def question_sum(questions, coefficients):
    if len(questions) != len(coefficients):
        raise ValueError("There must be one coefficient per question.")
    return sum(coefficient * question for coefficient, question in zip(coefficients, questions))


def state_from_answers(questions, coordinates):
    dimension = len(questions[0])
    return (np.eye(dimension) + question_sum(questions, coordinates)) / dimension


def pair_rotation(first, second, angle):
    return (math.cos(angle / 2) * np.eye(len(first))
            - math.sin(angle / 2) * (first @ second))


class QuestionAlgebraOrientationTests(unittest.TestCase):
    def test_uniform_sequential_neutrality_is_an_exact_projector_identity(self):
        for count in (2, 3, 4):
            questions = question_family(count)
            identity = np.eye(len(questions[0]), dtype=np.int64)
            for first_index, first in enumerate(questions):
                for second_index, second in enumerate(questions):
                    if first_index == second_index:
                        continue
                    for first_sign in (-1, 1):
                        unnormalized_projector = identity + first_sign * first
                        np.testing.assert_array_equal(unnormalized_projector @ second @ unnormalized_projector,
                                                      np.zeros_like(identity))
                        projector = unnormalized_projector / 2
                        for second_sign in (-1, 1):
                            effect = (identity + second_sign * second) / 2
                            np.testing.assert_array_equal(projector @ effect @ projector, projector / 2)

    def test_promised_preparations_do_not_certify_full_background_neutrality(self):
        first = np.kron(np.eye(2), REAL_X)
        second = np.kron(REAL_Z, (REAL_X + REAL_Z) / np.sqrt(2))
        projector = (np.eye(4) + first) / 2
        effect = (np.eye(4) + second) / 2
        promised = shared_memory_state(0, 1)
        hidden_flag_known = np.kron(np.diag([1., 0.]), (np.eye(2) + REAL_X) / 2)
        self.assertAlmostEqual(np.trace(effect @ projector @ promised @ projector), 0.5)
        biased = np.trace(effect @ projector @ hidden_flag_known @ projector)
        self.assertAlmostEqual(biased, (1 + 1 / np.sqrt(2)) / 2)
        self.assertGreater(np.linalg.norm(first @ second + second @ first), 1)

    def test_sharp_midpoint_identity_and_a_commuting_counterexample(self):
        for count in (2, 3, 4):
            questions = question_family(count)
            identity = np.eye(len(questions[0]), dtype=np.int64)
            for first_index, first in enumerate(questions):
                np.testing.assert_array_equal(first.T, first)
                np.testing.assert_array_equal(first @ first, identity)
                for second in questions[first_index + 1:]:
                    np.testing.assert_array_equal((first + second) @ (first + second), 2 * identity)
                    np.testing.assert_array_equal(first @ second + second @ first, 0 * identity)
        first = np.kron(REAL_Z, np.eye(2, dtype=np.int64))
        second = np.kron(np.eye(2, dtype=np.int64), REAL_Z)
        self.assertFalse(np.array_equal((first + second) @ (first + second), 2 * np.eye(4)))

    def test_three_question_product_is_a_central_orthogonal_complex_structure(self):
        for multiplicity in (1, 2, 3):
            questions = question_family(3, multiplicity)
            derived = ordered_product(questions)
            identity = np.eye(len(derived), dtype=np.int64)
            np.testing.assert_array_equal(derived.T, -derived)
            np.testing.assert_array_equal(derived @ derived, -identity)
            np.testing.assert_array_equal(derived.T @ derived, identity)
            for question in questions:
                np.testing.assert_array_equal(derived @ question, question @ derived)

    def test_canonical_block_certificate_requires_an_even_half_dimension(self):
        for multiplicity in (1, 2, 3):
            questions = question_family(3, multiplicity)
            half = len(questions[0]) // 2
            block = questions[2][:half, half:]
            np.testing.assert_array_equal(block.T, -block)
            np.testing.assert_array_equal(block @ block, -np.eye(half, dtype=np.int64))
            self.assertEqual(half, 2 * multiplicity)
            np.testing.assert_array_equal(ordered_product(questions), np.kron(np.eye(2, dtype=int), -block))

    def test_word_gram_matrices_certify_independence_and_centers_exactly(self):
        expected = {2: (4, [0], [0, 1, 2]),
                    3: (8, [0, 7], [0, 1, 2, 4]),
                    4: (16, [0], [0, 1, 2, 4, 8, 15])}
        for count, (size, center, symmetric) in expected.items():
            questions = question_family(count)
            words = word_basis(questions)
            dimension = len(questions[0])
            gram = np.array([[np.sum(left * right) for right in words] for left in words])
            np.testing.assert_array_equal(gram, dimension * np.eye(size, dtype=np.int64))
            self.assertEqual(central_word_indices(questions), center)
            self.assertEqual(symmetric_word_indices(questions), symmetric)
            for word in words:
                for question in questions:
                    product = word @ question
                    matches = [np.array_equal(product, candidate) or np.array_equal(product, -candidate)
                               for candidate in words]
                    self.assertEqual(sum(matches), 1)

    def test_order_only_changes_orientation_sign_and_real_basis_changes_are_covariant(self):
        questions = question_family(3)
        derived = ordered_product(questions)
        for order in permutations(range(3)):
            inversions = sum(order[first] > order[second] for first in range(3)
                             for second in range(first + 1, 3))
            np.testing.assert_array_equal(ordered_product(tuple(questions[index] for index in order)),
                                          (-1) ** inversions * derived)
        generator = np.random.default_rng(210)
        orthogonal, _ = np.linalg.qr(generator.normal(size=(4, 4)))
        changed = tuple(orthogonal @ question @ orthogonal.T for question in questions)
        np.testing.assert_allclose(ordered_product(changed), orthogonal @ derived @ orthogonal.T, atol=1e-15)

    def test_pair_controls_preserve_the_derived_structure_and_rotate_questions(self):
        questions = question_family(3)
        derived = ordered_product(questions)
        for first_index in range(3):
            for second_index in range(first_index + 1, 3):
                first, second = questions[first_index], questions[second_index]
                for angle in (-1.2, 0., 0.7, math.pi / 2):
                    gate = pair_rotation(first, second, angle)
                    np.testing.assert_allclose(gate.T @ gate, np.eye(4), atol=2e-16)
                    np.testing.assert_array_equal(gate @ derived, derived @ gate)
                    np.testing.assert_allclose(gate @ first @ gate.T,
                                               math.cos(angle) * first + math.sin(angle) * second, atol=3e-16)

    def test_answer_ball_positivity_and_projector_probabilities(self):
        questions = question_family(3)
        generator = np.random.default_rng(1210)
        for _ in range(12):
            coordinates = generator.normal(size=3)
            coordinates /= 1 + np.linalg.norm(coordinates)
            axis = generator.normal(size=3)
            axis /= np.linalg.norm(axis)
            density = state_from_answers(questions, coordinates)
            expected_spectrum = np.repeat([(1 - np.linalg.norm(coordinates)) / 4,
                                          (1 + np.linalg.norm(coordinates)) / 4], 2)
            np.testing.assert_allclose(np.linalg.eigvalsh(density), expected_spectrum, atol=2e-16)
            for sign in (-1, 1):
                projector = (np.eye(4) + sign * question_sum(questions, axis)) / 2
                probability = (1 + sign * np.dot(axis, coordinates)) / 2
                self.assertAlmostEqual(np.trace(projector @ density), probability)
                output = projector @ density @ projector
                np.testing.assert_allclose(output, probability * state_from_answers(questions, sign * axis),
                                           atol=2e-16)

    def test_arbitrary_real_density_has_a_valid_answer_ball_quotient(self):
        questions = question_family(3, 2)
        generator = np.random.default_rng(2210)
        for _ in range(12):
            raw = generator.normal(size=(8, 8))
            density = raw @ raw.T
            density /= np.trace(density)
            coordinates = np.array([np.trace(density @ question) for question in questions])
            self.assertLessEqual(np.linalg.norm(coordinates), 1)
            representative = state_from_answers(questions, coordinates)
            for question in questions:
                self.assertAlmostEqual(np.trace(density @ question), np.trace(representative @ question))

    def test_generated_instruments_preserve_the_quotient_for_adaptive_records(self):
        questions = question_family(3)
        generator = np.random.default_rng(3210)
        raw = generator.normal(size=(4, 4))
        density = raw @ raw.T
        density /= np.trace(density)
        coordinates = [np.trace(density @ question) for question in questions]
        representative = state_from_answers(questions, coordinates)
        self.assertGreater(np.linalg.norm(density - representative), 0.01)
        branches = [(density, representative, 0)]
        for depth in range(3):
            next_branches = []
            for actual, reduced, history in branches:
                index = (history + depth) % 3
                positive = (np.eye(4) + questions[index]) / 2
                negative = np.eye(4) - positive
                gate = pair_rotation(questions[index], questions[(index + 1) % 3], 0.37)
                for result in (0, 1):
                    coefficient = 0.3 if result == 0 else 0.7
                    operation = (math.sqrt(coefficient) * positive
                                 + math.sqrt(1 - coefficient) * negative) @ gate
                    new_actual = operation @ actual @ operation.T
                    new_reduced = operation @ reduced @ operation.T
                    self.assertAlmostEqual(np.trace(new_actual), np.trace(new_reduced))
                    for question in questions:
                        self.assertAlmostEqual(np.trace(new_actual @ question), np.trace(new_reduced @ question))
                    next_branches.append((new_actual, new_reduced, 2 * history + result))
            branches = next_branches
        self.assertAlmostEqual(sum(np.trace(actual) for actual, _, _ in branches), 1)

    def test_two_question_interface_hides_the_third_until_cross_interface_control(self):
        questions = question_family(3)
        positive = state_from_answers(questions, [0., 0., 1.])
        negative = state_from_answers(questions, [0., 0., -1.])
        pair_words = word_basis(questions[:2])
        for index in symmetric_word_indices(questions[:2]):
            self.assertEqual(np.trace(pair_words[index] @ (positive - negative)), 0)
        gate = pair_rotation(questions[0], questions[2], -math.pi / 2)
        effect = (np.eye(4) + questions[0]) / 2
        self.assertAlmostEqual(np.trace(effect @ gate @ positive @ gate.T), 1)
        self.assertAlmostEqual(np.trace(effect @ gate @ negative @ gate.T), 0)

    def test_fourth_question_breaks_intrinsic_orientation_but_keeps_pair_sharpness(self):
        questions = question_family(4)
        old_orientation = ordered_product(questions[:3])
        fourth = questions[3]
        np.testing.assert_array_equal(old_orientation @ fourth, -fourth @ old_orientation)
        self.assertEqual(np.linalg.norm(old_orientation @ fourth - fourth @ old_orientation, ord=2), 2)
        words = word_basis(questions)
        self.assertEqual(central_word_indices(questions), [0])
        self.assertTrue(all(np.array_equal(words[index].T, words[index])
                            for index in central_word_indices(questions)))

    def test_four_question_associative_closure_exposes_a_fifth_sharp_direction(self):
        questions = question_family(4)
        fifth = ordered_product(questions)
        np.testing.assert_array_equal(fifth.T, fifth)
        np.testing.assert_array_equal(fifth @ fifth, np.eye(8, dtype=int))
        for question in questions:
            np.testing.assert_array_equal(question @ fifth + fifth @ question, np.zeros((8, 8), dtype=int))
        coordinates = np.array([1., 2., -1., 3., 2.]) / np.sqrt(19)
        observable = question_sum(questions + (fifth,), coordinates)
        np.testing.assert_allclose(observable @ observable, np.eye(8), atol=3e-16)

    def test_global_real_controls_need_not_preserve_the_selected_orientation(self):
        questions = question_family(3)
        derived = ordered_product(questions)
        real_gate = np.kron(np.eye(2, dtype=int), REAL_Z)
        np.testing.assert_array_equal(real_gate @ derived @ real_gate.T, -derived)
        np.testing.assert_array_equal(real_gate @ questions[0] @ real_gate.T, questions[0])
        np.testing.assert_array_equal(real_gate @ questions[1] @ real_gate.T, questions[1])
        np.testing.assert_array_equal(real_gate @ questions[2] @ real_gate.T, -questions[2])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(QuestionAlgebraOrientationTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 210,
        "antecedent_rounds": [80, 85, 87, 88, 99, 183, 189, 208, 209],
        "scope": "Finite real matrix candidate; full-background sequential neutrality of sharp questions",
        "inputs": ["real symmetric effects, positive states and trace probabilities",
               "symmetric involutions Q_i with selective projector updates",
               "each other question unbiased after either result, for all real input density matrices",
                   "three directions for the central-J theorem",
                   "generated-algebra permission closure if claiming a complete complex-equivalent interface"],
        "neutrality_identity": "P_i,s Q_j P_i,s=0 for both s and i!=j",
        "anticommutation_and_sharp_circle_follow_from_neutrality": True,
        "promised_preparations_suffice_for_full_background_neutrality": False,
        "old_209_known_flag_followup_positive_probability": "(1+1/sqrt(2))/2",
        "common_complex_structure_assumed_in_theorem": False,
        "derived_J": "Q_1 Q_2 Q_3",
        "J_transpose": "-J",
        "J_squared": "-I",
        "J_commutes_with_generated_three_question_algebra": True,
        "nonzero_real_representation_dimension": "multiple of 4",
        "three_question_algebra": "M_2(C) viewed as a real algebra; statement proved via real block form",
        "three_question_symmetric_dimension": 4,
        "answer_state_region": "unit three-dimensional ball, within the declared full matrix candidate",
        "generated_algebra_Kraus_instruments_preserve_all_finite_quotient_records": True,
        "center_rows": [{"questions": count, "real_carrier_dimension": len(question_family(count)[0]),
                         "algebra_dimension": len(word_basis(question_family(count))),
                         "symmetric_dimension": len(symmetric_word_indices(question_family(count))),
                         "center_basis_masks": central_word_indices(question_family(count))}
                        for count in (2, 3, 4)],
        "fourth_question_anticommutes_with_old_J": True,
        "four_question_algebra_has_a_central_J": False,
        "four_question_algebra_has_no_commuting_J_outside_it": "not asserted",
        "four_question_symmetric_closure_extra_question": "Q_1 Q_2 Q_3 Q_4",
        "full_real_controls_automatically_preserve_selected_J": False,
        "three_directions_or_sharpness_derived_from_SoCA": False,
        "all_composite_systems_or_quantum_incompleteness_derived": False,
        "physical_time_and_resource_costs_derived": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("question_algebra_orientation_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()