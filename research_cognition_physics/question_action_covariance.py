"""Round 211: a conditional dimension selector for question-action coupling.

Classify linear SO(n)-equivariant maps from question vectors to skew
generators. Full symmetry and a nonzero map select n=3, not cognition alone.
Even sign flips give an exact exclusion certificate; plane quarter-turns
and a rank-two integer minor finish n=3. Reflections and a fixed-plane
four-dimensional counterexample expose the assumptions.
"""

import argparse
import json
import unittest
from itertools import combinations
from pathlib import Path

import numpy as np

from question_algebra_orientation import ordered_product, question_family, question_sum


def skew_unit(dimension, row, column):
    matrix = np.zeros((dimension, dimension), dtype=np.int64)
    matrix[row, column] = 1
    matrix[column, row] = -1
    return matrix


def even_sign_flips(dimension):
    signs = []
    for first, second in combinations(range(dimension), 2):
        diagonal = np.ones(dimension, dtype=np.int64)
        diagonal[first] = diagonal[second] = -1
        signs.append(diagonal)
    return tuple(signs)


def surviving_coefficients(dimension):
    signs = even_sign_flips(dimension)
    return [(source, row, column) for source in range(dimension)
            for row, column in combinations(range(dimension), 2)
            if all(diagonal[source] == diagonal[row] * diagonal[column] for diagonal in signs)]


def quarter_turn(dimension, row, column):
    rotation = np.eye(dimension, dtype=np.int64)
    rotation[row, row] = rotation[column, column] = 0
    rotation[row, column] = -1
    rotation[column, row] = 1
    return rotation


def covariance_constraints(dimension, coefficients):
    columns = []
    rotations = [quarter_turn(dimension, row, column)
                 for row, column in combinations(range(dimension), 2)]
    for source, row, column in coefficients:
        generator = skew_unit(dimension, row, column)
        blocks = []
        for rotation in rotations:
            for input_axis in range(dimension):
                transformed_input = rotation[source, input_axis] * generator
                transformed_output = int(source == input_axis) * rotation @ generator @ rotation.T
                blocks.append((transformed_input - transformed_output).reshape(-1))
        columns.append(np.concatenate(blocks))
    return np.array(columns, dtype=np.int64).T


def nonzero_minor_certificate(constraints):
    nonzero = np.unique(constraints[np.any(constraints != 0, axis=1)], axis=0)
    for first, second in combinations(range(len(nonzero)), 2):
        for left, right in combinations(range(nonzero.shape[1]), 2):
            block = nonzero[np.ix_([first, second], [left, right])]
            determinant = int(block[0, 0] * block[1, 1] - block[0, 1] * block[1, 0])
            if determinant:
                return {"matrix": block.tolist(), "determinant": determinant}
    raise AssertionError("Expected a nonzero two-by-two minor.")


def cross_generator(coordinates):
    first, second, third = coordinates
    return np.array([[0, -third, second], [third, 0, -first], [-second, first, 0]])


def fixed_three_plane_generator(coordinates):
    generator = np.zeros((len(coordinates), len(coordinates)))
    generator[:3, :3] = cross_generator(coordinates[:3])
    return generator


def lifted_question_generator(questions, coordinates, scale=1.0):
    derived = ordered_product(questions)
    return -scale * derived @ question_sum(questions, coordinates) / 2


class QuestionActionCovarianceTests(unittest.TestCase):
    def test_even_sign_flips_leave_only_three_dimensional_volume_coefficients(self):
        for dimension in range(1, 11):
            survivors = surviving_coefficients(dimension)
            expected = [(0, 1, 2), (1, 0, 2), (2, 0, 1)] if dimension == 3 else []
            self.assertEqual(survivors, expected)
            for diagonal in even_sign_flips(dimension):
                self.assertEqual(int(np.prod(diagonal)), 1)

    def test_every_eliminated_coefficient_has_a_concrete_even_flip_witness(self):
        for dimension in range(2, 9):
            survivors = set(surviving_coefficients(dimension))
            for source in range(dimension):
                for row, column in combinations(range(dimension), 2):
                    if (source, row, column) in survivors:
                        continue
                    witnesses = [diagonal for diagonal in even_sign_flips(dimension)
                                 if diagonal[source] != diagonal[row] * diagonal[column]]
                    self.assertTrue(witnesses)

    def test_three_dimensional_covariance_has_an_exact_one_dimensional_kernel(self):
        constraints = covariance_constraints(3, surviving_coefficients(3))
        coefficients = np.array([-1, 1, -1], dtype=np.int64)
        np.testing.assert_array_equal(constraints @ coefficients, np.zeros(len(constraints), dtype=int))
        certificate = nonzero_minor_certificate(constraints)
        self.assertNotEqual(certificate["determinant"], 0)
        self.assertEqual(constraints.shape[1], 3)

    def test_surviving_coefficients_reconstruct_cross_product_on_every_basis_axis(self):
        survivors = surviving_coefficients(3)
        coefficients = (-1, 1, -1)
        for axis in range(3):
            generator = sum(value * int(source == axis) * skew_unit(3, row, column)
                            for value, (source, row, column) in zip(coefficients, survivors))
            np.testing.assert_array_equal(generator, cross_generator(np.eye(3, dtype=int)[axis]))

    def test_cross_map_is_linear_skew_and_conserves_its_own_question(self):
        first = np.array([2, -1, 3], dtype=np.int64)
        second = np.array([-1, 4, 2], dtype=np.int64)
        generator = cross_generator(first)
        np.testing.assert_array_equal(generator.T, -generator)
        np.testing.assert_array_equal(generator @ first, np.zeros(3, dtype=int))
        np.testing.assert_array_equal(cross_generator(2 * first - 3 * second),
                                      2 * generator - 3 * cross_generator(second))
        for state in (first, second, np.array([1, 1, 1], dtype=int)):
            self.assertEqual(first @ generator @ state, 0)

    def test_full_so3_covariance_beyond_the_finite_certificate(self):
        random = np.random.default_rng(211)
        for _ in range(20):
            rotation, _ = np.linalg.qr(random.normal(size=(3, 3)))
            if np.linalg.det(rotation) < 0:
                rotation[:, 0] *= -1
            coordinates = random.normal(size=3)
            np.testing.assert_allclose(cross_generator(rotation @ coordinates),
                                       rotation @ cross_generator(coordinates) @ rotation.T, atol=3e-15)

    def test_reflection_covariance_would_kill_the_nonzero_solution(self):
        reflection = np.diag([-1, 1, 1])
        coordinates = np.array([1, 2, 3])
        transformed = cross_generator(reflection @ coordinates)
        conjugated = reflection @ cross_generator(coordinates) @ reflection.T
        np.testing.assert_array_equal(transformed, -conjugated)
        self.assertFalse(np.array_equal(transformed, conjugated))

    def test_a_four_dimensional_fixed_plane_preserves_answers_but_breaks_full_covariance(self):
        coordinates = np.array([1., 2., 3., 4.])
        generator = fixed_three_plane_generator(coordinates)
        np.testing.assert_array_equal(generator.T, -generator)
        np.testing.assert_array_equal(generator @ coordinates, np.zeros(4))
        basis_question = np.array([1., 0., 0., 0.])
        rotation = quarter_turn(4, 0, 3)
        lhs = fixed_three_plane_generator(rotation @ basis_question)
        rhs = rotation @ fixed_three_plane_generator(basis_question) @ rotation.T
        np.testing.assert_array_equal(lhs, np.zeros((4, 4)))
        self.assertAlmostEqual(np.linalg.norm(lhs - rhs), np.sqrt(2))

    def test_the_derived_J_lift_realizes_the_selected_question_action_map(self):
        questions = question_family(3)
        derived = ordered_product(questions)
        for coordinates in (np.array([1, 0, 0]), np.array([2, -1, 3])):
            generator = lifted_question_generator(questions, coordinates)
            observable = question_sum(questions, coordinates)
            np.testing.assert_array_equal(generator.T, -generator)
            np.testing.assert_array_equal(generator @ observable, observable @ generator)
            np.testing.assert_array_equal(generator @ derived, derived @ generator)
            induced = np.array([[np.trace(left @ (generator @ right - right @ generator)) / len(left)
                                 for right in questions] for left in questions])
            np.testing.assert_array_equal(induced, cross_generator(coordinates))

    def test_central_phase_changes_do_not_change_the_answer_generator(self):
        questions = question_family(3)
        derived = ordered_product(questions)
        generator = lifted_question_generator(questions, [1, 2, 3])
        altered = generator + 7 * derived
        self.assertGreater(np.linalg.norm(altered - generator), 1)
        for question in questions:
            np.testing.assert_array_equal(altered @ question - question @ altered,
                                          generator @ question - question @ generator)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(QuestionActionCovarianceTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    constraints = covariance_constraints(3, surviving_coefficients(3))
    report = {
        "round": 211,
        "antecedent_rounds": [183, 189, 204, 209, 210],
        "scope": "Linear SO(n)-equivariant question-to-skew-generator map without extra directional data",
        "candidate_conditions": ["linear dependence on real question vector",
                                 "reversible orthogonal answer-space action",
                                 "full proper-rotation covariance",
                                 "nonzero correspondence",
                                 "own-question expectation conserved for all states"],
        "nonzero_map_exists_iff": "n=3 for n>=1",
        "solution": "D(x) r = kappa * (x cross r), kappa nonzero",
        "scale_or_sign_fixed_by_axioms": False,
        "integer_minor_certificate": nonzero_minor_certificate(constraints),
        "three_survivor_null_vector": [-1, 1, -1],
        "dimension_rows": [{"n": dimension, "linear_coefficients": dimension * dimension * (dimension - 1) // 2,
                            "after_even_flips": len(surviving_coefficients(dimension)),
                            "full_covariance_dimension": int(dimension == 3)} for dimension in range(1, 9)],
        "including_reflections_allows_nonzero_correspondence": False,
        "fixed_three_plane_in_four_dimensions_satisfies_conservation": True,
        "fixed_three_plane_in_four_dimensions_satisfies_full_SO4_covariance": False,
        "real_question_algebra_lift": "A_x = -(kappa/2) J Q(x), modulo a central phase generator",
        "question_to_action_relation_derived_from_SoCA": False,
        "all_conditions_claimed_logically_independent": False,
        "state_dependent_or_reference_dependent_maps_excluded_by_theorem": False,
        "all_cognitive_models_or_all_composite_theories_classified": False,
        "physical_time_or_hbar_derived": False,
        "quantum_necessity_or_incompleteness_proved": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("question_action_covariance_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()