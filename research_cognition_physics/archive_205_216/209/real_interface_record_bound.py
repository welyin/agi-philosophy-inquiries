"""Round 209: a sharp finite-table bound and shared-memory counterexamples.

Six preparations and three binary measurements are compared at fixed real
Hilbert dimension two without source-detector shared side information.
The analytic probability error is visibility/6. One shared random bit or
a four-dimensional real carrier removes this finite-table obstruction.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from common_orientation_structure import encode_state
from complex_control_from_reference import real_lift
from two_real_form_control import flow, imaginary_generator, real_generator


IDENTITY = np.eye(2)
PAULI_X = np.array([[0., 1.], [1., 0.]])
PAULI_Y = np.array([[0., -1j], [1j, 0.]])
PAULI_Z = np.diag([1., -1.])
PAULI_AXES = (PAULI_X, PAULI_Y, PAULI_Z)
SIGNS = (1, -1)


def target_table(visibility=1.0):
    return np.array([[[0.5 * (1 + sign * visibility * int(measurement == preparation))
                       for sign in SIGNS] for preparation in range(3)]
                     for measurement in range(3)])


def contrast_matrix(table):
    return table[:, :, 0] - table[:, :, 1]


def real_density(vector, sign=1):
    return (IDENTITY + sign * (vector[0] * PAULI_X + vector[1] * PAULI_Z)) / 2


def real_table(preparation_vectors, measurement_vectors):
    return np.array([[[np.trace(real_density(measurement_vector)
                               @ real_density(preparation_vector, sign)).real
                       for sign in SIGNS] for preparation_vector in preparation_vectors]
                     for measurement_vector in measurement_vectors])


def optimal_real_table(visibility=1.0):
    directions = np.array([[1., 0.], [-0.5, np.sqrt(3) / 2], [-0.5, -np.sqrt(3) / 2]])
    vectors = np.sqrt(2 * visibility / 3) * directions
    return real_table(vectors, vectors)


def shared_bit_table(visibility=1.0):
    tables = []
    for shared_sign in SIGNS:
        vectors = np.array([[1., 0.], [0., 1.],
                            [shared_sign / np.sqrt(2), shared_sign / np.sqrt(2)]])
        tables.append(real_table(vectors, visibility * vectors))
    return sum(tables) / 2


def classical_shared_table():
    references = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
    return [[[sum(Fraction(int(sign * reference[preparation] * reference[measurement] == 1), 4)
                  for reference in references)
              for sign in SIGNS] for preparation in range(3)] for measurement in range(3)]


def shared_bit_luders_sequence():
    initial = real_density(np.array([1., 0.]))
    records = np.zeros((2, 2))
    for shared_sign in SIGNS:
        intermediate_axis = np.array([shared_sign, shared_sign]) / np.sqrt(2)
        for first_index, first_sign in enumerate(SIGNS):
            projector = real_density(intermediate_axis, first_sign)
            branch = projector @ initial @ projector
            for final_index, final_sign in enumerate(SIGNS):
                final_effect = real_density(np.array([1., 0.]), final_sign)
                records[first_index, final_index] += np.trace(final_effect @ branch).real / 2
    return records


def shared_memory_state(preparation, sign):
    joint = np.zeros((4, 4))
    for memory_index, shared_sign in enumerate(SIGNS):
        vectors = np.array([[1., 0.], [0., 1.],
                            [shared_sign / np.sqrt(2), shared_sign / np.sqrt(2)]])
        block = slice(2 * memory_index, 2 * memory_index + 2)
        joint[block, block] = real_density(vectors[preparation], sign) / 2
    return joint


def refreshed_memory_kraus(measurement, outcome):
    kraus = []
    for old_index, old_sign in enumerate(SIGNS):
        old_vectors = np.array([[1., 0.], [0., 1.], [old_sign, old_sign]])
        old_vectors[2] /= np.sqrt(2)
        _, old_eigenvectors = np.linalg.eigh(real_density(old_vectors[measurement], outcome))
        for new_index, new_sign in enumerate(SIGNS):
            new_vectors = np.array([[1., 0.], [0., 1.], [new_sign, new_sign]])
            new_vectors[2] /= np.sqrt(2)
            _, new_eigenvectors = np.linalg.eigh(real_density(new_vectors[measurement], outcome))
            operation = np.zeros((4, 4))
            operation[2 * new_index:2 * new_index + 2, 2 * old_index:2 * old_index + 2] = (
                np.outer(new_eigenvectors[:, -1], old_eigenvectors[:, -1]) / np.sqrt(2))
            kraus.append(operation)
    return kraus


def pulse_gates():
    real_axis = real_generator(2, 0, 1)
    second_axis = imaginary_generator(2, 0, 1)
    return ((flow(real_axis, -np.pi / 4), flow(real_axis, np.pi / 4)),
            (flow(second_axis, np.pi / 4), flow(second_axis, -np.pi / 4)),
            (IDENTITY, flow(real_axis, np.pi / 2)))


def sample_budget(visibility=1.0, risk=0.01):
    tolerance = visibility / 24
    repetitions = math.ceil(math.log(36 / risk) / (2 * tolerance * tolerance))
    return {"visibility": visibility, "risk": risk, "probability_tolerance": tolerance,
            "repetitions_per_cell": repetitions, "total_fresh_preparations": 18 * repetitions,
            "union_bound": 36 * math.exp(-2 * repetitions * tolerance * tolerance)}


class RealInterfaceRecordBoundTests(unittest.TestCase):
    def test_pulse_preparations_and_inverse_readout_reproduce_all_records(self):
        initial = np.diag([1., 0.])
        gates = pulse_gates()
        actual = np.zeros((3, 3, 2))
        for measurement in range(3):
            readout = gates[measurement][0].conj().T
            for preparation in range(3):
                for sign_index, sign in enumerate(SIGNS):
                    gate = gates[preparation][sign_index]
                    density = gate @ initial @ gate.conj().T
                    np.testing.assert_allclose(density, (IDENTITY + sign * PAULI_AXES[preparation]) / 2,
                                               atol=5e-16)
                    actual[measurement, preparation, sign_index] = np.trace(
                        initial @ readout @ density @ readout.conj().T).real
        np.testing.assert_allclose(actual, target_table(), atol=1e-15)

    def test_all_target_preparation_pairs_have_the_same_average(self):
        for axis in PAULI_AXES:
            pair = [(IDENTITY + sign * axis) / 2 for sign in SIGNS]
            np.testing.assert_array_equal(sum(pair) / 2, IDENTITY / 2)

    def test_arbitrary_biased_real_effects_factor_through_two_coordinates(self):
        generator = np.random.default_rng(209)
        for _ in range(12):
            vectors = generator.normal(size=(3, 2, 2))
            vectors /= 1 + np.linalg.norm(vectors, axis=2, keepdims=True)
            effects = generator.normal(size=(3, 2))
            effects *= 0.2 / (1 + np.linalg.norm(effects, axis=1, keepdims=True))
            biases = np.array([0.3, 0.5, 0.7])
            table = np.zeros((3, 3, 2))
            for measurement in range(3):
                effect = (biases[measurement] * IDENTITY + effects[measurement, 0] * PAULI_X
                          + effects[measurement, 1] * PAULI_Z)
                self.assertGreater(np.linalg.eigvalsh(effect).min(), 0)
                self.assertLess(np.linalg.eigvalsh(effect).max(), 1)
                for preparation in range(3):
                    for sign_index in range(2):
                        table[measurement, preparation, sign_index] = np.trace(
                            effect @ real_density(vectors[preparation, sign_index])).real
            expected = effects @ (vectors[:, 0, :] - vectors[:, 1, :]).T
            np.testing.assert_allclose(contrast_matrix(table), expected, atol=3e-16)
            self.assertLess(np.linalg.svd(expected, compute_uv=False)[-1], 1e-15)
            self.assertGreaterEqual(np.max(np.abs(table - target_table())), 1 / 6 - 1e-14)

    def test_optimal_real_construction_attains_the_exact_error(self):
        for visibility in (0., 0.07, 0.5, 0.9, 1.):
            actual = optimal_real_table(visibility)
            expected_contrast = visibility * (np.eye(3) - np.ones((3, 3)) / 3)
            np.testing.assert_allclose(contrast_matrix(actual), expected_contrast, atol=3e-16)
            np.testing.assert_allclose(np.abs(actual - target_table(visibility)),
                                       np.full((3, 3, 2), visibility / 6), atol=2e-16)
            self.assertGreaterEqual(actual.min(), 0)
            self.assertLessEqual(actual.max(), 1)

    def test_rational_probability_certificate_covers_every_table_cell(self):
        for visibility in (Fraction(0), Fraction(1, 7), Fraction(1, 2), Fraction(1)):
            for measurement in range(3):
                for preparation in range(3):
                    diagonal = int(measurement == preparation)
                    for sign in SIGNS:
                        ideal = (1 + sign * visibility * diagonal) / 2
                        achieved = (1 + sign * visibility * (diagonal - Fraction(1, 3))) / 2
                        self.assertEqual(abs(ideal - achieved), visibility / 6)

    def test_shared_random_bit_removes_the_obstruction(self):
        for visibility in (0., 0.2, 0.7, 1.):
            np.testing.assert_allclose(shared_bit_table(visibility), target_table(visibility), atol=2e-16)
        self.assertEqual(np.linalg.matrix_rank(contrast_matrix(shared_bit_table())), 3)

    def test_one_classical_message_bit_and_two_shared_bits_also_match(self):
        exact = classical_shared_table()
        for measurement in range(3):
            for preparation in range(3):
                for sign_index, sign in enumerate(SIGNS):
                    expected = Fraction(1 + sign * int(measurement == preparation), 2)
                    self.assertEqual(exact[measurement][preparation][sign_index], expected)

    def test_fixed_shared_bit_luders_updates_fail_a_two_measurement_record(self):
        actual = shared_bit_luders_sequence()
        np.testing.assert_allclose(actual, [[3 / 8, 1 / 8], [3 / 8, 1 / 8]], atol=2e-16)
        initial = (IDENTITY + PAULI_X) / 2
        expected = np.zeros((2, 2))
        for first_index, first_sign in enumerate(SIGNS):
            projector = (IDENTITY + first_sign * PAULI_Z) / 2
            branch = projector @ initial @ projector
            for final_index, final_sign in enumerate(SIGNS):
                effect = (IDENTITY + final_sign * PAULI_X) / 2
                expected[first_index, final_index] = np.trace(effect @ branch).real
        np.testing.assert_array_equal(expected, np.full((2, 2), 1 / 4))
        self.assertAlmostEqual(np.abs(actual - expected).sum() / 2, 1 / 4)

    def test_refreshed_internal_memory_restores_each_pauli_instrument_branch(self):
        for measurement in range(3):
            instruments = [refreshed_memory_kraus(measurement, outcome) for outcome in SIGNS]
            np.testing.assert_allclose(sum(operation.T @ operation for branch in instruments
                                           for operation in branch), np.eye(4), atol=1e-15)
            for preparation in range(3):
                for sign_index, sign in enumerate(SIGNS):
                    density = shared_memory_state(preparation, sign)
                    for outcome_index, outcome in enumerate(SIGNS):
                        output = sum(operation @ density @ operation.T
                                     for operation in instruments[outcome_index])
                        probability = (1 + sign * outcome * int(preparation == measurement)) / 2
                        expected = probability * shared_memory_state(measurement, outcome)
                        np.testing.assert_allclose(output, expected, atol=5e-16)

    def test_this_shared_memory_encoding_does_not_extend_positively_to_the_full_ball(self):
        joint = (shared_memory_state(0, 1) + shared_memory_state(2, 1)
                 - np.eye(4) / 2) / np.sqrt(2) + np.eye(4) / 4
        self.assertLess(np.linalg.eigvalsh(joint).min(), -0.07)
        density = (IDENTITY + (PAULI_X + PAULI_Z) / np.sqrt(2)) / 2
        self.assertGreaterEqual(np.linalg.eigvalsh(density).min(), -1e-15)

    def test_four_dimensional_real_carrier_matches_without_shared_classical_bits(self):
        for measurement in range(3):
            effect = real_lift((IDENTITY + PAULI_AXES[measurement]) / 2)
            np.testing.assert_array_equal(effect, effect.T)
            self.assertGreaterEqual(np.linalg.eigvalsh(effect).min(), 0)
            self.assertLessEqual(np.linalg.eigvalsh(effect).max(), 1)
            for preparation in range(3):
                for sign_index, sign in enumerate(SIGNS):
                    density = encode_state((IDENTITY + sign * PAULI_AXES[preparation]) / 2)
                    np.testing.assert_array_equal(density, density.T)
                    self.assertGreaterEqual(np.linalg.eigvalsh(density).min(), 0)
                    self.assertEqual(np.trace(density), 1)
                    self.assertAlmostEqual(np.trace(effect @ density).real,
                                           target_table()[measurement, preparation, sign_index])

    def test_nonorthogonal_frame_gram_has_the_predicted_small_singular_value(self):
        for angle in (0., 0.1, 0.5, np.pi / 2, 2.5, np.pi):
            gram = np.array([[1., 0., np.cos(angle)], [0., 1., 0.],
                             [np.cos(angle), 0., 1.]])
            self.assertAlmostEqual(np.linalg.det(gram), np.sin(angle) ** 2)
            self.assertAlmostEqual(np.linalg.svd(gram, compute_uv=False)[-1], 1 - abs(np.cos(angle)))

    def test_finite_sample_budget_covers_all_eighteen_cells(self):
        for visibility in (1., 0.5):
            budget = sample_budget(visibility)
            self.assertLessEqual(budget["union_bound"], budget["risk"])
            self.assertLess(2 * budget["probability_tolerance"], visibility / 6)
            self.assertEqual(budget["total_fresh_preparations"], 18 * budget["repetitions_per_cell"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(RealInterfaceRecordBoundTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 209,
        "antecedent_rounds": [85, 87, 171, 187, 188, 201, 203, 208],
        "scope": "Six preparations, three binary measurements, real dimension 2, no shared side information",
        "extra_shared_source_detector_randomness_or_schedule_allowed_in_bound": False,
        "private_preparation_and_measurement_randomness_allowed": True,
        "preparation_pair_equal_mixture_assumed_in_lower_bound": False,
        "arbitrary_biased_binary_effects_allowed": True,
        "ideal_probability_table_plus": target_table().tolist(),
        "contrast_real_rank_upper": 2,
        "target_contrast": "visibility * I_3",
        "minimum_maximum_probability_error_exact": "visibility/6 for 0<=visibility<=1",
        "attaining_real_contrast": "visibility * (I_3 - ones(3,3)/3)",
        "one_shared_random_bit_plus_rebit_matches_table": True,
        "two_shared_random_bits_plus_one_classical_message_bit_matches_table": True,
        "four_dimensional_real_carrier_matches_without_shared_classical_bits": True,
        "prepare_X_plus_then_Z_X_fixed_shared_bit_luders_record_TV_exact": "1/4",
        "refreshed_internal_bit_real_CP_instrument_matches_all_finite_Pauli_histories": True,
        "refreshed_model_scope": "six axis preparations and mixtures; only three Pauli projective instruments",
        "same_memory_encoding_is_positive_on_full_Bloch_ball": False,
        "fresh_memory_randomness_and_reset_cost_quantified": False,
        "all_possible_shared_memory_update_models_excluded": False,
        "general_frame_lower_bound": "visibility*(1-abs(cos(phi)))/6; sharpness not asserted",
        "sampling_assumptions": "fresh independent Bernoulli trials for each fixed preparation and setting",
        "sample_budgets": [sample_budget(1.), sample_budget(0.5)],
        "samples_collected": 0,
        "minimal_real_dimension_for_full_complex_theory_proved": False,
        "dimension_constraint_or_relative_frame_derived_from_cognition": False,
        "quantum_necessity_or_incompleteness_proved": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("real_interface_record_bound_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()