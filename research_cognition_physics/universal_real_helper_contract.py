"""Round 216: a universal retained real helper and its return contract.

For every real symmetric A, J_2 tensor A/2 is an analytic orthogonal-control
generator preserving I_2 tensor A. A maximally mixed helper keeps its marginal
but stores correlations. Requiring exact product return for every data input
instead forces a same-dimension orthogonal data channel; round 215 applies.
"""

import argparse
import json
import math
import unittest
from itertools import combinations, product
from pathlib import Path

import numpy as np

from common_orientation_structure import encode_state
from complex_control_from_reference import real_lift
from observable_action_compatibility import left_jordan_commutator
from partition_interface_closure import sector_bases
from question_algebra_orientation import REAL_J
from real_interface_record_bound import IDENTITY, PAULI_X, PAULI_Z
from real_spectral_action_bound import commutator


def generator_for(observable):
    return np.kron(REAL_J, observable) / 2


def embedded_observable(observable):
    return np.kron(IDENTITY, observable)


def natural_helper_generator(observable, data_scale, helper_scale):
    data_term = data_scale * observable + helper_scale * np.trace(observable) * np.eye(len(observable))
    return np.kron(REAL_J, data_term)


def matrix_trigonometry(observable, parameter):
    eigenvalues, eigenvectors = np.linalg.eigh(observable)
    cosine = (eigenvectors * np.cos(parameter * eigenvalues / 2)) @ eigenvectors.T
    sine = (eigenvectors * np.sin(parameter * eigenvalues / 2)) @ eigenvectors.T
    return cosine, sine


def joint_gate(observable, parameter):
    cosine, sine = matrix_trigonometry(observable, parameter)
    return np.kron(IDENTITY, cosine) + np.kron(REAL_J, sine)


def data_marginal(joint):
    dimension = len(joint) // 2
    return np.trace(joint.reshape(2, dimension, 2, dimension), axis1=0, axis2=2)


def helper_marginal(joint):
    dimension = len(joint) // 2
    return np.trace(joint.reshape(2, dimension, 2, dimension), axis1=1, axis2=3)


def reduced_channel(observable, parameter, density):
    cosine, sine = matrix_trigonometry(observable, parameter)
    return cosine @ density @ cosine + sine @ density @ sine


def entropy_bits(density):
    eigenvalues = np.linalg.eigvalsh(density)
    if eigenvalues.min() < -1e-12:
        raise ValueError("Entropy requires a positive state.")
    positive = eigenvalues[eigenvalues > 1e-14]
    return float(-np.sum(positive * np.log2(positive)))


def trace_distance(first, second):
    return float(np.abs(np.linalg.eigvalsh(first - second)).sum() / 2)


class UniversalRealHelperContractTests(unittest.TestCase):
    def test_universal_assignment_is_linear_skew_and_preserves_its_embedded_question(self):
        random = np.random.default_rng(216)
        for dimension in (2, 3, 4):
            raw_first = random.integers(-3, 4, size=(dimension, dimension))
            raw_second = random.integers(-3, 4, size=(dimension, dimension))
            first, second = raw_first + raw_first.T, raw_second + raw_second.T
            generator = generator_for(first)
            np.testing.assert_array_equal(generator.T, -generator)
            np.testing.assert_array_equal(commutator(embedded_observable(first), generator), np.zeros_like(generator))
            np.testing.assert_array_equal(generator_for(2 * first - 3 * second),
                                          2 * generator - 3 * generator_for(second))

    def test_joint_gate_is_real_orthogonal_reversible_and_preserves_all_question_probabilities(self):
        observable = np.array([[1., 2., 0.], [2., -1., 1.], [0., 1., 3.]])
        gate = joint_gate(observable, 0.37)
        np.testing.assert_allclose(gate.T @ gate, np.eye(6), rtol=0, atol=2e-15)
        np.testing.assert_allclose(joint_gate(observable, -0.37) @ gate, np.eye(6), rtol=0, atol=2e-15)
        np.testing.assert_allclose(gate.T @ embedded_observable(observable) @ gate,
                                   embedded_observable(observable), rtol=0, atol=8e-15)
        _, eigenvectors = np.linalg.eigh(observable)
        for index in range(3):
            projector = embedded_observable(np.outer(eigenvectors[:, index], eigenvectors[:, index]))
            np.testing.assert_allclose(gate.T @ projector @ gate, projector, rtol=0, atol=2e-15)

    def test_joint_control_is_covariant_under_actual_real_data_basis_changes(self):
        random = np.random.default_rng(1216)
        orthogonal, _ = np.linalg.qr(random.normal(size=(3, 3)))
        observable = np.diag([-1., 0., 2.])
        expanded = np.kron(IDENTITY, orthogonal)
        np.testing.assert_allclose(generator_for(orthogonal @ observable @ orthogonal.T),
                                   expanded @ generator_for(observable) @ expanded.T, rtol=0, atol=3e-16)

    def test_covariant_two_parameter_family_separates_data_control_and_helper_only_rotation(self):
        for dimension in (2, 3, 4):
            observable = np.diag(np.arange(dimension))
            basis_change = np.eye(dimension, dtype=int)[::-1]
            basis_change[0] *= -1
            expanded = np.kron(IDENTITY, basis_change)
            for data_scale, helper_scale in ((0.5, 0), (1, -1), (0, 2)):
                actual = natural_helper_generator(basis_change @ observable @ basis_change.T, data_scale, helper_scale)
                expected = expanded @ natural_helper_generator(observable, data_scale, helper_scale) @ expanded.T
                np.testing.assert_array_equal(actual, expected)
                np.testing.assert_array_equal(commutator(embedded_observable(observable),
                                                         natural_helper_generator(observable, data_scale, helper_scale)),
                                              np.zeros_like(actual))
            helper_only = natural_helper_generator(observable, 0, 1)
            symmetric, _ = sector_bases(dimension)
            for question in symmetric:
                np.testing.assert_array_equal(commutator(helper_only, embedded_observable(question)),
                                              np.zeros_like(helper_only))

    def test_data_reduction_matches_the_cosine_channel_for_every_real_helper_state(self):
        observable = np.diag([-1., 0., 2.])
        density = np.ones((3, 3)) / 3
        expected = reduced_channel(observable, 0.71, density)
        helpers = (IDENTITY / 2, np.diag([1., 0.]), (IDENTITY + PAULI_X) / 2,
                   (IDENTITY + 0.3 * PAULI_X + 0.4 * PAULI_Z) / 2)
        gate = joint_gate(observable, 0.71)
        for helper in helpers:
            actual = data_marginal(gate @ np.kron(helper, density) @ gate.T)
            np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-16)
        eigenvalues = np.diag(observable)
        multipliers = np.cos(0.71 * (eigenvalues[:, None] - eigenvalues[None, :]) / 2)
        np.testing.assert_allclose(expected, multipliers * density, rtol=0, atol=2e-16)

    def test_instantaneous_data_generator_vanishes_but_second_derivative_does_not(self):
        initial = (IDENTITY + PAULI_X) / 2
        joint = np.kron(IDENTITY / 2, initial)
        generator = generator_for(PAULI_Z)
        np.testing.assert_array_equal(data_marginal(commutator(generator, joint)), np.zeros((2, 2)))
        second = data_marginal(commutator(generator, commutator(generator, joint)))
        np.testing.assert_array_equal(second, -commutator(PAULI_Z, commutator(PAULI_Z, initial)) / 4)
        np.testing.assert_array_equal(second, -PAULI_X / 2)

    def test_fixed_helper_marginal_survives_alternating_controls_as_correlations_accumulate(self):
        initial = (IDENTITY + PAULI_X) / 2
        joint = np.kron(IDENTITY / 2, initial)
        complex_density = initial.astype(complex)
        for observable, parameter in ((PAULI_Z, 0.7), (PAULI_X, -0.3), (PAULI_Z, 0.4)):
            gate = joint_gate(observable, parameter)
            cosine, sine = matrix_trigonometry(observable, parameter)
            complex_gate = cosine - 1j * sine
            np.testing.assert_allclose(gate, real_lift(complex_gate), rtol=0, atol=0)
            joint = gate @ joint @ gate.T
            complex_density = complex_gate @ complex_density @ complex_gate.conj().T
            np.testing.assert_allclose(joint, encode_state(complex_density), rtol=0, atol=2e-15)
            np.testing.assert_allclose(helper_marginal(joint), IDENTITY / 2, rtol=0, atol=2e-15)

    def test_maximally_mixed_helper_return_is_not_product_return(self):
        initial = (IDENTITY + PAULI_X) / 2
        gate = joint_gate(PAULI_Z, math.pi / 2)
        output = gate @ np.kron(IDENTITY / 2, initial) @ gate.T
        np.testing.assert_allclose(data_marginal(output), IDENTITY / 2, rtol=0, atol=2e-16)
        np.testing.assert_allclose(helper_marginal(output), IDENTITY / 2, rtol=0, atol=2e-16)
        independent = np.kron(helper_marginal(output), data_marginal(output))
        self.assertAlmostEqual(trace_distance(output, independent), 0.5)
        mutual_information = entropy_bits(helper_marginal(output)) + entropy_bits(data_marginal(output)) - entropy_bits(output)
        self.assertAlmostEqual(mutual_information, 1)
        self.assertAlmostEqual(np.trace(output @ output), 0.5)
        self.assertAlmostEqual(np.trace(independent @ independent), 0.25)

    def test_retained_inverse_restores_data_but_fresh_helper_inverse_does_not(self):
        initial = (IDENTITY + PAULI_X) / 2
        parameter = math.pi / 2
        forward = joint_gate(PAULI_Z, parameter)
        joint = forward @ np.kron(IDENTITY / 2, initial) @ forward.T
        restored = data_marginal(forward.T @ joint @ forward)
        np.testing.assert_allclose(restored, initial, rtol=0, atol=2e-15)
        fresh_inverse = reduced_channel(PAULI_Z, -parameter, data_marginal(joint))
        self.assertAlmostEqual(trace_distance(fresh_inverse, initial), 0.5)
        self.assertAlmostEqual(np.trace(initial @ restored), 1)
        self.assertAlmostEqual(np.trace(initial @ fresh_inverse), 0.5)

    def test_fresh_helper_channels_do_not_follow_the_retained_group_law(self):
        initial = (IDENTITY + PAULI_X) / 2
        first_parameter, second_parameter = 0.7, 0.9
        fresh = reduced_channel(PAULI_Z, second_parameter,
                                reduced_channel(PAULI_Z, first_parameter, initial))
        retained = reduced_channel(PAULI_Z, first_parameter + second_parameter, initial)
        expected_distance = abs(math.sin(first_parameter) * math.sin(second_parameter)) / 2
        self.assertAlmostEqual(trace_distance(fresh, retained), expected_distance)

    def test_auxiliary_derivative_leaves_the_data_space_despite_double_commutator_matching(self):
        data_derivative = commutator(generator_for(PAULI_X), embedded_observable(PAULI_Z))
        self.assertGreater(np.linalg.norm(data_derivative), 1)
        np.testing.assert_array_equal(data_marginal(data_derivative), np.zeros((2, 2)))
        symmetric, _ = sector_bases(2)
        for first, second, third in product(symmetric, repeat=3):
            double_generator = commutator(generator_for(first), generator_for(second))
            actual = commutator(double_generator, embedded_observable(third))
            expected = -embedded_observable(left_jordan_commutator(first, second, third))
            np.testing.assert_array_equal(actual, expected)

    def test_minimal_closed_symmetric_space_contains_the_full_complex_equivalent_interface(self):
        for dimension in (2, 3, 4):
            symmetric, skew = sector_bases(dimension)
            basis = [embedded_observable(matrix) for matrix in symmetric]
            basis.extend(np.kron(REAL_J, matrix) for matrix in skew)
            self.assertEqual(len(basis), dimension * dimension)
            gram = np.array([[np.sum(first * second) for second in basis] for first in basis])
            np.testing.assert_array_equal(gram, np.diag(np.diag(gram)))
            self.assertTrue(np.all(np.diag(gram) > 0))
            orientation = np.kron(REAL_J, np.eye(dimension))
            for matrix in basis:
                np.testing.assert_array_equal(matrix, matrix.T)
                np.testing.assert_array_equal(commutator(matrix, orientation), np.zeros_like(matrix))
            for row, column in combinations(range(dimension), 2):
                diagonal = np.zeros((dimension, dimension))
                diagonal[row, row] = 1
                off_diagonal = np.zeros_like(diagonal)
                off_diagonal[row, column] = off_diagonal[column, row] = 1
                expected = np.kron(REAL_J, diagonal @ off_diagonal - off_diagonal @ diagonal) / 2
                np.testing.assert_array_equal(commutator(generator_for(diagonal), embedded_observable(off_diagonal)), expected)

    def test_product_return_has_the_purity_constraint_and_orthogonal_examples_saturate_it(self):
        helper = np.diag([0.3, 0.7])
        target_rotation = math.cos(0.31) * IDENTITY + math.sin(0.31) * REAL_J
        joint = np.kron(IDENTITY, target_rotation)
        for density in ((IDENTITY + PAULI_X) / 2, np.diag([0.2, 0.8])):
            output = joint @ np.kron(helper, density) @ joint.T
            target = target_rotation @ density @ target_rotation.T
            np.testing.assert_allclose(output, np.kron(helper, target), rtol=0, atol=2e-15)
            self.assertAlmostEqual(np.trace(output @ output), np.trace(helper @ helper) * np.trace(density @ density))
            self.assertAlmostEqual(np.trace(target @ target), np.trace(density @ density))

    def test_retained_control_reversal_preserves_arbitrary_external_relations(self):
        random = np.random.default_rng(2216)
        raw = random.normal(size=(12, 12))
        joint = raw @ raw.T
        joint /= np.trace(joint)
        observable = np.array([[1., 1., 0.], [1., 0., -1.], [0., -1., 2.]])
        gate = np.kron(joint_gate(observable, 0.27), np.eye(2))
        output = gate @ joint @ gate.T
        np.testing.assert_allclose(gate.T @ output @ gate, joint, rtol=0, atol=6e-16)
        self.assertGreater(np.linalg.eigvalsh(output).min(), 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(UniversalRealHelperContractTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 216,
        "antecedent_rounds": [2, 56, 85, 88, 122, 124, 172, 173, 214, 215],
        "scope": "Universal real data-observable control with a retained helper; marginal versus product return",
        "joint_generator": "K_A=J_2 tensor A/2 for every A in Sym_d(R)",
        "helper_dimension": 2,
        "data_observable_embedding": "I_2 tensor A; not all Sym_(2d)(R)",
        "linear_analytic_assignment_preserves_embedded_A_for_all_joint_inputs": True,
        "classification_with_one_helper_and_full_real_data_basis_covariance": "All linear equivariant K:Sym_d(R)->so(2d) have K_A=J_2 tensor (alpha*A+beta*Tr(A)*I), for d>=2",
        "basis_covariance_or_minimal_helper_derived_from_cognition": False,
        "reduced_data_channel": "C rho C+S rho S, C=cos(sA/2), S=sin(sA/2)",
        "reduced_eigenbasis_multiplier": "cos(s*(lambda_i-lambda_j)/2)",
        "reduced_first_derivative_at_zero": "zero",
        "reduced_second_derivative_at_zero": "-[A,[A,rho]]/4",
        "maximally_mixed_helper_marginal_preserved_on_encoded_histories": True,
        "helper_product_return_for_all_inputs": False,
        "qubit_witness": {
            "A": "Z",
            "parameter": "pi/2",
            "initial_data": "+X",
            "initial_helper": "I/2",
            "both_final_marginals": "I/2",
            "distance_from_product_of_marginals": "1/2",
            "mutual_information_bits": "1",
            "retained_inverse_X_plus_probability": "1",
            "fresh_helper_inverse_X_plus_probability": "1/2"
        },
        "embedded_observable_first_derivative_closes_inside_data_space": False,
        "double_commutator_matches_round214_on_embedded_data": True,
        "minimal_symmetric_derivative_closed_space": "I tensor Sym_d(R) plus J_2 tensor Skew_d(R), dimension d^2",
        "resulting_interface_is_common_J_complex_equivalent": True,
        "whole_ordinary_real_theory_has_been_restricted_to_that_interface": False,
        "exact_product_return_theorem": "Finite fixed input-independent helper restored uncorrelated for every data state by a real orthogonal joint operation forces an orthogonal data channel",
        "continuous_full_A_own_question_preserving_product_return_family": "identity data channel",
        "product_return_assumption_follows_from_internal_resource_accounting": False,
        "source_purifier_and_clock_costs_fully_computed": False,
        "cognitive_necessity_or_quantum_incompleteness_proved": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("universal_real_helper_contract_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()