"""Round 215: spectral obstruction for a continuous full real correspondence.

On all real symmetric observables, a continuous assignment of real skew
generators commuting with their own observable must vanish: simple spectra
are dense. Degenerate encoded interfaces and retained ancillas are explicit
exceptions to the full-domain premise, not violations of the theorem.
"""

import argparse
import json
import math
import unittest
from itertools import combinations
from pathlib import Path

import numpy as np

from common_orientation_structure import encode_state
from complex_control_from_reference import real_lift
from question_action_covariance import skew_unit
from real_interface_record_bound import IDENTITY, PAULI_X, PAULI_Z
from question_algebra_orientation import REAL_J
from observable_action_compatibility import jordan, left_jordan_commutator


def commutator(first, second):
    return first @ second - second @ first


def spectral_energy(eigenvalues, generator):
    differences = np.asarray(eigenvalues)[:, None] - np.asarray(eigenvalues)[None, :]
    return np.sum(differences * differences * generator * generator)


def skew_commutator_gram(eigenvalues):
    observable = np.diag(eigenvalues)
    basis = [skew_unit(len(eigenvalues), row, column)
             for row, column in combinations(range(len(eigenvalues)), 2)]
    images = [commutator(observable, generator) for generator in basis]
    return np.array([[np.sum(first * second) for second in images] for first in images])


def preserved_skew_dimension(eigenvalues):
    return sum(eigenvalues[row] == eigenvalues[column]
               for row, column in combinations(range(len(eigenvalues)), 2))


def complex_generator(observable):
    centered = observable - np.trace(observable) * np.eye(len(observable)) / len(observable)
    return -0.5j * centered


class RealSpectralActionBoundTests(unittest.TestCase):
    def test_spectral_commutator_energy_is_an_exact_integer_identity(self):
        random = np.random.default_rng(215)
        for dimension in range(2, 7):
            eigenvalues = np.arange(dimension) ** 2
            raw = random.integers(-3, 4, size=(dimension, dimension))
            generator = raw - raw.T
            image = commutator(np.diag(eigenvalues), generator)
            np.testing.assert_array_equal(image.T, image)
            self.assertEqual(np.sum(image * image), spectral_energy(eigenvalues, generator))
            gap = min(abs(eigenvalues[row] - eigenvalues[column])
                      for row, column in combinations(range(dimension), 2))
            self.assertGreaterEqual(np.sum(image * image), gap * gap * np.sum(generator * generator))

    def test_nondegenerate_observable_has_an_exact_full_rank_commutator_gram(self):
        for dimension in range(2, 7):
            eigenvalues = np.arange(dimension)
            expected = np.diag([2 * (eigenvalues[row] - eigenvalues[column]) ** 2
                                for row, column in combinations(range(dimension), 2)])
            np.testing.assert_array_equal(skew_commutator_gram(eigenvalues), expected)
            self.assertTrue(np.all(np.diag(expected) > 0))
            self.assertEqual(preserved_skew_dimension(eigenvalues), 0)

    def test_spectral_gap_bound_is_covariant_under_real_basis_changes(self):
        random = np.random.default_rng(1215)
        for dimension in (2, 3, 4, 5):
            orthogonal, _ = np.linalg.qr(random.normal(size=(dimension, dimension)))
            eigenvalues = np.arange(dimension, dtype=float)
            raw = random.normal(size=(dimension, dimension))
            original_generator = raw - raw.T
            generator = orthogonal @ original_generator @ orthogonal.T
            observable = orthogonal @ np.diag(eigenvalues) @ orthogonal.T
            image = commutator(observable, generator)
            self.assertAlmostEqual(np.linalg.norm(image) ** 2, spectral_energy(eigenvalues, original_generator))
            self.assertGreaterEqual(np.linalg.norm(image) + 1e-13, np.linalg.norm(generator))

    def test_worst_state_expectation_rate_is_the_commutator_operator_norm(self):
        observable = np.diag([0., 1., 3.])
        generator = np.array([[0., 1., -2.], [-1., 0., 0.5], [2., -0.5, 0.]])
        image = commutator(observable, generator)
        eigenvalues, eigenvectors = np.linalg.eigh(image)
        index = np.argmax(np.abs(eigenvalues))
        state = np.outer(eigenvectors[:, index], eigenvectors[:, index])
        rate = np.trace(observable @ (generator @ state - state @ generator))
        self.assertAlmostEqual(abs(rate), np.linalg.norm(image, ord=2))
        self.assertGreaterEqual(abs(rate), np.linalg.norm(generator) / np.sqrt(3))

    def test_degenerate_eigenspace_supports_real_rotation_but_splitting_removes_it(self):
        generator = skew_unit(3, 0, 1)
        degenerate = np.diag([0., 0., 1.])
        np.testing.assert_array_equal(commutator(degenerate, generator), np.zeros((3, 3)))
        self.assertEqual(preserved_skew_dimension([0, 0, 1]), 1)
        for split in (0.5, 0.25, 0.125):
            eigenvalues = [0, split, 1]
            self.assertEqual(preserved_skew_dimension(eigenvalues), 0)
            self.assertEqual(np.sum(commutator(np.diag(eigenvalues), generator) ** 2), 2 * split ** 2)

    def test_nonzero_degenerate_only_assignment_would_be_discontinuous(self):
        at_degeneracy = skew_unit(3, 0, 1)
        self.assertEqual(np.sum(at_degeneracy ** 2), 2)
        for split in (0.5, 0.25, 0.125):
            self.assertEqual(preserved_skew_dimension([0, split, 1]), 0)
            nearby_only_option = np.zeros((3, 3))
            self.assertEqual(np.sum((at_degeneracy - nearby_only_option) ** 2), 2)

    def test_complex_phases_preserve_a_simple_spectrum_but_change_other_answers(self):
        observable = np.diag([-1., 1.])
        generator = complex_generator(observable)
        np.testing.assert_array_equal(generator.conj().T, -generator)
        np.testing.assert_array_equal(commutator(observable, generator), np.zeros((2, 2)))
        self.assertGreater(np.linalg.norm(generator), 0)
        initial = (IDENTITY + PAULI_X) / 2
        phase = np.diag(np.exp(np.diag(generator) * np.pi / 2))
        output = phase @ initial @ phase.conj().T
        self.assertAlmostEqual(np.trace(observable @ output).real, np.trace(observable @ initial).real)
        self.assertAlmostEqual(np.trace(PAULI_X @ output).real, 0)
        self.assertGreater(np.linalg.norm(output - initial), 0.9)

    def test_common_J_real_encoding_has_exact_spectral_doubling_and_nonzero_control(self):
        observable = np.diag([-1., 1.])
        lifted_observable = real_lift(observable)
        lifted_generator = real_lift(complex_generator(observable))
        np.testing.assert_array_equal(np.linalg.eigvalsh(lifted_observable), [-1., -1., 1., 1.])
        self.assertEqual(preserved_skew_dimension([-1, -1, 1, 1]), 2)
        np.testing.assert_array_equal(lifted_generator.T, -lifted_generator)
        np.testing.assert_array_equal(commutator(lifted_observable, lifted_generator), np.zeros((4, 4)))
        initial = (IDENTITY + PAULI_X) / 2
        phase = np.diag(np.exp(np.diag(complex_generator(observable)) * np.pi / 2))
        real_gate = real_lift(phase)
        output = real_gate @ encode_state(initial) @ real_gate.T
        np.testing.assert_allclose(output, encode_state(phase @ initial @ phase.conj().T), atol=1e-16)
        self.assertAlmostEqual(np.trace(real_lift(PAULI_X) @ output), 0)

    def test_breaking_encoding_degeneracy_breaks_own_observable_conservation(self):
        observable = real_lift(np.diag([-1., 1.]))
        generator = real_lift(complex_generator(np.diag([-1., 1.])))
        for split in (0.5, 0.25, 0.125):
            perturbed = observable + split * np.diag([0., 0., 1., 1.])
            self.assertEqual(preserved_skew_dimension(np.diag(perturbed)), 0)
            self.assertEqual(np.sum(commutator(perturbed, generator) ** 2), split ** 2)
            self.assertEqual(np.sum(generator ** 2), 1)

    def test_retained_real_ancilla_allows_nontrivial_data_question_preserving_control(self):
        question = np.kron(PAULI_Z, IDENTITY)
        angle = 0.7
        rotation = math.cos(angle) * IDENTITY + math.sin(angle) * REAL_J
        gate = np.zeros((4, 4))
        gate[:2, :2] = IDENTITY
        gate[2:, 2:] = rotation
        initial = np.kron((IDENTITY + PAULI_X) / 2, np.diag([1., 0.]))
        output = gate @ initial @ gate.T
        np.testing.assert_allclose(gate.T @ question @ gate, question, rtol=0, atol=2e-16)
        np.testing.assert_allclose(gate.T @ output @ gate, initial, atol=2e-16)
        self.assertAlmostEqual(np.trace(np.kron(PAULI_X, IDENTITY) @ output), math.cos(angle))
        self.assertEqual(preserved_skew_dimension([1, 1, -1, -1]), 2)

    def test_one_preparation_conservation_is_weaker_than_all_state_conservation(self):
        question = PAULI_Z
        generator = REAL_J
        initial = IDENTITY / 2
        self.assertEqual(np.trace(question @ (generator @ initial - initial @ generator)), 0)
        self.assertGreater(np.linalg.norm(commutator(question, generator)), 1)

    def test_zero_real_correspondence_cannot_satisfy_nonassociative_jordan_compatibility(self):
        np.testing.assert_array_equal(left_jordan_commutator(PAULI_Z, PAULI_X, PAULI_X), PAULI_Z)
        np.testing.assert_array_equal(jordan(PAULI_X, PAULI_Z), np.zeros((2, 2)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(RealSpectralActionBoundTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 215,
        "antecedent_rounds": [80, 85, 88, 173, 210, 211, 212, 213, 214],
        "scope": "All real symmetric observables, continuous assignment to real orthogonal generators on the same full carrier",
        "own_question_preservation_for_all_states": "[A,K(A)]=0",
        "simple_real_spectrum_commuting_skew_generator": "zero",
        "continuous_full_domain_assignment": "identically zero by density of simple spectra",
        "linearity_or_SO_dimension_isotropy_required": False,
        "real_degenerate_commutant_dimension": "sum_j m_j*(m_j-1)/2",
        "complex_degenerate_commutant_dimension_mod_global_phase": "sum_j m_j^2-1",
        "spectral_energy_identity": "norm([A,K])_F^2=sum_ij (lambda_i-lambda_j)^2*K_ij^2",
        "simple_spectrum_gap_bound": "norm([A,K])_F >= delta*norm(K)_F",
        "worst_expectation_rate_bound": "sup_rho abs(Tr(A[K,rho])) >= delta*norm(K)_F/sqrt(d)",
        "all_real_Kraus_or_dilated_instruments_excluded": False,
        "common_J_encoded_observables_have_even_spectral_multiplicities": True,
        "retained_real_ancilla_can_preserve_a_restricted_data_question_nontrivially": True,
        "degenerate_only_discontinuous_assignment_excluded_without_continuity": False,
        "all_state_condition_follows_from_one_preparation": False,
        "ordinary_real_full_orthogonal_correspondence_meets_round214_nontrivial_compatibility": False,
        "complex_equivalent_real_representation_excluded": False,
        "continuous_full_domain_coupling_derived_from_cognition": False,
        "physical_time_or_general_resource_gap_derived": False,
        "quantum_necessity_or_incompleteness_proved": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("real_spectral_action_bound_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()