"""All real completely positive extensions of the calibrated fresh branches.

The local symmetric-sector map is fixed. The remaining scalar is constrained
by an exact interval, not chosen by a cognitive principle in this project.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from ancillary_operation_equivalence import (J, REAL_BASIS, apply_operation, apply_on_first,
                                             determinant_parameter, local_operation_summary)
from bipartite_composition import HIDDEN, interaction
from certified_intervals import Interval as I, SCALE, sin_interval
from quantum_interface_audit import (ALPHA, IDENTITY, PAULI_X, PAULI_Y, PAULI_Z,
                                     branch_kraus, effective_readout, rotation_unitary,
                                     two_rebit_states)


def pauli_weights(attenuation, gamma):
    """Order I, X, J=-iY, Z. Also accepts Fraction for exact witnesses."""
    return ((1 + 2*attenuation + gamma)/4, (1-gamma)/4,
            (gamma - (2*attenuation - 1))/4, (1-gamma)/4)


def extension_noise_kraus(attenuation, gamma):
    if not 0 <= attenuation <= 1 or not 2*attenuation - 1 <= gamma <= 1:
        raise ValueError("A real CP completion needs 0<=a<=1 and 2a-1<=gamma<=1.")
    return tuple(math.sqrt(float(p)) * matrix for p, matrix in
                 zip(pauli_weights(attenuation, gamma), (IDENTITY.real, PAULI_X.real, J, PAULI_Z.real)))


def measurement_filter(setting, outcome, contrast):
    if outcome not in (0, 1) or not 0 <= contrast <= 1:
        raise ValueError("Use outcome 0/1 and contrast in [0,1].")
    sign = 2*outcome - 1
    rotation = rotation_unitary(setting).real
    diagonal = np.diag(np.sqrt([(1 + sign*contrast)/2, (1 - sign*contrast)/2]))
    return rotation @ diagonal @ rotation.T


def extension_branch_kraus(setting, outcome, contrast, gamma, attenuation=ALPHA):
    filter_matrix = measurement_filter(setting, outcome, contrast)
    return tuple(filter_matrix @ k for k in extension_noise_kraus(attenuation, gamma))


def kappa_bounds(contrast, attenuation=ALPHA):
    if not 0 <= contrast <= 1:
        raise ValueError("Use contrast in [0,1].")
    scale = math.sqrt(1 - contrast**2)/2
    return scale*(2*attenuation - 1), scale


def normalized_choi_from_action(action):
    result = np.zeros((4, 4), dtype=complex)
    for i, j in product(range(2), repeat=2):
        matrix = np.zeros((2, 2))
        matrix[i, j] = 1
        result += np.kron(action(matrix), matrix)/2
    return result


def diagonal_action(matrix, attenuation, gamma):
    return (np.trace(matrix)*IDENTITY + attenuation*np.trace(PAULI_X @ matrix)*PAULI_X
            + gamma*np.trace(PAULI_Y @ matrix)*PAULI_Y
            + attenuation*np.trace(PAULI_Z @ matrix)*PAULI_Z)/2


def extension_interval_certificate(contrast=Fraction(1, 2)):
    eta = Fraction(contrast)
    if not 0 <= eta <= 1:
        raise ValueError("Use contrast in [0,1].")
    alpha = 4 * sin_interval(I.rational(1, 4))
    scale = I.exact(0) if eta == 1 else (1 - I.exact(eta)**2).sqrt()/2
    lower, upper = scale*(2*alpha - 1), scale
    width = 2*scale*(1-alpha)
    readout_distance = alpha*width
    return {"contrast_exact": str(eta), "kappa_lower_interval": lower.floats(),
            "kappa_upper_interval": upper.floats(), "kappa_width_interval": width.floats(),
            "width_strict_lower_exact": str(Fraction(width.lo, SCALE)),
            "full_instrument_extreme_probe_TV_interval": readout_distance.floats(),
            "unique_extension": eta == 1}


def instrument_probe_distribution(contrast, gammas):
    density = two_rebit_states()[0]
    gate = interaction(math.pi/2)
    result = {}
    for outcome, gamma in enumerate(gammas):
        output = apply_on_first(extension_branch_kraus(.37, outcome, contrast, gamma), density)
        output = gate @ output @ gate.conj().T
        for readout in (0, 1):
            effect = np.kron(IDENTITY, effective_readout(0, readout, 1.))
            result[outcome, readout] = float(np.trace(effect @ output).real)
    return result


class FreshExtensionClassificationTests(unittest.TestCase):
    def test_entire_noise_interval_has_real_kraus_constructions(self):
        for attenuation, fraction in product((0., .2, .6, ALPHA, 1.), (0., .3, 1.)):
            gamma = (2*attenuation-1) + fraction*(2-2*attenuation)
            gamma = min(1., max(2*attenuation-1, gamma))
            kraus = extension_noise_kraus(attenuation, gamma)
            for matrix in (*REAL_BASIS, PAULI_Y):
                np.testing.assert_allclose(apply_operation(kraus, matrix), diagonal_action(matrix, attenuation, gamma), atol=3e-16)
            np.testing.assert_allclose(sum(k.T @ k for k in kraus), np.eye(2), atol=3e-16)

    def test_bell_probe_eigenvalues_are_exactly_the_four_pauli_weights(self):
        for attenuation, gamma in ((.6, .4), (ALPHA, 1.), (ALPHA, ALPHA), (.6, .1), (.6, 1.1)):
            choi = normalized_choi_from_action(lambda matrix: diagonal_action(matrix, attenuation, gamma))
            expected = (np.eye(4) + attenuation*np.kron(PAULI_X, PAULI_X)
                        - gamma*HIDDEN + attenuation*np.kron(PAULI_Z, PAULI_Z))/4
            np.testing.assert_allclose(choi, expected, atol=1e-16)
            np.testing.assert_allclose(np.linalg.eigvalsh(choi), sorted(pauli_weights(attenuation, gamma)), atol=3e-16)

    def test_invalid_extensions_have_exact_negative_bell_eigenvalues(self):
        attenuation = Fraction(3, 5)
        for gamma in (Fraction(1, 10), Fraction(11, 10)):
            self.assertEqual(min(pauli_weights(attenuation, gamma)), Fraction(-1, 40))
            with self.assertRaises(ValueError):
                extension_noise_kraus(attenuation, gamma)

    def test_every_allowed_fresh_extension_keeps_the_entire_original_local_map(self):
        for setting, outcome, eta, gamma in product((0., .31, -1.2), (0, 1), (0., .5, .93, 1.), (2*ALPHA-1, ALPHA, 1.)):
            actual = extension_branch_kraus(setting, outcome, eta, gamma)
            original = branch_kraus(setting, outcome, eta)
            np.testing.assert_allclose(local_operation_summary(actual), local_operation_summary(original), atol=7e-16)
            self.assertAlmostEqual(determinant_parameter(actual), gamma*math.sqrt(1-eta**2)/2, places=14)

    def test_invertible_filter_reduces_the_complete_positive_constraint_to_noise(self):
        for setting, outcome, eta, gamma in product((0., .6), (0, 1), (.1, .8), (2*ALPHA-1, .99, 1.)):
            filter_matrix = measurement_filter(setting, outcome, eta)
            kraus = extension_branch_kraus(setting, outcome, eta, gamma)
            choi = normalized_choi_from_action(lambda matrix: apply_operation(kraus, matrix))
            inverse = np.kron(np.linalg.inv(filter_matrix), np.eye(2))
            recovered = inverse @ choi @ inverse.T
            expected = normalized_choi_from_action(lambda matrix: diagonal_action(matrix, ALPHA, gamma))
            np.testing.assert_allclose(recovered, expected, atol=1e-15)

    def test_unit_strength_branches_have_zero_kappa_and_unique_full_extension(self):
        for setting, outcome in product((0., .8), (0, 1)):
            maps = [extension_branch_kraus(setting, outcome, 1., gamma) for gamma in (2*ALPHA-1, ALPHA, 1.)]
            for matrix in (*REAL_BASIS, PAULI_Y):
                for kraus in maps[1:]:
                    np.testing.assert_allclose(apply_operation(kraus, matrix), apply_operation(maps[0], matrix), atol=4e-16)
            for kraus in maps:
                self.assertAlmostEqual(determinant_parameter(kraus), 0., places=15)
        self.assertEqual(kappa_bounds(1.), (0., 0.))
        self.assertTrue(extension_interval_certificate(Fraction(1))["unique_extension"])

    def test_the_two_outcome_extensions_can_be_chosen_independently(self):
        for eta, first, second in product((0., .5, 1.), (2*ALPHA-1, 1.), (2*ALPHA-1, 1.)):
            branches = [extension_branch_kraus(.3, outcome, eta, gamma) for outcome, gamma in enumerate((first, second))]
            np.testing.assert_allclose(sum(k.T @ k for branch in branches for k in branch), np.eye(2), atol=7e-16)
            probabilities = instrument_probe_distribution(eta, (first, second))
            self.assertAlmostEqual(sum(probabilities.values()), 1., places=14)
            self.assertGreaterEqual(min(probabilities.values()), -1e-15)

    def test_exact_all_ancilla_distance_bound_is_saturated_by_one_real_auxiliary(self):
        rng = np.random.default_rng(53)
        first, second = [extension_noise_kraus(ALPHA, gamma) for gamma in (2*ALPHA-1, 1.)]
        bound = 1-ALPHA
        for density in two_rebit_states():
            difference = apply_on_first(first, density) - apply_on_first(second, density)
            self.assertAlmostEqual(np.linalg.norm(difference, ord="nuc")/2, bound, places=14)
        for dimension in (2, 3, 5):
            for _ in range(8):
                matrix = rng.normal(size=(2*dimension, 2*dimension))
                density = matrix @ matrix.T
                density /= np.trace(density)
                difference = apply_on_first(first, density) - apply_on_first(second, density)
                self.assertLessEqual(np.linalg.norm(difference, ord="nuc")/2, bound + 3e-16)

    def test_complete_instrument_record_distance_matches_the_extreme_formula(self):
        for eta in (0., .5, .85, 1.):
            probabilities = [instrument_probe_distribution(eta, (gamma, gamma)) for gamma in (2*ALPHA-1, 1.)]
            distance = sum(abs(probabilities[0][key] - probabilities[1][key]) for key in probabilities[0])/2
            self.assertAlmostEqual(distance, ALPHA*(1-ALPHA)*math.sqrt(1-eta**2), places=14)

    def test_individual_record_probabilities_include_mass_and_kappa(self):
        eta = .7
        gammas = (2*ALPHA-1, ALPHA)
        probabilities = instrument_probe_distribution(eta, gammas)
        for (outcome, readout), value in probabilities.items():
            kappa = gammas[outcome]*math.sqrt(1-eta**2)/2
            self.assertAlmostEqual(value, .25 + (2*readout-1)*ALPHA*kappa/2, places=14)

    def test_sequential_noise_completion_stays_inside_the_derived_interval(self):
        for first, second in product((.1, .7, ALPHA), repeat=2):
            for gamma, delta in product((2*first-1, 1.), (2*second-1, 1.)):
                self.assertGreaterEqual(gamma*delta, 2*first*second-1 - 3e-16)
                self.assertLessEqual(gamma*delta, 1.)
                a = extension_noise_kraus(first, gamma)
                b = extension_noise_kraus(second, delta)
                composed = tuple(kb @ ka for ka in a for kb in b)
                for matrix in (*REAL_BASIS, PAULI_Y):
                    np.testing.assert_allclose(apply_operation(composed, matrix), diagonal_action(matrix, first*second, gamma*delta), atol=5e-16)

    def test_interval_certificates_have_positive_width_except_at_unit_strength(self):
        for eta in (Fraction(0), Fraction(1, 2), Fraction(17, 20), Fraction(999, 1000)):
            report = extension_interval_certificate(eta)
            self.assertGreater(Fraction(report["width_strict_lower_exact"]), 0)
            lower, upper = kappa_bounds(float(eta))
            self.assertAlmostEqual(report["kappa_lower_interval"][0], lower, places=14)
            self.assertAlmostEqual(report["kappa_upper_interval"][1], upper, places=14)
        with self.assertRaises(ValueError):
            kappa_bounds(1.1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(FreshExtensionClassificationTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"scope": "All real 2x2 Kraus extensions of each fixed original fresh branch; not all generalized probability theories",
              "noise_completion_interval_exact": "2*alpha-1 <= gamma <= 1",
              "fresh_branch_extension_interval_exact": "(2*alpha-1)*sqrt(1-eta^2)/2 <= kappa_s <= sqrt(1-eta^2)/2",
              "outcome_intervals_independent": True,
              "all_directions_same_interval": True,
              "unique_full_branch_extension_iff_eta_equals_one": True,
              "equal_local_maps_all_real_ancillas_max_output_trace_distance": "abs(kappa_Phi-kappa_Psi)/2; one real auxiliary saturates",
              "original_completion": "upper endpoint gamma=1",
              "isotropic_completion": "gamma=alpha",
              "certificates": [extension_interval_certificate(eta) for eta in (Fraction(0), Fraction(1, 2), Fraction(17, 20), Fraction(1))],
              "cognition_selects_kappa": False,
              "quantum_theory_derived_from_cognition": False,
              "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("fresh_extension_classification_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
