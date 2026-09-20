"""Bipartite disk interfaces: hidden joint memory and reversible interactions.

The real-quantum composite and its joint gates are explicit candidate additions,
not assumptions derived from cognition. See research_note_47.md for scope.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import permutations, product
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I, SCALE, sin_interval
from outcome_updates import WINDOW_ATTENUATION as ALPHA
from quantum_interface_audit import (IDENTITY, PAULI_X, PAULI_Y, PAULI_Z,
                                     branch_kraus, density_from_summary,
                                     effective_readout, two_rebit_states)
from retention_tradeoff import fresh_branch_matrix


LOCAL_BASIS = (IDENTITY, PAULI_Z, PAULI_X)
VISIBLE_BASIS = tuple(np.kron(a, b) for a in LOCAL_BASIS for b in LOCAL_BASIS)
HIDDEN = np.kron(PAULI_Y, PAULI_Y)
FULL_BASIS = VISIBLE_BASIS + (HIDDEN,)
GENERATOR_LABELS = ("YI", "IY", "YX", "YZ", "XY", "ZY")
GENERATORS = tuple(np.kron(a, b) for a, b in (
    (PAULI_Y, IDENTITY), (IDENTITY, PAULI_Y), (PAULI_Y, PAULI_X),
    (PAULI_Y, PAULI_Z), (PAULI_X, PAULI_Y), (PAULI_Z, PAULI_Y)))


def visible_summary(density):
    return np.array([np.trace(operator @ density).real for operator in VISIBLE_BASIS]).reshape(3, 3)


def full_summary(density):
    return np.array([np.trace(operator @ density).real for operator in FULL_BASIS])


def reconstruct(visible, hidden_coordinate):
    coefficients = np.r_[np.asarray(visible).reshape(9), hidden_coordinate]
    return sum(value * operator for value, operator in zip(coefficients, FULL_BASIS)) / 4


def local_effect_vector(setting, outcome, contrast):
    sign, beta = 2 * outcome - 1, ALPHA * contrast
    return np.array([1., sign * beta * math.cos(setting), sign * beta * math.sin(setting)]) / 2


def joint_probability(visible, setting_a, outcome_a, setting_b, outcome_b, contrast=1.):
    return float(local_effect_vector(setting_a, outcome_a, contrast) @ visible
                 @ local_effect_vector(setting_b, outcome_b, contrast))


def disk_branch(setting, outcome, contrast):
    scale = np.diag([1., ALPHA, ALPHA])
    return np.linalg.solve(scale, fresh_branch_matrix(setting, outcome, contrast) @ scale)


def local_update(density, party, setting, outcome, contrast):
    if party not in (0, 1):
        raise ValueError("Use party 0 or 1.")
    local = branch_kraus(setting, outcome, contrast)
    operators = [np.kron(k, IDENTITY) if party == 0 else np.kron(IDENTITY, k) for k in local]
    return sum(k @ density @ k.conj().T for k in operators)


def interaction(angle, generator="YX"):
    if generator not in GENERATOR_LABELS:
        raise ValueError("Unknown real-orthogonal generator.")
    h = GENERATORS[GENERATOR_LABELS.index(generator)]
    # Every listed Hermitian Pauli string squares to the identity.
    return math.cos(angle / 2) * np.eye(4) - 1j * math.sin(angle / 2) * h


def joint_update(density, angle, generator="YX"):
    u = interaction(angle, generator)
    return u @ density @ u.conj().T


def full_generator(h):
    return np.array([[np.trace(a @ (-.5j * (h @ b - b @ h))).real / 4
                      for b in FULL_BASIS] for a in FULL_BASIS])


def leakage_matrix():
    return np.column_stack([full_generator(h)[:9, 9] for h in GENERATORS])


def integer_determinant(matrix):
    size = len(matrix)
    total = 0
    for order in permutations(range(size)):
        inversions = sum(order[i] > order[j] for i in range(size) for j in range(i + 1, size))
        term = (-1)**inversions
        for i, j in enumerate(order):
            term *= matrix[i][j]
        total += term
    return total


def exact_generator_certificate():
    # Analytic Pauli commutators give this integer table. The tests below also
    # compare it entry by entry with the matrix implementation.
    table = [[0] * 6 for _ in range(9)]
    for row, column, sign in ((1, 2, 1), (2, 3, -1), (3, 4, 1), (6, 5, -1)):
        table[row][column] = sign
    minor = [[table[row][column] for column in (2, 3, 4, 5)] for row in (1, 2, 3, 6)]
    determinant = integer_determinant(minor)
    return {"verified": determinant == 1, "generator_order": list(GENERATOR_LABELS),
            "visible_row_order": ["II", "IZ", "IX", "ZI", "ZZ", "ZX", "XI", "XZ", "XX"],
            "hidden_to_visible_integer_matrix": table, "nonzero_minor_determinant": determinant,
            "rank_exact": 4, "kernel_dimension_exact": 2,
            "autonomous_quotient_generators": ["YI", "IY"],
            "scope": "Classification within the full real-orthogonal two-rebit dynamics. It does not classify all possible disk composites."}


def bell_state():
    # The candidate joint gate prepares Phi+ from an allowed product preparation.
    ket00 = np.array([1, 0, 0, 0], dtype=complex)
    return joint_update(np.outer(ket00, ket00.conj()), math.pi / 2)


def prepare_hidden_pair(sign):
    ket00, ket01 = np.eye(4, dtype=complex)[:2]
    product00, product01 = np.outer(ket00, ket00.conj()), np.outer(ket01, ket01.conj())
    if sign not in (-1, 1):
        raise ValueError("Use sign -1 or +1.")
    return (joint_update(product00, -sign * math.pi / 2)
            + joint_update(product01, sign * math.pi / 2)) / 2


def correlation(density, angle_a, angle_b, contrast=1.):
    visible = visible_summary(density)
    return sum((2 * a - 1) * (2 * b - 1) * joint_probability(visible, angle_a, a, angle_b, b, contrast)
               for a, b in product((0, 1), repeat=2))


def chsh_value(density, contrast=1.):
    a0, a1, b0, b1 = 0., math.pi / 2, math.pi / 4, -math.pi / 4
    return (correlation(density, a0, b0, contrast) + correlation(density, a0, b1, contrast)
            + correlation(density, a1, b0, contrast) - correlation(density, a1, b1, contrast))


def bell_certificate(contrast=Fraction(17, 20)):
    eta = Fraction(contrast)
    if not 0 <= eta <= 1:
        raise ValueError("Use contrast in [0,1].")
    alpha = 4 * sin_interval(I.rational(1, 4))
    value = 2 * I.exact(2).sqrt() * (alpha * eta)**2
    margin = value - 2
    return {"contrast_exact": str(eta), "value_interval": value.floats(),
            "margin_above_general_local_bound_exact": str(Fraction(margin.lo, SCALE)),
            "violates_general_local_bound": margin.lo > 0,
            "scope": "Prediction of the explicitly added joint-state/gate candidate, with independent settings and no communication during the separated readouts; not a prediction already derived from cognition."}


class BipartiteCompositionTests(unittest.TestCase):
    def test_ten_coordinates_reconstruct_real_joint_states_and_zero_filling_can_fail(self):
        gram = np.array([[np.trace(a @ b).real / 4 for b in FULL_BASIS] for a in FULL_BASIS])
        np.testing.assert_array_equal(gram, np.eye(10))
        rng = np.random.default_rng(47)
        for _ in range(12):
            matrix = rng.normal(size=(4, 4))
            density = matrix @ matrix.T
            density /= np.trace(density)
            coordinates = full_summary(density)
            np.testing.assert_allclose(reconstruct(coordinates[:9].reshape(3, 3), coordinates[9]), density, atol=2e-16)
        projected = reconstruct(visible_summary(bell_state()), 0)
        np.testing.assert_allclose(np.linalg.eigvalsh(projected), [-.25, .25, .25, .75], atol=3e-16)

    def test_product_probabilities_and_branch_updates_use_the_same_nine_coordinates(self):
        states = (bell_state(), *two_rebit_states(), np.eye(4) / 4)
        for density, party, outcome in product(states, (0, 1), (0, 1)):
            setting, eta = .31, .73
            visible = visible_summary(density)
            branch = disk_branch(setting, outcome, eta)
            actual = visible_summary(local_update(density, party, setting, outcome, eta))
            expected = branch @ visible if party == 0 else visible @ branch.T
            np.testing.assert_allclose(actual, expected, atol=4e-16)
        for density, a, b in product(states, (0, 1), (0, 1)):
            effect = np.kron(effective_readout(.23, a, .8), effective_readout(-.41, b, .8))
            actual = joint_probability(visible_summary(density), .23, a, -.41, b, .8)
            self.assertAlmostEqual(actual, np.trace(effect @ density).real, places=14)

    def test_local_measurement_stage_is_normalized_positive_and_nonsignaling(self):
        for density in (bell_state(), *two_rebit_states()):
            v = visible_summary(density)
            for aa, bb, eta in product((0., .4, 1.1), (-.3, .7), (.15, .85, 1.)):
                probabilities = np.array([[joint_probability(v, aa, a, bb, b, eta) for b in (0, 1)] for a in (0, 1)])
                self.assertGreaterEqual(probabilities.min(), 0)
                self.assertAlmostEqual(probabilities.sum(), 1)
                np.testing.assert_allclose(probabilities.sum(axis=0),
                                           [local_effect_vector(bb, b, eta) @ v[0] for b in (0, 1)], atol=3e-16)
                np.testing.assert_allclose(probabilities.sum(axis=1),
                                           [local_effect_vector(aa, a, eta) @ v[:, 0] for a in (0, 1)], atol=3e-16)

    def test_local_adaptive_protocols_with_classical_messages_cannot_reveal_hidden_coordinate(self):
        totals = np.zeros(2)
        for record in product((0, 1), repeat=5):
            outputs = list(two_rebit_states())
            for step, outcome in enumerate(record):
                # Settings may depend on the complete earlier public record.
                setting = .17 * step + .29 * sum(record[:step])
                outputs = [local_update(density, step % 2, setting, outcome, (.3, 1., .7, .5, .9)[step]) for density in outputs]
                np.testing.assert_allclose(visible_summary(outputs[0]), visible_summary(outputs[1]), atol=3e-16)
            traces = np.array([np.trace(density).real for density in outputs])
            self.assertAlmostEqual(traces[0], traces[1], places=14)
            totals += traces
        np.testing.assert_allclose(totals, 1, atol=5e-16)

    def test_single_reversible_gate_reveals_hidden_memory_with_existing_noisy_readout(self):
        for angle in (0., .13, .8, math.pi / 2):
            u = interaction(angle)
            np.testing.assert_allclose(u.imag, 0, atol=0)
            np.testing.assert_allclose(interaction(-angle) @ u, np.eye(4), atol=3e-16)
            for sign, density in zip((1, -1), two_rebit_states()):
                output = joint_update(density, angle)
                expected = (np.eye(4) + sign * (math.cos(angle) * HIDDEN + math.sin(angle) * np.kron(IDENTITY, PAULI_Z))) / 4
                np.testing.assert_allclose(output, expected, atol=3e-16)
                for eta in (.15, 1.):
                    effect = np.kron(IDENTITY, effective_readout(0, 1, eta))
                    self.assertAlmostEqual(np.trace(output @ effect).real, (1 + sign * ALPHA * eta * math.sin(angle)) / 2, places=14)

    def test_hidden_states_are_preparable_from_products_with_the_declared_gate(self):
        for sign, target in zip((1, -1), two_rebit_states()):
            np.testing.assert_allclose(prepare_hidden_pair(sign), target, atol=3e-16)
        self.assertAlmostEqual(np.trace(bell_state() @ HIDDEN).real, -1)

    def test_generator_classification_has_an_exact_nonzero_integer_minor(self):
        report = exact_generator_certificate()
        self.assertTrue(report["verified"])
        np.testing.assert_array_equal(leakage_matrix(), report["hidden_to_visible_integer_matrix"])
        self.assertEqual(report["nonzero_minor_determinant"], 1)
        self.assertEqual(report["autonomous_quotient_generators"], ["YI", "IY"])

    def test_ten_coordinate_dynamics_closes_but_interaction_derivatives_need_the_tenth(self):
        density = bell_state()
        for h in GENERATORS:
            generator = full_generator(h)
            np.testing.assert_array_equal(generator + generator.T, np.zeros((10, 10)))
            derivative = -.5j * (h @ density - density @ h)
            np.testing.assert_allclose(generator @ full_summary(density), full_summary(derivative), atol=3e-16)
        pulled = interaction(math.pi / 2).conj().T @ np.kron(IDENTITY, PAULI_Z) @ interaction(math.pi / 2)
        np.testing.assert_allclose(pulled, HIDDEN, atol=3e-16)

    def test_separable_disk_witness_and_general_bell_bound_are_distinct(self):
        # This separable preparation saturates Tr(T)<=1; it is not the class of
        # every possible Bell-local hidden-variable model.
        z = density_from_summary([1, ALPHA, 0])
        x = density_from_summary([1, 0, ALPHA])
        separable = (np.kron(z, z) + np.kron(x, x)) / 2
        self.assertAlmostEqual(np.trace(visible_summary(separable)[1:, 1:]), 1)
        self.assertAlmostEqual(np.trace(visible_summary(bell_state())[1:, 1:]), 2)
        scores = [a0*b0+a0*b1+a1*b0-a1*b1 for a0,a1,b0,b1 in product((-1,1), repeat=4)]
        self.assertEqual(max(scores), 2)
        self.assertEqual(min(scores), -2)

    def test_bell_candidate_violates_the_local_bound_with_the_actual_readout_noise(self):
        for eta in (.15, .85, 1.):
            self.assertAlmostEqual(chsh_value(bell_state(), eta), 2 * math.sqrt(2) * (ALPHA * eta)**2, places=14)
        self.assertTrue(bell_certificate()["violates_general_local_bound"])
        self.assertGreater(Fraction(bell_certificate()["margin_above_general_local_bound_exact"]), 0)
        self.assertTrue(bell_certificate(Fraction(1))["violates_general_local_bound"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(BipartiteCompositionTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"generator_certificate": exact_generator_certificate(),
              "homogeneous_predictive_dimensions": {"local_protocol_quotient": 9, "with_declared_interaction": 10},
              "new_joint_coordinate": "Tr(rho sigma_y tensor sigma_y)",
              "hidden_pair_total_variation_after_gate": "alpha*eta*abs(sin(theta))",
              "hidden_pair_single_read_success_at_eta_1_theta_pi_over_2": (1 + ALPHA) / 2,
              "zero_fill_bell_projection_eigenvalues": [-.25, .25, .25, .75],
              "bell_certificates": [bell_certificate(Fraction(17,20)), bell_certificate(Fraction(1))],
              "bell_contrast_threshold_diagnostic": 2**(-.25) / ALPHA,
              "quantum_theory_derived_from_cognition": False,
              "universal_disk_no_interaction_theorem_proven_here": False,
              "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("bipartite_composition_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
