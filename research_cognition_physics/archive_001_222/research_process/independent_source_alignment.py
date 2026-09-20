"""Round 62: align two independent real references through one middle party.

Order before the middle measurement is A, B_left, B_right, C. After discarding
B_right it is A, B, C. Source tensor independence and feedforward are explicit.
"""

import argparse
import json
import math
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from bilocal_record_tomography import embed_operator
from bipartite_composition import interaction
from quantum_interface_audit import ALPHA, PAULI_X, branch_kraus
from role_symmetry_and_swap import pauli_word
from shared_real_reference_tomography import shared_reference


def keep_systems(density, keep, systems):
    keep = tuple(keep)
    if len(set(keep)) != len(keep) or any(i < 0 or i >= systems for i in keep):
        raise ValueError("Keep distinct system indices in range.")
    rest = tuple(i for i in range(systems) if i not in keep)
    order = keep+rest+tuple(systems+i for i in keep+rest)
    grouped = density.reshape((2,)*(2*systems)).transpose(order)
    return np.einsum("abcb->ac", grouped.reshape(2**len(keep), 2**len(rest), 2**len(keep), 2**len(rest)))


def pair_reference(correlation=1.):
    if not -1 <= correlation <= 1:
        raise ValueError("Use correlation in [-1,1].")
    return (np.eye(4)+correlation*pauli_word("YY"))/4


def independent_sources(left=1., right=1.):
    return np.kron(pair_reference(left), pair_reference(right))


def after_middle_gate(left=1., right=1.):
    gate = embed_operator(interaction(math.pi/2), (1, 2), 4)
    return gate @ independent_sources(left, right) @ gate.conj().T


def middle_branch(sign, contrast=1., left=1., right=1.):
    if sign not in (-1, 1):
        raise ValueError("Use sign -1 or +1.")
    density = after_middle_gate(left, right)
    output = np.zeros((16, 16), dtype=complex)
    for k in branch_kraus(0., (sign+1)//2, contrast):
        extended = embed_operator(k, (2,), 4)
        output += extended @ density @ extended.conj().T
    return keep_systems(output, (0, 1, 3), 4)


def conditional_reference(sign, beta, left=1., right=1.):
    if sign not in (-1, 1) or not 0 <= beta <= 1:
        raise ValueError("Use a binary sign and beta in [0,1].")
    pair_reference(left)
    pair_reference(right)
    return (np.eye(8)+left*pauli_word("YYI")+sign*beta*right*pauli_word("IYY")
            +sign*beta*left*right*pauli_word("YIY"))/8


def real_x_with_reset(density, site, systems):
    """An old YX(pi) flow on reset helper/target gives X on the target."""
    joint = np.kron(np.diag([1., 0]), density)
    gate = embed_operator(interaction(math.pi), (0, site+1), systems+1)
    joint = gate @ joint @ gate.conj().T
    return keep_systems(joint, tuple(range(1, systems+1)), systems+1)


def feedforward_reference(contrast=1., left=1., right=1.):
    return middle_branch(1, contrast, left, right)+real_x_with_reset(middle_branch(-1, contrast, left, right), 2, 3)


class IndependentSourceAlignmentTests(unittest.TestCase):
    def test_sources_are_an_actual_tensor_product_and_outer_marginal_is_independent(self):
        for left, right in product((0., .4, 1.), repeat=2):
            rho = independent_sources(left, right)
            np.testing.assert_allclose(keep_systems(rho, (0, 1), 4), pair_reference(left), atol=0)
            np.testing.assert_allclose(keep_systems(rho, (2, 3), 4), pair_reference(right), atol=0)
            np.testing.assert_allclose(keep_systems(rho, (0, 3), 4), np.eye(4)/4, atol=0)

    def test_actual_middle_gate_and_full_old_instrument_match_all_coefficients(self):
        for sign, eta, left, right in product((-1, 1), (0., .4, 1.), (.3, 1.), (.2, 1.)):
            actual = middle_branch(sign, eta, left, right)
            expected = conditional_reference(sign, ALPHA*eta, left, right)/2
            np.testing.assert_allclose(actual, expected, atol=3e-16)

    def test_each_middle_flag_has_probability_half(self):
        for sign, eta in product((-1, 1), (.1, .7, 1.)):
            self.assertAlmostEqual(np.trace(middle_branch(sign, eta)).real, .5, places=14)

    def test_local_x_correction_is_compiled_from_a_real_reset_and_an_old_flow(self):
        rng = np.random.default_rng(62)
        matrix = rng.normal(size=(8, 8))+1j*rng.normal(size=(8, 8))
        rho = matrix @ matrix.conj().T
        rho /= np.trace(rho)
        for site in range(3):
            target = embed_operator(PAULI_X, (site,), 3)
            np.testing.assert_allclose(real_x_with_reset(rho, site, 3), target @ rho @ target, atol=2e-16)

    def test_feedforward_aligns_both_flag_branches_to_the_same_normalized_state(self):
        for eta in (.2, .85, 1.):
            positive = 2*middle_branch(1, eta)
            corrected_negative = 2*real_x_with_reset(middle_branch(-1, eta), 2, 3)
            np.testing.assert_allclose(positive, corrected_negative, atol=3e-16)
            np.testing.assert_allclose(feedforward_reference(eta), conditional_reference(1, ALPHA*eta), atol=3e-16)

    def test_forgetting_the_flag_without_feedback_removes_outer_alignment(self):
        rho = middle_branch(1)+middle_branch(-1)
        np.testing.assert_allclose(rho, np.kron(pair_reference(), np.eye(2)/2), atol=3e-16)
        np.testing.assert_allclose(keep_systems(rho, (0, 2), 3), np.eye(4)/4, atol=2e-16)

    def test_corrected_reference_has_the_claimed_exact_distance_and_local_marginals(self):
        for beta in (0., .2, ALPHA, 1.):
            rho = conditional_reference(1, beta)
            distance = np.abs(np.linalg.eigvalsh(rho-shared_reference(3))).sum()/2
            self.assertAlmostEqual(distance, (1-beta)/2, places=14)
            for site in range(3):
                np.testing.assert_allclose(keep_systems(rho, (site,), 3), np.eye(2)/2, atol=0)

    def test_weak_source_errors_and_detector_error_are_all_retained(self):
        for left, right, beta in product((.2, .7, 1.), repeat=3):
            rho = conditional_reference(1, beta, left, right)
            distance = np.abs(np.linalg.eigvalsh(rho-shared_reference(3))).sum()/2
            self.assertAlmostEqual(distance, 1-(1+left)*(1+beta*right)/4, places=14)

    def test_all_conditional_states_remain_real_positive_and_normalized(self):
        for sign, beta, left, right in product((-1, 1), (0., .3, 1.), (-1., .2, 1.), (-1., .6, 1.)):
            rho = conditional_reference(sign, beta, left, right)
            np.testing.assert_array_equal(rho.imag, np.zeros((8, 8)))
            self.assertGreaterEqual(np.linalg.eigvalsh(rho).min(), -2e-16)
            self.assertAlmostEqual(np.trace(rho).real, 1, places=14)

    def test_invalid_resources_and_flags_are_rejected(self):
        with self.assertRaises(ValueError):
            pair_reference(1.1)
        with self.assertRaises(ValueError):
            middle_branch(0)
        with self.assertRaises(ValueError):
            conditional_reference(1, -.1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(IndependentSourceAlignmentTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"round": 62, "preparation": "tau_A_BL tensor tau_BR_C, independent sources",
              "middle_operation": "Old YX(pi/2) on BL,BR; old fresh Z read on BR; discard BR",
              "flag_probabilities": [0.5, 0.5], "conditional_state": "[III+YYI+r beta(IYY+YIY)]/8",
              "feedforward": "Send r to C; apply a locally compiled real X when r=-1",
              "postselection_required": False, "one_read_reference_distance": (1-ALPHA)/2,
              "ignored_flag_outer_state": "I_A/2 tensor I_C/2",
              "source_independence_alone_forbids_alignment": False,
              "quantum_theory_derived_from_cognition": False,
              "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("independent_source_alignment_results.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
