"""Round 59: local record completeness is an extra tomography assumption.

Unassisted protocols mean one unknown joint preparation, independent local
ancillas and classical communication. A preshared joint reference is different.
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
from complex_bipartite_closure import local_effect
from quantum_interface_audit import ALPHA, two_rebit_states
from role_symmetry_and_swap import pauli_word


REAL_WORDS = tuple(a+b for a, b in product("IXYZ", repeat=2) if (a+b).count("Y") % 2 == 0)
COMPLEX_WORDS = tuple(a+b for a, b in product("IXYZ", repeat=2))
PRODUCT_REAL_WORDS = tuple(a+b for a, b in product("IXZ", repeat=2))


def pairing_matrix(effect_words, state_words):
    return np.array([[np.trace(pauli_word(e) @ pauli_word(s)).real/4 for s in state_words] for e in effect_words])


def tensor_power(matrix, count):
    result = np.ones((1, 1), dtype=complex)
    for _ in range(count):
        result = np.kron(result, matrix)
    return result


def partial_transpose(density, sites, systems):
    order = list(range(2*systems))
    for site in sites:
        order[site], order[systems+site] = order[systems+site], order[site]
    return density.reshape((2,)*(2*systems)).transpose(order).reshape(density.shape)


def separated_pair_probabilities(density, a, b):
    return np.array([np.trace(density @ np.kron(local_effect(a, s), local_effect(b, t))).real
                     for s, t in product((-1, 1), repeat=2)])


class LocalRecordAxiomTests(unittest.TestCase):
    def test_exact_real_product_pairing_has_only_the_yy_kernel(self):
        matrix = pairing_matrix(PRODUCT_REAL_WORDS, REAL_WORDS)
        nonzero = [i for i in range(10) if np.any(matrix[:, i])]
        np.testing.assert_array_equal(matrix[:, nonzero], np.eye(9))
        self.assertEqual([REAL_WORDS[i] for i in range(10) if i not in nonzero], ["YY"])

    def test_complex_product_pairing_separates_the_full_hermitian_space(self):
        np.testing.assert_array_equal(pairing_matrix(COMPLEX_WORDS, COMPLEX_WORDS), np.eye(16))

    def test_physical_real_counterstates_are_positive_distinct_and_locally_identical(self):
        plus, minus = two_rebit_states()
        for rho in (plus, minus):
            np.testing.assert_allclose(np.linalg.eigvalsh(rho), [0, 0, .5, .5], atol=1e-16)
        for a, b in product("XZ", repeat=2):
            np.testing.assert_allclose(separated_pair_probabilities(plus, a, b), separated_pair_probabilities(minus, a, b), atol=1e-16, rtol=0)
        gate = interaction(math.pi/2)
        effect = gate.conj().T @ np.kron(np.eye(2), local_effect("Z", 1)) @ gate
        self.assertAlmostEqual(np.trace((plus-minus) @ effect).real, ALPHA, places=14)

    def test_arbitrary_real_local_filters_still_cannot_reveal_the_yy_difference(self):
        rng = np.random.default_rng(59)
        delta = two_rebit_states()[0]-two_rebit_states()[1]
        for _ in range(20):
            # Every fine-grained LOCC leaf has this form; random filters are
            # cross-checks of the proof Tr(Y K^T K)=0, not its replacement.
            a = np.eye(2)
            b = np.eye(2)
            for _ in range(4):
                a = rng.normal(size=(2, 2)) @ a/3
                b = rng.normal(size=(2, 2)) @ b/3
            effect = np.kron(a.T @ a, b.T @ b)
            self.assertAlmostEqual(np.trace(delta @ effect).real, 0, places=14)

    def test_marginal_summaries_do_not_determine_even_a_classical_joint_state(self):
        correlated = np.diag([.5, 0, 0, .5])
        anticorrelated = np.diag([0, .5, .5, 0])
        for rho in (correlated, anticorrelated):
            np.testing.assert_array_equal(np.einsum("abcb->ac", rho.reshape(2, 2, 2, 2)), np.eye(2)/2)
            np.testing.assert_array_equal(np.einsum("abad->bd", rho.reshape(2, 2, 2, 2)), np.eye(2)/2)
        distance = np.abs(separated_pair_probabilities(correlated, "Z", "Z")-separated_pair_probabilities(anticorrelated, "Z", "Z")).sum()/2
        self.assertAlmostEqual(distance, ALPHA**2, places=14)

    def test_locally_factorized_statistics_do_not_certify_independent_preparation(self):
        rho = (np.eye(4)+.7*pauli_word("YY"))/4
        for a, b in product("XZ", repeat=2):
            np.testing.assert_allclose(separated_pair_probabilities(rho, a, b), [.25]*4, atol=1e-16)
        self.assertGreater(np.linalg.norm(rho-np.eye(4)/4), .3)

    def test_opposite_q_states_remain_partial_transposes_for_arbitrary_sampled_copy_counts(self):
        plus, minus = two_rebit_states()
        for n in range(1, 4):
            np.testing.assert_array_equal(partial_transpose(tensor_power(plus, n), tuple(range(0, 2*n, 2)), 2*n), tensor_power(minus, n))

    def test_collective_real_local_copy_measurements_cannot_determine_the_sign(self):
        rng = np.random.default_rng(159)
        plus, minus = two_rebit_states()
        for n in range(1, 4):
            a = rng.normal(size=(2**n, 2**n))
            b = rng.normal(size=(2**n, 2**n))
            ea, eb = a.T @ a, b.T @ b
            ea /= np.trace(ea)
            eb /= np.trace(eb)
            effect = embed_operator(ea, tuple(range(0, 2*n, 2)), 2*n) @ embed_operator(eb, tuple(range(1, 2*n, 2)), 2*n)
            self.assertAlmostEqual(np.trace((tensor_power(plus, n)-tensor_power(minus, n)) @ effect).real, 0, places=15)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(LocalRecordAxiomTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"round": 59, "unassisted_real_product_rank": 9, "real_joint_dimension": 10, "kernel": ["YY"],
              "complex_product_rank": 16, "complex_joint_dimension": 16,
              "finite_locc_completeness_equivalent_to_product_separation": True,
              "scope": "Fixed-system operational state space; independent local ancillary preparations; no preshared joint resource",
              "marginal_state_sufficiency_fails_even_classically": True,
              "local_factorization_implies_independent_preparation": False,
              "opposite_q_sign_hidden_under_unassisted_collective_real_local_copies": True,
              "quantum_theory_derived_from_cognition": False,
              "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("local_record_axiom_results.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
