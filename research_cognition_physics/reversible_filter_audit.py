"""Round 182: homogeneous cones and probabilistically reversible filters.

Forward and reverse success branches compose to a scalar identity, including
external references. Failure branches remain part of the physical instrument.
"""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import unittest

import numpy as np


def positive_root(matrix, inverse=False):
    values, vectors = np.linalg.eigh(matrix)
    if inverse and values.min() <= 0:
        raise ValueError("Inverse square root requires a positive definite matrix.")
    powers = 1 / np.sqrt(values) if inverse else np.sqrt(np.maximum(values, 0))
    return (vectors * powers) @ vectors.conj().T


def reversible_filters(rho, sigma):
    operator = positive_root(sigma) @ positive_root(rho, inverse=True)
    inverse = np.linalg.inv(operator)
    c = min(1., 1 / np.linalg.norm(operator, 2) ** 2)
    d = min(1., 1 / np.linalg.norm(inverse, 2) ** 2)
    forward, backward = np.sqrt(c) * operator, np.sqrt(d) * inverse
    failure = positive_root(np.eye(len(rho)) - forward.conj().T @ forward)
    return forward, backward, failure, c, d


class ReversibleFilterAuditTests(unittest.TestCase):
    def test_classical_full_support_conversion_and_reversal_exactly(self):
        F = Fraction
        p, q = [F(1, 2), F(1, 2)], [F(3, 4), F(1, 4)]
        ratios = [q_i / p_i for p_i, q_i in zip(p, q)]
        c = min(p_i / q_i for p_i, q_i in zip(p, q))
        d = min(q_i / p_i for p_i, q_i in zip(p, q))
        self.assertEqual((c, d, c * d), (F(2, 3), F(1, 2), F(1, 3)))
        self.assertEqual([c * t * value for t, value in zip(ratios, p)], [c * value for value in q])
        for t in ratios:
            self.assertEqual(c * t * d / t, F(1, 3))
            self.assertTrue(0 <= c * t <= 1)
            self.assertTrue(0 <= d / t <= 1)

    def test_deterministic_classical_reversibility_does_not_connect_mixed_states(self):
        p, q = np.array([.5, .5]), np.array([.75, .25])
        for permutation in (np.eye(2), np.array([[0, 1], [1, 0]])):
            np.testing.assert_array_equal(permutation @ p, p)
            self.assertFalse(np.array_equal(permutation @ p, q))

    def test_real_and_complex_known_state_filters_preserve_their_contract(self):
        rng = np.random.default_rng(182)
        for complex_model in (False, True):
            for dimension in (2, 3):
                states = []
                for _ in range(2):
                    raw = rng.normal(size=(dimension, dimension))
                    if complex_model:
                        raw = raw + 1j * rng.normal(size=raw.shape)
                    state = raw @ raw.conj().T + np.eye(dimension)
                    states.append(state / np.trace(state))
                rho, sigma = states
                a, b, failure, c, d = reversible_filters(rho, sigma)
                np.testing.assert_allclose(a @ rho @ a.conj().T, c * sigma, atol=2e-15)
                np.testing.assert_allclose(b @ sigma @ b.conj().T, d * rho, atol=2e-15)
                np.testing.assert_allclose(b @ a, np.sqrt(c * d) * np.eye(dimension), atol=2e-15)
                np.testing.assert_allclose(a.conj().T @ a + failure.conj().T @ failure,
                                           np.eye(dimension), atol=2e-15)
                self.assertLessEqual(np.linalg.norm(b, 2), 1 + 1e-15)
                if not complex_model:
                    self.assertFalse(np.iscomplexobj(a))

    def test_successful_round_trip_preserves_unknown_external_entanglement(self):
        rho, sigma = np.eye(2) / 2, np.array([[.6, .2j], [-.2j, .4]])
        a, b, _, c, d = reversible_filters(rho, sigma)
        psi = np.array([1, 1j, -2, .5, 0, 3j], dtype=complex)
        psi /= np.linalg.norm(psi)
        joint = np.outer(psi, psi.conj())
        aa, bb = np.kron(a, np.eye(3)), np.kron(b, np.eye(3))
        result = bb @ aa @ joint @ aa.conj().T @ bb.conj().T
        np.testing.assert_allclose(result, c * d * joint, atol=4e-16)

    def test_failure_branch_cannot_be_omitted_from_unconditional_statistics(self):
        rho, sigma = np.eye(2) / 2, np.diag([.75, .25])
        a, _, failure, c, _ = reversible_filters(rho, sigma)
        success = a @ rho @ a.T
        failed = failure @ rho @ failure.T
        self.assertAlmostEqual(np.trace(success), 2 / 3)
        self.assertAlmostEqual(np.trace(failed), 1 / 3)
        np.testing.assert_allclose(success + failed, rho, atol=3e-16)
        self.assertGreater(np.linalg.norm(success / c - rho), .3)

    def test_conditioned_normalization_is_not_an_affine_deterministic_map(self):
        a, _, _, _, _ = reversible_filters(np.eye(2) / 2, np.diag([.75, .25]))
        def conditioned(state):
            result = a @ state @ a.T
            return result / np.trace(result)
        p0, p1 = np.diag([1., 0.]), np.diag([0., 1.])
        mixture_after = (conditioned(p0) + conditioned(p1)) / 2
        conditioned_mixture = conditioned((p0 + p1) / 2)
        self.assertGreater(np.linalg.norm(mixture_after - conditioned_mixture), .3)

    def test_boundary_ranks_and_vanishing_success_are_explicit(self):
        with self.assertRaises(ValueError):
            reversible_filters(np.diag([1., 0.]), np.eye(2) / 2)
        for epsilon in (1e-1, 1e-2, 1e-4):
            rho, sigma = np.diag([1 - epsilon, epsilon]), np.eye(2) / 2
            _, _, _, c, _ = reversible_filters(rho, sigma)
            self.assertAlmostEqual(c, 2 * epsilon)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(ReversibleFilterAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = {
        "round": 182,
        "quantum_order_automorphism": "X -> A X A^dagger, A=sqrt(sigma)*rho^(-1/2)",
        "classical_example_forward_success": "2/3",
        "classical_example_reverse_conditional_success": "1/2",
        "all_input_round_trip_success": "1/3",
        "success_branch_preserves_external_relations_on_reversal": True,
        "homogeneity_means_deterministic_normalized_reversibility": False,
        "arbitrary_cone_automorphisms_are_physically_allowed_by_definition": False,
        "fixed_protocol_interior_to_boundary_uniform_success_bound": False,
        "scope": "Positive definite known design states; selected two-success branch, not unconditional recovery",
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("reversible_filter_audit_results.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
