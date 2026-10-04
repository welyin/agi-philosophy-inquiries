"""Round 195: a fixed Kahler carrier supports exact classical memory processes.

This is an embedding with declared instruments, not a derivation of quantum
theory or an identification of mixed density matrices with pure-state rays.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
import json
from pathlib import Path
import unittest

import numpy as np
import history_process_audit as classical


def permutation(action, noise_bit=0):
    matrix = np.zeros((4, 4), dtype=complex)
    for column, state in enumerate(classical.STATES):
        successor, _ = classical.step(state, action, noise_bit)
        matrix[classical.STATES.index(successor), column] = 1
    return matrix


def instruments(action, noise=F(0)):
    if action not in classical.ACTIONS or not 0 <= noise <= 1:
        raise ValueError("Invalid command or noise.")
    if action == "read":
        return {("read", bit): [np.diag([int(s[0] == bit) for s in classical.STATES])]
                for bit in (0, 1)}
    weights = ((0, 1-noise), (1, noise)) if action == "tick" else ((0, F(1)),)
    return {("done", action): [np.sqrt(float(weight))*permutation(action, bit)
                               for bit, weight in weights if weight]}


def branch_density(rho, operators):
    return sum((operator @ rho @ operator.conj().T for operator in operators),
               np.zeros_like(rho, dtype=complex))


def realification(matrix):
    return np.block([[matrix.real, -matrix.imag], [matrix.imag, matrix.real]])


def encoded_distribution(initial):
    return np.diag([float(initial.get(state, 0)) for state in classical.STATES]).astype(complex)


def quantum_transcripts(initial, policy, horizon, noise=F(0)):
    active = {(): encoded_distribution(initial)}
    for _ in range(horizon):
        following = {}
        for history, rho in active.items():
            for record, operators in instruments(policy(history), noise).items():
                updated = branch_density(rho, operators)
                if np.trace(updated).real > 0:
                    following[history+(record,)] = updated
        active = following
    return {history: float(np.trace(rho).real) for history, rho in active.items()}


def visible_marginal(rho):
    return np.einsum("imjm->ij", rho.reshape(2, 2, 2, 2))


class GeometryProcessBridgeTests(unittest.TestCase):
    def test_same_geometry_supports_distinct_hamiltonians(self):
        identity = np.eye(2)
        sigma_x = np.array([[0, 1], [1, 0]])
        j = np.block([[np.zeros((2, 2)), -identity], [identity, np.zeros((2, 2))]])
        omega = -j
        psi = np.array([1, 0], dtype=complex)
        for time in (0.0, 0.7, np.pi):
            moving = np.cos(time/2)*identity-1j*np.sin(time/2)*sigma_x
            real_u = realification(moving)
            np.testing.assert_allclose(real_u.T @ real_u, np.eye(4), atol=1e-14)
            np.testing.assert_allclose(real_u.T @ omega @ real_u, omega, atol=1e-14)
            np.testing.assert_allclose(real_u @ j, j @ real_u, atol=1e-14)
            self.assertAlmostEqual(abs((moving @ psi)[1])**2, np.sin(time/2)**2)
        self.assertEqual(abs(psi[1])**2, 0)
        self.assertAlmostEqual(abs(((-1j*sigma_x) @ psi)[1])**2, 1)

    def test_memory_permutations_preserve_common_linear_kahler_structure(self):
        identity, zero = np.eye(4), np.zeros((4, 4))
        j = np.block([[zero, -identity], [identity, zero]])
        for action in classical.ACTIONS:
            for bit in (0, 1):
                u = permutation(action, bit)
                np.testing.assert_array_equal(u.conj().T @ u, identity)
                real_u = realification(u)
                np.testing.assert_array_equal(real_u.T @ real_u, np.eye(8))
                np.testing.assert_array_equal(real_u @ j, j @ real_u)
                np.testing.assert_array_equal(real_u.T @ (-j) @ real_u, -j)
        tick = permutation("tick")
        np.testing.assert_array_equal(tick @ tick @ tick, identity)

    def test_every_outcome_and_unnormalized_successor_matches_classical_kernel(self):
        for action in classical.ACTIONS:
            for noise in (F(0), F(1, 5), F(1, 2)):
                operations = instruments(action, noise)
                completeness = sum((k.conj().T @ k for ks in operations.values() for k in ks),
                                   np.zeros((4, 4), dtype=complex))
                np.testing.assert_allclose(completeness, np.eye(4), atol=1e-14)
                for state in classical.STATES:
                    rho = encoded_distribution({state: F(1)})
                    expected = defaultdict(lambda: np.zeros(4))
                    for (successor, record), weight in classical.branches(state, action, noise).items():
                        expected[record][classical.STATES.index(successor)] += float(weight)
                    for record, operators in operations.items():
                        np.testing.assert_allclose(branch_density(rho, operators),
                                                   np.diag(expected[record]), atol=1e-14)

    def test_adaptive_full_records_match_with_hidden_noise(self):
        initial = {state: F(i+1, 10) for i, state in enumerate(classical.STATES)}
        def policy(history):
            if len(history) % 2 == 0:
                return "read"
            return "tick" if history[-1][-1] else "flip_memory"
        for noise in (F(0), F(1, 5), F(1, 2)):
            expected = classical.transcript_distribution(initial, policy, 7, noise)
            actual = quantum_transcripts(initial, policy, 7, noise)
            self.assertEqual(expected.keys(), actual.keys())
            for record in expected:
                self.assertAlmostEqual(actual[record], float(expected[record]))

    def test_embedding_preserves_declared_mixtures_and_independent_composition(self):
        p, q = np.array([.1, .2, .3, .4]), np.array([.4, .3, .2, .1])
        np.testing.assert_allclose(np.diag(.3*p+.7*q), .3*np.diag(p)+.7*np.diag(q))
        np.testing.assert_array_equal(np.diag(np.kron(p, q)), np.kron(np.diag(p), np.diag(q)))
        joint = np.diag(np.kron(p, q))
        operation = np.kron(permutation("tick"), permutation("swap"))
        evolved = operation @ joint @ operation.conj().T
        p_new = permutation("tick") @ p
        q_new = permutation("swap") @ q
        np.testing.assert_array_equal(evolved, np.diag(np.kron(p_new, q_new)))

    def test_discarding_memory_breaks_update_despite_identical_current_marginals(self):
        first = encoded_distribution({(0, 0): F(1)})
        second = encoded_distribution({(0, 1): F(1)})
        np.testing.assert_array_equal(visible_marginal(first), visible_marginal(second))
        u = permutation("tick")
        first_next = visible_marginal(u @ first @ u.conj().T)
        second_next = visible_marginal(u @ second @ u.conj().T)
        np.testing.assert_array_equal(first_next, np.diag([1, 0]))
        np.testing.assert_array_equal(second_next, np.diag([0, 1]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checked = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(GeometryProcessBridgeTests))
    if not checked.wasSuccessful():
        raise SystemExit(1)
    result = {
        "round": 195,
        "scope": "Exact affine classical process embedding, not derivation of quantum structure",
        "linear_kahler_carrier": "C^4 = R^8 with g=Re inner product, omega=Im inner product, J=i",
        "pure_ray_space_of_carrier": "CP^3; only four basis rays needed by the deterministic classical labels",
        "encoded_ensemble_states": "diag(p); mixed matrices are not single CP^3 points",
        "same_geometry_distinct_processes": True,
        "complete_selected_interface_history_preserved": True,
        "noise_source_publicly_readable": False,
        "memoryless_local_projection_closed_for_all_preparations": False,
        "noiseless_process_embedding_requires_equilibrium": False,
        "full_soca_or_quantum_necessity_proved": False,
        "automated_tests": checked.testsRun,
    }
    if args.write_results:
        Path(__file__).with_name("geometry_process_bridge_results.json").write_text(
            json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"tests": checked.testsRun, "round": 195}))


if __name__ == "__main__":
    main()
