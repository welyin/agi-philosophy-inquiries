"""Round 176: marginal-only autonomous prediction is stronger than cut consistency."""

import argparse
from itertools import permutations
import json
from pathlib import Path
import unittest

import numpy as np

from cut_consistency_audit import trace_keep


def classical_marginal(dims, side):
    reduction = np.zeros((dims[side], np.prod(dims)), dtype=int)
    for flat, pair in enumerate(np.ndindex(*dims)):
        reduction[pair[side], flat] = 1
    return reduction


def autonomous_classical_map(kernel, dims=(2, 2), side=0):
    """Return L iff R K = L R on the full simplex (strong lumpability here)."""
    reduction = classical_marginal(dims, side)
    output = reduction @ kernel
    reduced = np.empty((dims[side], dims[side]), dtype=output.dtype)
    for label in range(dims[side]):
        columns = np.flatnonzero(reduction[label])
        reduced[:, label] = output[:, columns[0]]
        if not np.allclose(output[:, columns], reduced[:, label, None], atol=0, rtol=0):
            return None
    return reduced


def quantum_heisenberg_residual(unitary, side=0):
    """Checks complete qubit Hermitian basis, including Y for complex theory."""
    basis = (np.eye(2), np.array([[0., 1.], [1., 0.]]),
             np.diag([1., -1.]), np.array([[0., -1j], [1j, 0.]]))
    largest = 0.
    for a in basis:
        local = np.kron(a, np.eye(2)) if side == 0 else np.kron(np.eye(2), a)
        pulled = unitary.conj().T @ local @ unitary
        marginal = trace_keep(pulled, (2, 2), (side,)) / 2
        candidate = np.kron(marginal, np.eye(2)) if side == 0 else np.kron(np.eye(2), marginal)
        largest = max(largest, float(np.linalg.norm(pulled - candidate)))
    return largest


class MarginalDynamicsAuditTests(unittest.TestCase):
    def test_all_two_bit_permutations_have_only_four_two_sided_autonomous_cases(self):
        counts = [0, 0]
        for destinations in permutations(range(4)):
            matrix = np.eye(4, dtype=int)[:, destinations]
            first = autonomous_classical_map(matrix, side=0)
            second = autonomous_classical_map(matrix, side=1)
            counts[0] += first is not None
            counts[1] += first is not None and second is not None
        self.assertEqual(counts, [8, 4])

    def test_classical_reduction_factorization_and_a_hidden_direction_witness(self):
        cnot = np.eye(4, dtype=int)[[0, 1, 3, 2]]
        reduction = classical_marginal((2, 2), 0)
        local = autonomous_classical_map(cnot, side=0)
        np.testing.assert_array_equal(reduction @ cnot, local @ reduction)
        self.assertIsNone(autonomous_classical_map(cnot, side=1))
        hidden = np.array([1, -1, -1, 1])
        np.testing.assert_array_equal(reduction @ hidden, [0, 0])
        target_reduction = classical_marginal((2, 2), 1)
        self.assertGreater(np.linalg.norm(target_reduction @ cnot @ hidden), 0)

    def test_identical_both_marginals_give_opposite_future_outputs(self):
        cnot = np.eye(4, dtype=int)[[0, 1, 3, 2]]
        first = np.diag([.5, 0., 0., .5])
        second = np.diag([0., .5, .5, 0.])
        for side in (0, 1):
            np.testing.assert_array_equal(trace_keep(first, (2, 2), (side,)), np.eye(2) / 2)
            np.testing.assert_array_equal(trace_keep(second, (2, 2), (side,)), np.eye(2) / 2)
        out_a = trace_keep(cnot @ first @ cnot.T, (2, 2), (1,))
        out_b = trace_keep(cnot @ second @ cnot.T, (2, 2), (1,))
        np.testing.assert_array_equal(out_a, np.diag([1, 0]))
        np.testing.assert_array_equal(out_b, np.diag([0, 1]))
        self.assertEqual(np.linalg.norm(out_a - out_b, ord='nuc') / 2, 1)
        midpoint = (out_a + out_b) / 2
        self.assertEqual(np.linalg.norm(out_a - midpoint, ord='nuc') / 2, .5)

    def test_product_unitaries_are_autonomous_but_cnot_and_swap_are_not(self):
        rng = np.random.default_rng(176)
        for _ in range(5):
            a, _ = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))
            b, _ = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))
            for side in (0, 1):
                self.assertLess(quantum_heisenberg_residual(np.kron(a, b), side), 3e-15)
        for matrix in (np.eye(4)[[0, 1, 3, 2]], np.eye(4)[[0, 2, 1, 3]]):
            for side in (0, 1):
                self.assertGreater(quantum_heisenberg_residual(matrix, side), 1)

    def test_control_qubit_coherence_is_affected_by_target_even_for_real_inputs(self):
        cnot = np.eye(4)[[0, 1, 3, 2]]
        plus = np.ones((2, 2)) / 2
        outputs = []
        for target in (np.diag([1., 0.]), plus):
            out = cnot @ np.kron(plus, target) @ cnot.T
            outputs.append(trace_keep(out, (2, 2), (0,)))
        np.testing.assert_array_equal(outputs[0], np.eye(2) / 2)
        np.testing.assert_array_equal(outputs[1], plus)

    def test_fixed_independent_environment_does_allow_a_reduced_channel(self):
        cnot = np.eye(4)[[0, 1, 3, 2]]
        environment = np.diag([.3, .7])
        apply = lambda rho: trace_keep(cnot @ np.kron(rho, environment) @ cnot.T, (2, 2), (0,))
        first, second = np.ones((2, 2)) / 2, np.diag([.8, .2])
        np.testing.assert_allclose(apply(.4 * first + .6 * second), .4 * apply(first) + .6 * apply(second), atol=2e-16)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(MarginalDynamicsAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 176,
        "classical_condition": "ker(R) subset ker(RK), equivalently RK=LR",
        "two_bit_permutations": 24,
        "one_side_classically_autonomous": 8,
        "both_sides_classically_autonomous": 4,
        "same_marginals_distinct_relations_output_trace_distance": "1",
        "two_preparation_state_only_prediction_minimax_error": "1/2",
        "all_input_marginal_autonomy_for_finite_closed_complex_unitary_implies_product": True,
        "fixed_independent_environment_may_have_reduced_channel": True,
        "marginal_autonomy_is_derived_from_cut_consistency": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("marginal_dynamics_audit_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
