"""Round 125: lossless flagged intake of independently encoded subjects.

Ideal real orthogonal control and reference-sector readout are declared here.
The retained flag preserves old interfaces; it is not unconditional access to
the full common complex interface. No old noisy-gate compilation is claimed.
"""

import argparse
import json
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from common_orientation_structure import encode_state
from complex_control_from_reference import real_lift
from encoded_composition_audit import (
    aligned_isometry, independent_encoding, product_state,
)
from quantum_interface_audit import PAULI_X, PAULI_Y


def orientation_flags(count):
    if count < 1:
        raise ValueError("At least one independently encoded subject is needed.")
    return tuple((0,) + bits for bits in product((0, 1), repeat=count - 1))


def joining_kraus(dimensions):
    """Input: all private references, then targets. Output: one reference, targets."""
    if not dimensions or any(d < 1 for d in dimensions):
        raise ValueError("Positive target dimensions are required.")
    base = aligned_isometry(len(dimensions))
    target_identity = np.eye(int(np.prod(dimensions)))
    result = {}
    for flag in orientation_flags(len(dimensions)):
        flip = product_state([PAULI_X.real if bit else np.eye(2) for bit in flag]).real
        result[flag] = np.kron((flip @ base).T, target_identity)
    return result


def intake_matrix(dimensions):
    """Orthogonal change of encoding; output order is flag, reference, targets."""
    return np.vstack(tuple(joining_kraus(dimensions).values()))


def flag_blocks(matrix, block_size):
    """An actual unread/read classical flag: keep its diagonal blocks."""
    result = np.zeros_like(matrix)
    for start in range(0, len(matrix), block_size):
        result[start:start + block_size, start:start + block_size] = (
            matrix[start:start + block_size, start:start + block_size]
        )
    return result


def branch_target(states, flag):
    return product_state([rho.conj() if bit else rho for rho, bit in zip(states, flag)])


def partial_transpose(matrix, dimensions, positions):
    count = len(dimensions)
    order = list(range(2 * count))
    for pos in positions:
        order[pos], order[pos + count] = order[pos + count], order[pos]
    return matrix.reshape(tuple(dimensions) * 2).transpose(order).reshape(matrix.shape)


def random_state(rng, dimension):
    raw = rng.normal(size=(dimension, dimension)) + 1j * rng.normal(size=(dimension, dimension))
    density = raw @ raw.conj().T
    return density / np.trace(density).real


def old_history_operator(a, b):
    """Two sequential local instruments, with Bob's setting depending on Alice."""
    phase = np.diag([1, np.exp(.37j)])
    alice = (np.diag([np.sqrt(.8), np.sqrt(.3)]),
             np.diag([np.sqrt(.2), np.sqrt(.7)]))[a] @ phase
    theta = .31 + .23 * a
    unitary = np.array([[np.cos(theta), -1j * np.sin(theta)],
                       [-1j * np.sin(theta), np.cos(theta)]])
    bob = (np.diag([np.sqrt(.6), np.sqrt(.1)]),
           np.diag([np.sqrt(.4), np.sqrt(.9)]))[b] @ unitary
    return alice, bob


class FlaggedSubjectJoiningTests(unittest.TestCase):
    def test_intake_is_orthogonal_without_extra_quantum_capacity(self):
        for dims in ((2, 2), (2, 3), (2, 2, 2)):
            w = intake_matrix(dims)
            expected = 2 ** len(dims) * int(np.prod(dims))
            self.assertEqual(w.shape, (expected, expected))
            np.testing.assert_allclose(w.T @ w, np.eye(expected), atol=1e-15)
            np.testing.assert_allclose(w @ w.T, np.eye(expected), atol=1e-15)

    def test_all_flags_have_exact_branch_formula_on_unknown_mixed_states(self):
        rng = np.random.default_rng(125)
        for dims in ((2, 2), (2, 3), (3, 2), (2, 2, 2)):
            states = [random_state(rng, d) for d in dims]
            source = independent_encoding(states)
            weight = float(Fraction(1, 2 ** (len(dims) - 1)))
            for flag, k in joining_kraus(dims).items():
                branch = k @ source @ k.T
                self.assertAlmostEqual(np.trace(branch), weight)
                np.testing.assert_allclose(
                    branch, weight * encode_state(branch_target(states, flag)), atol=2e-16,
                )

    def test_reading_flag_preserves_entire_promised_density_under_inverse(self):
        rng = np.random.default_rng(225)
        for dims in ((2, 3), (2, 2, 2)):
            source = independent_encoding([random_state(rng, d) for d in dims])
            w = intake_matrix(dims)
            transformed = w @ source @ w.T
            recorded = flag_blocks(transformed, 2 * int(np.prod(dims)))
            np.testing.assert_allclose(recorded, transformed, atol=2e-17)
            np.testing.assert_allclose(w.T @ recorded @ w, source, atol=2e-16)

    def test_coherent_intake_is_reversible_but_classical_readout_needs_input_promise(self):
        rng = np.random.default_rng(325)
        vector = rng.normal(size=16)
        vector /= np.linalg.norm(vector)
        source = np.outer(vector, vector)
        w = intake_matrix((2, 2))
        transformed = w @ source @ w.T
        np.testing.assert_allclose(w.T @ transformed @ w, source, atol=5e-16)
        recovered = w.T @ flag_blocks(transformed, 8) @ w
        self.assertGreater(np.linalg.norm(recovered - source), .1)

    def test_adaptive_old_local_instruments_preserve_all_records_and_conditional_states(self):
        rng = np.random.default_rng(425)
        states = [random_state(rng, 2), random_state(rng, 2)]
        source = independent_encoding(states)
        total = 0.
        for a, b in product((0, 1), repeat=2):
            operators = old_history_operator(a, b)
            updated = [k @ rho @ k.conj().T for k, rho in zip(operators, states)]
            expected_probability = np.trace(product_state(updated)).real
            total += expected_probability
            actual_probability = 0.
            for flag, k in joining_kraus((2, 2)).items():
                old = k @ source @ k.T
                controlled = real_lift(product_state([
                    op.conj() if bit else op for op, bit in zip(operators, flag)
                ]))
                actual = controlled @ old @ controlled.T
                np.testing.assert_allclose(
                    actual, encode_state(branch_target(updated, flag)) / 2, atol=1e-16,
                )
                actual_probability += np.trace(actual)
            self.assertAlmostEqual(actual_probability, expected_probability)
        self.assertAlmostEqual(total, 1.)

    def test_flagged_ppt_joint_effect_works_but_bell_effect_does_not(self):
        rng = np.random.default_rng(525)
        states = [random_state(rng, 2), random_state(rng, 2)]
        source = independent_encoding(states)
        effect = (np.eye(4) + .6 * np.kron(PAULI_Y, PAULI_Y)) / 2
        actual = 0.
        for flag, k in joining_kraus((2, 2)).items():
            adjusted = partial_transpose(effect, (2, 2), [i for i, bit in enumerate(flag) if bit])
            self.assertGreaterEqual(np.linalg.eigvalsh(adjusted).min(), 0.)
            self.assertLessEqual(np.linalg.eigvalsh(adjusted).max(), 1.)
            actual += np.trace((k @ source @ k.T) @ real_lift(adjusted))
        self.assertAlmostEqual(actual, np.trace(product_state(states) @ effect).real)
        bell = np.outer([1., 0., 0., 1.], [1., 0., 0., 1.]) / 2
        self.assertAlmostEqual(np.linalg.eigvalsh(partial_transpose(bell, (2, 2), (1,))).min(), -.5)

    def test_forgetting_flag_loses_an_old_y_distinction(self):
        outputs = []
        for sign in (-1, 1):
            source = independent_encoding((np.eye(2) / 2, (np.eye(2) + sign * PAULI_Y) / 2))
            outputs.append(sum(k @ source @ k.T for k in joining_kraus((2, 2)).values()))
        np.testing.assert_allclose(*outputs, atol=1e-16)

    def test_rollback_after_failure_does_not_supply_fresh_independent_retry(self):
        rng = np.random.default_rng(625)
        source = independent_encoding([random_state(rng, 2), random_state(rng, 2)])
        ks = joining_kraus((2, 2))
        failure = ks[(0, 1)]
        branch = failure @ source @ failure.T
        recovered_conditional = failure.T @ branch @ failure / np.trace(branch)
        self.assertAlmostEqual(np.trace(recovered_conditional), 1.)
        self.assertAlmostEqual(np.trace(ks[(0, 0)] @ recovered_conditional @ ks[(0, 0)].T), 0.)
        self.assertAlmostEqual(np.trace(failure @ recovered_conditional @ failure.T), 1.)
        self.assertGreater(np.linalg.norm(recovered_conditional - source), .1)

    def test_known_real_new_input_is_a_deterministic_exception(self):
        rng = np.random.default_rng(725)
        host = random_state(rng, 3)
        newcomer = np.array([[.6, .2], [.2, .4]])
        source = independent_encoding((host, newcomer))
        output = sum(k @ source @ k.T for k in joining_kraus((3, 2)).values())
        np.testing.assert_allclose(output, encode_state(np.kron(host, newcomer)), atol=2e-16)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(FlaggedSubjectJoiningTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 125,
        "input": "One independent standard real encoding E(rho_i) per subject",
        "intake": "A real orthogonal change of encoding, with all orientation flags retained",
        "n_subjects_flag_bits": "n-1",
        "flag_probability_exact": "2^(1-n)",
        "branch_formula": "2^(1-n) E(tensor_i conjugate^flag_i(rho_i))",
        "all_old_local_complex_instruments_and_classical_feedback_preserved": True,
        "unconditional_inverse_recovers_full_promised_input_density": True,
        "arbitrary_external_purifications_preserved_after_classical_flag_readout": False,
        "coherent_flag_version_is_invertible_on_arbitrary_real_input": True,
        "same_unknown_state_after_observed_failure_is_a_fresh_retry": False,
        "retry_this_filter_success_after_failure_exact": "0",
        "simple_reinterpretation_supplies_all_common_effects": False,
        "bell_partial_transpose_smallest_eigenvalue_exact": "-1/2",
        "conjugation_invariant_new_input_admits_deterministic_join": True,
        "ideal_reference_sector_readout": True,
        "compiled_into_old_noisy_primitives": False,
        "physical_gate_time_and_distributed_communication_costs_certified": False,
        "all_records_and_unused_systems_retained_in_total_accounting": True,
        "quantum_theory_derived_from_cognition": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("flagged_subject_joining_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
