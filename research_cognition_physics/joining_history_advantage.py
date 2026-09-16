"""Round 136: a sharp advantage from reopening retained joining history.

First A,B use the balanced optimal channel. With its environment inaccessible,
every second-stage CPTP decoder has worst final error >= 11/36. Reopening the
coherent history permits 1/4 while preserving the same complete AB marginal.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import unittest

import numpy as np

from approximate_joining_optimality import six_probes
from approximate_subject_joining import apply_kraus, deterministic_joining_kraus, trace_distance
from common_orientation_structure import encode_state
from compiled_joining_frontier import compiled_sharing, sharing_circuit
from compiled_subject_joining import circuit_matrix
from complex_control_from_reference import real_lift
from encoded_composition_audit import independent_encoding, product_state
from flagged_subject_joining import random_state
from joining_disturbance_sharing import predicted_logical_output
from multisubject_joining_optimum import (
    local_spectral_certificate, moment_matrices, multisubject_kraus,
    predicted_logical, reorder_qubits, y_projector,
)
from quantum_interface_audit import PAULI_Y


def balanced_pair(first, second):
    return predicted_logical_output(first, second, F(1, 2))


def restricted_input(states):
    return independent_encoding((balanced_pair(*states[:2]), states[2]))


def restricted_output(states):
    return apply_kraus(deterministic_joining_kraus(4, 2), restricted_input(states))


def common_ab_marginal(encoded):
    return np.trace(encoded.reshape((2,)*8), axis1=3, axis2=7).reshape(8, 8)


def restricted_certificate():
    local_spectral_certificate()
    a, b = F(1, 2), F(1, 3)
    lam = F(2, 3)
    a_noisy, b_noisy = lam*a+(1-lam)/4, lam*b+(1-lam)/4
    ab_match, ab_mismatch = a*a_noisy, b*b_noisy
    same = max(ab_match*a, ab_mismatch*b)
    different = max(ab_match*b, ab_mismatch*a)
    score = 4*(same+different)
    return {
        "single_matched_noisy_moment_norm_exact": str(a_noisy),
        "single_mismatched_noisy_moment_norm_exact": str(b_noisy),
        "ab_matched_moment_norm_exact": str(ab_match),
        "ab_mismatched_moment_norm_exact": str(ab_mismatch),
        "aligned_unweighted_block_norm_exact": str(same),
        "opposed_unweighted_block_norm_exact": str(different),
        "all_second_stage_channels_support_bound_exact": str(score),
        "optimal_restricted_worst_error_exact": str(1-score),
        "old_ab_preservation_needed_to_prove_lower_bound": False,
        "proof_uses_float_optimizer": False,
    }


@lru_cache(maxsize=1)
def restricted_score_matrix():
    a, b = (m/6 for m in moment_matrices())
    an, bn = (2*a/3+np.eye(4)/12), (2*b/3+np.eye(4)/12)
    ca = (np.kron(a, an)+np.kron(an, a))/2
    cb = (np.kron(b, bn)+np.kron(bn, b))/2
    grouped = np.zeros((512, 512), dtype=complex)
    for first, second, out in product((0, 1), repeat=3):
        refs = product_state((y_projector(first).T, y_projector(second).T, y_projector(out)))
        targets = np.kron(ca if first == out else cb, a if second == out else b)
        grouped += np.kron(refs, targets)/4
    result = reorder_qubits(grouped, (0, 1, 3, 5, 7, 2, 4, 6, 8))
    if np.max(abs(result.imag)) > 1e-14:
        raise ArithmeticError("Expected a real score matrix.")
    return result.real


def restricted_dual():
    aligned = (np.eye(4)+np.kron(PAULI_Y, PAULI_Y).real)/2
    return np.kron(5*aligned/192+5*(np.eye(4)-aligned)/288, np.eye(8))


class JoiningHistoryAdvantageTests(unittest.TestCase):
    def test_exact_tensor_moment_certificate(self):
        cert = restricted_certificate()
        self.assertEqual(cert["ab_matched_moment_norm_exact"], "5/24")
        self.assertEqual(cert["ab_mismatched_moment_norm_exact"], "11/108")
        self.assertEqual(cert["all_second_stage_channels_support_bound_exact"], "25/36")
        self.assertEqual(cert["optimal_restricted_worst_error_exact"], "11/36")
        self.assertFalse(cert["old_ab_preservation_needed_to_prove_lower_bound"])

    def test_factorized_score_equals_actual_216_fresh_product_preparations(self):
        direct = np.zeros((512, 512))
        for states in product(six_probes(), repeat=3):
            direct += np.kron(restricted_input(states).T, real_lift(product_state(states)))/216
        np.testing.assert_allclose(restricted_score_matrix(), direct, atol=5e-17)

    def test_all_channel_dual_is_saturated_by_host_preserving_second_stage(self):
        r, y = restricted_score_matrix(), restricted_dual()
        self.assertAlmostEqual(np.trace(y), 25/36)
        slack = np.kron(y, np.eye(16))-r
        # Numerical audit of the full matrix; the exact proof is the local
        # projector certificate plus the tensor-moment inequalities above.
        self.assertGreaterEqual(np.linalg.eigvalsh(slack).min(), -2e-16)
        ks = deterministic_joining_kraus(4, 2)
        for k in ks:
            np.testing.assert_allclose(slack@k.T.ravel(), 0., atol=2e-16)
        score = sum(k.T.ravel()@r@k.T.ravel() for k in ks)
        self.assertAlmostEqual(score, 25/36)

    def test_all_pure_probe_errors_have_strict_history_advantage(self):
        for states in product(six_probes(), repeat=3):
            ideal = encode_state(product_state(states))
            restricted = restricted_output(states)
            whole_access = apply_kraus(multisubject_kraus(3), independent_encoding(states))
            self.assertAlmostEqual(trace_distance(restricted, ideal), 11/36)
            self.assertAlmostEqual(trace_distance(whole_access, ideal), 1/4)

    def test_entire_old_ab_joint_state_is_the_same_in_both_protocols(self):
        rng = np.random.default_rng(136)
        for _ in range(12):
            states = [random_state(rng, 2) for _ in range(3)]
            before = encode_state(balanced_pair(*states[:2]))
            joint = apply_kraus(multisubject_kraus(3), independent_encoding(states))
            for result in (joint, restricted_output(states)):
                np.testing.assert_allclose(common_ab_marginal(result), before, atol=5e-16)
            # Equality covers arbitrary later AB-only joint effects, not only locals.
            raw = rng.normal(size=(8, 8))
            effect = raw@raw.T
            effect /= np.linalg.eigvalsh(effect).max()
            self.assertAlmostEqual(np.trace(common_ab_marginal(joint)@effect), np.trace(before@effect))

    def test_coherent_old_history_reopens_without_unknown_copies_or_source_purification(self):
        rng = np.random.default_rng(236)
        states = [random_state(rng, 2) for _ in range(3)]
        source_ab = independent_encoding(states[:2])
        v, _ = compiled_sharing(F(1, 2))
        gates, systems = sharing_circuit(F(1, 2))
        old_whole = v@source_ab@v.T
        u = circuit_matrix(gates, systems)
        recovered = u.T@old_whole@u
        expected = np.kron(source_ab, np.diag([1.]+[0.]*7))
        np.testing.assert_allclose(recovered, expected, atol=1e-15)
        original_ab = recovered[::8, ::8]  # Unused ancillary |000> remains in the whole.
        with_new = reorder_qubits(np.kron(original_ab, encode_state(states[2])), (0,1,4,2,3,5))
        np.testing.assert_allclose(with_new, independent_encoding(states), atol=1e-15)
        np.testing.assert_allclose(apply_kraus(multisubject_kraus(3), with_new),
                                   encode_state(predicted_logical(states)), atol=1e-15)

    def test_errors_are_anticorrelated_in_the_three_party_joint_protocol(self):
        zero = np.diag([1., 0.])
        joint = predicted_logical((zero, zero, zero))
        expected = np.zeros(8)
        expected[0] = .75
        expected[[1, 2, 4]] = 1/12
        np.testing.assert_allclose(joint, np.diag(expected), atol=2e-16)
        old = balanced_pair(zero, zero)
        from approximate_subject_joining import logical_newcomer
        sequential = np.kron(old, logical_newcomer(zero))
        self.assertAlmostEqual(sequential[0, 0], 25/36)
        self.assertGreater(sequential[3, 3], 0.)
        self.assertAlmostEqual(joint[3, 3], 0.)

    def test_abstract_real_dilation_and_exact_gain_accounting(self):
        ks = multisubject_kraus(3)
        gram = np.array([[np.sum(a*b) for b in ks] for a in ks])
        self.assertEqual(np.linalg.matrix_rank(gram), 10)
        padded = np.zeros((256, 64))
        padded[:160] = np.vstack(ks)
        np.testing.assert_allclose(padded.T@padded, np.eye(64), atol=2e-15)
        self.assertEqual(F(11,36)-F(1,4), F(1,18))
        self.assertEqual(F(1,18)/F(11,36), F(2,11))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(JoiningHistoryAdvantageTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"round": 136,
        "first_stage": "Fixed optimal balanced A,B joining, with coherent history retained",
        "restricted_second_stage": "Any CPTP map accessing only common R,A,B and the fresh encoded C",
        "restricted_certificate": restricted_certificate(),
        "whole_history_optimal_error_exact": "1/4",
        "error_improvement_exact": "1/18",
        "relative_error_reduction_exact": "2/11",
        "complete_old_ab_marginal_preserved_on_promised_inputs": True,
        "arbitrary_old_external_correlations_preserved_at_active_interface_claimed": False,
        "unobserved_environment_is_deleted": False,
        "old_record_measurements_or_unmodelled_intervening_actions_allowed": False,
        "unknown_copies_or_source_purifications_required": False,
        "three_party_kraus_rank": 10,
        "abstract_additional_pure_rebits_from_fresh_six_rebit_input": 2,
        "three_party_native_gate_compilation_completed_in_this_round": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("joining_history_advantage_results.json").write_text(
            json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
