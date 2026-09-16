"""Round 129: exact two-qubit frontier for sharing deterministic joining error.

For marginal worst trace errors delta_A and delta_B in the common interface,
delta_A + delta_B >= 1/6. A real CPTP family attains the entire lower edge.
This concerns the fixed common output, not arbitrary flag-dependent old access.
"""

import argparse
import json
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from approximate_joining_optimality import dual_certificate, six_probes
from approximate_subject_joining import (
    apply_kraus, encoded_host_marginal, logical_newcomer, trace_distance,
    transpose_kraus,
)
from common_orientation_structure import decode_state, encode_state, normalizer_flip
from complex_control_from_reference import real_lift
from encoded_composition_audit import independent_encoding
from flagged_subject_joining import joining_kraus, partial_transpose, random_state
from quantum_interface_audit import PAULI_X, PAULI_Y, PAULI_Z


def sharing_kraus(newcomer_share):
    theta = float(newcomer_share)
    if not np.isfinite(theta) or not 0 <= theta <= 1:
        raise ValueError("The newcomer's share must lie between zero and one.")
    intake = joining_kraus((2, 2))
    result = [intake[(0, 0)]]
    for k in transpose_kraus(2):
        if theta:
            result.append(np.sqrt(theta) * np.kron(np.eye(4), k) @ intake[(0, 1)])
        if theta < 1:
            repair_host = np.kron(np.eye(2), np.kron(k, np.eye(2))) @ normalizer_flip(4)
            result.append(np.sqrt(1 - theta) * repair_host @ intake[(0, 1)])
    return tuple(result)


def sharing_output(host, newcomer, newcomer_share):
    return apply_kraus(sharing_kraus(newcomer_share), independent_encoding((host, newcomer)))


def predicted_logical_output(host, newcomer, newcomer_share):
    theta = float(newcomer_share)
    return (theta * np.kron(host, logical_newcomer(newcomer)) +
            (1 - theta) * np.kron(logical_newcomer(host), newcomer))


def newcomer_marginal(encoded):
    return np.einsum("abcdbf->acdf", encoded.reshape(2, 2, 2, 2, 2, 2)).reshape(4, 4)


def sharing_row(newcomer_share):
    theta = Fraction(newcomer_share)
    if not 0 <= theta <= 1:
        raise ValueError("Invalid share.")
    return {
        "newcomer_share_exact": str(theta),
        "host_worst_marginal_error_exact": str((1 - theta) / 6),
        "newcomer_worst_marginal_error_exact": str(theta / 6),
        "host_distinguishability_factor_exact": str(1 - (1 - theta) / 3),
        "newcomer_distinguishability_factor_exact": str(1 - theta / 3),
        "joint_worst_trace_distance_exact": "1/6",
    }


def bell_task(newcomer_share):
    host = (np.eye(2) + PAULI_X) / 2
    new = np.diag([1., 0.])
    cnot = np.array([[1, 0, 0, 0], [0, 1, 0, 0],
                     [0, 0, 0, 1], [0, 0, 1, 0]], dtype=float)
    operation = real_lift(cnot)
    encoded = operation @ sharing_output(host, new, newcomer_share) @ operation.T
    logical = decode_state(encoded)
    bell = np.outer([1., 0., 0., 1.], [1., 0., 0., 1.]) / 2
    return logical, float(np.trace(logical @ bell).real), float(
        np.linalg.eigvalsh(partial_transpose(logical, (2, 2), (1,))).min())


class JoiningDisturbanceSharingTests(unittest.TestCase):
    def test_complete_real_family_preserves_trace_for_all_tested_shares(self):
        for theta in (0, Fraction(1, 4), Fraction(1, 2), Fraction(3, 4), 1):
            ks = sharing_kraus(theta)
            np.testing.assert_allclose(sum(k.T @ k for k in ks), np.eye(16), atol=1e-15)
            self.assertTrue(all(np.isrealobj(k) for k in ks))

    def test_full_mixed_state_and_both_marginal_formulas(self):
        rng = np.random.default_rng(129)
        for theta in (0., .3, .5, 1.):
            for _ in range(3):
                host, new = random_state(rng, 2), random_state(rng, 2)
                actual = sharing_output(host, new, theta)
                predicted = predicted_logical_output(host, new, theta)
                np.testing.assert_allclose(actual, encode_state(predicted), atol=3e-16)
                host_state = theta * host + (1 - theta) * logical_newcomer(host)
                new_state = (1 - theta) * new + theta * logical_newcomer(new)
                np.testing.assert_allclose(encoded_host_marginal(actual, 2, 2), encode_state(host_state), atol=3e-16)
                np.testing.assert_allclose(newcomer_marginal(actual), encode_state(new_state), atol=3e-16)

    def test_pure_inputs_attain_entire_sharp_marginal_frontier_and_joint_bound(self):
        for theta in (Fraction(0), Fraction(1, 2), Fraction(1)):
            for host, new in zip(six_probes(), reversed(six_probes())):
                output = sharing_output(host, new, theta)
                self.assertAlmostEqual(trace_distance(encoded_host_marginal(output, 2, 2), encode_state(host)),
                                       float((1 - theta) / 6))
                self.assertAlmostEqual(trace_distance(newcomer_marginal(output), encode_state(new)),
                                       float(theta / 6))
                self.assertAlmostEqual(trace_distance(output, encode_state(np.kron(host, new))), 1 / 6)

    def test_commuting_support_union_bound_underlying_marginal_lower_bound(self):
        for first in six_probes():
            for second in six_probes():
                a = real_lift(np.kron(first, np.eye(2)))
                b = real_lift(np.kron(np.eye(2), second))
                joint = real_lift(np.kron(first, second))
                np.testing.assert_array_equal(a @ b, joint)
                slack = np.eye(8) - a - b + joint
                np.testing.assert_array_equal(slack @ slack, slack)
                self.assertGreaterEqual(np.linalg.eigvalsh(slack).min(), -1e-15)
        self.assertEqual(dual_certificate()["average_support_acceptance_upper_bound_exact"], "5/6")

    def test_balanced_interface_retains_five_sixths_of_each_local_distinction(self):
        fixed = np.array([[.7, .1j], [-.1j, .3]])
        for axis in (PAULI_X, PAULI_Y, PAULI_Z):
            minus, plus = [(np.eye(2) + sign * axis) / 2 for sign in (-1, 1)]
            first, second = [newcomer_marginal(sharing_output(fixed, new, .5)) for new in (minus, plus)]
            self.assertAlmostEqual(trace_distance(first, second), 5 / 6)
            first, second = [encoded_host_marginal(sharing_output(host, fixed, .5), 2, 2)
                             for host in (minus, plus)]
            self.assertAlmostEqual(trace_distance(first, second), 5 / 6)

    def test_common_entangling_task_remains_available_with_quantified_loss(self):
        for theta in (0., .25, .5, .75, 1.):
            logical, acceptance, minimum_pt = bell_task(theta)
            self.assertAlmostEqual(acceptance, 5 / 6)
            self.assertAlmostEqual(minimum_pt, -1 / 3)
            self.assertAlmostEqual(np.trace(logical), 1.)

    def test_whole_state_can_keep_all_environment_labels_without_postselection(self):
        ks = sharing_kraus(.5)
        self.assertEqual(len(ks), 7)
        isometry = np.vstack(ks)
        np.testing.assert_allclose(isometry.T @ isometry, np.eye(16), atol=1e-15)
        source = independent_encoding((six_probes()[1], six_probes()[2]))
        whole = isometry @ source @ isometry.T
        np.testing.assert_allclose(isometry.T @ whole @ isometry, source, atol=3e-16)

    def test_exact_shares_and_invalid_inputs(self):
        row = sharing_row(Fraction(1, 2))
        self.assertEqual(row["host_worst_marginal_error_exact"], "1/12")
        self.assertEqual(row["newcomer_worst_marginal_error_exact"], "1/12")
        for theta in (-.1, 1.1, float("nan")):
            with self.assertRaises(ValueError):
                sharing_kraus(theta)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(JoiningDisturbanceSharingTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 129,
        "scope": "Two independent standard unknown qubit encodings, deterministic common output",
        "universal_marginal_error_lower_bound": "delta_A + delta_B >= 1/6",
        "entire_lower_frontier_attained": True,
        "optimal_equal_marginal_error_exact": "1/12",
        "joint_worst_trace_distance_for_all_frontier_channels_exact": "1/6",
        "rows": [sharing_row(theta) for theta in (Fraction(0), Fraction(1, 4), Fraction(1, 2),
                                                Fraction(3, 4), Fraction(1))],
        "bell_task_acceptance_exact": "5/6",
        "bell_task_partial_transpose_minimum_exact": "-1/3",
        "balanced_channel_kraus_count": 7,
        "balanced_dilation_environment_dimension_sufficient": 7,
        "balanced_additional_pure_rebits_sufficient_for_binary_orthogonal_extension": 2,
        "all_records_and_coherent_environment_retained": True,
        "individual_conditional_fine_records_have_same_error_guarantee": False,
        "claimed_advantage_over_all_original_real_joint_protocols": False,
        "compiled_into_old_noisy_primitives": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("joining_disturbance_sharing_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
