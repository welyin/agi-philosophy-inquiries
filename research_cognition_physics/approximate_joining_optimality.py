"""Round 128: exact finite dual certificate for deterministic qubit joining.

36 product preparations bound the average target-support acceptance by 5/6
for every real CPTP map from independent encodings to the common output.
Consequently its worst trace-distance error is at least 1/6. Round 127 attains
that lower bound for all unknown input states.
"""

import argparse
import json
import unittest
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

from approximate_subject_joining import deterministic_joining_kraus, join_output, trace_distance
from common_orientation_structure import encode_state
from complex_control_from_reference import real_lift
from encoded_composition_audit import independent_encoding
from quantum_interface_audit import PAULI_X, PAULI_Y, PAULI_Z
from universal_joining_bound import integer_kernel_basis


def six_probes():
    return tuple((np.eye(2) + sign * axis) / 2
                 for axis in (PAULI_X, PAULI_Y, PAULI_Z) for sign in (-1, 1))


@lru_cache(maxsize=1)
def integer_dual():
    # C is integer, C.T C = 2 P_aligned.
    c = integer_kernel_basis()[0]
    y_scaled = 24 * np.eye(16, dtype=np.int64) + 6 * c.T @ c
    # All small dyadic arithmetic is exact in binary64; no division by 36 here.
    average_scaled = 16 * sum(
        np.kron(independent_encoding((a, b)).T, real_lift(np.kron(a, b)))
        for a in six_probes() for b in six_probes())
    if not np.array_equal(average_scaled, np.rint(average_scaled)):
        raise ArithmeticError("The score matrix did not have the expected exact integer scale.")
    r_scaled = average_scaled.astype(np.int64)
    slack = np.kron(y_scaled, np.eye(8, dtype=np.int64)) - r_scaled
    return y_scaled, r_scaled, slack


@lru_cache(maxsize=1)
def dual_certificate():
    y, _, slack = integer_dual()
    if not np.array_equal(slack, slack.T):
        raise ArithmeticError("Dual slack must be symmetric.")
    eigenvalue_candidates = (0, 16, 20, 24, 32, 36)
    polynomial = np.eye(128, dtype=np.int64)
    arithmetic_bound = 1
    for root in eigenvalue_candidates:
        factor = slack - root * np.eye(128, dtype=np.int64)
        # The product of absolute row-sum bounds also bounds all partial sums
        # in each matrix multiplication, preventing int64 overflow.
        arithmetic_bound *= int(np.abs(factor).sum(axis=1).max())
        if arithmetic_bound >= 2 ** 63:
            raise ArithmeticError("Insufficient integer range for this certificate.")
        polynomial = polynomial @ factor
    if np.any(polynomial):
        raise ArithmeticError("The exact nonnegative spectral polynomial failed.")
    score = Fraction(int(np.trace(y)), 576)
    return {
        "product_preparations": 36,
        "input_dimension": 16,
        "output_dimension": 8,
        "choi_dimension": 128,
        "common_integer_scale": 576,
        "dual_input_trace_numerator": int(np.trace(y)),
        "average_support_acceptance_upper_bound_exact": str(score),
        "worst_trace_distance_lower_bound_exact": str(1 - score),
        "symmetric_integer_slack": True,
        "nonnegative_polynomial_roots_for_scaled_slack": list(eigenvalue_candidates),
        "exact_matrix_polynomial_is_zero": True,
        "integer_partial_sum_absolute_bound": arithmetic_bound,
        "uses_float_eigensolver_or_optimizer_to_certify_psd": False,
    }


def choi_from_kraus(kraus):
    # Input then output order; vec(K) means K.T.ravel() in this convention.
    vectors = [k.T.ravel() for k in kraus]
    return sum(np.outer(vector, vector.conj()) for vector in vectors)


def average_support_score(kraus):
    result = 0.
    for first in six_probes():
        for second in six_probes():
            source = independent_encoding((first, second))
            target_effect = real_lift(np.kron(first, second))
            output = sum(k @ source @ k.T for k in kraus)
            result += np.trace(output @ target_effect).real / 36
    return result


class ApproximateJoiningOptimalityTests(unittest.TestCase):
    def test_exact_dual_certificate_proves_one_sixth_lower_bound(self):
        report = dual_certificate()
        self.assertEqual(report["average_support_acceptance_upper_bound_exact"], "5/6")
        self.assertEqual(report["worst_trace_distance_lower_bound_exact"], "1/6")
        self.assertEqual(report["integer_partial_sum_absolute_bound"], 2997485568)

    def test_score_matrix_agrees_with_direct_channel_probabilities(self):
        ks = deterministic_joining_kraus()
        _, r, _ = integer_dual()
        choi = choi_from_kraus(ks)
        self.assertAlmostEqual(np.trace(choi @ r).real / 576, average_support_score(ks))
        self.assertAlmostEqual(average_support_score(ks), 5 / 6)

    def test_actual_choi_is_trace_preserving_and_saturates_dual(self):
        choi = choi_from_kraus(deterministic_joining_kraus())
        partial = np.einsum("abcb->ac", choi.reshape(16, 8, 16, 8))
        np.testing.assert_allclose(partial, np.eye(16), atol=1e-15)
        slack = integer_dual()[2]
        np.testing.assert_allclose(slack @ choi, 0., atol=3e-14)

    def test_all_36_probes_attain_lower_bound_in_trace_distance(self):
        for first in six_probes():
            for second in six_probes():
                output = join_output(first, second)
                ideal = encode_state(np.kron(first, second))
                self.assertAlmostEqual(trace_distance(output, ideal), 1 / 6)
                support = real_lift(np.kron(first, second))
                self.assertAlmostEqual(np.trace(output @ support), 5 / 6)

    def test_random_unrestricted_real_channels_obey_certificate(self):
        rng = np.random.default_rng(128)
        for _ in range(4):
            # General isometry 16 -> 8 x 3, without imposed orientation symmetry.
            raw = rng.normal(size=(24, 16))
            q, _ = np.linalg.qr(raw)
            ks = tuple(q[i:i + 8] for i in range(0, 24, 8))
            score = average_support_score(ks)
            self.assertLessEqual(score, 5 / 6 + 1e-14)
            choi = choi_from_kraus(ks)
            y, r, slack = integer_dual()
            self.assertAlmostEqual(np.trace(choi @ slack) / 576,
                                   np.trace(y) / 576 - np.trace(choi @ r) / 576, places=13)

    def test_target_support_test_is_a_valid_effect_not_a_pure_real_density(self):
        for first, second in ((six_probes()[0], six_probes()[3]),
                              (six_probes()[2], six_probes()[5])):
            effect = real_lift(np.kron(first, second))
            np.testing.assert_array_equal(effect @ effect, effect)
            self.assertEqual(np.trace(effect), 2.)
            self.assertEqual(np.trace(encode_state(np.kron(first, second)) @ effect), 1.)

    def test_integer_dual_trace_and_orientation_weights(self):
        y, _, _ = integer_dual()
        c = integer_kernel_basis()[0]
        aligned = c.T @ c / 2
        expected_scaled = 36 * aligned + 24 * (np.eye(16) - aligned)
        np.testing.assert_array_equal(y, expected_scaled)
        self.assertEqual(int(np.trace(y)), 480)

    def test_six_state_second_moment_matches_analytic_identity(self):
        swap = np.eye(4).reshape(2, 2, 2, 2).transpose(0, 1, 3, 2).reshape(4, 4)
        moment = sum(np.kron(rho, rho) for rho in six_probes())
        np.testing.assert_array_equal(moment, np.eye(4) + swap)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(ApproximateJoiningOptimalityTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 128,
        "task": "Deterministic common encoding of two independent unknown standard qubit encodings",
        "optimization_class": "All real completely positive trace-preserving maps, including fixed independent auxiliaries",
        "optimal_worst_trace_distance_exact": "1/6",
        "optimal_average_target_support_acceptance_exact": "5/6",
        "finite_dual_certificate": dual_certificate(),
        "attainment": "Round 127, with the old host preserved exactly",
        "any_host_dimension_at_least_two_and_qubit_newcomer_same_optimum": True,
        "higher_newcomer_dimension_optimum_claimed": False,
        "retained_records_deleted_physically": False,
        "arbitrary_record_dependent_old_interfaces_equated_to_common_encoding": False,
        "compiled_into_old_noisy_primitives": False,
        "quantum_theory_derived_from_cognition": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("approximate_joining_optimality_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
