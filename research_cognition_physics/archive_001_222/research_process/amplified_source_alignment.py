"""Round 64: amplify reference alignment using one resettable noisy pointer.

The middle YX gate is applied once. Its Z pointer is copied to fresh/reset
helpers; directly rereading the original destructive instrument is different.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from bilocal_record_tomography import embed_operator
from independent_source_alignment import after_middle_gate, conditional_reference, keep_systems, real_x_with_reset
from operational_effect_closure import compiled_copy_gate, majority_certificate, majority_error
from quantum_interface_audit import ALPHA, apply_branch, branch_kraus
from shared_real_reference_tomography import shared_reference


def pointer_branch(density, wire, systems, outcome, contrast=1.):
    joint = np.kron(np.diag([1., 0]), density)
    copy = embed_operator(compiled_copy_gate(), (0, wire+1), systems+1)
    joint = copy @ joint @ copy.conj().T
    output = np.zeros_like(joint)
    for k in branch_kraus(0., outcome, contrast):
        extended = embed_operator(k, (0,), systems+1)
        output += extended @ joint @ extended.conj().T
    return keep_systems(output, tuple(range(1, systems+1)), systems+1)


def alignment_record_states(count, contrast=1.):
    majority_error(count, .1)  # Check positive odd count.
    states = {(): after_middle_gate()}
    for _ in range(count):
        states = {record+(2*outcome-1,): pointer_branch(rho, 2, 4, outcome, contrast)
                  for record, rho in states.items() for outcome in (0, 1)}
    return {record: keep_systems(rho, (0, 1, 3), 4) for record, rho in states.items()}


def amplified_reference_by_circuit(count, contrast=1.):
    states = alignment_record_states(count, contrast)
    return sum((rho if sum(record) > 0 else real_x_with_reset(rho, 2, 3)) for record, rho in states.items())


def direct_read_records(z_sign, count):
    states = {(): np.diag([(1+z_sign)/2, (1-z_sign)/2]).astype(complex)}
    for _ in range(count):
        states = {record+(2*outcome-1,): apply_branch(rho, 0., outcome, 1.)
                  for record, rho in states.items() for outcome in (0, 1)}
    return {record: float(np.trace(rho).real) for record, rho in states.items()}


def reference_certificate(count):
    certificate = majority_certificate(count)
    upper = Fraction(certificate["error_upper_exact"])
    return {"pointer_reads": count, "reference_trace_distance_upper_exact": str(upper),
            "reference_trace_distance_diagnostic": certificate["error_diagnostic"],
            "strictly_positive_finite_error": certificate["strictly_positive_error"],
            "two_independent_pair_sources_per_run": 2,
            "middle_initial_alignment_flows": 1, "middle_copy_flows": count,
            "middle_reusable_extra_pointers": 1,
            "maximum_remote_correction_flows": 1, "maximum_feedforward_bits": 1,
            "postselection_success_probability": 1,
            "below_1e_minus_9": upper < Fraction(1, 10**9),
            "below_1e_minus_12": upper < Fraction(1, 10**12)}


class AmplifiedSourceAlignmentTests(unittest.TestCase):
    def test_pointer_records_follow_the_full_conditional_binary_channel_likelihood(self):
        states = alignment_record_states(3)
        for record, state in states.items():
            likelihoods = [math.prod((1+sign*ALPHA*outcome)/2 for outcome in record) for sign in (1, -1)]
            expected = sum(likelihood*conditional_reference(sign, 1.)/2 for sign, likelihood in zip((1, -1), likelihoods))
            np.testing.assert_allclose(state, expected, atol=4e-16)

    def test_all_record_probabilities_sum_to_one_and_no_trials_are_discarded(self):
        states = alignment_record_states(3, .7)
        self.assertAlmostEqual(sum(np.trace(rho).real for rho in states.values()), 1., places=14)
        for flag in (-1, 1):
            mass = sum(np.trace(rho).real for record, rho in states.items() if (1 if sum(record) > 0 else -1) == flag)
            self.assertAlmostEqual(mass, .5, places=14)

    def test_actual_reset_pointer_circuit_has_the_predicted_amplified_reference(self):
        for count, eta in ((1, .7), (3, .7), (3, 1.)):
            error = majority_error(count, (1-ALPHA*eta)/2)
            np.testing.assert_allclose(amplified_reference_by_circuit(count, eta), conditional_reference(1, 1-2*error), atol=6e-16)

    def test_pointer_copy_does_not_destroy_the_original_z_value(self):
        for sign, outcome in product((-1, 1), (0, 1)):
            rho = np.diag([(1+sign)/2, (1-sign)/2]).astype(complex)
            updated = pointer_branch(rho, 0, 1, outcome)
            np.testing.assert_allclose(updated, (1+sign*ALPHA*(2*outcome-1))*rho/2, atol=3e-16)

    def test_direct_repeated_old_reading_retains_only_the_first_read_information(self):
        for count in (1, 3, 5):
            plus, minus = direct_read_records(1, count), direct_read_records(-1, count)
            distance = sum(abs(plus[key]-minus[key]) for key in plus)/2
            self.assertAlmostEqual(distance, ALPHA, places=14)

    def test_majority_alignment_is_optimal_for_the_specified_binary_pointer_records(self):
        count = 5
        error = 0.
        for record in product((-1, 1), repeat=count):
            plus = math.prod((1+ALPHA*s)/2 for s in record)
            minus = math.prod((1-ALPHA*s)/2 for s in record)
            error += min(plus, minus)/2
        self.assertAlmostEqual(error, majority_error(count, (1-ALPHA)/2), places=18)

    def test_alignment_distance_equals_majority_error_without_omitted_failure_mass(self):
        for count in (1, 3):
            actual = amplified_reference_by_circuit(count)
            error = np.abs(np.linalg.eigvalsh(actual-shared_reference(3))).sum()/2
            self.assertAlmostEqual(error, majority_error(count, (1-ALPHA)/2), places=14)

    def test_nine_and_fifteen_read_certificates_are_strict_and_positive(self):
        for count, threshold in ((9, Fraction(4673, 10**13)), (15, Fraction(33, 10**16))):
            certificate = reference_certificate(count)
            upper = Fraction(certificate["reference_trace_distance_upper_exact"])
            self.assertGreater(upper, 0)
            self.assertLess(upper, threshold)
            self.assertTrue(certificate["strictly_positive_finite_error"])

    def test_even_or_empty_pointer_counts_are_rejected(self):
        for count in (0, 2):
            with self.assertRaises(ValueError):
                alignment_record_states(count)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(AmplifiedSourceAlignmentTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"round": 64, "scope": "Two perfect independent real pair references; old middle flow; independent reset pointer cycles; optional one-bit physical feedforward",
              "alignment_distance": "delta_m = binomial majority error for p=(1-alpha)/2",
              "finite_protocol_has_no_postselection": True,
              "direct_destructive_rereading_amplifies_alignment": False,
              "certificates": [reference_certificate(count) for count in (1, 3, 9, 15)],
              "global_optimal_resource_cost_claimed": False,
              "quantum_theory_derived_from_cognition": False,
              "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("amplified_source_alignment_results.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
