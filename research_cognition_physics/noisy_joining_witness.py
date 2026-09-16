"""Round 132: finite old-readout witness for the declared encoded interface.

Each logical Y read uses a real YX reference-target gate and the original noisy
Z instrument. Sources are prepared with pure resets, real gates, and retained
purifying systems. Preparation labels and purifications are hidden from the
joining device. This tests a conversion promise, not the field of nature.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path

import numpy as np

from bilocal_record_tomography import embed_operator
from bipartite_composition import interaction
from certified_intervals import Interval as I, SCALE, sin_interval, ceil_div
from common_orientation_structure import encode_state
from complex_control_from_reference import real_lift
from independent_source_alignment import keep_systems
from joining_disturbance_sharing import sharing_output
from joint_relation_records import interval_fields
from mixed_state_joining_frontier import (
    allowed_support_score, ideal_support_score, optimum, validate_radius,
)
from operational_effect_closure import compiled_copy_gate
from quantum_interface_audit import (
    ALPHA, PAULI_X, PAULI_Y, PAULI_Z, branch_kraus, rotation_unitary,
)
from reference_information_cost import log_interval


AXES = (PAULI_X, PAULI_Y, PAULI_Z)


def labelled_state(axis, sign, radius=F(1)):
    if axis not in (0, 1, 2) or sign not in (-1, 1):
        raise ValueError("Use an X/Y/Z axis and a sign +/-1.")
    return (np.eye(2) + float(validate_radius(radius)) * sign * AXES[axis]) / 2


def prepared_source_whole(axis, sign, radius=F(1)):
    """Four physical rebits in order R,E,T,F; E,F remain in the whole."""
    labelled_state(axis, sign, radius)
    r = float(validate_radius(radius))
    vector = np.zeros(16)
    vector[0] = 1.
    operations = (
        embed_operator(rotation_unitary(math.pi/2), (0,), 4),
        embed_operator(compiled_copy_gate(), (1, 0), 4),
        embed_operator(rotation_unitary(math.acos(r)), (2,), 4),
        embed_operator(compiled_copy_gate(), (3, 2), 4),
    )
    for gate in operations:
        vector = gate @ vector
    if axis == 0:
        vector = embed_operator(rotation_unitary(sign * math.pi/2), (2,), 4) @ vector
    elif axis == 1:
        vector = embed_operator(interaction(-sign * math.pi/2, "YX"), (0, 2), 4) @ vector
    elif sign == -1:
        vector = embed_operator(rotation_unitary(math.pi), (2,), 4) @ vector
    if np.max(np.abs(vector.imag)) > 1e-15:
        raise ArithmeticError("The declared source circuit must remain real.")
    return np.outer(vector.real, vector.real)


def measurement_rotation(axis, target):
    if target not in (1, 2) or axis not in (0, 1, 2):
        raise ValueError("Use a target A/B and an X/Y/Z axis.")
    if axis == 0:
        return embed_operator(rotation_unitary(-math.pi/2), (target,), 3).real
    if axis == 1:
        return embed_operator(interaction(math.pi/2, "YX"), (0, target), 3).real
    return np.eye(8)


@lru_cache(maxsize=64)
def actual_read_kraus(axis, target, sign, contrast=.5):
    if sign not in (-1, 1):
        raise ValueError("Use an outcome sign +/-1.")
    rotate = measurement_rotation(axis, target)
    return tuple(embed_operator(k, (target,), 3).real @ rotate
                 for k in branch_kraus(0., (sign+1)//2, contrast))


def actual_record_probability(encoded, axis_a, axis_b, sign_a, sign_b, contrast=.5):
    after_a = sum(k @ encoded @ k.T for k in actual_read_kraus(axis_a, 1, sign_a, contrast))
    return float(sum(np.trace(k @ after_a @ k.T)
                     for k in actual_read_kraus(axis_b, 2, sign_b, contrast)))


def corrected_record_score(prep_a, prep_b, result_a, result_b, visibility):
    if not 0 < visibility <= 1:
        raise ValueError("A positive calibrated visibility is required.")
    return (1 + prep_a*result_a/visibility) * (1 + prep_b*result_b/visibility) / 4


def mean_witness(radius, actual_join=True, contrast=.5):
    result = 0.
    for axis_a, prep_a, axis_b, prep_b in product(range(3), (-1, 1), range(3), (-1, 1)):
        a, b = labelled_state(axis_a, prep_a, radius), labelled_state(axis_b, prep_b, radius)
        output = sharing_output(a, b, .5) if actual_join else encode_state(np.kron(a, b))
        for sign_a, sign_b in product((-1, 1), repeat=2):
            p = actual_record_probability(output, axis_a, axis_b, sign_a, sign_b, contrast)
            result += p * corrected_record_score(prep_a, prep_b, sign_a, sign_b, ALPHA*contrast) / 36
    return result


def sample_certificate(radius=F(1), contrast=F(1, 2), total_error=F(1, 100)):
    r, eta, risk = validate_radius(radius), F(contrast), F(total_error)
    if r <= 0 or not 0 < eta <= 1 or not 0 < risk < 1:
        raise ValueError("Positive radius, contrast and a risk in (0,1) are required.")
    visibility = 4 * sin_interval(F(1, 4)) * I.exact(eta)
    score_range = (1 + visibility) / (2 * visibility**2)
    gap = I.exact(optimum(r))
    threshold = (ideal_support_score(r) + allowed_support_score(r)) / 2
    count_bound = 2 * score_range**2 * log_interval(2/risk) / gap**2
    count = ceil_div(count_bound.hi, SCALE)
    if count != ceil_div(count_bound.lo, SCALE):
        raise ArithmeticError("The interval does not determine the integer sample count.")
    # Two one-sided errors are each <= risk/2 by the bounded conditional
    # exponential-moment inequality. Count is sufficient for this bound,
    # not claimed statistically optimal.
    return {
        "radius_exact": str(r),
        "old_readout_contrast_exact": str(eta),
        "total_two_sided_error_budget_exact": str(risk),
        "visibility_interval": interval_fields(visibility),
        "corrected_score_range_interval": interval_fields(score_range),
        "ideal_mean_exact": str(ideal_support_score(r)),
        "all_channel_mean_upper_bound_exact": str(allowed_support_score(r)),
        "mean_gap_exact": str(optimum(r)),
        "decision_threshold_exact": str(threshold),
        "sufficient_independent_fresh_inputs": count,
        "sample_count_expression_interval": interval_fields(count_bound),
        "target_original_noisy_reads": 2*count,
        "pure_rebit_resets_for_explicit_source_circuit": 8*count,
        "source_and_verifier_pair_rotations_upper_bound": 8*count,
        "source_and_verifier_local_rotations_upper_bound": 20*count,
        "source_purifying_rebits_retained": 4*count,
        "preparation_label_fixed_length_bits_if_revealed_after_commitment": 6*count,
        "measurement_outcome_bits": 2*count,
        "unbiased_random_bits_expected_for_rejection_sampled_labels": 8*count,
        "randomness_generation_and_joining_device_costs_additional": True,
    }


class NoisyJoiningWitnessTests(unittest.TestCase):
    def test_actual_purified_real_source_prepares_all_six_encoded_states(self):
        for radius in (F(0), F(1, 2), F(1)):
            for axis, sign in product(range(3), (-1, 1)):
                whole = prepared_source_whole(axis, sign, radius)
                np.testing.assert_allclose(whole @ whole, whole, atol=6e-16)
                np.testing.assert_allclose(keep_systems(whole, (0, 2), 4),
                                           encode_state(labelled_state(axis, sign, radius)), atol=4e-16)

    def test_all_logical_axes_are_read_by_real_gates_and_old_noisy_z(self):
        for axis, target, sign in product(range(3), (1, 2), (-1, 1)):
            ks = actual_read_kraus(axis, target, sign)
            effect = sum(k.T @ k for k in ks)
            local = (np.eye(2) + sign*(ALPHA/2)*AXES[axis]) / 2
            logical = np.kron(local, np.eye(2)) if target == 1 else np.kron(np.eye(2), local)
            np.testing.assert_allclose(effect, real_lift(logical), atol=3e-16)
            self.assertTrue(all(np.isrealobj(k) for k in ks))

    def test_sequential_shared_reference_measurements_have_correct_product_effects(self):
        from flagged_subject_joining import random_state
        rng = np.random.default_rng(132)
        encoded = encode_state(random_state(rng, 4))
        for a, b in product(range(3), repeat=2):
            total = 0.
            for sa, sb in product((-1, 1), repeat=2):
                actual = actual_record_probability(encoded, a, b, sa, sb)
                ea = (np.eye(2) + sa*(ALPHA/2)*AXES[a]) / 2
                eb = (np.eye(2) + sb*(ALPHA/2)*AXES[b]) / 2
                expected = np.trace(encoded @ real_lift(np.kron(ea, eb)))
                self.assertAlmostEqual(actual, expected)
                total += actual
            self.assertAlmostEqual(total, 1.)

    def test_complete_noisy_records_recover_the_certified_mean_scores(self):
        for radius in (F(1, 2), F(1)):
            self.assertAlmostEqual(mean_witness(radius, True), float(allowed_support_score(radius)))
            self.assertAlmostEqual(mean_witness(radius, False), float(ideal_support_score(radius)))

    def test_score_range_and_signed_weights_are_explicit(self):
        v = ALPHA / 2
        scores = [corrected_record_score(1, 1, a, b, v) for a, b in product((-1, 1), repeat=2)]
        self.assertLess(min(scores), 0.)
        self.assertGreater(max(scores), 1.)
        self.assertAlmostEqual(max(scores)-min(scores), (1+v)/(2*v*v))
        self.assertEqual(len(scores), 4)

    def test_sample_counts_are_certified_by_intervals_and_grow_with_weaker_inputs(self):
        rows = [sample_certificate(r) for r in (F(1), F(1, 2), F(1, 10))]
        counts = [row["sufficient_independent_fresh_inputs"] for row in rows]
        self.assertTrue(counts[0] < counts[1] < counts[2])
        for row in rows:
            n = row["sufficient_independent_fresh_inputs"]
            self.assertEqual(row["target_original_noisy_reads"], 2*n)
            self.assertEqual(row["source_purifying_rebits_retained"], 4*n)

    def test_revealed_preparation_labels_change_the_task(self):
        # Preparing the requested ideal target from an early supplied label
        # beats the bound, but is not one fixed channel on the supplied input.
        for radius in (F(1, 2), F(1)):
            self.assertGreater(mean_witness(radius, False), float(allowed_support_score(radius)))

    def test_zero_visibility_or_radius_does_not_provide_a_finite_separation(self):
        for args in ((F(0), F(1, 2)), (F(1), F(0)), (F(1), F(2))):
            with self.assertRaises(ValueError):
                sample_certificate(*args)
        with self.assertRaises(ValueError):
            corrected_record_score(1, 1, 1, 1, 0.)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(NoisyJoiningWitnessTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 132,
        "logical_y_readout": "One old real YX reference-target rotation followed by the old noisy Z instrument",
        "extra_reference_and_joint_control_permissions_explicit": True,
        "new_native_single_rebit_y_measurement_assumed": False,
        "input_preparations_compiled_with_purifications_retained": True,
        "all_outcome_records_used": True,
        "preparation_labels_hidden_until_common_output_committed": True,
        "source_purifications_inaccessible_to_joining_device": True,
        "fresh_secret_labels_required_for_history_adaptive_device_bound": True,
        "calibrated_readout_visibility_assumed_known": True,
        "score_is_an_unbounded_probability": False,
        "score_is_a_bounded_signed_statistic": True,
        "actual_experimental_samples_collected": False,
        "statistically_optimal_sample_count_claimed": False,
        "joining_channel_itself_compiled_into_old_noisy_primitives": False,
        "test_claims_to_identify_real_vs_complex_field_of_nature": False,
        "rows": [sample_certificate(r) for r in (F(1), F(1, 2), F(1, 10))],
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("noisy_joining_witness_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
