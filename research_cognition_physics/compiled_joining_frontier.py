"""Round 134: compile the full disturbance frontier and close the resource ledger.

The selector is a retained coherent rebit, not a free random classical coin.
The displayed gate counts are constructive upper bounds with ideal joint access.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
import math
from pathlib import Path
import unittest

import numpy as np

from approximate_joining_optimality import choi_from_kraus
from approximate_subject_joining import apply_kraus, trace_distance
from common_orientation_structure import encode_state
from compiled_subject_joining import (
    Gate, circuit_matrix, controlled_ry, controlled_yx, controlled_z,
    endpoint_circuit, intake_circuit, isometry_and_kraus, repair_circuit,
    resource_counts, ry,
)
from encoded_composition_audit import independent_encoding
from flagged_subject_joining import random_state
from joining_disturbance_sharing import predicted_logical_output, sharing_kraus
from mixed_state_joining_frontier import allowed_support_score, optimum
from noisy_joining_witness import (
    actual_record_probability, corrected_record_score, labelled_state, sample_certificate,
)
from quantum_interface_audit import ALPHA


def control_circuit(control, gates):
    result = ()
    for gate in gates:
        if control in gate.sites:
            raise ValueError("Selector must be separate from the controlled circuit.")
        if gate.kind == "Ry":
            result += controlled_ry(control, gate.sites[0], gate.angle)
        elif gate.kind == "YX":
            result += controlled_yx(control, *gate.sites, gate.angle)
        else:
            raise ValueError("Uncompiled gate.")
    return result


def sharing_circuit(newcomer_share):
    theta = float(newcomer_share)
    if not np.isfinite(theta) or not 0 <= theta <= 1:
        raise ValueError("The newcomer's share must be in [0,1].")
    if theta in (0., 1.):
        return endpoint_circuit(repair_host=theta == 0.), 6
    selector = 6
    host_repair = controlled_z(0, 1, 4) + repair_circuit(2)
    # Selector |1> has weight theta. Conjugation by Ry(pi) controls the other arm on |0>.
    gates = (intake_circuit() + (ry(selector, 2*math.asin(math.sqrt(theta))),)
             + control_circuit(selector, repair_circuit(3))
             + (ry(selector, math.pi),) + control_circuit(selector, host_repair)
             + (ry(selector, -math.pi),))
    return gates, 7


@lru_cache(maxsize=8)
def compiled_sharing(newcomer_share):
    gates, systems = sharing_circuit(newcomer_share)
    return isometry_and_kraus(gates, systems)


def compiled_witness_mean(radius, newcomer_share=F(1, 2), contrast=.5):
    _, ks = compiled_sharing(newcomer_share)
    result = 0.
    for aa, sa, ab, sb in product(range(3), (-1, 1), range(3), (-1, 1)):
        states = labelled_state(aa, sa, radius), labelled_state(ab, sb, radius)
        output = apply_kraus(ks, independent_encoding(states))
        for ra, rb in product((-1, 1), repeat=2):
            probability = actual_record_probability(output, aa, ab, ra, rb, contrast)
            result += probability*corrected_record_score(sa, sb, ra, rb, ALPHA*contrast)/36
    return result


def workflow_resources(radius, newcomer_share):
    old = sample_certificate(radius)
    gates, systems = sharing_circuit(newcomer_share)
    join = resource_counts(gates, systems-4)
    count = old["sufficient_independent_fresh_inputs"]
    pair = 8 + join["joint_yx_rotations"]
    local = 20 + join["local_ry_rotations"]
    resets = 8 + systems - 4
    return {
        "radius_exact": str(F(radius)),
        "newcomer_share_exact": str(F(newcomer_share)),
        "input_pairs": count,
        "per_trial_pure_rebit_initializations": resets,
        "per_trial_pair_rotation_upper_bound": pair,
        "per_trial_local_rotation_upper_bound": local,
        "per_trial_original_noisy_reads": 2,
        "total_pure_rebit_initializations": resets*count,
        "total_pair_rotation_upper_bound": pair*count,
        "total_local_rotation_upper_bound": local*count,
        "total_original_noisy_reads": 2*count,
        "retained_source_purification_rebits_per_trial": 4,
        "retained_joining_environment_rebits_per_trial": systems-3,
        "joining_resources": join,
        "label_bits_if_revealed_after_commitment": 6*count,
        "outcome_bits": 2*count,
        "fair_label_random_bits_expected": 8*count,
        "ideal_gate_angles_and_calibrated_readout_assumed": True,
        "source_purifications_or_labels_used_by_joiner": False,
        "microscopic_readout_bath_and_label_rng_costs_included": False,
        "discarded_systems_or_erasure_assumed": False,
    }


def per_gate_angle_tolerance(newcomer_share, extra_error=F(1, 1000)):
    epsilon = F(extra_error)
    if not 0 < epsilon <= 1:
        raise ValueError("An additional trace-distance budget in (0,1] is required.")
    gates, _ = sharing_circuit(newcomer_share)
    return 2*epsilon / len(gates)


class CompiledJoiningFrontierTests(unittest.TestCase):
    def test_selector_compilation_matches_control_of_a_complete_noncommuting_circuit(self):
        from compiled_subject_joining import yx
        old = (ry(1, .3), yx(1, 2, -.8), ry(2, -.4))
        actual = circuit_matrix(control_circuit(0, old), 3)
        target = circuit_matrix((ry(0, .3), yx(0, 1, -.8), ry(1, -.4)), 2)
        expected = np.zeros((8, 8))
        expected[:4, :4] = np.eye(4)
        expected[4:, 4:] = target
        np.testing.assert_allclose(actual, expected, atol=1e-15)

    def test_entire_channel_matches_all_frontier_shares(self):
        for theta in (F(0), F(1, 4), F(1, 2), F(3, 4), F(1)):
            v, ks = compiled_sharing(theta)
            np.testing.assert_allclose(v.T@v, np.eye(16), atol=2e-14)
            np.testing.assert_allclose(choi_from_kraus(ks), choi_from_kraus(sharing_kraus(theta)), atol=1e-14)

    def test_both_unknown_states_follow_the_claimed_logical_formula(self):
        rng = np.random.default_rng(134)
        for theta in (F(0), F(1, 4), F(1, 2), F(1)):
            _, ks = compiled_sharing(theta)
            for _ in range(3):
                host, new = random_state(rng, 2), random_state(rng, 2)
                actual = apply_kraus(ks, independent_encoding((host, new)))
                np.testing.assert_allclose(actual,
                    encode_state(predicted_logical_output(host, new, theta)), atol=3e-15)

    def test_balanced_circuit_attains_mixed_state_optimum_in_all_six_directions(self):
        _, ks = compiled_sharing(F(1, 2))
        for radius in (F(1, 10), F(1, 2), F(1)):
            for aa, ab, sa, sb in product(range(3), range(3), (-1, 1), (-1, 1)):
                host, new = labelled_state(aa, sa, radius), labelled_state(ab, sb, radius)
                actual = apply_kraus(ks, independent_encoding((host, new)))
                self.assertAlmostEqual(trace_distance(actual, encode_state(np.kron(host, new))),
                                       float(optimum(radius)), places=13)

    def test_entire_seven_rebit_circuit_retains_external_correlations(self):
        gates, systems = sharing_circuit(F(1, 2))
        whole = circuit_matrix(gates, systems)
        np.testing.assert_allclose(whole.T@whole, np.eye(128), atol=2e-14)
        rng = np.random.default_rng(234)
        original = rng.normal(size=(128, 2))
        original /= np.linalg.norm(original)
        changed = circuit_matrix(gates, systems, original)
        inverse = tuple(Gate(g.kind, g.sites, -g.angle) for g in reversed(gates))
        np.testing.assert_allclose(circuit_matrix(inverse, systems, changed), original, atol=2e-15)

    def test_original_noisy_records_close_the_preparation_joining_detection_chain(self):
        for radius in (F(1, 10), F(1, 2), F(1)):
            self.assertAlmostEqual(compiled_witness_mean(radius), float(allowed_support_score(radius)), places=13)

    def test_complete_resource_ledger_and_endpoints(self):
        for theta, counts in ((F(0), (23, 12, 2)), (F(1, 2), (131, 35, 3)), (F(1), (18, 7, 2))):
            gates, systems = sharing_circuit(theta)
            report = resource_counts(gates, systems-4)
            self.assertEqual((report["joint_yx_rotations"], report["local_ry_rotations"],
                              report["extra_pure_rebit_initializations"]), counts)
        row = workflow_resources(F(1), F(1, 2))
        self.assertEqual(row["input_pairs"], 3555)
        self.assertEqual(row["total_pure_rebit_initializations"], 39105)
        self.assertEqual(row["total_pair_rotation_upper_bound"], 494145)
        self.assertEqual(row["total_local_rotation_upper_bound"], 195525)
        self.assertEqual(row["total_original_noisy_reads"], 7110)

    def test_angle_precision_bound_and_invalid_inputs(self):
        self.assertEqual(per_gate_angle_tolerance(F(1)), F(1, 12500))
        self.assertEqual(per_gate_angle_tolerance(F(1, 2)), F(1, 83000))
        theta = F(1, 2)
        gates, systems = sharing_circuit(theta)
        errors = [((-1)**i)*1e-4 for i in range(len(gates))]
        perturbed = tuple(Gate(g.kind, g.sites, g.angle+e) for g, e in zip(gates, errors))
        _, exact_ks = compiled_sharing(theta)
        _, changed_ks = isometry_and_kraus(perturbed, systems)
        rng = np.random.default_rng(334)
        source = independent_encoding((random_state(rng, 2), random_state(rng, 2)))
        error = trace_distance(apply_kraus(exact_ks, source), apply_kraus(changed_ks, source))
        self.assertLessEqual(error, sum(abs(e) for e in errors)/2 + 1e-14)
        self.assertGreater(error, 1e-8)
        for bad in (-.1, 1.1, float("nan")):
            with self.assertRaises(ValueError):
                sharing_circuit(bad)
        with self.assertRaises(ValueError):
            per_gate_angle_tolerance(theta, F(0))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(CompiledJoiningFrontierTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 134,
        "all_frontier_shares_compiled_with_retained_coherent_selector": True,
        "extra_selector_pure_initialization_counted": True,
        "balanced_common_error_exact": "r(1+r)/12",
        "balanced_each_marginal_error_exact": "r/12",
        "old_host_exact_endpoint_common_error_exact": "r/6",
        "native_primitives": ["Ry", "oriented YX"],
        "joining_readouts": 0,
        "all_environments_retained": True,
        "gate_noise_model": "Only additive primitive rotation angle errors in the robustness bound",
        "additional_trace_distance_bound": "min(1, sum_j abs(delta_angle_j)/2)",
        "per_gate_angle_tolerance_for_extra_error_1_over_1000": {
            "host_exact_endpoint_radians": str(per_gate_angle_tolerance(F(1))),
            "balanced_radians": str(per_gate_angle_tolerance(F(1, 2))),
        },
        "sample_counts_assume_ideal_source_and_calibrated_readout": True,
        "counts_claimed_gate_optimal": False,
        "experimental_data_collected": False,
        "workflow_rows": [workflow_resources(r, theta)
                          for r in (F(1), F(1, 2), F(1, 10)) for theta in (F(1), F(1, 2))],
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("compiled_joining_frontier_results.json").write_text(
            json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
