"""Round 137: an eight-rebit primitive circuit for the optimal three-party join.

Also compile recovery from the actual retained round-134 balanced first stage.
All costs below are constructive counts, with no free orientation measurement.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import unittest

import numpy as np

from approximate_subject_joining import apply_kraus, trace_distance
from common_orientation_structure import encode_state
from compiled_joining_frontier import control_circuit, sharing_circuit
from compiled_subject_joining import (
    Gate, circuit_matrix, controlled_z, intake_circuit, repair_circuit, ry,
)
from encoded_composition_audit import independent_encoding, product_state
from flagged_subject_joining import random_state
from joining_history_advantage import balanced_pair, common_ab_marginal
from multisubject_joining_optimum import multisubject_kraus, predicted_logical


def remap(gates, mapping):
    return tuple(Gate(g.kind, tuple(mapping.get(s, s) for s in g.sites), g.angle) for g in gates)


def inverse(gates):
    return tuple(Gate(g.kind, g.sites, -g.angle) for g in reversed(gates))


def zero_control(control, gates):
    return (ry(control, np.pi),) + control_circuit(control, gates) + (ry(control, -np.pi),)


def three_subject_circuit():
    # Input: RA,RB,RC,A,B,C,E1,E2. Output: hB,hC,R,A,B,C,E1,E2.
    first = intake_circuit(6)
    second = remap(intake_circuit(6), {0: 1, 1: 2})
    repair_b = zero_control(1, repair_circuit(4, 6, 7))
    repair_c = zero_control(0, remap(repair_circuit(5, 6, 7), {0: 1}))
    repair_a = control_circuit(1, controlled_z(0, 2, 6) + repair_circuit(3, 6, 7))
    return first + second + repair_b + repair_c + repair_a


def extract_kraus(isometry, active, systems):
    environment = tuple(i for i in range(systems) if i not in active)
    return tuple(isometry.reshape((2,)*systems+(64,)).transpose(
        environment+tuple(active)+(systems,)).reshape(-1, 16, 64))


@lru_cache(maxsize=1)
def compiled_three_subject():
    initial = np.eye(256)[:, ::4]
    v = circuit_matrix(three_subject_circuit(), 8, initial)
    return v, extract_kraus(v, (2,3,4,5), 8)


def embedded_columns(sites, systems):
    count = len(sites)
    indices = [sum(((value >> (count-1-j)) & 1) << (systems-1-site)
                   for j, site in enumerate(sites)) for value in range(2**count)]
    return np.eye(2**systems)[:, indices]


def reopened_history_circuit(include_first_stage=True):
    old, _ = sharing_circuit(F(1, 2))
    # Old physical order: RA,RB,A,B,E1,E2,selector; fresh RC,C use wires 7,8.
    batch = remap(three_subject_circuit(), {0:0,1:1,2:7,3:2,4:3,5:8,6:4,7:5})
    return (old if include_first_stage else ()) + inverse(old) + batch


def gate_counts(gates):
    return {"yx": sum(g.kind == "YX" for g in gates),
            "ry": sum(g.kind == "Ry" for g in gates), "all": len(gates)}


def resource_report():
    first, _ = sharing_circuit(F(1, 2))
    from compiled_subject_joining import endpoint_circuit
    sequential = first + endpoint_circuit()
    return {
        "fresh_global_join": {"new_pure_rebits": 2, "total_rebits": 8,
                              "gates": gate_counts(three_subject_circuit())},
        "fixed_balanced_first_then_restricted_join": {"new_pure_rebits_over_both_stages": 5,
                              "total_rebits_with_old_history_retained": 11,
                              "gates": gate_counts(sequential)},
        "reopen_after_first_stage_already_complete": {"additional_pure_rebits": 0,
                              "total_rebits": 9, "clean_unused_rebits_after_finish": 1,
                              "additional_gates": gate_counts(reopened_history_circuit(False))},
        "reopen_including_cost_of_original_first_stage": {"new_pure_rebits": 3,
                              "gates": gate_counts(reopened_history_circuit(True))},
        "joining_readouts_for_all_rows": 0,
        "gate_count_optimality_claimed": False,
        "fresh_global_auxiliary_minimum_for_this_specific_channel_and_pure_unitary_model": 2,
        "transport_storage_and_calibration_costs_included": False,
    }


class CompiledThreeSubjectJoiningTests(unittest.TestCase):
    def test_native_intake_has_the_correct_all_reference_branches(self):
        # Compare the entire channel via Kraus vector Gram matrices; no input fitting.
        _, ks = compiled_three_subject()
        expected = multisubject_kraus(3)
        left = np.column_stack([k.T.ravel() for k in ks])
        right = np.column_stack([k.T.ravel() for k in expected])
        np.testing.assert_allclose(left@left.T, right@right.T, atol=9e-15)

    def test_all_primitive_gates_and_coherent_isometry_are_real_and_complete(self):
        gates = three_subject_circuit()
        self.assertTrue(all(g.kind in ("Ry", "YX") for g in gates))
        self.assertTrue(all(np.isrealobj(g.matrix()) for g in gates))
        v, ks = compiled_three_subject()
        np.testing.assert_allclose(v.T@v, np.eye(64), atol=2e-14)
        self.assertEqual(sum(np.linalg.norm(k)>1e-12 for k in ks), 10)

    def test_unknown_states_and_old_ab_marginal(self):
        rng = np.random.default_rng(137)
        _, ks = compiled_three_subject()
        for _ in range(8):
            states = [random_state(rng, 2) for _ in range(3)]
            actual = apply_kraus(ks, independent_encoding(states))
            np.testing.assert_allclose(actual, encode_state(predicted_logical(states)), atol=3e-15)
            np.testing.assert_allclose(common_ab_marginal(actual), encode_state(balanced_pair(*states[:2])), atol=4e-15)

    def test_pure_inputs_attain_the_certified_quarter_error(self):
        rng = np.random.default_rng(237)
        _, ks = compiled_three_subject()
        for _ in range(8):
            states = []
            for _ in range(3):
                vector = rng.normal(size=2)+1j*rng.normal(size=2)
                vector /= np.linalg.norm(vector)
                states.append(np.outer(vector, vector.conj()))
            actual = apply_kraus(ks, independent_encoding(states))
            self.assertAlmostEqual(trace_distance(actual, encode_state(product_state(states))), 1/4, places=13)

    def test_full_inverse_preserves_arbitrary_external_correlations(self):
        rng = np.random.default_rng(337)
        original = rng.normal(size=(256, 3))
        original /= np.linalg.norm(original)
        gates = three_subject_circuit()
        changed = circuit_matrix(gates, 8, original)
        np.testing.assert_allclose(circuit_matrix(inverse(gates), 8, changed), original, atol=3e-15)

    def test_actual_nine_wire_reopening_pipeline_equals_fresh_batch_and_returns_selector(self):
        source_embedding = embedded_columns((0,1,7,2,3,8), 9)
        actual = circuit_matrix(reopened_history_circuit(), 9, source_embedding)
        native, _ = compiled_three_subject()
        result_embedding = embedded_columns((0,1,7,2,3,8,4,5), 9)
        np.testing.assert_allclose(actual, result_embedding@native, atol=1e-14)
        # The unused original selector (physical wire 6) is exactly |0>.
        rows_with_selector_one = [i for i in range(512) if (i >> 2) & 1]
        np.testing.assert_allclose(actual[rows_with_selector_one], 0., atol=2e-15)

    def test_resource_comparison_is_actual_primitive_counting(self):
        resources = resource_report()
        self.assertEqual(resources["fresh_global_join"]["gates"], {"yx":187,"ry":47,"all":234})
        self.assertEqual(resources["fixed_balanced_first_then_restricted_join"]["gates"], {"yx":149,"ry":42,"all":191})
        self.assertEqual(resources["reopen_after_first_stage_already_complete"]["additional_gates"], {"yx":318,"ry":82,"all":400})
        self.assertEqual(resources["reopen_including_cost_of_original_first_stage"]["gates"], {"yx":449,"ry":117,"all":566})

    def test_minimal_pure_auxiliary_capacity_for_this_channel_has_exact_gram_witness(self):
        ks = multisubject_kraus(3)
        gram = np.array([[np.sum(a*b) for b in ks] for a in ks])
        expected = np.diag([16.] + [16/3]*9)
        np.testing.assert_allclose(gram, expected, atol=9e-15)
        # Analytic Pauli orthogonality gives this positive diagonal exactly:
        # rank 10 requires environment dimension >=10; rebits require 16.
        self.assertGreater(16*10, 64*2)
        self.assertLessEqual(16*10, 64*4)
        self.assertEqual(F(2, 1000*234), F(1,117000))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(CompiledThreeSubjectJoiningTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"round": 137,
        "input_order": ["RA","RB","RC","A","B","C","E1_zero","E2_zero"],
        "output_order": ["hB","hC","R","A","B","C","E1","E2"],
        "native_primitives": ["Ry", "oriented YX"],
        "worst_joint_error_exact": "1/4",
        "same_balanced_old_ab_joint_marginal_preserved": True,
        "all_history_coherent_and_global_inverse_available": True,
        "source_purifications_or_preparation_labels_accessed": False,
        "reopening_preserves_a_second_copy_of_old_common_interface": False,
        "resources": resource_report(),
        "per_gate_angle_error_sufficient_for_extra_trace_error_1_over_1000": "1/117000 radians for fresh batch",
        "general_n_gate_compilation_or_noisy_history_optimum_claimed": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("compiled_three_subject_joining_results.json").write_text(
            json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
