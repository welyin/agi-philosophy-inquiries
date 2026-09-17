"""Round 151: physically compile coarse history and read only two label wires.

The third label wire remains in the inaccessible residual environment. It is
neither read for free nor erased. Gate and read counts are constructive bounds.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import unittest

import numpy as np

from centered_history_recovery import GROUPS, centered_decoder, centered_kraus, exact_error_interval, purification
from certified_intervals import Interval as I, SCALE
from classical_joining_history import independent_encoding
from coarsened_history_bound import pure_extension_metrics
from compiled_joining_frontier import compiled_sharing
from compiled_subject_joining import circuit_matrix
from compiled_three_subject_joining import remap, gate_counts
from external_correlation_recovery_bound import channel
from flagged_subject_joining import random_state
from independent_source_alignment import keep_systems
from joint_relation_records import interval_fields
from minimal_coherent_joining_history import conditional_unitary, environment_support
from noisy_classical_history import majority_error, physical_read_branch
from quantum_interface_audit import ALPHA
from walsh_coherent_history_compiler import walsh_pattern_ry, walsh_two_level_rotation, walsh_certificate


# Fine labels 0..5 -> (two reported bits, one unobserved residual bit).
FINE_CODES = (1,4,0,3,6,2)
NEGATIVE_PERMUTATION = (4,6,1,3,0,2,5,7)
NEGATIVE_ROTATIONS = ((12,14,F(-1)),(13,14,F(1,2)),(10,14,F(-1,2)),
                      (9,14,F(1,2)),(8,12,F(1,2)))


@lru_cache(maxsize=1)
def coarse_environment_circuit():
    gates = walsh_pattern_ry(1,((0,0),),2*np.arcsin(1/np.sqrt(3)))
    gates += walsh_pattern_ry(2,((0,0),),np.pi/2)
    gates += walsh_pattern_ry(3,((0,0),(1,1)),-np.pi/2)
    for first, second, angle in NEGATIVE_ROTATIONS:
        gates += walsh_two_level_rotation(first,second,float(angle)*np.pi)
    return gates


@lru_cache(maxsize=1)
def coarse_physical_isometry():
    return circuit_matrix(remap(coarse_environment_circuit(),{0:0,1:4,2:5,3:6}),
                          7,compiled_sharing(F(1,2))[0])


def confusion_matrix(error):
    if isinstance(error,I): valid=0<=error.lo<=error.hi<=SCALE//2
    else: valid=0<=error<=F(1,2)
    if not valid: raise ValueError("A symmetric bit error in [0,1/2] is required.")
    return tuple(tuple(error**((g^h).bit_count())*(1-error)**(2-(g^h).bit_count())
                       for h in range(4)) for g in range(4))


def noisy_kraus(error, report=None):
    matrix=confusion_matrix(error)
    if report is not None and (type(report) is not int or report not in range(4)):
        raise ValueError("Use a two-bit report 0..3.")
    return tuple(np.sqrt(float(matrix[g][h])/6)*centered_decoder(h) @ conditional_unitary(r+1)
                 for g, labels in enumerate(GROUPS) for r in labels
                 for h in (range(4) if report is None else (report,)))


def noisy_error_upper(error):
    delta=exact_error_interval()
    return delta+(1-delta)*(2*error-error**2)


def actual_two_read_recovery(source, contrast=.5, external_dimension=1):
    external_bits=external_dimension.bit_length()-1
    if 2**external_bits!=external_dimension: raise ValueError("Power-of-two external dimension required.")
    v=np.kron(coarse_physical_isometry(),np.eye(external_dimension))
    whole=v @ source @ v.T; result=np.zeros_like(source)
    for guessed in product((0,1),repeat=2):
        branch=whole
        for site, bit in zip((4,5),guessed): branch=physical_read_branch(branch,site,bit,contrast)
        active=keep_systems(branch,(0,1,2,3)+tuple(range(7,7+external_bits)),7+external_bits)
        decoder=np.kron(centered_decoder(2*guessed[0]+guessed[1]),np.eye(external_dimension))
        result+=decoder @ active @ decoder.T
    return result


def certificate(count,contrast=F(1,2)):
    error=majority_error(count,contrast); bound=noisy_error_upper(error)
    cost=gate_counts(coarse_environment_circuit()); reads=2*count
    return {"reads_per_label_bit":count,"contrast_exact":str(F(contrast)),
        "majority_bit_error":interval_fields(error),
        "promised_external_error_upper":interval_fields(bound),
        "strictly_beats_classical_only_half_error":bound.hi<SCALE//2,
        "original_noisy_reads":reads,"fresh_pointer_initializations":reads,
        "decoded_message_bits":2,"accessible_coherent_history_rebits":1,
        "unread_residual_history_rebits":1,"raw_readout_bits_retained":reads,
        "new_yx_after_first_join":cost["yx"]+reads+16,
        "new_ry_after_first_join":cost["ry"]+3*reads+17,
        "new_total_native_gates_after_first_join":cost["all"]+4*reads+33,
        "all_readout_results_accepted":True,
        "source_prior_join_transport_calibration_and_microscopic_bath_costs_included":False,
        "gate_angles_and_quantum_history_storage_ideal_in_this_table":True,
        "average_channel_bound_not_a_per_report_guarantee":True,
        "globally_optimal_error_or_read_count_claimed":False}


class NoisyTwoBitHistoryTests(unittest.TestCase):
    def test_environment_transform_is_real_orthogonal_and_maps_all_seven_occupied_vectors(self):
        matrix=circuit_matrix(coarse_environment_circuit(),4)
        np.testing.assert_allclose(matrix.T @ matrix,np.eye(16),atol=2e-14)
        np.testing.assert_allclose(matrix[8:,8:],np.eye(8)[:,NEGATIVE_PERMUTATION],atol=1e-14)
        target=np.zeros((16,7)); target[list(FINE_CODES),0]=1/np.sqrt(6)
        for r, code in enumerate(FINE_CODES): target[8+code,r+1]=1
        np.testing.assert_allclose(matrix @ environment_support(),target,atol=1e-14)

    def test_actual_seven_wire_circuit_leaves_exactly_the_required_coarse_branches(self):
        v=coarse_physical_isometry()
        branches=v.reshape((2,)*7+(16,)).transpose(4,5,6,0,1,2,3,7).reshape(8,16,16)
        for r, code in enumerate(FINE_CODES):
            np.testing.assert_allclose(branches[code],conditional_unitary(r+1)/np.sqrt(6),atol=1e-14)
        np.testing.assert_allclose(branches[[5,7]],np.zeros((2,16,16)),atol=1e-14)
        for g, labels in enumerate(GROUPS): self.assertEqual(set(FINE_CODES[r]//2 for r in labels),{g})

    def test_physical_two_pointer_reads_match_the_full_external_channel_without_reading_residual(self):
        rng=np.random.default_rng(151); raw=rng.normal(size=(32,3))+1j*rng.normal(size=(32,3))
        source=raw @ raw.conj().T; source/=np.trace(source)
        error=(1-ALPHA/2)/2
        expected=channel(tuple(np.kron(k,np.eye(2)) for k in noisy_kraus(error)),source)
        np.testing.assert_allclose(actual_two_read_recovery(source,.5,2),expected,atol=5e-15)

    def test_noisy_reports_remain_state_independent_with_the_nonuniform_group_prior(self):
        error=F(1,5); matrix=confusion_matrix(error)
        for row in matrix: self.assertEqual(sum(row),1)
        for report in range(4):
            prob=sum(F(len(GROUPS[g]),6)*matrix[g][report] for g in range(4))
            ks=noisy_kraus(error,report)
            np.testing.assert_allclose(sum(k.T @ k for k in ks),float(prob)*np.eye(16),atol=4e-16)
        np.testing.assert_allclose(sum(k.T @ k for k in noisy_kraus(error)),np.eye(16),atol=1e-15)

    def test_correct_group_component_is_the_exact_centered_channel(self):
        error=F(1,5); matrix=confusion_matrix(error)
        correct=tuple(np.sqrt(float(matrix[g][g])/6)*centered_decoder(g) @ conditional_unitary(r+1)
                      for g, labels in enumerate(GROUPS) for r in labels)
        rng=np.random.default_rng(251); raw=rng.normal(size=(16,3)); state=raw @ raw.T
        np.testing.assert_allclose(channel(correct,state),float((1-error)**2)*channel(centered_kraus(),state),atol=3e-15)

    def test_unknown_promised_external_purifications_obey_the_noisy_bound(self):
        rng=np.random.default_rng(351)
        for count in (1,3,15,25):
            error=sum(majority_error(count).floats())/2
            bound=sum(noisy_error_upper(majority_error(count)).floats())/2
            for _ in range(3):
                omega=independent_encoding((random_state(rng,2),random_state(rng,2)))
                self.assertLessEqual(pure_extension_metrics(noisy_kraus(error),purification(omega))["trace_error"],bound+1e-14)

    def test_interval_certificate_and_comparison_use_actual_read_and_gate_counts(self):
        self.assertEqual(gate_counts(coarse_environment_circuit()),{"yx":276,"ry":96,"all":372})
        result=certificate(3); six=walsh_certificate(3)
        self.assertEqual(result["original_noisy_reads"],6)
        self.assertEqual(result["fresh_pointer_initializations"],6)
        self.assertEqual(result["new_total_native_gates_after_first_join"],429)
        self.assertTrue(result["strictly_beats_classical_only_half_error"])
        self.assertFalse(certificate(1)["strictly_beats_classical_only_half_error"])
        self.assertEqual(six["new_total_native_gates_after_first_join"]-result["new_total_native_gates_after_first_join"],8)

    def test_stored_common_interface_is_unchanged_and_invalid_parameters_are_rejected(self):
        source=np.eye(16)/16; v=coarse_physical_isometry()
        np.testing.assert_allclose(keep_systems(v @ source @ v.T,(1,2,3),7),
                                   channel(compiled_sharing(F(1,2))[1],source),atol=4e-15)
        with self.assertRaises(ValueError): certificate(2)
        with self.assertRaises(ValueError): confusion_matrix(F(3,4))
        with self.assertRaises(ValueError): noisy_kraus(F(1,4),4)


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(NoisyTwoBitHistoryTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":151,"fine_label_codes":FINE_CODES,"negative_permutation":NEGATIVE_PERMUTATION,
        "negative_planar_rotations_pi_units":[[a,b,str(t)] for a,b,t in NEGATIVE_ROTATIONS],
        "environment_gate_counts":gate_counts(coarse_environment_circuit()),
        "new_pure_environment_helpers":0,"error_upper_formula":"delta4+(1-delta4)(2e-e^2)",
        "ideal_label_error_floor":interval_fields(exact_error_interval()),
        "closed_residual_source_purifications_pointers_and_baths_retained":True,
        "certificates":[certificate(m) for m in (1,3,15,25)],
        "comparison_scope":"specified equal odd repetitions per bit; unequal or adaptive reads not optimized",
        "automated_checks":{"run":checks.testsRun,"failures":len(checks.failures),"errors":len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name("noisy_two_bit_history_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    summary={k:v for k,v in report.items() if k!="certificates"}
    summary["certificates_summary"]=[{k:v for k,v in row.items() if k in
        ("original_noisy_reads","new_total_native_gates_after_first_join","promised_external_error_upper")} for row in report["certificates"]]
    print(json.dumps(summary,indent=2))


if __name__=="__main__": main()
