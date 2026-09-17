"""Round 147: finite native label reads while retaining the coherent history bit.

Every input, including arbitrary external correlations, has the same uniform
true-label probabilities. A classical confusion risk therefore bounds the
complete recovery trace error. The risk bound is not claimed to be tight for
the quantum channel or globally optimal in read count.
"""

import argparse
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import unittest

import numpy as np

from approximate_subject_joining import trace_distance
from certified_intervals import Interval as I, SCALE
from compiled_coherent_history import compressed_physical_isometry, resources
from coherent_history_noise_bound import phase_kraus
from external_correlation_recovery_bound import channel, max_entangled_metrics
from independent_source_alignment import keep_systems
from joint_relation_records import interval_fields
from minimal_coherent_joining_history import conditional_unitary
from noisy_classical_history import majority_error, physical_read_branch
from quantum_interface_audit import ALPHA


def decode_label(raw):
    if not isinstance(raw,int) or isinstance(raw,bool) or not 0<=raw<8:
        raise ValueError("Use an integer three-bit majority label.")
    return raw if raw<6 else raw-4  # 110 -> 010, 111 -> 011: nearest valid labels.


def confusion_matrix(bit_error):
    e=bit_error
    if isinstance(e,I): valid=0<=e.lo<=e.hi<=SCALE//2
    else: valid=0<=e<=F(1,2)
    if not valid: raise ValueError("A symmetric majority error in [0,1/2] is required.")
    result=[[0 for _ in range(6)] for _ in range(6)]
    for true,raw in product(range(6),range(8)):
        errors=(true ^ raw).bit_count()
        result[true][decode_label(raw)]+=e**errors*(1-e)**(3-errors)
    return result


def label_risk(bit_error):
    e=bit_error
    return (8*e-7*e**2+2*e**3)/3


def conditional_label_risks(bit_error):
    matrix=confusion_matrix(bit_error)
    return tuple(1-matrix[g][g]/sum(matrix[r][g] for r in range(6)) for g in range(6))


def recovery_kraus(bit_error,reported_label=None):
    matrix=confusion_matrix(bit_error)
    if reported_label is not None and reported_label not in range(6): raise ValueError("Use a valid decoded label.")
    result=[]
    for r in range(6):
        for g in range(6) if reported_label is None else (reported_label,):
            result.append(np.sqrt(float(matrix[r][g])/6)*conditional_unitary(g+1).T @ conditional_unitary(r+1))
    return tuple(result)


def actual_one_read_recovery(source,contrast=.5,external_dimension=1):
    v=np.kron(compressed_physical_isometry(),np.eye(external_dimension))
    whole=v @ source @ v.T
    external_bits=external_dimension.bit_length()-1
    if 2**external_bits!=external_dimension: raise ValueError("Use a power-of-two external dimension.")
    result=np.zeros_like(source)
    for guessed in product((0,1),repeat=3):
        branch=whole
        for site,bit in zip((4,5,6),guessed):
            branch=physical_read_branch(branch,site,bit,contrast)
        active=keep_systems(branch,(0,1,2,3)+tuple(range(7,7+external_bits)),7+external_bits)
        label=decode_label(4*guessed[0]+2*guessed[1]+guessed[2])
        decoder=np.kron(conditional_unitary(label+1).T,np.eye(external_dimension))
        result+=decoder @ active @ decoder.T
    return result


def certificate(count,contrast=F(1,2)):
    error=majority_error(count,contrast); risk=label_risk(error)
    base=resources(); reads=3*count
    return {"reads_per_label_bit":count,"contrast_exact":str(F(contrast)),
        "majority_bit_error":interval_fields(error),
        "uniform_external_recovery_error_upper":interval_fields(risk),
        "conditional_error_upper_by_report":[interval_fields(v) for v in conditional_label_risks(error)],
        "strictly_beats_classical_only_half_error":risk.hi<SCALE//2,
        "original_noisy_reads":reads,"fresh_pointer_initializations":reads,
        "retained_accessible_coherent_history_rebits":1,"decoded_message_bits":3,
        "raw_readout_bits_retained":reads,
        "new_yx_after_first_join":base["yx"]+reads+16,
        "new_ry_after_first_join":base["ry"]+3*reads+17,
        "new_total_native_gates_after_first_join":base["all"]+4*reads+33,
        "old_common_marginal_preserved_until_recovery_starts":True,
        "all_readout_results_accepted":True,
        "quantum_history_read_in_computational_basis":False,
        "gate_angles_and_quantum_history_storage_ideal_in_this_table":True,
        "source_prior_join_transport_calibration_and_microscopic_bath_costs_included":False,
        "quantum_error_bound_or_read_count_globally_optimal":False}


class NoisyCoherentHistoryRecoveryTests(unittest.TestCase):
    def test_exact_confusion_risk_polynomial_includes_invalid_outcomes(self):
        for error in (F(0),F(1,20),F(1,4),F(1,2)):
            matrix=confusion_matrix(error)
            self.assertTrue(all(sum(row)==1 for row in matrix))
            self.assertEqual(1-sum(matrix[r][r] for r in range(6))/6,label_risk(error))
        self.assertEqual(decode_label(6),2); self.assertEqual(decode_label(7),3)

    def test_decoding_is_maximum_likelihood_for_the_three_hard_majority_bits(self):
        for error in (F(1,20),F(1,4)):
            for raw in range(8):
                probabilities=[error**((r^raw).bit_count())*(1-error)**(3-(r^raw).bit_count()) for r in range(6)]
                self.assertEqual(probabilities[decode_label(raw)],max(probabilities))

    def test_actual_native_pointer_tree_preserves_quantum_flag_and_matches_external_map(self):
        rng=np.random.default_rng(147)
        raw=rng.normal(size=(32,3))+1j*rng.normal(size=(32,3)); source=raw @ raw.conj().T; source/=np.trace(source)
        error=(1-ALPHA/2)/2
        expected=channel(tuple(np.kron(k,np.eye(2)) for k in recovery_kraus(error)),source)
        actual=actual_one_read_recovery(source,.5,2)
        np.testing.assert_allclose(actual,expected,atol=4e-15)

    def test_correct_label_component_bounds_arbitrary_external_recovery(self):
        for error in (F(0),F(1,20),F(1,4)):
            ks=recovery_kraus(error)
            np.testing.assert_allclose(sum(k.T @ k for k in ks),np.eye(16),atol=1e-15)
            self.assertLessEqual(max_entangled_metrics(ks)["trace_error"],float(label_risk(error))+2e-15)
        rng=np.random.default_rng(247); vector=rng.normal(size=48)+1j*rng.normal(size=48); vector/=np.linalg.norm(vector)
        source=np.outer(vector,vector.conj()); error=F(1,20)
        output=channel(tuple(np.kron(k,np.eye(3)) for k in recovery_kraus(error)),source)
        self.assertLessEqual(trace_distance(output,source),float(label_risk(error))+1e-14)

    def test_report_conditioning_has_state_independent_probability_and_honest_error(self):
        error=F(1,10); matrix=confusion_matrix(error); risks=conditional_label_risks(error)
        source=np.eye(16)/16
        for g in range(6):
            ks=recovery_kraus(error,g); probability=float(sum(matrix[r][g] for r in range(6))/6)
            np.testing.assert_allclose(sum(k.T @ k for k in ks),probability*np.eye(16),atol=2e-16)
            normalized=tuple(k/np.sqrt(probability) for k in ks)
            self.assertLessEqual(max_entangled_metrics(normalized)["trace_error"],float(risks[g])+2e-15)
            self.assertAlmostEqual(np.trace(channel(ks,source)),probability)

    def test_finite_interval_certificates_beat_the_classical_only_bound(self):
        self.assertGreater(label_risk(majority_error(1)).lo,SCALE//2)
        self.assertLess(label_risk(majority_error(3)).hi,SCALE//2)
        self.assertLess(label_risk(majority_error(15)).hi,SCALE//20)
        self.assertLess(label_risk(majority_error(25)).hi,SCALE//100)
        self.assertTrue(certificate(3)["strictly_beats_classical_only_half_error"])

    def test_combined_independent_quantum_flag_noise_has_a_finite_upper_bound(self):
        error=F(1,20); coherence=F(9,10); phase_error=(1-coherence)/2
        bound=phase_error+(1-phase_error)*label_risk(error)
        combined=tuple(a @ b for a in recovery_kraus(error) for b in phase_kraus((coherence,)))
        self.assertLessEqual(max_entangled_metrics(combined)["trace_error"],float(bound)+1e-14)

    def test_native_resource_ledger_and_guards(self):
        report=certificate(3)
        self.assertEqual((report["original_noisy_reads"],report["fresh_pointer_initializations"]),(9,9))
        self.assertEqual((report["new_yx_after_first_join"],report["new_ry_after_first_join"]),(440,379))
        self.assertEqual(report["new_total_native_gates_after_first_join"],819)
        for value in (-1,F(3,4)):
            with self.assertRaises(ValueError): confusion_matrix(value)
        with self.assertRaises(ValueError): decode_label(8)
        with self.assertRaises(ValueError): certificate(2)


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(NoisyCoherentHistoryRecoveryTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":147,"label_risk_formula":"(8e - 7e^2 + 2e^3)/3",
        "hard_majority_invalid_label_decoding":{"110":"010","111":"011"},
        "uniform_error_includes_arbitrary_unknown_external_correlations":True,
        "combined_independent_flag_noise_upper":"q + (1-q) label_risk, q=(1-lambda)/2",
        "noise_bath_and_source_purification_inaccessible":True,
        "native_readout_fresh_pointers_and_all_records_retained":True,
        "certificates":[certificate(m) for m in (1,3,15,25)],
        "automated_checks":{"run":checks.testsRun,"failures":len(checks.failures),"errors":len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name("noisy_coherent_history_recovery_results.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    summary={key:value for key,value in report.items() if key!="certificates"}
    summary["certificates_summary"]=[{k:v for k,v in row.items() if k in ("reads_per_label_bit","original_noisy_reads","uniform_external_recovery_error_upper","new_total_native_gates_after_first_join")} for row in report["certificates"]]
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=="__main__": main()
