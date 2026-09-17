"""Round 140: finite original noisy reads of classical joining history.

Fresh copied pointers protect each history bit. An explicit variable-length
policy accepts every outcome. The uniform error bound follows from the exact
label confusion risk, not from a claim of optimal recovery or sample count.
"""

import argparse
from fractions import Fraction as F
from itertools import product
import json
import math
from pathlib import Path
import unittest

import numpy as np

from approximate_subject_joining import apply_kraus, trace_distance
from certified_intervals import Interval as I, SCALE, sin_interval
from classical_joining_history import (
    WEIGHTS, history_kraus, physical_record_class, recovery_isometry,
)
from common_orientation_structure import encode_state
from compiled_joining_frontier import compiled_sharing
from encoded_composition_audit import independent_encoding, product_state
from flagged_subject_joining import random_state
from independent_source_alignment import keep_systems
from joining_history_advantage import balanced_pair, common_ab_marginal
from joint_relation_records import interval_fields
from multisubject_joining_optimum import multisubject_kraus, reorder_qubits
from operational_effect_closure import compiled_copy_gate
from quantum_interface_audit import ALPHA, branch_kraus


def fine_records():
    result=[]
    for bits in product((0,1),repeat=4):
        try: r=physical_record_class(bits)
        except ValueError: continue
        result.append((bits,r,F(1,4) if r==0 else F(1,12)))
    return tuple(result)


def decode_majorities(bits):
    h,e1,e2,c=bits
    if any(bit not in (0,1) for bit in bits): raise ValueError("Binary estimates required.")
    if h==0 or (e1,e2)==(1,1): return 0
    return physical_record_class(bits)


def confusion_matrix(bit_error):
    e=bit_error
    result=[[0 for _ in range(7)] for _ in range(7)]
    for true,r,weight in fine_records():
        conditional_weight=weight/WEIGHTS[r]
        for guessed in product((0,1),repeat=4):
            errors=sum(a!=b for a,b in zip(true,guessed))
            result[r][decode_majorities(guessed)]+=conditional_weight*e**errors*(1-e)**(4-errors)
    return result


def label_risk(bit_error):
    e=bit_error
    return (5*e-6*e**2+3*e**3-e**4)/2


def majority_error(count,contrast=F(1,2)):
    if not isinstance(count,int) or isinstance(count,bool) or count<1 or count%2==0:
        raise ValueError("A positive odd number of fresh pointer reads is required.")
    eta=F(contrast)
    if not 0<eta<=1: raise ValueError("Use 0 < contrast <= 1.")
    visibility=4*sin_interval(F(1,4))*I.exact(eta)
    correct=(1+visibility)/2
    return sum((math.comb(count,k)*correct**k*(1-correct)**(count-k)
                for k in range((count+1)//2)),I.exact(0))


def noisy_recovered_input(source,bit_error):
    matrix=confusion_matrix(bit_error)
    result=np.zeros((16,16),dtype=source.dtype)
    for r,k in enumerate(history_kraus()):
        branch=k@source@k.T
        for guessed in range(7):
            restore=recovery_isometry(guessed)
            result+=float(matrix[r][guessed])*(restore@branch@restore.T)
    return result


def noisy_three_party_output(states,bit_error):
    restored=noisy_recovered_input(independent_encoding(states[:2]),bit_error)
    source=reorder_qubits(np.kron(restored,encode_state(states[2])),(0,1,4,2,3,5))
    return apply_kraus(multisubject_kraus(3),source)


def protected_read_kraus(guessed_bit,contrast=.5):
    if guessed_bit not in (0,1): raise ValueError("Use a binary result.")
    # Old outcome 1 is positive Z and therefore estimates computational bit 0.
    result=[]
    copy=compiled_copy_gate().real  # Pointer first, retained history bit second.
    for k in branch_kraus(0.,1-guessed_bit,contrast):
        joint=np.kron(k.real,np.eye(2))@copy
        for final_pointer in (0,1):
            result.append(joint[2*final_pointer:2*final_pointer+2,:2])
    return tuple(result)


def physical_read_branch(whole,site,guessed_bit,contrast=.5):
    systems=(len(whole)).bit_length()-1
    bits=(np.arange(len(whole)) >> (systems-1-site)) & 1
    result=np.zeros_like(whole)
    for k in protected_read_kraus(guessed_bit,contrast):
        if np.max(abs(k-np.diag(np.diag(k))))>1e-15:
            raise ArithmeticError("Expected a diagonal nondemolition history instrument.")
        values=np.diag(k)[bits]
        result+=values[:,None]*whole*values[None,:]
    return result


def full_one_read_recovery(source,contrast=.5):
    v,_=compiled_sharing(F(1,2)); whole=v@source@v.T
    result=np.zeros((16,16))
    for h in (0,1):
        first=physical_read_branch(whole,0,h,contrast)
        if h==0:
            active=keep_systems(first,(1,2,3),7)
            restore=recovery_isometry(0)
            result+=restore@active@restore.T
        else:
            for e1,e2,c in product((0,1),repeat=3):
                branch=first
                for site,value in ((4,e1),(5,e2),(6,c)):
                    branch=physical_read_branch(branch,site,value,contrast)
                active=keep_systems(branch,(1,2,3),7)
                restore=recovery_isometry(decode_majorities((h,e1,e2,c)))
                result+=restore@active@restore.T
    return result


def certificate(count=15,contrast=F(1,2)):
    error=majority_error(count,contrast); risk=label_risk(error)
    upper=I.exact(F(1,4))+risk
    margin=I.exact(F(11,36))-upper
    return {
        "reads_per_queried_history_bit":count,"contrast_exact":str(F(contrast)),
        "majority_bit_error":interval_fields(error),"canonical_label_error":interval_fields(risk),
        "uniform_final_joint_trace_error_upper":interval_fields(upper),
        "uniform_old_ab_deviation_from_previous_interface_upper":interval_fields(risk),
        "margin_below_no_history_optimum":interval_fields(margin),
        "strict_improvement_certified":margin.lo>0,
        "expected_original_reads_exact":str(F(5*count,2)),"maximum_original_reads":4*count,
        "expected_fresh_pointer_initializations_exact":str(F(5*count,2)),
        "maximum_fresh_pointer_initializations":4*count,
        "additional_fresh_rebits_for_recovery_and_batch_join":3,
        "maximum_new_pure_initializations_after_first_join":4*count+3,
        "maximum_yx_gates_after_first_join":4*count+195,
        "maximum_ry_gates_after_first_join":12*count+59,
        "decoded_message_bits_if_sent_to_a_separate_controller":3,
        "maximum_raw_outcome_bits_retained":4*count,
        "source_preparation_prior_first_join_storage_transport_and_readout_bath_costs_included":False,
        "all_outcomes_accepted":True,"old_ab_preserved_exactly_at_finite_count":False,
        "bound_claimed_tight_or_read_count_globally_minimal":False,
    }


class NoisyClassicalHistoryTests(unittest.TestCase):
    def test_actual_copied_pointer_instrument_matches_noisy_classical_record(self):
        for contrast in (.5,1.):
            v=ALPHA*contrast
            for guessed in (0,1):
                ks=protected_read_kraus(guessed,contrast)
                for i,j in product(range(2),repeat=2):
                    unit=np.zeros((2,2)); unit[i,j]=1
                    actual=apply_kraus(ks,unit)
                    expected=np.zeros((2,2))
                    if i==j: expected[i,i]=(1+v*(1 if i==guessed else -1))/2
                    np.testing.assert_allclose(actual,expected,atol=3e-16)

    def test_repeated_fresh_pointers_are_independent_only_conditioned_on_the_kept_bit(self):
        correct=(1+ALPHA/2)/2
        for true in (0,1):
            for record in product((0,1),repeat=3):
                original=np.diag([1.,0.]) if true==0 else np.diag([0.,1.])
                branch=original
                expected=1.
                for value in record:
                    branch=apply_kraus(protected_read_kraus(value),branch)
                    expected*=correct if value==true else 1-correct
                np.testing.assert_allclose(branch,expected*original,atol=4e-16)

    def test_exact_confusion_risk_matches_all_true_and_reported_codes(self):
        for e in (F(0),F(1,100),F(1,10),F(1,2)):
            table=confusion_matrix(e)
            self.assertTrue(all(sum(row)==1 for row in table))
            risk=1-sum(WEIGHTS[r]*table[r][r] for r in range(7))
            self.assertEqual(risk,label_risk(e))
        self.assertEqual(sum(weight for _,_,weight in fine_records()),F(1))

    def test_actual_seven_rebit_history_tree_matches_confusion_channel(self):
        rng=np.random.default_rng(140)
        source=independent_encoding([random_state(rng,2) for _ in range(2)])
        actual=full_one_read_recovery(source)
        expected=noisy_recovered_input(source,(1-ALPHA/2)/2)
        np.testing.assert_allclose(actual,expected,atol=2e-15)

    def test_state_and_three_party_bounds_hold_without_postselection(self):
        rng=np.random.default_rng(240)
        for e in (F(0),F(1,100),F(1,10)):
            risk=float(label_risk(e))
            for _ in range(4):
                states=[random_state(rng,2) for _ in range(3)]
                original=independent_encoding(states[:2])
                restored=noisy_recovered_input(original,e)
                self.assertAlmostEqual(np.trace(restored),1.)
                self.assertLessEqual(trace_distance(restored,original),risk+2e-15)
                output=noisy_three_party_output(states,e)
                ideal=encode_state(product_state(states))
                self.assertLessEqual(trace_distance(output,ideal),.25+risk+2e-15)
                old=encode_state(balanced_pair(*states[:2]))
                self.assertLessEqual(trace_distance(common_ab_marginal(output),old),risk+2e-15)

    def test_certified_threshold_for_this_bound_and_policy_is_fifteen(self):
        for m in range(1,15,2): self.assertFalse(certificate(m)["strict_improvement_certified"])
        row=certificate(15)
        self.assertTrue(row["strict_improvement_certified"])
        upper=F(row["uniform_final_joint_trace_error_upper"]["upper_exact"])
        self.assertLess(upper,F(295031,1000000))
        self.assertGreater(upper,F(295030,1000000))
        self.assertTrue(certificate(1,F(1))["strict_improvement_certified"])

    def test_expected_read_cost_and_all_historical_resources_are_explicit(self):
        e=F(1,10); expected=F(0)
        for bits,_,weight in fine_records():
            h=bits[0]
            for guessed in (0,1):
                probability=1-e if guessed==h else e
                expected+=weight*probability*(1+3*guessed)
        self.assertEqual(expected,F(5,2))
        row=certificate(15)
        self.assertEqual(row["expected_original_reads_exact"],"75/2")
        self.assertEqual(row["maximum_original_reads"],60)
        self.assertEqual(row["maximum_new_pure_initializations_after_first_join"],63)
        self.assertEqual(row["maximum_yx_gates_after_first_join"],255)
        self.assertEqual(row["maximum_ry_gates_after_first_join"],239)

    def test_invalid_read_budgets_and_contrasts_are_rejected(self):
        for m in (0,2,-1,2.5,True):
            with self.assertRaises(ValueError): majority_error(m)
        for eta in (F(0),F(-1),F(2)):
            with self.assertRaises(ValueError): majority_error(1,eta)
        with self.assertRaises(ValueError): decode_majorities((0,0,0,2))


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(NoisyClassicalHistoryTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":140,"history_policy":"Read h first; if estimated h=1 read E1,E2,selector; invalid E1E2=11 returns none",
        "all_reads_use_fresh_copied_pointers_and_original_noisy_z":True,
        "ideal_projective_history_readout_assumed":False,
        "history_classicalization_requires_arbitrary_purification_preservation":False,
        "risk_polynomial_exact":"(5e-6e^2+3e^3-e^4)/2",
        "uniform_joint_error_upper_formula":"1/4 + label_risk(e)",
        "physical_history_registers_and_readout_environments_retained":True,
        "bounds_apply_after_averaging_all_records_not_each_conditioned_record":True,
        "rows":[certificate(m,eta) for eta,m in ((F(1,2),1),(F(1,2),13),(F(1,2),15),(F(1,2),25),(F(1),1))],
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("noisy_classical_history_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
