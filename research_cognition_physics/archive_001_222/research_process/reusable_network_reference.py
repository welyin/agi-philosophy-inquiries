"""Round 122: retain the distributed reference through actual old instruments.

Exact orientation sectors are preserved. A mixed shared reference correlates
different trials; unchanged average marginal is not an independent fresh reset.
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

from bell_network_statistics import SETTINGS, endpoint_observables, bob_effects, network_table
from finite_network_separation import compiled_axis_gate
from independent_source_alignment import pair_reference, keep_systems
from joint_relation_records import interval_fields
from one_bit_real_network import target_assemblage, records_for_reference, visibility_interval
from operational_effect_closure import majority_certificate, reusable_pointer_branch
from quantum_interface_audit import IDENTITY, PAULI_Y, PAULI_Z
from reference_information_cost import target_intervals


def direction_vectors(x,z):
    axes=np.eye(3)
    i,j,sign=SETTINGS[z]
    return axes[x],(axes[i]+sign*axes[j])/math.sqrt(2)


@lru_cache(maxsize=90)
def local_record_kraus(direction,sign,count):
    if sign not in (-1,1): raise ValueError("Binary sign required.")
    visibility=majority_certificate(count)["effective_visibility_diagnostic"]
    gate=compiled_axis_gate(direction)
    result=[]
    for eigenvalue in (-1,1):
        projection=(IDENTITY+eigenvalue*PAULI_Z)/2
        weight=(1+sign*visibility*eigenvalue)/2
        result.append(np.kron(IDENTITY,math.sqrt(weight)*projection)@gate)
    return tuple(result)


@lru_cache(maxsize=4)
def assemblage(count):
    visibility=majority_certificate(count)["effective_visibility_diagnostic"]
    return target_assemblage(visibility)


def retained_branch(reference,x,z,a,b,c,count=5):
    """References R_A,R_C retained; target outputs are only traced for this view."""
    left,right=direction_vectors(x,z)
    # Initial order R_A,R_C,T_A,T_C -> R_A,T_A,R_C,T_C, a bookkeeping permutation.
    initial=keep_systems(np.kron(reference,assemblage(count)[b]),(0,2,1,3),4)
    output=np.zeros((16,16),dtype=complex)
    for ka,kc in product(local_record_kraus(tuple(left),a,count),
                         local_record_kraus(tuple(right),c,count)):
        joint=np.kron(ka,kc)
        output+=joint@initial@joint.conj().T
    return keep_systems(output,(0,2),4)


@lru_cache(maxsize=4)
def sector_tables(count=5):
    visibility=majority_certificate(count)["effective_visibility_diagnostic"]
    return {sign:records_for_reference(pair_reference(sign),visibility,visibility,visibility)
            for sign in (-1,1)}


def predicted_branch(g,key,count=5):
    tables=sector_tables(count)
    plus=(1+g)*tables[1][key]/2
    minus=(1-g)*tables[-1][key]/2
    return plus*pair_reference(1)+minus*pair_reference(-1)


def posterior_correlation(g,key,count=5):
    tables=sector_tables(count)
    plus=(1+g)*tables[1][key]/2
    minus=(1-g)*tables[-1][key]/2
    if plus+minus<=0: raise ValueError("Impossible branch.")
    return (plus-minus)/(plus+minus)


def resources(trials):
    if not isinstance(trials,int) or trials<1: raise ValueError("Positive trial count required.")
    return {"task_trials":trials,"fresh_target_Bell_pairs":2*trials,
            "reference_link_pairs_once":2,"retained_outer_reference_rebits":2,
            "alignment_reads_once":5,"target_reads":20*trials,
            "total_old_reads":20*trials+5,
            "premeasurement_alignment_message_bits":1,
            "pair_rotations_after_sources_prepared_upper":23*trials+7,
            "fresh_reference_policy_total_old_reads":25*trials,
            "fresh_reference_policy_alignment_message_bits":trials,
            "fresh_reference_policy_reference_link_pairs":2*trials,
            "preparation_storage_maintenance_and_final_record_exchange_are_separate":True}


def reuse_certificate():
    g=visibility_interval(5)
    error=(1-g)/2
    _,_,matched_g,_=target_intervals()
    return {"initial_correlation":interval_fields(g),
            "whole_record_distance_to_perfect_shared_orientation_upper":interval_fields(error),
            "matched_reused_complex_private_resources_correlation":interval_fields(matched_g),
            "matched_reuse_whole_record_distance_upper":interval_fields((1-matched_g)/2),
            "bound_independent_of_finite_trial_count":True,
            "conditional_reference_can_change":True,
            "different_trials_are_not_iid_when_reference_is_imperfect":True,
            "ideal_old_control_and_no_storage_noise":True,
            "resources":[resources(n) for n in (1,10,1000)],
            "does_not_preserve_independence_of_entire_cross_trial_lab_resources":True}


class ReusableNetworkReferenceTests(unittest.TestCase):
    def test_coarse_pointer_instrument_is_the_actual_old_record_sum(self):
        rho=(IDENTITY+.31*PAULI_Y+.27*PAULI_Z)/2
        for count in (1,3,5):
            visibility=majority_certificate(count)["effective_visibility_diagnostic"]
            for sign in (-1,1):
                actual=np.zeros((2,2),dtype=complex)
                for record in product((0,1),repeat=count):
                    if (1 if sum(record)>count//2 else -1)!=sign: continue
                    state=rho.copy()
                    for outcome in record: state=reusable_pointer_branch(state,outcome)
                    actual+=state
                expected=sum((1+sign*visibility*t)/2*((IDENTITY+t*PAULI_Z)/2)@rho@((IDENTITY+t*PAULI_Z)/2)
                             for t in (-1,1))
                np.testing.assert_allclose(actual,expected,atol=3e-15,rtol=0)

    def test_each_actual_local_kraus_preserves_the_reference_orientation(self):
        generator=np.kron(PAULI_Y,IDENTITY)
        for x,z,sign in product(range(3),range(6),(-1,1)):
            for direction in direction_vectors(x,z):
                for k in local_record_kraus(tuple(direction),sign,5):
                    np.testing.assert_allclose(k.imag,0,atol=0)
                    np.testing.assert_allclose(k@generator,generator@k,atol=5e-16,rtol=0)

    def test_all_288_retained_branches_equal_the_orientation_filter(self):
        g=.61
        for key in sector_tables(5)[1]:
            actual=retained_branch(pair_reference(g),*key)
            np.testing.assert_allclose(actual,predicted_branch(g,key),atol=4e-16,rtol=0)

    def test_average_reference_is_unchanged_for_every_setting(self):
        rho=pair_reference(.73)
        for x,z in product(range(3),range(6)):
            output=sum(retained_branch(rho,x,z,a,b,c) for a,b,c in product((-1,1),range(4),(-1,1)))
            np.testing.assert_allclose(output,rho,atol=3e-15,rtol=0)

    def test_conditioned_reference_is_not_generally_the_initial_marginal(self):
        g=.4
        key=(1,4,1,0,1)
        actual=retained_branch(pair_reference(g),*key)
        probability=np.trace(actual).real
        updated=posterior_correlation(g,key)
        self.assertGreater(abs(updated-g),.2)
        np.testing.assert_allclose(actual/probability,pair_reference(updated),atol=3e-15,rtol=0)

    def test_two_complete_trials_match_a_shared_sector_not_fresh_resets(self):
        g=.61
        first_keys=[(1,4,a,b,c) for a,b,c in product((-1,1),range(4),(-1,1))]
        second_keys=[(1,0,a,b,c) for a,b,c in product((-1,1),range(4),(-1,1))]
        tables=sector_tables()
        distance_from_fresh=0.
        for first in first_keys:
            mid=retained_branch(pair_reference(g),*first)
            for second in second_keys:
                actual=retained_branch(mid,*second)
                expected=((1+g)/2*tables[1][first]*tables[1][second]*pair_reference(1)
                          +(1-g)/2*tables[-1][first]*tables[-1][second]*pair_reference(-1))
                np.testing.assert_allclose(actual,expected,atol=2e-16,rtol=0)
                p1=np.trace(predicted_branch(g,first)).real
                p2=np.trace(predicted_branch(g,second)).real
                distance_from_fresh+=abs(np.trace(actual).real-p1*p2)/2
        self.assertGreater(distance_from_fresh,.1)

    def test_perfect_shared_orientation_returns_the_same_visible_reference_in_each_branch(self):
        for key in sector_tables()[1]:
            actual=retained_branch(pair_reference(1),*key)
            np.testing.assert_allclose(actual,sector_tables()[1][key]*pair_reference(1),atol=4e-16,rtol=0)

    def test_setup_resources_are_charged_once_but_task_resources_continue(self):
        single=resources(1)
        self.assertEqual(single["total_old_reads"],25)
        self.assertEqual(single["pair_rotations_after_sources_prepared_upper"],30)
        batch=resources(1000)
        self.assertEqual(batch["total_old_reads"],20005)
        self.assertEqual(batch["premeasurement_alignment_message_bits"],1)
        cert=reuse_certificate()
        self.assertLess(F(cert["whole_record_distance_to_perfect_shared_orientation_upper"]["upper_exact"]),F(14,10**7))
        with self.assertRaises(ValueError): resources(0)

    def test_reusing_both_complex_private_references_matches_the_same_shared_real_sector(self):
        bias=.7
        vis=majority_certificate(5)["effective_visibility_diagnostic"]
        tables={}
        for left,right in product((-1,1),repeat=2):
            alice=endpoint_observables(vis,left)[0]
            charlie=endpoint_observables(vis,right)[1]
            tables[left,right]=network_table(alice,charlie,bob_effects(vis))
        pairs=[((1,4,1,0,1),(1,0,-1,2,1)),((0,2,-1,3,-1),(1,5,1,1,-1))]
        sectors=sector_tables()
        for first,second in pairs:
            actual=sum((1+left*bias)*(1+right*bias)/4*tables[left,right][first]*tables[left,right][second]
                       for left,right in product((-1,1),repeat=2))
            expected=sum((1+sign*bias*bias)/2*sectors[sign][first]*sectors[sign][second] for sign in (-1,1))
            self.assertAlmostEqual(actual,expected,places=15)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ReusableNetworkReferenceTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":122,"reuse_certificate":reuse_certificate(),
            "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0},
            "quantum_theory_derived_from_cognition":False}
    if args.write_results:
        Path(__file__).with_name("reusable_network_reference_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
