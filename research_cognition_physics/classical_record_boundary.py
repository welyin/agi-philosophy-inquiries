"""Round 101: retaining classical records does not restore the erased query."""

import argparse
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I
from complex_whole_real_interfaces import random_state
from deferred_relation_queries import memory_channel,memory_blocks
from incompatible_relation_writes import G,H,calibration_state,written_state,forget_memories,old_channel
from joint_relation_records import (memory_axis,read_memory_tree,validated_count,interval_fields)
from memory_dephasing_threshold import dephase_memory
from one_bit_real_network import visibility_interval
from persistent_relation_memory import record_tree
from quantum_interface_audit import ALPHA,rotation_unitary
from weak_relation_tradeoff import read_axis,trace_distance,forget_front_pointer


def protected_G_records(state,c,count):
    """Copy the old G read axis to reset pointers; keep M in its original frame."""
    validated_count(count)
    turn=np.kron(rotation_unitary(-read_axis(2*math.acos(c))),np.eye(len(state)//2))
    return {history:turn.conj().T@branch@turn
            for history,branch in record_tree(turn@state@turn.conj().T,count).items()}


def analytic_protected_branch(state,c,history):
    a=np.kron(memory_axis(c),np.eye(len(state)//2))
    result=np.zeros_like(state,dtype=complex)
    for z in (-1,1):
        p=(np.eye(len(state))+z*a)/2
        likelihood=math.prod((1+(2*r-1)*z*ALPHA)/2 for r in history)
        result+=likelihood*p@state@p
    return result


def block_record_distance(minus,plus,keep="MN"):
    if minus.keys()!=plus.keys(): raise ValueError("Record alphabets must agree.")
    total=0.
    for history in minus:
        a,b=minus[history],plus[history]
        if keep=="N": a,b=forget_front_pointer(a),forget_front_pointer(b)
        elif keep=="none": a,b=np.array([[np.trace(a)]]),np.array([[np.trace(b)]])
        elif keep!="MN": raise ValueError("Keep MN, N, or none.")
        total+=trace_distance(a,b)
    return total


def boundary_certificate(c,d,count):
    c,d=Fraction(c),Fraction(d)
    if not 0<=c<=1 or not 0<=d<=1: raise ValueError("Use overlaps in [0,1].")
    gamma=visibility_interval(count)
    v,w=(1-I.exact(c*c)).sqrt(),(1-I.exact(d*d)).sqrt()
    return {"c_exact":str(c),"d_exact":str(d),"G_classical_readouts":count,
        "G_record_alone_trace_distance":interval_fields(gamma*v),
        "G_record_plus_MN_trace_distance":interval_fields(v),
        "H_record_plus_MN_trace_distance":interval_fields(I.exact(c)*w),
        "H_record_plus_N_trace_distance":interval_fields(I.exact(c)*w),
        "H_record_alone_trace_distance_exact":"0",
        "copy_pair_rotations":count,"copy_and_frame_local_rotations":3*count+2,
        "additional_pointer_slots":1,"pointer_resets_after_first":count-1,
        "all_raw_record_bits_kept":count,
        "record_bits_can_be_processed_jointly_with_remaining_memories":True}


class ClassicalRecordBoundaryTests(unittest.TestCase):
    def test_actual_protected_instrument_matches_every_projected_branch(self):
        state=random_state(np.random.default_rng(101),4)
        for count in (1,3,5):
            leaves=protected_G_records(state,.8,count)
            for history,branch in leaves.items():
                np.testing.assert_allclose(branch,analytic_protected_branch(state,.8,history),atol=8e-16)

    def test_even_one_protected_read_completely_dephases_the_G_axis(self):
        state=random_state(np.random.default_rng(201),8)
        for count in (1,3):
            np.testing.assert_allclose(sum(protected_G_records(state,.8,count).values()),
                                       dephase_memory(state,.8,0),atol=1e-15)

    def test_retaining_all_records_and_both_memories_does_not_restore_H(self):
        for c,d in ((0,0),(.8,.8),(.3,.7),(1,.4)):
            for count in (1,3):
                minus,plus=[protected_G_records(memory_channel(calibration_state(H,s),c,d),c,count)
                            for s in (-1,1)]
                for keep in ("MN","N"):
                    self.assertAlmostEqual(block_record_distance(minus,plus,keep),
                                           c*math.sqrt(1-d*d),places=13)
                self.assertAlmostEqual(block_record_distance(minus,plus,"none"),0,places=14)

    def test_every_H_calibration_branch_factors_into_label_independent_M_and_N(self):
        c,d=.8,.6
        ad,_,dd=memory_blocks(d)
        for s in (-1,1):
            leaves=protected_G_records(memory_channel(calibration_state(H,s),c,d),c,3)
            n_state=(ad+s*c*dd)/2
            for history,branch in leaves.items():
                q=np.zeros((2,2),dtype=complex)
                for z in (-1,1):
                    probability=math.prod((1+(2*r-1)*z*ALPHA)/2 for r in history)
                    q+=probability*(np.eye(2)+z*memory_axis(c))/4
                np.testing.assert_allclose(branch,np.kron(q,n_state),atol=1e-15)

    def test_the_retained_G_memory_is_more_informative_than_the_noisy_classical_record(self):
        c,d=.8,.8
        for count in (1,3):
            minus,plus=[protected_G_records(memory_channel(calibration_state(G,s),c,d),c,count)
                        for s in (-1,1)]
            gamma=sum(visibility_interval(count).floats())/2
            self.assertAlmostEqual(block_record_distance(minus,plus,"MN"),.6,places=13)
            self.assertAlmostEqual(block_record_distance(minus,plus,"none"),gamma*.6,places=13)
            self.assertAlmostEqual(block_record_distance(minus,plus,"N"),gamma*.6,places=13)

    def test_direct_and_protected_reads_have_same_first_record_but_different_remaining_G_information(self):
        c,d=.8,.8
        direct=[read_memory_tree(memory_channel(calibration_state(G,s),c,d),0,c,1) for s in (-1,1)]
        protected=[protected_G_records(memory_channel(calibration_state(G,s),c,d),c,1) for s in (-1,1)]
        for left,right in zip(direct,protected):
            for history in left:
                self.assertAlmostEqual(np.trace(left[history]).real,np.trace(right[history]).real,places=14)
        self.assertAlmostEqual(block_record_distance(*direct),ALPHA*.6,places=14)
        self.assertAlmostEqual(block_record_distance(*protected),.6,places=14)

    def test_old_system_and_complex_reference_unconditional_state_stays_the_same(self):
        state=random_state(np.random.default_rng(301),8)
        leaves=protected_G_records(written_state(state,.8,.6),.8,3)
        np.testing.assert_allclose(sum(forget_memories(b) for b in leaves.values()),
                                   old_channel(state,.8,.6),atol=2e-15)

    def test_finite_certificates_do_not_turn_more_classical_repetitions_into_new_H_information(self):
        earlier=Fraction(0)
        for count in (1,3,5):
            row=boundary_certificate(Fraction(4,5),Fraction(4,5),count)
            self.assertAlmostEqual(row["H_record_plus_MN_trace_distance"]["diagnostic"],.48)
            self.assertGreater(Fraction(row["G_record_alone_trace_distance"]["lower_exact"]),earlier)
            earlier=Fraction(row["G_record_alone_trace_distance"]["upper_exact"])
            self.assertEqual(row["pointer_resets_after_first"],count-1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(ClassicalRecordBoundaryTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":101,"instrument":"B_history(omega)=sum_z likelihood(history|z) P_z omega P_z",
        "one_protected_record_already_erases_G_axis_offdiagonal_blocks":True,
        "ideal_H_distinguishability_with_all_classical_records_and_MN":"c*sqrt(1-d^2)",
        "ideal_G_distinguishability_with_all_classical_records_and_MN":"sqrt(1-c^2)",
        "direct_noisy_read_and_protected_read_have_same_first_effect_but_different_retained_state":True,
        "examples":[boundary_certificate(Fraction(4,5),Fraction(4,5),m) for m in (1,3,5)],
        "environment_quantum_access_included":False,
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("classical_record_boundary_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
