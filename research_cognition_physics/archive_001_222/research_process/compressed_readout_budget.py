"""Round 95: certified noisy readout after compression and explicit resource costs."""

import argparse
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I, SCALE
from coherent_memory_compression import stream_isometry
from complex_whole_real_interfaces import flagged_complex_whole, random_state
from one_bit_real_network import visibility_interval
from persistent_relation_memory import record_tree
from quantum_interface_audit import branch_kraus, rotation_unitary
from repeated_weak_memory import fresh_certificate
from weak_relation_tradeoff import coherence_channel, forget_front_pointer, read_axis


def compressed_record_tree(state, overlaps, count):
    if not isinstance(count,int) or count<1 or count%2!=1:
        raise ValueError("Use a positive odd read count.")
    # Round 94 proves that scratch is a deterministic |0> factor, so its removal
    # is a partial trace, never postselection. Its physical slot can now be read pointer.
    w = stream_isometry(overlaps).reshape(2,2,4,4)[:,0].reshape(8,4)
    w = np.kron(w,np.eye(len(state)//4))
    written = w@state@w.conj().T
    angle = read_axis(2*math.acos(math.prod(overlaps)))
    if count==1:
        leaves={}
        for outcome in (0,1):
            ks = [np.kron(k,np.eye(len(state))) for k in branch_kraus(angle,outcome,1.)]
            leaves[(outcome,)] = sum(k@written@k.conj().T for k in ks)
        return leaves
    turn = np.kron(rotation_unitary(-angle),np.eye(len(state)))
    return record_tree(turn@written@turn.conj().T,count)


def resource_cost(write_count, read_count, compression=True):
    if (not isinstance(write_count,int) or write_count<1 or
        not isinstance(read_count,int) or read_count<1 or read_count%2!=1):
        raise ValueError("Positive write count and positive odd read count required.")
    if not compression and write_count!=1:
        raise ValueError("The direct stronger-write alternative writes exactly once.")
    merge_count = write_count-1 if compression else 0
    protected = read_count>1
    copies = read_count if protected else 0
    slots = 2 if merge_count or protected else 1
    return {"original_system_writes":write_count,"compression_merges":merge_count,
        "pair_rotations":3*write_count+3*merge_count+copies,
        "local_rotations_excluding_settings_internal_to_old_readout":
            3*write_count+3*merge_count+3*copies+int(protected),
        "old_readouts":read_count,"peak_extra_quantum_slots":slots,
        "persistent_accumulator_qubits":1,
        "initial_pure_memory_preparations":slots,
        "read_pointer_resets_after_first_read":read_count-1 if protected else 0,
        "compression_resets":0,"outcome_bits_produced":read_count,
        "classical_counter_bits":math.ceil(math.log2(read_count+1)),
        "adaptive_message_rounds":0,
        "communication_scope":"All quantum gates co-located; exporting final decision costs one bit"}


def certificate(overlaps, count):
    cs = tuple(Fraction(c) for c in overlaps)
    if not cs or any(not 0<=c<=1 for c in cs):
        raise ValueError("Nonempty exact overlaps in [0,1] required.")
    costs = resource_cost(len(cs),count)
    c = math.prod(cs)
    v = (1-I.exact(c*c)).sqrt()
    ideal = (1+v)/2
    success = (1+visibility_interval(count)*v)/2
    gap = ideal-success
    return {"write_overlaps_exact":[str(x) for x in cs],"total_overlap_exact":str(c),
        "worst_old_state_trace_distance_exact":str((1-c)/2),
        "success_lower_exact":str(Fraction(success.lo,SCALE)),
        "success_upper_exact":str(Fraction(success.hi,SCALE)),
        "success_diagnostic":sum(success.floats())/2,
        "ideal_success_lower_exact":str(Fraction(ideal.lo,SCALE)),
        "ideal_success_upper_exact":str(Fraction(ideal.hi,SCALE)),
        "ideal_success_diagnostic":sum(ideal.floats())/2,
        "ideal_gap_upper_exact":str(Fraction(gap.hi,SCALE)),
        "stream_cost":costs,"equivalent_single_write_cost":resource_cost(1,count,False)}


class CompressedReadoutBudgetTests(unittest.TestCase):
    def test_actual_compiler_and_all_noisy_branches_match_success_certificate(self):
        for cs in ((Fraction(4,5),)*3,(Fraction(1,3),Fraction(9,10))):
            for m in (1,3,5):
                success=0.
                for sign in (-1,1):
                    leaves=compressed_record_tree(flagged_complex_whole(sign),tuple(map(float,cs)),m)
                    self.assertAlmostEqual(sum(np.trace(x).real for x in leaves.values()),1,places=13)
                    for history,branch in leaves.items():
                        if (1 if sum(history)>m//2 else -1)==sign:
                            success+=np.trace(branch).real/2
                self.assertAlmostEqual(success,certificate(cs,m)["success_diagnostic"],places=13)

    def test_unconditional_old_channel_is_identical_for_arbitrary_complex_input(self):
        state=random_state(np.random.default_rng(95),8)
        cs=(.8,.7,.9)
        for m in (1,3):
            output=sum(forget_front_pointer(b) for b in compressed_record_tree(state,cs,m).values())
            np.testing.assert_allclose(output,coherence_channel(state,math.prod(cs)),atol=1e-15)

    def test_every_record_keeps_a_definite_sector_old_whole(self):
        for sign in (-1,1):
            state=flagged_complex_whole(sign)
            for branch in compressed_record_tree(state,(.8,)*3,3).values():
                reduced=forget_front_pointer(branch)
                np.testing.assert_allclose(reduced,np.trace(branch).real*state,atol=5e-16)

    def test_three_writes_one_final_read_strictly_beat_three_separate_reads(self):
        baseline=fresh_certificate(Fraction(4,5),Fraction(3,5),3)
        compressed=certificate((Fraction(4,5),)*3,1)
        self.assertEqual(baseline["worst_old_state_trace_distance_exact"],
                         compressed["worst_old_state_trace_distance_exact"])
        self.assertGreater(Fraction(compressed["success_lower_exact"]),
                            Fraction(baseline["independent_direct_read_success_upper_exact"]))

    def test_finite_readouts_approach_the_same_disturbance_bound_from_below(self):
        previous=Fraction(1,2)
        for m in (1,3,5,11):
            row=certificate((Fraction(4,5),)*3,m)
            lower,upper=Fraction(row["success_lower_exact"]),Fraction(row["success_upper_exact"])
            self.assertGreater(lower,previous)
            self.assertLess(upper,Fraction(row["ideal_success_lower_exact"]))
            previous=upper
        c=Fraction(64,125); delta=(1-c)/2
        self.assertEqual(1-c*c,4*delta*(1-delta))

    def test_equivalent_single_write_reduces_compiled_resources(self):
        for n in (2,3,9):
            for m in (1,3):
                streaming=resource_cost(n,m)
                direct=resource_cost(1,m,False)
                self.assertEqual(streaming["pair_rotations"]-direct["pair_rotations"],6*(n-1))
                self.assertEqual(streaming["local_rotations_excluding_settings_internal_to_old_readout"]
                    -direct["local_rotations_excluding_settings_internal_to_old_readout"],6*(n-1))
                self.assertGreaterEqual(streaming["peak_extra_quantum_slots"],direct["peak_extra_quantum_slots"])

    def test_cost_bookkeeping_counts_pointer_resets_and_avoids_unneeded_copy_for_one_read(self):
        self.assertEqual(resource_cost(3,1)["pair_rotations"],15)
        self.assertEqual(resource_cost(3,3)["pair_rotations"],18)
        self.assertEqual(resource_cost(1,1,False)["peak_extra_quantum_slots"],1)
        self.assertEqual(resource_cost(3,3)["read_pointer_resets_after_first_read"],2)
        self.assertEqual(resource_cost(3,3)["compression_resets"],0)

    def test_invalid_budgets_and_trivial_overlap_endpoints(self):
        for m in (0,2,-1):
            with self.assertRaises(ValueError): certificate((Fraction(4,5),),m)
        self.assertAlmostEqual(certificate((Fraction(1),),3)["success_diagnostic"],.5)
        self.assertGreater(certificate((Fraction(0),),3)["success_diagnostic"],.9999)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(CompressedReadoutBudgetTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":95,"success":"(1+gamma_m*sqrt(1-C^2))/2; C=product(c_i)",
        "same_original_whole_not_new_unknown_copies":True,
        "single_write_dominates_this_compiler_if_strength_is_freely_tunable":True,
        "all_finite_protocols_at_equal_resources_optimized":False,
        "exact_gates_and_reset_preparations_assumed":True,
        "separate_three_read_baseline":fresh_certificate(Fraction(4,5),Fraction(3,5),3),
        "compressed_examples":[certificate((Fraction(4,5),)*3,m) for m in (1,3,5)],
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("compressed_readout_budget_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
