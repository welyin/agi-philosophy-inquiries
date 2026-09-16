"""Round 92: rereading one weak memory versus writing fresh memories."""

import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I, SCALE
from complex_whole_real_interfaces import flagged_complex_whole, random_state
from operational_effect_closure import majority_certificate
from persistent_relation_memory import record_tree
from quantum_interface_audit import ALPHA, branch_kraus, rotation_unitary
from weak_relation_tradeoff import (write_state, forget_front_pointer, coherence_channel,
    memory_state, memory_vectors, read_axis, trace_distance, tradeoff_certificate)


def same_memory_record_probability(sign,angle,history):
    v=math.sin(angle/2)
    return sum((1+sign*z*v)/2*np.prod([(1+(2*r-1)*z*ALPHA)/2 for r in history])
               for z in (-1,1))


def fresh_write_read_branch(state,angle,outcome):
    written=write_state(state,angle)
    output=np.zeros_like(written)
    for k in branch_kraus(read_axis(angle),outcome,1.):
        full=np.kron(k,np.eye(len(state)))
        output+=full@written@full.conj().T
    return forget_front_pointer(output)


def fresh_record_tree(state,angle,count):
    leaves={():state}
    for _ in range(count):
        leaves={history+(r,):fresh_write_read_branch(branch,angle,r)
                for history,branch in leaves.items() for r in (0,1)}
    return leaves


def fresh_certificate(c,v,count):
    c,v=Fraction(c),Fraction(v)
    if c*c+v*v!=1 or not 0 <= c <= 1 or not 0 < v <= 1:
        raise ValueError("Use nonnegative exact amplitudes with c^2+v^2=1 and v>0.")
    row=majority_certificate(count,v)
    ideal=(1+(1-I.exact(c**(2*count))).sqrt())/2
    return {"write_count":count,"old_readouts":count,
            "reset_memory_preparations":count,"write_pair_interactions":3*count,
            "memory_overlap_after_writes_exact":str(c**count),
            "worst_old_state_trace_distance_exact":str((1-c**count)/2),
            "independent_direct_read_success_lower_exact":str(1-Fraction(row["error_upper_exact"])),
            "independent_direct_read_success_upper_exact":str(1-Fraction(row["error_lower_exact"])),
            "independent_direct_read_success_diagnostic":1-float(Fraction(row["error_upper_exact"])),
            "ideal_joint_memory_success_upper_exact":str(Fraction(ideal.hi,SCALE)),
            "ideal_joint_memory_success_diagnostic":sum(ideal.floats())/2}


class RepeatedWeakMemoryTests(unittest.TestCase):
    def test_same_memory_full_records_match_a_latent_bit_mixture(self):
        for angle in (.3,2*math.atan2(3,4),2.7):
            turn=rotation_unitary(-read_axis(angle))
            for sign in (-1,1):
                state=turn@memory_state(sign,angle)@turn.conj().T
                leaves=record_tree(state,3)
                for history,branch in leaves.items():
                    self.assertAlmostEqual(np.trace(branch).real,
                        same_memory_record_probability(sign,angle,history),places=14)

    def test_repeated_same_memory_readouts_are_not_independent_given_the_original_label(self):
        angle=2*math.atan2(3,4)
        v=3/5
        for sign in (-1,1):
            rows={h:same_memory_record_probability(sign,angle,h) for h in product((0,1),repeat=2)}
            mean=sum((2*h[0]-1)*p for h,p in rows.items())
            cross=sum((2*h[0]-1)*(2*h[1]-1)*p for h,p in rows.items())
            self.assertAlmostEqual(mean,sign*v*ALPHA)
            self.assertAlmostEqual(cross,ALPHA**2)
            self.assertAlmostEqual(cross-mean*mean,ALPHA**2*(1-v*v))
            self.assertGreater(abs(rows[(1,1)]-((1+sign*v*ALPHA)/2)**2),.1)

    def test_same_memory_many_reads_approach_eighty_percent_not_certainty(self):
        c,v=Fraction(4,5),Fraction(3,5)
        previous=Fraction(1,2)
        for count in (1,3,5,11,31):
            row=tradeoff_certificate(c,v,count)
            self.assertLess(Fraction(row["success_upper_exact"]),Fraction(4,5))
            self.assertGreater(Fraction(row["success_lower_exact"]),previous)
            previous=Fraction(row["success_lower_exact"])

    def test_fresh_writes_on_the_same_original_whole_give_independent_label_records(self):
        angle=2*math.atan2(3,4)
        for sign in (-1,1):
            whole=flagged_complex_whole(sign)
            leaves=fresh_record_tree(whole,angle,3)
            for history,branch in leaves.items():
                probability=np.prod([(1+(2*r-1)*sign*ALPHA*.6)/2 for r in history])
                np.testing.assert_allclose(branch,probability*whole,atol=5e-16)
            self.assertAlmostEqual(sum(np.trace(b).real for b in leaves.values()),1,places=14)

    def test_fresh_recording_accumulates_old_coherence_loss(self):
        angle=2*math.atan2(3,4)
        state=random_state(np.random.default_rng(92),8)
        output=sum(fresh_record_tree(state,angle,3).values())
        np.testing.assert_allclose(output,coherence_channel(state,(4/5)**3),atol=2e-16)
        witness=np.diag([1.,0,0,0])
        self.assertAlmostEqual(trace_distance(witness,coherence_channel(witness,.8**3)),.244)

    def test_joint_memory_overlap_sets_an_independent_upper_bound(self):
        for angle in (.4,2*math.atan2(3,4),2.5):
            for count in (1,2,3,4):
                vectors=[]
                for start in memory_vectors(angle):
                    vector=np.array([1.])
                    for _ in range(count): vector=np.kron(vector,start)
                    vectors.append(vector)
                states=[np.outer(v,v) for v in vectors]
                expected=math.sqrt(1-math.cos(angle/2)**(2*count))
                self.assertAlmostEqual(trace_distance(*states),expected,places=13)

    def test_finite_fresh_certificate_matches_every_actual_record_and_exceeds_single_memory_ceiling(self):
        angle=2*math.atan2(3,4)
        rows=[fresh_record_tree(flagged_complex_whole(s),angle,3) for s in (-1,1)]
        success=sum(max(np.trace(rows[0][h]).real,np.trace(rows[1][h]).real) for h in rows[0])/2
        certificate=fresh_certificate(Fraction(4,5),Fraction(3,5),3)
        self.assertAlmostEqual(success,certificate["independent_direct_read_success_diagnostic"],places=14)
        self.assertGreater(Fraction(certificate["independent_direct_read_success_lower_exact"]),Fraction(4,5))
        self.assertEqual(Fraction(certificate["worst_old_state_trace_distance_exact"]),Fraction(61,250))
        self.assertLess(success,certificate["ideal_joint_memory_success_diagnostic"])

    def test_infinitely_small_writes_with_nonzero_total_strength_leave_finite_disturbance(self):
        for strength in (.3,1.,2.):
            limit=math.exp(-strength/2)
            previous=1.
            for count in (101,501,2001):
                overlap=math.cos(math.sqrt(strength/count))**count
                self.assertLess(overlap,limit)
                error=abs(overlap-limit)
                self.assertLess(error,previous)
                previous=error
            self.assertGreater((1-limit)/2,0)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(RepeatedWeakMemoryTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":92,
        "same_memory_records_are_conditionally_independent_given_original_label":False,
        "same_memory_readout_covariance":"alpha^2 c^2",
        "any_memory_only_protocol_success_ceiling":"(1+sqrt(1-c^2))/2",
        "same_memory_rereads_change_its_quantum_state":True,
        "fresh_writes_require_reinteraction_with_original_whole":True,
        "fresh_writes_require_new_unknown_input_copies":False,
        "fresh_overlap":"c^n","fresh_worst_old_state_distance":"(1-c^n)/2",
        "fresh_memory_ideal_joint_success":"(1+sqrt(1-c^(2n)))/2",
        "ideal_joint_memory_measurement_compiled_here":False,
        "finite_c_four_fifths_comparison":[
            {"same_memory":tradeoff_certificate(Fraction(4,5),Fraction(3,5),n),
             "fresh_writes":fresh_certificate(Fraction(4,5),Fraction(3,5),n)}
            for n in (1,3,5,9)],
        "weak_write_scaling":"theta_n=2 sqrt(lambda/n)",
        "limiting_overlap":"exp(-lambda/2)",
        "limiting_worst_old_state_distance":"(1-exp(-lambda/2))/2",
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("repeated_weak_memory_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
