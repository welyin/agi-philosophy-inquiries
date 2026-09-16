"""Round 104: finite real compression and one-bit recovery from accessible fragments."""
import argparse
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from bilocal_record_tomography import embed_operator
from certified_intervals import Interval as I
from coherent_memory_compression import compressor
from complex_whole_real_interfaces import random_state
from deferred_relation_queries import memory_channel,h_query_decoder
from environment_query_recovery import (branch_moments,analytic_environment_branch,
    interval_abs,interval_max)
from incompatible_relation_writes import H,calibration_state,written_state,forget_memories,old_channel
from independent_source_alignment import keep_systems
from joint_relation_records import (axis_read_memory_tree,read_memory_tree,interval_fields,validated_count)
from memory_dephasing_threshold import dephase_memory,dilate_memory
from one_bit_real_network import visibility_interval
from partial_environment_access import accessible_state,partition_overlaps
from weak_relation_tradeoff import forget_front_pointer,trace_distance


def compress_accessible(state,overlaps):
    """State starts with accessible E wires. Joint XY and YX control is required."""
    if not overlaps: raise ValueError("At least one accessible fragment is needed.")
    systems=int(round(math.log2(len(state))))
    total=np.eye(len(state),dtype=complex)
    a=overlaps[0]
    for site,lam in enumerate(overlaps[1:],1):
        total=embed_operator(compressor(a,lam),(0,site),systems)@total
        a*=lam
    return total@state@total.conj().T


def compressed_environment_branches(state,c,overlaps,accessible,count):
    accessible=tuple(accessible)
    if not accessible: raise ValueError("Use a nonempty accessible subset.")
    a,_=partition_overlaps(overlaps,accessible)
    joint=accessible_state(state,c,overlaps,accessible)
    compressed=compress_accessible(joint,[overlaps[i] for i in accessible])
    k=len(accessible)
    systems=int(round(math.log2(len(compressed))))
    effective=keep_systems(compressed,(0,)+tuple(range(k,systems)),systems)
    result={r:np.zeros_like(state,dtype=complex) for r in (-1,1)}
    for history,branch in axis_read_memory_tree(effective,0,math.acos(a),count).items():
        r=1 if sum(history)>count//2 else -1
        result[r]+=forget_front_pointer(branch)
    return result


def partial_factor(c,a,b,gamma):
    return sum(max(c*p,b*abs(t)) for r in (-1,1) for p,t in [branch_moments(a,gamma,r)])


def partial_strategy(c,a,b,gamma,r):
    p,t=branch_moments(a,gamma,r)
    return ("joint",1 if t>=0 else -1) if b*abs(t)>c*p else ("second_only",1)


def query_branches(branch,c,d,count,strategy,sign):
    if strategy=="joint":
        u=np.kron(h_query_decoder(c,d),np.eye(len(branch)//4))
        leaves=axis_read_memory_tree(u@branch@u.conj().T,0,0.,count)
    else:
        leaves=read_memory_tree(branch,1,d,count)
    return {h:(sign*(1 if sum(h)>count//2 else -1),rho) for h,rho in leaves.items()}


def partial_recovery_records(old,c,d,overlaps,accessible,environment_count,query_count):
    a,b=partition_overlaps(overlaps,accessible)
    gamma=sum(visibility_interval(environment_count).floats())/2
    records=compressed_environment_branches(written_state(old,c,d),c,overlaps,accessible,environment_count)
    return {(r,h):(guess,rho) for r,branch in records.items()
            for h,(guess,rho) in query_branches(branch,c,d,query_count,
                *partial_strategy(c,a,b,gamma,r)).items()}


def compiler_budget(fragment_count,accessible_count,environment_count,query_count):
    if not 1<=accessible_count<=fragment_count: raise ValueError("Nonempty accessible subset required.")
    validated_count(environment_count);validated_count(query_count)
    n,k,m,q=fragment_count,accessible_count,environment_count,query_count
    em=m if m>1 else 0
    qm=q if q>1 else 0
    reused=bool(m>1 and k>1)
    slots=n+2+int(m>1 and not reused)+int(q>1)
    return {"scope":"Declared batch compiler, worst query branch, no angle simplification; original writes included",
        "created_fragments":n,"accessed_fragments":k,
        "pair_rotations":4+n+3*(k-1)+em+1+qm,
        "extra_local_rotations":6+5*n+3*(k-1)+(3*m+1 if m>1 else 0)+5+3*qm,
        "environment_readouts":m,"final_readouts":q,"total_readouts":m+q,
        "global_extra_quantum_slots_peak":slots,"initial_pure_preparations":slots,
        "pointer_resets":(m-1 if m>1 else 0)+(q-1 if q>1 else 0),
        "environment_pointer_reuses_purified_compression_scratch":reused,
        "return_classical_message_bits":1,
        "requires_coherent_gates_between_accessible_fragments":k>1,
        "requires_coherent_return_gate_between_environment_and_memory":False,
        "fragment_transport_or_routing_cost_included":False}


def partial_certificate(c,d,overlaps,accessible,environment_count,query_count):
    c,d=Fraction(c),Fraction(d)
    if not 0<=c<=1 or not 0<=d<=1: raise ValueError("Use c,d in [0,1].")
    overlaps=tuple(map(Fraction,overlaps));accessible=tuple(accessible)
    a,b=partition_overlaps(overlaps,accessible)
    budget=compiler_budget(len(overlaps),len(accessible),environment_count,query_count)
    gamma=visibility_interval(environment_count);eta=visibility_interval(query_count)
    total=I.exact(0);choices=[]
    for r in (-1,1):
        p=(1+r*I.exact(a)*gamma)/2
        t=I.exact(b)*(I.exact(a)+r*gamma)/2
        local=I.exact(c)*p;joint=interval_abs(t)
        total+=interval_max(local,joint)
        choice=("joint_positive" if t.lo>=0 else "joint_negative") if joint.lo>local.hi else (
            "second_only" if local.lo>joint.hi else "tie_or_interval_overlap")
        choices.append({"sign":r,"certified_strategy":choice})
    w=(1-I.exact(d*d)).sqrt()
    return {"c_exact":str(c),"d_exact":str(d),"accessible_count":len(accessible),
        "accessible_overlap_exact":str(a),"inaccessible_overlap_exact":str(b),
        "recovery_factor":interval_fields(total),
        "finite_H_success":interval_fields((1+eta*w*total)/2),
        "without_environment_success":interval_fields((1+eta*w*I.exact(max(c,a*b)))/2),
        "accessible_ideal_H_distance":interval_fields(w*I.exact(max(c,b))),
        "strategies":choices,"resource_budget":budget,
        "finite_resource_global_optimality_proved":False}


class FragmentQueryCompilerTests(unittest.TestCase):
    def test_actual_compression_leaves_scratch_pure_with_complex_reference(self):
        rho=random_state(np.random.default_rng(104),8)
        overlaps=(.4,.7,.9);subset=(2,0)
        a,b=partition_overlaps(overlaps,subset)
        state=compress_accessible(accessible_state(rho,.8,overlaps,subset),[overlaps[i] for i in subset])
        effective=dilate_memory(dephase_memory(rho,.8,b),.8,a)
        # Expected E_eff, zero scratch, M,N,reference. Index change is bookkeeping.
        expected=np.kron(effective,np.diag([1.,0.]))
        expected=keep_systems(expected,(0,4,1,2,3),5)
        np.testing.assert_allclose(state,expected,atol=2e-15)

    def test_all_actual_compressed_environment_messages_equal_residual_dephasing_formula(self):
        rho=random_state(np.random.default_rng(204),4)
        for overlaps,subset in (((.2,.7,.95),(2,0)),((0.,.8),(0,1)),((.2,0.),(0,))):
            a,b=partition_overlaps(overlaps,subset)
            for count in (1,3):
                gamma=sum(visibility_interval(count).floats())/2
                for r,branch in compressed_environment_branches(rho,.8,overlaps,subset,count).items():
                    expected=dephase_memory(analytic_environment_branch(rho,.8,a,gamma,r),.8,b)
                    np.testing.assert_allclose(branch,expected,atol=2e-15)

    def test_message_block_H_trace_norm_attains_formula(self):
        for overlaps,subset in (((.7,.95),(0,)),((.9,.9),(0,1)),((0.,0.),(0,))):
            a,b=partition_overlaps(overlaps,subset)
            pair=[compressed_environment_branches(memory_channel(calibration_state(H,s),.8,.8),
                    .8,overlaps,subset,1) for s in (-1,1)]
            actual=sum(trace_distance(pair[0][r],pair[1][r]) for r in (-1,1))
            self.assertAlmostEqual(actual,.6*partial_factor(.8,a,b,sum(visibility_interval(1).floats())/2),places=13)

    def test_finite_actual_feedback_queries_attain_strict_interval_predictions(self):
        f=Fraction
        for overlaps,subset,m,q in (((f(7,10),f(19,20)),(0,),1,1),
                ((f(9,10),f(9,10)),(0,1),3,1),((f(0),f(0)),(0,),1,3)):
            cert=partial_certificate(f(4,5),f(4,5),overlaps,subset,m,q)
            success=0.
            for s in (-1,1):
                leaves=partial_recovery_records(calibration_state(H,s),.8,.8,list(map(float,overlaps)),subset,m,q)
                self.assertAlmostEqual(sum(np.trace(rho).real for _,rho in leaves.values()),1,places=13)
                success+=sum(np.trace(rho).real for guess,rho in leaves.values() if guess==s)/2
            self.assertAlmostEqual(success,cert["finite_H_success"]["diagnostic"],places=13)

    def test_perfect_and_uninformative_messages_obey_universal_access_bound(self):
        for c in (0.,.3,.8,1.):
            for a in (0.,.4,1.):
                for b in (0.,.6,.95,1.):
                    self.assertAlmostEqual(partial_factor(c,a,b,1),max(c,b))
                    self.assertAlmostEqual(partial_factor(c,a,b,0),max(c,a*b))
                    value=partial_factor(c,a,b,.97)
                    self.assertLessEqual(value,max(c,b)+1e-15)
                    self.assertGreaterEqual(value,max(c,a*b)-1e-15)

    def test_feedback_keeps_old_unconditional_complex_reference(self):
        old=random_state(np.random.default_rng(304),8)
        records=partial_recovery_records(old,.8,.6,(.7,.95),(0,),1,1)
        np.testing.assert_allclose(sum(forget_memories(rho) for _,rho in records.values()),
                                   old_channel(old,.8,.6),atol=2e-15)

    def test_batch_budget_matches_round102_and_reuses_only_local_pure_scratch(self):
        self.assertEqual([compiler_budget(1,1,1,1)[k] for k in
                          ("pair_rotations","extra_local_rotations","global_extra_quantum_slots_peak")],[6,16,3])
        self.assertEqual([compiler_budget(1,1,3,3)[k] for k in
                          ("pair_rotations","extra_local_rotations","global_extra_quantum_slots_peak")],[12,35,5])
        row=compiler_budget(3,2,3,3)
        self.assertEqual(row["pair_rotations"],17)
        self.assertEqual(row["global_extra_quantum_slots_peak"],6)
        self.assertTrue(row["environment_pointer_reuses_purified_compression_scratch"])
        self.assertFalse(row["fragment_transport_or_routing_cost_included"])

    def test_large_fragment_count_has_exact_scalar_certificate_without_tensor_simulation(self):
        f=Fraction
        row=partial_certificate(f(4,5),f(4,5),[f(19,20)]*100,range(96),3,1)
        self.assertGreater(Fraction(row["finite_H_success"]["lower_exact"]),
                           Fraction(row["without_environment_success"]["upper_exact"]))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(FragmentQueryCompilerTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    f=Fraction
    report={"round":104,"formula":"T=sum_r max(c*(1+r*a*gamma)/2,b*abs(a+r*gamma)/2)",
        "a":"product of accessible overlaps","b":"product of inaccessible overlaps",
        "ideal_environment_reading_attains_universal_accessible_H_bound":True,
        "finite_examples":[partial_certificate(f(4,5),f(4,5),[f(19,20)]*100,range(k),3,1)
                           for k in (95,96,99,100)],
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("fragment_query_compiler_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
