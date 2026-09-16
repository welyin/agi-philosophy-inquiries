"""Round 99: recover a choice of incompatible queries from coherent memory correlations."""

import argparse
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I
from complex_whole_real_interfaces import random_state
from incompatible_relation_writes import (G,H,K,memory_vector,sequence_isometry,
    written_state,calibration_state,forget_memories,old_channel)
from joint_relation_records import (memory_axis,pointer_effect,read_memory_tree,
    axis_read_memory_tree,interval_fields,record_certificate)
from one_bit_real_network import visibility_interval
from operational_effect_closure import compiled_copy_gate,majority_effect_by_circuit
from quantum_interface_audit import PAULI_X,PAULI_Y,PAULI_Z,rotation_unitary,effective_readout
from weak_relation_tradeoff import read_axis,trace_distance


def memory_channel(state,c,d):
    joint=written_state(state,c,d)
    dim=len(state)
    return np.einsum("abcb->ac",joint.reshape(4,dim,4,dim))


def memory_blocks(c):
    minus,plus=memory_vector(-1,c),memory_vector(1,c)
    return (np.outer(minus,minus)+np.outer(plus,plus),
            np.outer(minus,plus)+np.outer(plus,minus),
            np.outer(plus,plus)-np.outer(minus,minus))


def analytic_calibration_memories(observable,sign,c,d):
    ac,bc,_=memory_blocks(c)
    ad,_,dd=memory_blocks(d)
    if observable=="G":
        m=memory_vector(sign,c)
        return np.kron(np.outer(m,m),ad)/2
    if observable=="H":
        return (np.kron(ac,ad)+sign*np.kron(bc,dd))/4
    raise ValueError("Choose G or H.")


def query_observables(c,d):
    complement=c*PAULI_Z+math.sqrt((1-c)*(1+c))*PAULI_X
    return np.kron(memory_axis(c),np.eye(2)),np.kron(complement,memory_axis(d))


def h_query_decoder(c,d):
    """Two local axis rotations, then one compiled controlled-Ry(pi) parity gate."""
    rotate=np.kron(rotation_unitary(-math.acos(c)),
                   rotation_unitary(-read_axis(2*math.acos(d))))
    return compiled_copy_gate()@rotate


def chosen_memory_effect(c,d,query,count):
    if query=="G": return np.kron(pointer_effect(c,count,1),np.eye(2))
    if query!="H": raise ValueError("Choose G or H.")
    # This decoder transfers the product observable to Z on M.
    positive=effective_readout(0,1,1.) if count==1 else majority_effect_by_circuit(count)
    u=h_query_decoder(c,d)
    return u.conj().T@np.kron(positive,np.eye(2))@u


def chosen_record_tree(state,c,d,query,count):
    joint=written_state(state,c,d)
    if query=="G": return read_memory_tree(joint,0,c,count)
    if query!="H": raise ValueError("Choose G or H.")
    u=np.kron(h_query_decoder(c,d),np.eye(len(state)))
    return axis_read_memory_tree(u@joint@u.conj().T,0,0.,count)


def query_certificate(c,d,count):
    c,d=Fraction(c),Fraction(d)
    if not 0<=c<=1 or not 0<=d<=1: raise ValueError("Use overlaps in [0,1].")
    gamma=visibility_interval(count)
    v,w=(1-I.exact(c*c)).sqrt(),(1-I.exact(d*d)).sqrt()
    return {"c_exact":str(c),"d_exact":str(d),"old_reads_for_chosen_query":count,
        "G_chosen_success":interval_fields((1+gamma*v)/2),
        "H_chosen_success":interval_fields((1+gamma*w)/2),
        "H_from_second_memory_alone_success":interval_fields((1+gamma*I.exact(c)*w)/2),
        "H_joint_improvement":interval_fields(gamma*(1-I.exact(c))*w/2),
        "H_local_product_success_with_count_reads_per_memory":interval_fields((1+gamma*gamma*w)/2),
        "local_product_total_old_reads":2*count,"local_product_message_bits":1,
        "H_decoder_pair_rotations":1,"H_decoder_local_rotations":5,
        "two_answers_produced_in_one_run":False}


class DeferredRelationQueriesTests(unittest.TestCase):
    def test_joint_memory_calibration_states_match_analytic_blocks(self):
        for c,d in ((0,0),(.8,.8),(.4,.9),(1,.5),(.5,1)):
            for label,op in (("G",G),("H",H)):
                for sign in (-1,1):
                    np.testing.assert_allclose(memory_channel(calibration_state(op,sign),c,d),
                        analytic_calibration_memories(label,sign,c,d),atol=7e-16)

    def test_joint_memories_keep_each_binary_trace_distance_at_its_own_write_strength(self):
        for c,d in ((0,0),(.8,.8),(.3,.7),(1,.4),(.4,1)):
            for op,expected in ((G,math.sqrt(1-c*c)),(H,math.sqrt(1-d*d))):
                minus=memory_channel(calibration_state(op,-1),c,d)
                plus=memory_channel(calibration_state(op,1),c,d)
                self.assertAlmostEqual(trace_distance(minus,plus),expected,places=14)

    def test_chosen_observables_pull_back_to_original_unattenuated_by_other_write(self):
        for c,d in ((0,0),(.8,.8),(.3,.7),(1,.4),(.4,1)):
            v=sequence_isometry(c,d)
            og,oh=query_observables(c,d)
            np.testing.assert_allclose(v.conj().T@np.kron(og,np.eye(4))@v,
                                       math.sqrt(1-c*c)*G,atol=1e-15)
            np.testing.assert_allclose(v.conj().T@np.kron(oh,np.eye(4))@v,
                                       math.sqrt(1-d*d)*H,atol=1e-15)

    def test_actual_old_gate_compiles_the_H_product_observable_without_SWAPPING(self):
        for c,d in ((0,0),(.8,.8),(1,.4),(.4,1)):
            u=h_query_decoder(c,d)
            np.testing.assert_allclose(u.conj().T@np.kron(PAULI_Z,np.eye(2))@u,
                                       query_observables(c,d)[1],atol=8e-16)
            self.assertAlmostEqual(np.linalg.det(u).real,1)
            positive=np.zeros((4,4),dtype=complex)
            for result in (0,1):
                positive+=np.kron(effective_readout(math.acos(c),result,1.),
                    effective_readout(read_axis(2*math.acos(d)),result,1.))
            gamma=sum(visibility_interval(1).floats())/2
            np.testing.assert_allclose(positive,
                (np.eye(4)+gamma**2*query_observables(c,d)[1])/2,atol=4e-16)

    def test_full_Kraus_records_match_finite_certificates_and_operator_effects(self):
        for count in (1,3):
            row=query_certificate(Fraction(4,5),Fraction(4,5),count)
            for query,op in (("G",G),("H",H)):
                success=0.
                for sign in (-1,1):
                    leaves=chosen_record_tree(calibration_state(op,sign),.8,.8,query,count)
                    self.assertAlmostEqual(sum(np.trace(b).real for b in leaves.values()),1,places=13)
                    success+=sum(np.trace(b).real for h,b in leaves.items()
                        if (1 if sum(h)>count//2 else -1)==sign)/2
                self.assertAlmostEqual(success,row[query+"_chosen_success"]["diagnostic"],places=13)
                v=sequence_isometry(.8,.8)
                actual=v.conj().T@np.kron(chosen_memory_effect(.8,.8,query,count),np.eye(4))@v
                gamma=sum(visibility_interval(count).floats())/2
                np.testing.assert_allclose(actual,(np.eye(4)+gamma*.6*op)/2,atol=1e-15)
        # Reading either query acts only on the memories, even for complex old
        # states correlated with an unoperated reference.
        state=random_state(np.random.default_rng(199),8)
        for query in ("G","H"):
            leaves=chosen_record_tree(state,.8,.6,query,3)
            np.testing.assert_allclose(sum(forget_memories(b) for b in leaves.values()),
                                       old_channel(state,.8,.6),atol=2e-15)

    def test_ideal_optimal_pointer_queries_anticommute_and_have_unique_Helstrom_projectors(self):
        c,d=.8,.6
        og,oh=query_observables(c,d)
        np.testing.assert_allclose(og@oh+oh@og,0,atol=2e-16)
        for label,observable in (("G",og),("H",oh)):
            diff=analytic_calibration_memories(label,1,c,d)-analytic_calibration_memories(label,-1,c,d)
            eigen,basis=np.linalg.eigh(diff)
            self.assertGreater(np.min(np.abs(eigen)),.01)
            sign=(basis*np.sign(eigen))@basis.conj().T
            np.testing.assert_allclose(sign,observable,atol=1e-15)
        self.assertGreater(np.linalg.norm(og@oh-oh@og),3.9)

    def test_strong_writes_encode_one_logical_qubit_algebra_in_the_memories(self):
        og,oh=query_observables(0,0)
        ok=-1j*og@oh
        np.testing.assert_allclose(ok,np.kron(PAULI_Y,PAULI_Z),atol=0)
        rng=np.random.default_rng(99)
        for _ in range(12):
            state=random_state(rng,4)
            g,h,k=(np.trace(op@state).real for op in (G,H,K))
            expected=(np.eye(4)+g*og+h*oh+k*ok)/4
            np.testing.assert_allclose(memory_channel(state,0,0),expected,atol=6e-16)

    def test_certified_recovery_exceeds_second_memory_alone_but_is_a_choice(self):
        row=query_certificate(Fraction(4,5),Fraction(4,5),1)
        self.assertGreater(Fraction(row["H_joint_improvement"]["lower_exact"]),Fraction(59,1000))
        baseline=record_certificate(Fraction(4,5),Fraction(4,5),1,1)
        self.assertAlmostEqual(row["H_from_second_memory_alone_success"]["diagnostic"],
                               baseline["H_calibration_success"]["diagnostic"])
        self.assertFalse(row["two_answers_produced_in_one_run"])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(DeferredRelationQueriesTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":99,"G_query":"(c*X-sqrt(1-c^2)*Z)_M tensor I_N",
        "H_query":"(c*Z+sqrt(1-c^2)*X)_M tensor (d*X-sqrt(1-d^2)*Z)_N",
        "ideal_joint_memory_G_trace_distance":"sqrt(1-c^2)",
        "ideal_joint_memory_H_trace_distance":"sqrt(1-d^2)",
        "optimal_ideal_queries_anticommute":True,
        "one_run_attains_both_individual_ideal_optima_for_nontrivial_writes":False,
        "strong_write_memory_state":"(I-g*ZI-h*XZ+k*YZ)/4; g=<YY>,h=<IZ>,k=<YX>",
        "arbitrary_old_four_dimensional_state_fully_copied_to_memories":False,
        "examples":[query_certificate(Fraction(4,5),Fraction(4,5),m) for m in (1,3,5)],
        "strong_write_chosen_query":query_certificate(Fraction(0),Fraction(0),1),
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("deferred_relation_queries_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
