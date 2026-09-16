"""Round 97: sequential coherent writes of G=YY and H=IZ on one old system."""

import argparse
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from bilocal_record_tomography import embed_operator
from bipartite_composition import interaction
from complex_whole_real_interfaces import flagged_complex_whole, random_state
from independent_source_alignment import keep_systems
from quantum_interface_audit import IDENTITY, rotation_unitary
from role_symmetry_and_swap import pauli_word
from weak_relation_tradeoff import weak_write_unitary, trace_distance


G, H = pauli_word("YY"), pauli_word("IZ")
K = -1j*G@H


def memory_vector(sign,c):
    if sign not in (-1,1) or not 0<=c<=1:
        raise ValueError("Use sign +/-1 and overlap in [0,1].")
    return np.array([1.,0.]) if sign==-1 else np.array([c,math.sqrt((1-c)*(1+c))])


def projector(observable,sign):
    return (np.eye(4)+sign*observable)/2


def h_write_unitary(c):
    """H=+1 (old B=0) rotates memory, H=-1 leaves it at |0>.

    Control-on-zero copy: same one YX interaction and three local rotations
    as round 48, with the interaction angle sign reversed.
    """
    memory_vector(1,c)
    angle=2*math.acos(c)
    axis=np.kron(IDENTITY,rotation_unitary(-math.pi/2))
    pair=np.kron(rotation_unitary(angle/2),IDENTITY)@axis@interaction(angle/2,"YX")@axis.conj().T
    return embed_operator(pair,(0,2),3)


def sequence_unitary(c,d,reverse=False):
    """Order M(G record), N(H record), old A,B. Only write order reverses."""
    memory_vector(1,c); memory_vector(1,d)
    ug=embed_operator(weak_write_unitary(2*math.acos(c)),(0,2,3),4)
    uh=embed_operator(h_write_unitary(d),(1,2,3),4)
    return ug@uh if reverse else uh@ug


def sequence_isometry(c,d,reverse=False):
    return sequence_unitary(c,d,reverse)@np.kron(np.eye(4)[:,:1],np.eye(4))


def analytic_isometry(c,d,reverse=False):
    result=np.zeros((16,4),dtype=complex)
    for s in (-1,1):
        for t in (-1,1):
            p,q=projector(G,s),projector(H,t)
            memories=np.kron(memory_vector(s,c),memory_vector(t,d))
            result+=np.kron(memories[:,None],p@q if reverse else q@p)
    return result


def written_state(state,c,d,reverse=False):
    if len(state)%4: raise ValueError("Old dimension must be a multiple of four.")
    v=np.kron(sequence_isometry(c,d,reverse),np.eye(len(state)//4))
    return v@state@v.conj().T


def forget_memories(state):
    dim=len(state)//4
    return np.einsum("abad->bd",state.reshape(4,dim,4,dim))


def old_channel(state,c,d):
    extension=np.eye(len(state)//4)
    ops=(np.eye(len(state)),np.kron(G,extension),np.kron(H,extension),np.kron(H@G,extension))
    weights=((1+c)*(1+d),(1-c)*(1+d),(1+c)*(1-d),(1-c)*(1-d))
    return sum(weight*op@state@op.conj().T/4 for weight,op in zip(weights,ops))


def disturbance(c,d):
    return 1-(1+c)*(1+d)/4


def calibration_state(observable,sign):
    return (np.eye(4)+sign*observable)/4


class IncompatibleRelationWritesTests(unittest.TestCase):
    def test_control_zero_compiler_matches_analytic_H_isometry(self):
        for c in (0,.3,.8,1):
            v=h_write_unitary(c)@np.kron(np.eye(2)[:,:1],np.eye(4))
            expected=sum(np.kron(memory_vector(s,c)[:,None],projector(H,s)) for s in (-1,1))
            np.testing.assert_allclose(v,expected,atol=5e-16)

    def test_full_two_write_circuit_matches_both_time_ordered_isometries(self):
        np.testing.assert_array_equal(G@H+H@G,np.zeros((4,4)))
        for c,d in ((0,0),(.8,.8),(.3,.9),(1,1)):
            for reverse in (False,True):
                np.testing.assert_allclose(sequence_isometry(c,d,reverse),
                    analytic_isometry(c,d,reverse),atol=9e-16)
                u=sequence_unitary(c,d,reverse)
                np.testing.assert_allclose(u.conj().T@u,np.eye(16),atol=2e-15)

    def test_complete_old_channel_including_correlated_complex_spectators(self):
        rng=np.random.default_rng(97)
        for size in (4,8,12):
            state=random_state(rng,size)
            for reverse in (False,True):
                np.testing.assert_allclose(forget_memories(written_state(state,.8,.6,reverse)),
                    old_channel(state,.8,.6),atol=6e-16)

    def test_first_memory_and_untouched_reference_joint_state_survive_second_write(self):
        state=random_state(np.random.default_rng(197),8)
        before=written_state(state,.7,1)
        for d in (0,.5,.9):
            after=written_state(state,.7,d)
            # Five qubits M,N,A,B,C. The first memory jointly with reference C.
            np.testing.assert_allclose(keep_systems(after,(0,4),5),
                                       keep_systems(before,(0,4),5),atol=4e-16)

    def test_historical_record_survives_while_old_G_relation_changes(self):
        c,d=.8,.6
        for sign in (-1,1):
            state=flagged_complex_whole(sign)
            after=written_state(state,c,d)
            memory=keep_systems(after,(0,),5)
            m=memory_vector(sign,c)
            np.testing.assert_allclose(memory,np.outer(m,m),atol=7e-16)
            old=forget_memories(after)
            self.assertAlmostEqual(np.trace(old@np.kron(G,np.eye(2))).real,sign*d)
            self.assertAlmostEqual(trace_distance(old,state),(1-d)/2)

    def test_anticommuting_coordinates_have_different_old_retention_factors(self):
        for op,factor in ((G,.6),(H,.8),(K,.48)):
            state=calibration_state(op,1)
            self.assertAlmostEqual(np.trace(op@old_channel(state,.8,.6)).real,factor)

    def test_product_witness_attains_the_global_disturbance_bound(self):
        psi=np.array([1.,1.,0,0])/math.sqrt(2)  # old |0,+>
        frame=np.column_stack([op@psi for op in (np.eye(4),G,H,H@G)])
        np.testing.assert_allclose(frame.conj().T@frame,np.eye(4),atol=3e-16)
        state=np.outer(psi,psi)
        rng=np.random.default_rng(297)
        for c,d in ((0,0),(.8,.8),(.5,.3),(1,.7),(1,1)):
            self.assertAlmostEqual(trace_distance(old_channel(state,c,d),state),disturbance(c,d))
            rho=random_state(rng,12)
            self.assertLessEqual(trace_distance(old_channel(rho,c,d),rho),disturbance(c,d)+1e-14)

    def test_reversing_order_preserves_summed_old_channel_but_changes_memory(self):
        state=calibration_state(H,1)
        forward=written_state(state,.8,.6)
        reverse=written_state(state,.8,.6,True)
        np.testing.assert_allclose(forget_memories(forward),forget_memories(reverse),atol=5e-16)
        self.assertGreater(trace_distance(keep_systems(forward,(1,),4),
                                          keep_systems(reverse,(1,),4)),.05)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(IncompatibleRelationWritesTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    c=d=Fraction(4,5)
    report={"round":97,"observables":["YY","IZ"],"anticommute":True,
        "write_pair_rotations":4,"write_local_rotations":6,"new_memory_qubits":2,
        "first_memory_and_external_reference_marginal_unchanged_by_second_write":True,
        "worst_old_trace_distance_including_spectators":"1-(1+c)*(1+d)/4",
        "four_fifths_example":{"c":str(c),"d":str(d),"worst_old_disturbance_exact":str(disturbance(c,d)),
            "old_G_retention_exact":str(d),"old_H_retention_exact":str(c),
            "old_K_retention_exact":str(c*d),
            "old_disturbance_for_initial_G_sector_exact":str((1-d)/2)},
        "summed_old_channel_independent_of_write_order":True,
        "memory_record_independent_of_write_order":False,
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("incompatible_relation_writes_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
