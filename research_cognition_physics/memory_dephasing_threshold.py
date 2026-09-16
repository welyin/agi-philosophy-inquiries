"""Round 100: preserved historical answers and a threshold for future queries."""

import argparse
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I
from complex_whole_real_interfaces import random_state
from deferred_relation_queries import memory_channel, memory_blocks, h_query_decoder
from incompatible_relation_writes import (G,H,written_state,calibration_state,
    forget_memories,old_channel)
from joint_relation_records import memory_axis,read_memory_tree,axis_read_memory_tree,interval_fields
from one_bit_real_network import visibility_interval
from operational_effect_closure import compiled_copy_gate
from quantum_interface_audit import rotation_unitary
from weak_relation_tradeoff import read_axis,forget_front_pointer,trace_distance


def valid_overlap(value):
    if not 0<=value<=1: raise ValueError("Use a parameter in [0,1].")


def dephase_memory(state,c,lam):
    valid_overlap(c); valid_overlap(lam)
    a=np.kron(memory_axis(c),np.eye(len(state)//2))
    return (1+lam)*state/2+(1-lam)*a@state@a/2


def dephasing_unitary(c,lam):
    """Environment E first, memory M second. One pair gate, five local rotations."""
    valid_overlap(c); valid_overlap(lam)
    angle=read_axis(2*math.acos(c))
    rotate=np.kron(np.eye(2),rotation_unitary(angle))
    return rotate@compiled_copy_gate(2*math.acos(lam))@rotate.conj().T


def dilate_memory(state,c,lam):
    u=np.kron(dephasing_unitary(c,lam),np.eye(len(state)//2))
    return u@np.kron(np.diag([1.,0]),state)@u.conj().T


def h_strategy(c,lam):
    valid_overlap(c); valid_overlap(lam)
    return "joint" if lam>c else "second_only"


def dephased_query_records(state,c,d,lam,count,query="H"):
    noisy=dephase_memory(written_state(state,c,d),c,lam)
    if query=="G": return read_memory_tree(noisy,0,c,count)
    if query!="H": raise ValueError("Choose G or H.")
    if h_strategy(c,lam)=="second_only":
        return read_memory_tree(noisy,1,d,count)
    decode=np.kron(h_query_decoder(c,d),np.eye(len(state)))
    return axis_read_memory_tree(decode@noisy@decode.conj().T,0,0.,count)


def threshold_certificate(c,d,lam,count):
    c,d,lam=map(Fraction,(c,d,lam))
    for x in (c,d,lam): valid_overlap(x)
    gamma=visibility_interval(count)
    v,w=(1-I.exact(c*c)).sqrt(),(1-I.exact(d*d)).sqrt()
    best=I.exact(max(c,lam))*w
    return {"c_exact":str(c),"d_exact":str(d),"retained_coherence_exact":str(lam),
        "selected_H_strategy":h_strategy(c,lam),"final_old_readouts":count,
        "G_ideal_trace_distance":interval_fields(v),
        "H_ideal_trace_distance":interval_fields(best),
        "G_finite_success":interval_fields((1+gamma*v)/2),
        "H_finite_success":interval_fields((1+gamma*best)/2),
        "H_fixed_old_joint_axis_success":interval_fields((1+gamma*I.exact(lam)*w)/2),
        "old_system_disturbance_exact":str(1-(1+c)*(1+d)/4)}


def first_threshold_crossing(c,per_step):
    c,per_step=Fraction(c),Fraction(per_step)
    if not 0<c<1 or not 0<per_step<1:
        raise ValueError("Use nontrivial exact parameters.")
    n,coherence=0,Fraction(1)
    while coherence>c:
        n+=1; coherence*=per_step
    return {"steps":n,"previous_coherence_exact":str(per_step**(n-1)),
            "crossing_coherence_exact":str(coherence),"threshold_exact":str(c)}


class MemoryDephasingThresholdTests(unittest.TestCase):
    def test_actual_real_dilation_matches_dephasing_on_arbitrary_complex_correlations(self):
        rng=np.random.default_rng(100)
        for size in (4,8,12):
            state=random_state(rng,size)
            for c,lam in ((0,0),(.8,0),(.8,.7),(1,1)):
                np.testing.assert_allclose(forget_front_pointer(dilate_memory(state,c,lam)),
                                           dephase_memory(state,c,lam),atol=7e-16)

    def test_full_dilation_is_unitary_and_can_be_reversed_with_environment_access(self):
        state=random_state(np.random.default_rng(200),8)
        u=np.kron(dephasing_unitary(.8,.3),np.eye(4))
        initial=np.kron(np.diag([1.,0]),state)
        np.testing.assert_allclose(u.conj().T@u,np.eye(16),atol=2e-15)
        np.testing.assert_allclose(u.conj().T@dilate_memory(state,.8,.3)@u,initial,atol=1e-15)

    def test_exact_calibration_trace_distances_and_switching_point(self):
        for c,d in ((0,0),(.8,.8),(.3,.7),(1,.4),(.4,1)):
            for lam in (0,.2,c,.9,1):
                for op,target in ((G,math.sqrt(1-c*c)),(H,math.sqrt(1-d*d)*max(c,lam))):
                    outputs=[dephase_memory(memory_channel(calibration_state(op,s),c,d),c,lam)
                             for s in (-1,1)]
                    self.assertAlmostEqual(trace_distance(*outputs),target,places=14)

    def test_H_difference_spectrum_proves_the_threshold(self):
        c,d=.8,.6
        ac,bc,_=memory_blocks(c)
        _,_,dd=memory_blocks(d)
        for lam in (0,.4,.8,1):
            expected=np.kron(c*np.eye(2)+lam*(bc-c*np.eye(2)),dd)/2
            outputs=[dephase_memory(memory_channel(calibration_state(H,s),c,d),c,lam)
                     for s in (-1,1)]
            np.testing.assert_allclose(outputs[1]-outputs[0],expected,atol=8e-16)
            values=np.linalg.eigvalsh(c*np.eye(2)+lam*(bc-c*np.eye(2)))
            np.testing.assert_allclose(values,[c-lam,c+lam],atol=3e-16)

    def test_noise_changes_the_memory_but_not_G_axis_statistics(self):
        state=memory_channel(calibration_state(G,1),.8,.8)
        noisy=dephase_memory(state,.8,0)
        self.assertGreater(trace_distance(state,noisy),.3)
        a=np.kron(memory_axis(.8),np.eye(2))
        self.assertAlmostEqual(np.trace(a@state).real,np.trace(a@noisy).real,places=14)

    def test_complete_finite_records_attain_each_selected_strategy(self):
        for lam in (.4,.8,.9,1):
            for count in (1,3):
                row=threshold_certificate(Fraction(4,5),Fraction(4,5),Fraction(str(lam)),count)
                success=0.
                for s in (-1,1):
                    leaves=dephased_query_records(calibration_state(H,s),.8,.8,lam,count)
                    self.assertAlmostEqual(sum(np.trace(x).real for x in leaves.values()),1,places=13)
                    success+=sum(np.trace(b).real for h,b in leaves.items()
                        if (1 if sum(h)>count//2 else -1)==s)/2
                self.assertAlmostEqual(success,row["H_finite_success"]["diagnostic"],places=13)

    def test_old_system_channel_and_complex_reference_do_not_change_under_memory_noise(self):
        state=random_state(np.random.default_rng(300),8)
        for lam in (0,.5,1):
            noisy=dephase_memory(written_state(state,.8,.6),.8,lam)
            np.testing.assert_allclose(forget_memories(noisy),old_channel(state,.8,.6),atol=6e-16)

    def test_repeated_noise_multiplication_and_exact_finite_threshold_certificate(self):
        state=random_state(np.random.default_rng(400),4)
        actual=dephase_memory(dephase_memory(state,.8,.7),.8,.9)
        np.testing.assert_allclose(actual,dephase_memory(state,.8,.63),atol=3e-16)
        row=first_threshold_crossing(Fraction(4,5),Fraction(19,20))
        self.assertEqual(row["steps"],5)
        self.assertGreater(Fraction(row["previous_coherence_exact"]),Fraction(4,5))
        self.assertLess(Fraction(row["crossing_coherence_exact"]),Fraction(4,5))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(MemoryDephasingThresholdTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":100,"dephasing":"Lambda_lam(rho)=(1+lam)rho/2+(1-lam)A rho A/2; A=N_c on M",
        "ideal_H_distinguishability":"sqrt(1-d^2)*max(c,lam)",
        "ideal_G_distinguishability":"sqrt(1-c^2)",
        "dilation_pair_rotations":1,"dilation_local_rotations":5,"initialized_environment_qubits":1,
        "environment_access_required_for_inverting_dilation":True,
        "examples":[threshold_certificate(Fraction(4,5),Fraction(4,5),lam,1)
                    for lam in (Fraction(1),Fraction(9,10),Fraction(4,5),Fraction(2,5),Fraction(0))],
        "repeated_noise":first_threshold_crossing(Fraction(4,5),Fraction(19,20)),
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("memory_dephasing_threshold_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
