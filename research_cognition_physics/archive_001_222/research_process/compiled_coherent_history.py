"""Round 146: compile coherent history compression into original real gates.

Use the freedom on unoccupied environment states: a structured completion
replaces round 144's arbitrary orthogonal completion but acts identically on
its entire seven-dimensional support. No active-system helper is used.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
import math
from pathlib import Path
import unittest

import numpy as np

from compiled_joining_frontier import control_circuit, compiled_sharing
from compiled_subject_joining import circuit_matrix, ry
from compiled_three_subject_joining import inverse, remap, gate_counts
from external_correlation_recovery_bound import channel
from minimal_coherent_joining_history import (
    environment_support, environment_transform, conditional_unitary, coherent_decoder_circuit,
)


NEGATIVE_PERMUTATION=(1,4,0,3,2,5,7,6)
# Chronological planar rotations G(i,j,angle), on physical environment basis states.
NEGATIVE_ROTATIONS=((12,15,-1),(14,15,F(1,2)),(10,12,F(1,2)),
                    (9,12,F(1,2)),(8,9,F(1,2)))


def pattern_ry(target,controls,angle):
    """controls is a tuple of (site,required bit); all lowering is explicit."""
    sites=[site for site,_ in controls]
    if target in sites or len(set(sites))!=len(sites) or any(bit not in (0,1) for _,bit in controls):
        raise ValueError("Distinct control wires and binary conditions required.")
    gates=(ry(target,angle),)
    for site,_ in controls: gates=control_circuit(site,gates)
    flips=tuple(ry(site,np.pi) for site,bit in controls if bit==0)
    return flips+gates+inverse(flips)


def adjacent_ry(first,second,angle,systems=4):
    delta=first ^ second
    if first<0 or second<0 or max(first,second)>=2**systems or delta==0 or delta & (delta-1):
        raise ValueError("Use two in-range basis indices differing at exactly one bit.")
    target=systems-delta.bit_length()
    controls=tuple((site,(first>>(systems-1-site))&1) for site in range(systems) if site!=target)
    return pattern_ry(target,controls,angle)


def two_level_rotation(first,second,angle,systems=4):
    """G(first,second,angle) with block [[cos,-sin],[sin,cos]].

    A signed Gray path moves first next to second; the sign and the target-bit
    order are explicitly included in the middle Ry angle. No free permutation.
    """
    if not 0<=first<second<2**systems: raise ValueError("Ordered distinct basis indices required.")
    differing=[site for site in range(systems) if (first ^ second)&(1<<(systems-1-site))]
    current=first; sign=1; path=()
    for site in differing[:-1]:
        next_value=current ^ (1<<(systems-1-site))
        if (current>>(systems-1-site))&1: sign*=-1
        path+=adjacent_ry(current,next_value,np.pi,systems)
        current=next_value
    last=differing[-1]
    orientation=1 if ((current>>(systems-1-last))&1)==0 else -1
    middle=adjacent_ry(current,second,2*sign*orientation*angle,systems)
    return path+middle+inverse(path)


def planar_matrix(first,second,angle,dimension=16):
    result=np.eye(dimension); c,s=math.cos(angle),math.sin(angle)
    result[np.ix_((first,second),(first,second))]=((c,-s),(s,c))
    return result


@lru_cache(maxsize=1)
def environment_circuit():
    gates=pattern_ry(1,((0,0),),2*math.asin(1/math.sqrt(3)))
    gates+=pattern_ry(2,((0,0),(1,0)),np.pi/2)
    for first,second,multiple in NEGATIVE_ROTATIONS:
        gates+=two_level_rotation(first,second,float(multiple)*np.pi)
    return gates


@lru_cache(maxsize=1)
def compiled_environment():
    return circuit_matrix(environment_circuit(),4)


@lru_cache(maxsize=1)
def compressed_physical_isometry():
    # Original physical order h,R,A,B,E1,E2,C -> Q,R,A,B,M1,M2,M3.
    gates=remap(environment_circuit(),{0:0,1:4,2:5,3:6})
    return circuit_matrix(gates,7,compiled_sharing(F(1,2))[0])


def compressed_branches():
    v=compressed_physical_isometry()
    return tuple(v.reshape((2,)*7+(16,)).transpose(4,5,6,0,1,2,3,7).reshape(8,16,16))


def resources():
    counts=gate_counts(environment_circuit())
    counts.update({"new_pure_auxiliaries":0,"accessed_old_environment_wires":4,
        "retained_accessible_coherent_history_rebits":1,"classical_label_wires":3,
        "conditional_decoder_yx_upper":16,"conditional_decoder_ry_upper":17,
        "gate_minimality_claimed":False})
    return counts


class CompiledCoherentHistoryTests(unittest.TestCase):
    def test_multiple_controls_are_original_gates_with_correct_zero_patterns(self):
        for controls in (((0,0),),((0,0),(1,1)),((0,1),(1,0),(2,1))):
            target=len(controls); systems=target+1; angle=.43
            expected=np.eye(2**systems)
            index=sum(bit<<(systems-1-site) for site,bit in controls)
            expected[np.ix_((index,index+1),(index,index+1))]=((math.cos(angle/2),-math.sin(angle/2)),(math.sin(angle/2),math.cos(angle/2)))
            np.testing.assert_allclose(circuit_matrix(pattern_ry(target,controls,angle),systems),expected,atol=2e-15)

    def test_signed_gray_paths_implement_nonadjacent_planar_rotations(self):
        for first,second in ((0,15),(3,12),(9,12),(10,12),(12,15)):
            for angle in (.37,-.81):
                actual=circuit_matrix(two_level_rotation(first,second,angle),4)
                np.testing.assert_allclose(actual,planar_matrix(first,second,angle),atol=7e-15)

    def test_five_exact_planar_rotations_realize_the_negative_label_permutation(self):
        actual=np.eye(16)
        for i,j,multiple in NEGATIVE_ROTATIONS:
            actual=planar_matrix(i,j,float(multiple)*np.pi) @ actual
        expected=np.eye(16); expected[8:,8:]=np.eye(8)[:,NEGATIVE_PERMUTATION]
        np.testing.assert_allclose(actual,expected,atol=3e-16)

    def test_compilation_matches_old_transform_on_every_occupied_environment_vector(self):
        actual=compiled_environment(); support=environment_support()
        np.testing.assert_allclose(actual @ support,environment_transform() @ support,atol=2e-14)
        np.testing.assert_allclose(actual.T @ actual,np.eye(16),atol=3e-14)
        self.assertAlmostEqual(np.linalg.det(actual),1.,places=12)

    def test_actual_seven_wire_circuit_keeps_six_reversible_branches(self):
        branches=compressed_branches()
        for label in range(6):
            np.testing.assert_allclose(branches[label],conditional_unitary(label+1)/np.sqrt(6),atol=2e-14)
            decoder=circuit_matrix(coherent_decoder_circuit(label+1),4)
            np.testing.assert_allclose(decoder @ branches[label],np.eye(16)/np.sqrt(6),atol=2e-14)
        np.testing.assert_allclose(branches[6:],np.zeros((2,16,16)),atol=2e-14)

    def test_native_environment_unitary_and_inverse_preserve_external_columns(self):
        rng=np.random.default_rng(146); columns=rng.normal(size=(16,5))+1j*rng.normal(size=(16,5))
        columns/=np.linalg.norm(columns)
        recovered=circuit_matrix(inverse(environment_circuit()),4,circuit_matrix(environment_circuit(),4,columns))
        np.testing.assert_allclose(recovered,columns,atol=2e-14)

    def test_current_common_marginal_remains_the_actual_old_channel(self):
        rng=np.random.default_rng(246); raw=rng.normal(size=(16,4))+1j*rng.normal(size=(16,4))
        source=raw @ raw.conj().T; source/=np.trace(source)
        joint=channel(compressed_branches(),source)
        marginal=np.trace(joint.reshape(2,8,2,8),axis1=0,axis2=2)
        np.testing.assert_allclose(marginal,channel(compiled_sharing(F(1,2))[1],source),atol=5e-15)

    def test_all_counts_use_only_native_environment_operations_and_guards(self):
        gates=environment_circuit()
        self.assertEqual(sum(g.kind=="YX" for g in gates),415)
        self.assertTrue(all(g.kind in ("Ry","YX") and all(0<=s<4 for s in g.sites) for g in gates))
        self.assertEqual(resources()["new_pure_auxiliaries"],0)
        with self.assertRaises(ValueError): adjacent_ry(0,3,.1)
        with self.assertRaises(ValueError): two_level_rotation(3,3,.1)
        with self.assertRaises(ValueError): pattern_ry(0,((0,1),),.1)


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CompiledCoherentHistoryTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":146,"resources":resources(),
        "environment_gate_list":[{"kind":g.kind,"sites":list(g.sites),"angle_radians":g.angle} for g in environment_circuit()],
        "equivalent_to_round_144_on_full_occupied_support":True,
        "same_action_on_unoccupied_complement_required":False,
        "original_common_marginal_preserved":True,"arbitrary_external_recovery_with_exact_labels":True,
        "finite_noisy_label_readout_completed_in_this_round":False,
        "gates_have_ideal_controllable_angles":True,
        "gate_angle_trace_error_bound":"sum_j abs(delta_angle_j)/2",
        "per_gate_tolerance_for_extra_error_one_thousandth_radians":2e-3/(len(environment_circuit())+33),
        "automated_checks":{"run":checks.testsRun,"failures":len(checks.failures),"errors":len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name("compiled_coherent_history_results.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    summary={key:value for key,value in report.items() if key!="environment_gate_list"}
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=="__main__": main()
