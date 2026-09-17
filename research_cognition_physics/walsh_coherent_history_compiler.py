"""Round 148: commuting Pauli expansion reduces the native compression cost.

Replace nested control lowering by the exact expansion of a bit-pattern
projector. Every Pauli exponential is built by native YX conjugations of Ry.
The entire structured environment operation, including its unused complement,
is the same as round 146. This is an improved upper bound, not optimal synthesis.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
import math
from pathlib import Path
import unittest

import numpy as np

from compiled_coherent_history import (
    NEGATIVE_ROTATIONS, compiled_environment, compressed_physical_isometry,
    pattern_ry, planar_matrix, resources as old_resources,
)
from compiled_joining_frontier import compiled_sharing
from compiled_subject_joining import circuit_matrix, ry, yx
from compiled_three_subject_joining import inverse, remap, gate_counts
from external_correlation_recovery_bound import channel
from independent_source_alignment import keep_systems
from minimal_coherent_joining_history import conditional_unitary
from noisy_classical_history import physical_read_branch
from noisy_coherent_history_recovery import certificate, decode_label, recovery_kraus
from quantum_interface_audit import ALPHA, PAULI_Y, PAULI_Z


def walsh_pattern_ry(target,controls,angle):
    sites=[site for site,_ in controls]
    if target in sites or len(set(sites))!=len(sites) or any(bit not in (0,1) for _,bit in controls):
        raise ValueError("Distinct control wires and binary conditions required.")
    result=()
    for mask in range(2**len(controls)):
        selected=[controls[j] for j in range(len(controls)) if (mask>>j)&1]
        sign=(-1)**sum(bit for _,bit in selected)
        if not selected:
            result+=(ry(target,angle/2**len(controls)),)
            continue
        start=selected[0][0]; current=start; chain=()
        for next_site in [site for site,_ in selected[1:]]+[target]:
            chain+=(yx(next_site,current,np.pi/2),); current=next_site
        result+=inverse(chain)+(ry(start,sign*angle/2**len(controls)),)+chain
    return result


def adjacent_walsh_ry(first,second,angle,systems=4):
    delta=first ^ second
    if first<0 or second<0 or max(first,second)>=2**systems or delta==0 or delta & (delta-1):
        raise ValueError("Use adjacent in-range binary basis indices.")
    target=systems-delta.bit_length()
    controls=tuple((site,(first>>(systems-1-site))&1) for site in range(systems) if site!=target)
    return walsh_pattern_ry(target,controls,angle)


def walsh_two_level_rotation(first,second,angle,systems=4):
    if not 0<=first<second<2**systems: raise ValueError("Ordered distinct basis indices required.")
    differing=[site for site in range(systems) if (first ^ second)&(1<<(systems-1-site))]
    current=first; sign=1; path=()
    for site in differing[:-1]:
        next_value=current ^ (1<<(systems-1-site))
        if (current>>(systems-1-site))&1: sign*=-1
        path+=adjacent_walsh_ry(current,next_value,np.pi,systems); current=next_value
    last=differing[-1]
    orientation=1 if ((current>>(systems-1-last))&1)==0 else -1
    return path+adjacent_walsh_ry(current,second,2*sign*orientation*angle,systems)+inverse(path)


@lru_cache(maxsize=1)
def walsh_environment_circuit():
    gates=walsh_pattern_ry(1,((0,0),),2*math.asin(1/math.sqrt(3)))
    gates+=walsh_pattern_ry(2,((0,0),(1,0)),np.pi/2)
    for first,second,multiple in NEGATIVE_ROTATIONS:
        gates+=walsh_two_level_rotation(first,second,float(multiple)*np.pi)
    return gates


@lru_cache(maxsize=1)
def walsh_environment():
    return circuit_matrix(walsh_environment_circuit(),4)


@lru_cache(maxsize=1)
def walsh_physical_isometry():
    gates=remap(walsh_environment_circuit(),{0:0,1:4,2:5,3:6})
    return circuit_matrix(gates,7,compiled_sharing(F(1,2))[0])


def walsh_certificate(count,contrast=F(1,2)):
    result=certificate(count,contrast); counts=gate_counts(walsh_environment_circuit())
    reads=result["original_noisy_reads"]
    result.update({"environment_compiler":"commuting pattern-projector expansion",
        "new_yx_after_first_join":counts["yx"]+reads+16,
        "new_ry_after_first_join":counts["ry"]+3*reads+17,
        "new_total_native_gates_after_first_join":counts["all"]+4*reads+33,
        "native_gates_saved_vs_round_147":old_resources()["all"]-counts["all"]})
    return result


class WalshCoherentHistoryCompilerTests(unittest.TestCase):
    def test_native_conjugation_migrates_y_and_leaves_z_on_the_previous_site(self):
        coupling=circuit_matrix((yx(1,0,np.pi/2),),2)
        np.testing.assert_allclose(coupling @ np.kron(PAULI_Y,np.eye(2)) @ coupling.T,
                                   np.kron(PAULI_Z,PAULI_Y),atol=3e-16)

    def test_all_control_patterns_up_to_three_controls_match_the_old_native_compiler(self):
        for count in range(4):
            for mask in range(2**count):
                controls=tuple((j,(mask>>j)&1) for j in range(count))
                for angle in (.37,-.83):
                    actual=circuit_matrix(walsh_pattern_ry(count,controls,angle),count+1)
                    expected=circuit_matrix(pattern_ry(count,controls,angle),count+1)
                    np.testing.assert_allclose(actual,expected,atol=4e-15)

    def test_controlled_pattern_gate_counts_follow_the_exact_subset_sum(self):
        for count in range(1,4):
            counts=gate_counts(walsh_pattern_ry(count,tuple((j,0) for j in range(count)),.3))
            self.assertEqual(counts["yx"],count*2**count)
            self.assertEqual(counts["ry"],2**count)

    def test_signed_gray_rotations_still_match_full_two_level_matrices(self):
        for first,second in ((0,15),(3,12),(9,12),(12,15)):
            actual=circuit_matrix(walsh_two_level_rotation(first,second,-.47),4)
            np.testing.assert_allclose(actual,planar_matrix(first,second,-.47),atol=7e-15)

    def test_entire_environment_and_actual_seven_wire_isometry_are_unchanged(self):
        np.testing.assert_allclose(walsh_environment(),compiled_environment(),atol=2e-14)
        np.testing.assert_allclose(walsh_physical_isometry(),compressed_physical_isometry(),atol=2e-14)

    def test_finite_pointer_recovery_with_new_circuit_has_the_same_channel(self):
        rng=np.random.default_rng(148); raw=rng.normal(size=(16,4))+1j*rng.normal(size=(16,4))
        source=raw @ raw.conj().T; source/=np.trace(source)
        v=walsh_physical_isometry(); whole=v @ source @ v.T; actual=np.zeros_like(source)
        for raw_label in range(8):
            branch=whole
            for site,shift in ((4,2),(5,1),(6,0)):
                branch=physical_read_branch(branch,site,(raw_label>>shift)&1,.5)
            active=keep_systems(branch,(0,1,2,3),7)
            decoder=conditional_unitary(decode_label(raw_label)+1).T
            actual+=decoder @ active @ decoder.T
        expected=channel(recovery_kraus((1-ALPHA/2)/2),source)
        np.testing.assert_allclose(actual,expected,atol=4e-15)

    def test_external_columns_survive_native_compression_and_inverse(self):
        rng=np.random.default_rng(248); columns=rng.normal(size=(16,3))+1j*rng.normal(size=(16,3)); columns/=np.linalg.norm(columns)
        gates=walsh_environment_circuit()
        restored=circuit_matrix(inverse(gates),4,circuit_matrix(gates,4,columns))
        np.testing.assert_allclose(restored,columns,atol=2e-14)

    def test_complete_resource_reduction_and_same_finite_error_certificate(self):
        counts=gate_counts(walsh_environment_circuit())
        self.assertEqual((counts["yx"],counts["ry"],counts["all"]),(274,94,368))
        for m in (3,15,25):
            new,old=walsh_certificate(m),certificate(m)
            self.assertEqual(new["uniform_external_recovery_error_upper"],old["uniform_external_recovery_error_upper"])
            self.assertEqual(old["new_total_native_gates_after_first_join"]-new["new_total_native_gates_after_first_join"],382)
        self.assertEqual(walsh_certificate(3)["new_total_native_gates_after_first_join"],437)
        with self.assertRaises(ValueError): walsh_pattern_ry(0,((0,1),),.1)
        with self.assertRaises(ValueError): adjacent_walsh_ry(0,3,.1)


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(WalshCoherentHistoryCompilerTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":148,"environment_native_counts":gate_counts(walsh_environment_circuit()),
        "previous_environment_native_counts":{k:old_resources()[k] for k in ("yx","ry","all")},
        "environment_native_gates_saved":382,"extra_pure_ancillas_for_compression":0,
        "full_environment_matrix_unchanged":True,"finite_external_error_guarantees_unchanged":True,
        "gate_count_globally_minimal":False,
        "conditional_decoder_counts_unchanged":{"yx_upper":16,"ry_upper":17},
        "environment_gate_list":[{"kind":g.kind,"sites":list(g.sites),"angle_radians":g.angle} for g in walsh_environment_circuit()],
        "certificates":[walsh_certificate(m) for m in (3,15,25)],
        "automated_checks":{"run":checks.testsRun,"failures":len(checks.failures),"errors":len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name("walsh_coherent_history_compiler_results.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    summary={k:v for k,v in report.items() if k not in ("environment_gate_list","certificates")}
    summary["total_native_gates_by_read_count"]={str(row["original_noisy_reads"]):row["new_total_native_gates_after_first_join"] for row in report["certificates"]}
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=="__main__": main()
