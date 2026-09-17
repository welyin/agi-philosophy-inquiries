"""Round 154: execute adaptive label policies with the original native instrument.

All raw results are accepted and retained. Only label wires are queried; the
coherent history bit and closed residual/source systems are never read. The
classical risk certificate bounds the quantum channel but need not be tight.
"""

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import unittest

import numpy as np

from adaptive_history_readout import adaptive_solution, adaptive_certificate, updated
from centered_history_recovery import (GROUPS, centered_decoder, centered_decoder_circuit,
                                      exact_error_interval, purification)
from classical_joining_history import independent_encoding
from coarsened_history_bound import pure_extension_metrics
from compiled_subject_joining import circuit_matrix
from compiled_three_subject_joining import gate_counts
from external_correlation_recovery_bound import channel
from flagged_subject_joining import random_state
from independent_source_alignment import keep_systems
from joint_relation_records import interval_fields
from minimal_coherent_joining_history import conditional_unitary, coherent_decoder_circuit
from noisy_classical_history import physical_read_branch, majority_error
from noisy_two_bit_history import coarse_physical_isometry, coarse_environment_circuit, noisy_error_upper
from quantum_interface_audit import ALPHA
from walsh_coherent_history_compiler import walsh_physical_isometry, walsh_environment_circuit


def policy_confusion(solution,correct_probability):
    model=solution.model; matrix=[[0 for _ in range(model.labels)] for _ in range(model.labels)]
    def visit(depth,signed,true,probability):
        choice=solution.nodes[(depth,signed)][1]
        if depth==model.total:
            matrix[true][choice]+=probability; return
        true_bit=(true>>(model.width-1-choice))&1
        for bit in (0,1):
            factor=correct_probability if bit==true_bit else 1-correct_probability
            visit(depth+1,updated(signed,choice,bit),true,probability*factor)
    for true in range(model.labels): visit(0,(0,)*model.width,true,1)
    return matrix


def policy_kraus(solution,correct_probability,report=None):
    model=solution.model; matrix=policy_confusion(solution,correct_probability)
    if report is not None and (type(report) is not int or report not in range(model.labels)):
        raise ValueError("Report outside the selected label alphabet.")
    result=[]
    for true in range(model.labels):
        fine=GROUPS[true] if model.coarse else (true,)
        for r in fine:
            for g in (range(model.labels) if report is None else (report,)):
                decoder=centered_decoder(g) if model.coarse else conditional_unitary(g+1).T
                result.append(np.sqrt(float(matrix[true][g])/6)*decoder @ conditional_unitary(r+1))
    return tuple(result)


def actual_recovery(solution,source,external_dimension=1):
    model=solution.model
    if type(external_dimension) is not int or external_dimension<1 or external_dimension & (external_dimension-1):
        raise ValueError("Power-of-two external dimension required.")
    if source.shape!=(16*external_dimension,16*external_dimension): raise ValueError("Wrong input dimension.")
    extra=external_dimension.bit_length()-1
    isometry=coarse_physical_isometry() if model.coarse else walsh_physical_isometry()
    v=np.kron(isometry,np.eye(external_dimension)); whole=v @ source @ v.T
    result=np.zeros_like(source)
    def visit(depth,signed,branch):
        nonlocal result
        choice=solution.nodes[(depth,signed)][1]
        if depth==model.total:
            active=keep_systems(branch,(0,1,2,3)+tuple(range(7,7+extra)),7+extra)
            gates=centered_decoder_circuit(choice) if model.coarse else coherent_decoder_circuit(choice+1)
            decoder=np.kron(circuit_matrix(gates,4),np.eye(external_dimension))
            result+=decoder @ active @ decoder.T
            return
        for bit in (0,1):
            next_branch=physical_read_branch(branch,4+choice,bit,float(model.contrast))
            visit(depth+1,updated(signed,choice,bit),next_branch)
    visit(0,(0,)*model.width,whole)
    return result


def report_certificate(solution):
    model=solution.model; matrix=policy_confusion(solution,model.p)
    priors=(F(1,3),F(1,3),F(1,6),F(1,6)) if model.coarse else (F(1,6),)*6
    gains=((1-F(3,2)*exact_error_interval()),)*2+(1,1) if model.coarse else (1,)*6
    reports=[]
    for g in range(model.labels):
        probability=sum(priors[r]*matrix[r][g] for r in range(model.labels))
        if probability.lo<=0: raise ArithmeticError("This certificate requires a report of positive probability.")
        bound=1-priors[g]*gains[g]*matrix[g][g]/probability
        reports.append({"report":g,"probability":interval_fields(probability),
                        "conditional_external_error_upper":interval_fields(bound)})
    return reports


def physical_certificate(reads,coarse=False):
    solution=adaptive_solution(reads,coarse); result=adaptive_certificate(reads,coarse)
    env=gate_counts(coarse_environment_circuit() if coarse else walsh_environment_circuit())
    result.update({"environment_gate_counts":env,"new_pure_environment_helpers":0,
        "new_yx_after_first_join":env["yx"]+reads+16,
        "new_ry_after_first_join":env["ry"]+3*reads+17,
        "total_native_gates_including_the_166_gate_first_join":env["all"]+4*reads+33+166,
        "reports":report_certificate(solution),
        "reachable_policy_states":len(solution.policy_rows()),
        "accessible_coherent_history_rebits":1,"decoded_message_bits":2 if coarse else 3,
        "closed_residual_history_rebits":1 if coarse else 0,
        "all_source_purifications_readout_records_pointers_and_noise_baths_retained":True,
        "ideal_gate_angles_and_coherent_storage_assumed":True,
        "source_preparation_transport_calibration_and_classical_controller_time_included":False})
    return result


class PhysicalAdaptiveHistoryTests(unittest.TestCase):
    def test_complete_policy_confusion_is_normalized_and_matches_the_gain_polynomial(self):
        for coarse,n in ((False,5),(True,4)):
            solution=adaptive_solution(n,coarse); model=solution.model
            rational=policy_confusion(solution,F(3,4))
            self.assertTrue(all(sum(row)==1 for row in rational))
            matrix=policy_confusion(solution,model.p)
            weights=((2-3*exact_error_interval())/6,)*2+(F(1,6),)*2 if coarse else (F(1,6),)*6
            direct=1-sum(weights[g]*matrix[g][g] for g in range(model.labels))
            self.assertLessEqual(direct.lo,solution.error_interval.hi)
            self.assertGreaterEqual(direct.hi,solution.error_interval.lo)

    def test_five_actual_native_reads_and_feedback_match_the_six_label_external_channel(self):
        solution=adaptive_solution(5); rng=np.random.default_rng(154)
        raw=rng.normal(size=(32,3))+1j*rng.normal(size=(32,3)); state=raw @ raw.conj().T; state/=np.trace(state)
        expected=channel(tuple(np.kron(k,np.eye(2)) for k in policy_kraus(solution,(1+ALPHA/2)/2)),state)
        np.testing.assert_allclose(actual_recovery(solution,state,2),expected,atol=6e-15)

    def test_four_actual_native_reads_leave_the_coarse_residual_unread(self):
        solution=adaptive_solution(4,True); rng=np.random.default_rng(254)
        raw=rng.normal(size=(32,4))+1j*rng.normal(size=(32,4)); state=raw @ raw.conj().T; state/=np.trace(state)
        expected=channel(tuple(np.kron(k,np.eye(2)) for k in policy_kraus(solution,(1+ALPHA/2)/2)),state)
        np.testing.assert_allclose(actual_recovery(solution,state,2),expected,atol=6e-15)
        self.assertTrue(all(row.get("physical_wire",4) in (4,5) for row in solution.policy_rows()))

    def test_report_probabilities_are_input_independent_and_conditional_bounds_are_honest(self):
        for coarse,n in ((False,5),(True,4)):
            solution=adaptive_solution(n,coarse); p=(1+ALPHA/2)/2
            for report in report_certificate(solution):
                ks=policy_kraus(solution,p,report["report"]); probability=report["probability"]["diagnostic"]
                np.testing.assert_allclose(sum(k.T @ k for k in ks),probability*np.eye(16),atol=6e-16)
                metric=pure_extension_metrics(tuple(k/np.sqrt(probability) for k in ks),np.eye(16)/4)
                self.assertLessEqual(metric["trace_error"],report["conditional_external_error_upper"]["diagnostic"]+1e-14)

    def test_average_bounds_protect_unknown_external_states_under_the_declared_promises(self):
        rng=np.random.default_rng(354)
        for coarse,n in ((False,5),(True,4)):
            solution=adaptive_solution(n,coarse); ks=policy_kraus(solution,(1+ALPHA/2)/2)
            np.testing.assert_allclose(sum(k.T @ k for k in ks),np.eye(16),atol=2e-15)
            for _ in range(4):
                if coarse:
                    omega=independent_encoding((random_state(rng,2),random_state(rng,2))); c=purification(omega)
                else:
                    c=rng.normal(size=(16,5))+1j*rng.normal(size=(16,5)); c/=np.linalg.norm(c)
                self.assertLessEqual(pure_extension_metrics(ks,c)["trace_error"],sum(solution.error_interval.floats())/2+1e-14)

    def test_five_read_six_label_policy_improves_the_previous_six_read_four_label_certificate(self):
        result=physical_certificate(5)
        self.assertEqual((result["new_yx_after_first_join"],result["new_ry_after_first_join"]),(295,126))
        self.assertEqual(result["new_native_gates_after_first_join"],421)
        self.assertEqual(result["total_native_gates_including_the_166_gate_first_join"],587)
        self.assertLess(adaptive_solution(5).error_interval.hi,noisy_error_upper(majority_error(3)).lo)
        self.assertEqual(physical_certificate(4,True)["new_native_gates_after_first_join"],421)

    def test_current_equal_budget_comparison_and_raw_record_storage_are_explicit(self):
        six=physical_certificate(7); four=physical_certificate(6,True)
        self.assertEqual((six["new_native_gates_after_first_join"],four["new_native_gates_after_first_join"]),(429,429))
        self.assertEqual((six["fresh_pointer_initializations"],four["fresh_pointer_initializations"]),(7,6))
        self.assertLess(adaptive_solution(7).error_interval.hi,adaptive_solution(6,True).error_interval.lo)
        self.assertFalse(six["source_preparation_transport_calibration_and_classical_controller_time_included"])

    def test_dimension_and_report_guards(self):
        solution=adaptive_solution(5)
        with self.assertRaises(ValueError): actual_recovery(solution,np.eye(16)/16,3)
        with self.assertRaises(ValueError): actual_recovery(solution,np.eye(8)/8)
        with self.assertRaises(ValueError): policy_kraus(solution,.75,6)


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(PhysicalAdaptiveHistoryTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":154,"certificates":[physical_certificate(n,c) for n,c in ((5,False),(4,True),(7,False),(6,True))],
        "read_policy_artifact":"adaptive_history_readout_results.json",
        "actual_raw_outcome_trees_checked_with_unknown_external_system":True,
        "quantum_trace_distance_minimax_or_minimum_read_count_proved":False,
        "next_question":"replace the label-risk surrogate by recovery-sensitive quantum error costs",
        "automated_checks":{"run":checks.testsRun,"failures":len(checks.failures),"errors":len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name("physical_adaptive_history_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"round":154,"checks":report["automated_checks"],"selected":[
        {"labels":r["labels"],"reads":r["reads"],"bound":r["external_error_upper"]["diagnostic"],
         "gates":r["new_native_gates_after_first_join"],"policy_states":r["reachable_policy_states"],
         "worst_report_bound":max(x["conditional_external_error_upper"]["diagnostic"] for x in r["reports"])}
        for r in report["certificates"]]},indent=2))


if __name__=="__main__": main()
