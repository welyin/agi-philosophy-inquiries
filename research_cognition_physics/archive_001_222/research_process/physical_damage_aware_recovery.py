"""Round 157: three original reads already have a promised error below one half.

The three-read cost-optimal protocol equals the old once-per-bit protocol.
The improved guarantee comes from sharp conditional costs on promised inputs.
The sum-of-costs bound still has demonstrable slack for mixtures of errors.
"""

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import unittest

import numpy as np

from adaptive_history_readout import updated
from centered_history_recovery import purification
from certified_intervals import Interval as I, SCALE
from classical_joining_history import independent_encoding
from coarsened_history_bound import pure_extension_metrics
from damage_aware_history_readout import damage_solution, damage_certificate
from external_correlation_recovery_bound import channel
from flagged_subject_joining import random_state
from history_recovery_loss_matrix import cost_interval, symbols, conditional_kraus, attaining_states
from joint_relation_records import interval_fields
from noisy_coherent_history_recovery import confusion_matrix, decode_label, recovery_kraus
from physical_adaptive_history import policy_confusion, policy_kraus, actual_recovery
from quantum_interface_audit import ALPHA


def severe_and_reduced_probabilities(error):
    return ((4*error-5*error**2+4*error**3)/6,(4*error-3*error**2)/2)


def once_per_bit_bound(error):
    severe,reduced=severe_and_reduced_probabilities(error)
    return severe+cost_interval("s")*reduced


def conditional_report_certificates(solution):
    model=solution.model; matrix=policy_confusion(solution,model.p)
    priors=(F(1,3),F(1,3),F(1,6),F(1,6)) if model.coarse else (F(1,6),)*6
    result=[]
    for report in range(model.labels):
        probability=sum(priors[true]*matrix[true][report] for true in range(model.labels))
        if not isinstance(probability,I): probability=I.exact(probability)
        if probability.hi==0: continue
        if probability.lo<=0: raise ArithmeticError("Unresolved zero-probability report.")
        loss=sum(priors[true]*matrix[true][report]*cost_interval(symbols(model.coarse)[true][report])
                 for true in range(model.labels))
        result.append({"report":report,"probability":interval_fields(probability),
                       "conditional_promised_error_upper":interval_fields(loss/probability)})
    return result


def physical_damage_certificate(reads,coarse=False):
    solution=damage_solution(reads,coarse); result=damage_certificate(reads,coarse)
    result.update({"new_yx_after_first_join":(292 if coarse else 290)+reads,
        "new_ry_after_first_join":(113 if coarse else 111)+3*reads,
        "total_gates_including_166_gate_first_join":(571 if coarse else 567)+4*reads,
        "reports":conditional_report_certificates(solution),
        "accessible_coherent_history_rebits":1,"decoded_message_bits":2 if coarse else 3,
        "fresh_pointer_initializations":reads,"raw_readout_bits_retained":reads,
        "all_results_accepted":True,"source_purification_and_noise_baths_closed":True,
        "all_environment_and_readout_systems_retained":True,
        "ideal_gate_angles_and_coherent_storage_assumed":True,
        "classical_compute_transport_calibration_and_source_costs_included":False,
        "only_sum_of_sharp_conditional_costs_optimized":True})
    return result


class PhysicalDamageAwareRecoveryTests(unittest.TestCase):
    def test_three_read_policy_is_exactly_the_existing_once_per_bit_protocol(self):
        solution=damage_solution(3)
        for row in solution.policy_rows():
            if row["depth"]<3: self.assertEqual(row["read_label_wire"],row["depth"])
        for raw in range(8):
            signed=(0,0,0)
            for depth in range(3): signed=updated(signed,depth,(raw>>(2-depth))&1)
            self.assertEqual(solution.nodes[(3,signed)][1],decode_label(raw))

    def test_confusion_and_cost_probability_polynomials_match_exact_rational_enumeration(self):
        for error in (F(0),F(1,10),F(1,4),F(1,2)):
            matrix=policy_confusion(damage_solution(3),1-error)
            self.assertEqual(matrix,confusion_matrix(error))
            severe=sum(matrix[r][g]/6 for r in range(6) for g in range(6) if symbols()[r][g]=="1")
            reduced=sum(matrix[r][g]/6 for r in range(6) for g in range(6) if symbols()[r][g]=="s")
            self.assertEqual((severe,reduced),severe_and_reduced_probabilities(error))

    def test_actual_original_instrument_and_external_channel_are_unchanged_at_three_reads(self):
        rng=np.random.default_rng(157); raw=rng.normal(size=(32,3))+1j*rng.normal(size=(32,3))
        source=raw @ raw.conj().T; source/=np.trace(source)
        error=(1-ALPHA/2)/2
        expected=channel(tuple(np.kron(k,np.eye(2)) for k in recovery_kraus(error)),source)
        np.testing.assert_allclose(actual_recovery(damage_solution(3),source,2),expected,atol=5e-15)

    def test_new_five_read_and_four_label_policies_also_match_native_feedback_trees(self):
        rng=np.random.default_rng(257)
        for reads,coarse in ((5,False),(4,True)):
            solution=damage_solution(reads,coarse)
            raw=rng.normal(size=(32,4))+1j*rng.normal(size=(32,4)); source=raw @ raw.conj().T; source/=np.trace(source)
            expected=channel(tuple(np.kron(k,np.eye(2)) for k in policy_kraus(solution,(1+ALPHA/2)/2)),source)
            np.testing.assert_allclose(actual_recovery(solution,source,2),expected,atol=6e-15)

    def test_unknown_promised_purifications_and_every_report_obey_the_new_bounds(self):
        rng=np.random.default_rng(357)
        for reads,coarse in ((3,False),(5,False),(4,True)):
            solution=damage_solution(reads,coarse); p=(1+ALPHA/2)/2; ks=policy_kraus(solution,p)
            for _ in range(3):
                states=(random_state(rng,2),random_state(rng,2)); c=purification(independent_encoding(states))
                self.assertLessEqual(pure_extension_metrics(ks,c)["trace_error"],sum(solution.error_interval.floats())/2+1e-14)
                for row in conditional_report_certificates(solution):
                    report_ks=policy_kraus(solution,p,row["report"]); probability=row["probability"]["diagnostic"]
                    np.testing.assert_allclose(sum(k.T @ k for k in report_ks),probability*np.eye(16),atol=7e-16)
                    actual=pure_extension_metrics(tuple(k/np.sqrt(probability) for k in report_ks),c)["trace_error"]
                    self.assertLessEqual(actual,row["conditional_promised_error_upper"]["diagnostic"]+1e-14)

    def test_three_read_interval_and_full_native_resource_account(self):
        solution=damage_solution(3); formula=once_per_bit_bound(1-solution.model.p)
        self.assertLessEqual(formula.lo,solution.error_interval.hi); self.assertGreaterEqual(formula.hi,solution.error_interval.lo)
        self.assertLess(formula.hi,SCALE//2)
        report=physical_damage_certificate(3)
        self.assertEqual((report["new_yx_after_first_join"],report["new_ry_after_first_join"]),(293,120))
        self.assertEqual(report["new_native_gates_after_first_join"],413)
        self.assertEqual(report["total_gates_including_166_gate_first_join"],579)
        self.assertEqual(report["fresh_pointer_initializations"],3)

    def test_separately_sharp_costs_can_be_strictly_loose_after_mixing_errors(self):
        # True A_Y; choose A_X or A_Z with equal probability: each edge has cost 1.
        mixture=tuple(conditional_kraus(1,g)[0]/np.sqrt(2) for g in (0,2))
        known_pair=conditional_kraus(0,2,True)
        for actual,expected in zip(mixture,known_pair): np.testing.assert_allclose(actual,expected,atol=3e-16)
        c=purification(independent_encoding(attaining_states(0,2,True)))
        self.assertAlmostEqual(pure_extension_metrics(mixture,c)["trace_error"],sum(cost_interval("b").floats())/2)
        self.assertLess(cost_interval("b").hi,SCALE)

    def test_zero_probability_reports_are_excluded_without_postselection(self):
        reports=conditional_report_certificates(damage_solution(0))
        self.assertEqual(len(reports),1)
        self.assertEqual(reports[0]["probability"]["diagnostic"],1.)
        with self.assertRaises(ValueError): damage_solution(-1)


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(PhysicalDamageAwareRecoveryTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    solution=damage_solution(3)
    report={"round":157,"three_read_protocol_identical_to_original_once_per_label_bit_protocol":True,
        "three_read_bound_formula":"(4e-5e^2+4e^3)/6 + (sqrt(3)/2)*(4e-3e^2)/2",
        "three_read_legal_maximally_mixed_input_purification_diagnostic":
            pure_extension_metrics(policy_kraus(solution,(1+ALPHA/2)/2),np.eye(16)/4),
        "diagnostic_is_not_a_worst_case_certificate":True,
        "certificates":[physical_damage_certificate(n,c) for n,c in ((3,False),(5,False),(4,True),(7,False))],
        "strict_slack_example":{"true_label":"A_Y","equal_reports":["A_X","A_Z"],
            "sum_of_sharp_costs":1,"actual_sharp_mixture_cost":interval_fields(cost_interval("b"))},
        "next_question":"optimize full averaged recovery error with one common unknown source state",
        "automated_checks":{"run":checks.testsRun,"failures":len(checks.failures),"errors":len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name("physical_damage_aware_recovery_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"round":157,"checks":report["automated_checks"],"selected":[
        {"labels":r["labels"],"reads":r["reads"],"bound":r["new_policy_with_sharp_costs"]["diagnostic"],
         "gates":r["new_native_gates_after_first_join"],
         "worst_report_bound":max(x["conditional_promised_error_upper"]["diagnostic"] for x in r["reports"])}
        for r in report["certificates"]]},indent=2))


if __name__=="__main__": main()
