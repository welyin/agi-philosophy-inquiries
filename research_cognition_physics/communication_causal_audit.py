"""Round 73: identical declared records can have different message dependence.

Match round 69's 288 records, plus Bob's retained flag, with calibrated real
feedforward. The matched interface omits trusted transport logs and internal
pointer/randomness histories. Blocking the message distinguishes the models.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from bell_network_statistics import network_score
from certified_intervals import Interval as I, SCALE, sin_interval
from finite_network_separation import finite_bias, finite_protocol_table
from independent_source_alignment import real_x_with_reset
from one_bit_real_network import forget_flag, records_for_reference, reference_flags_by_circuit, visibility_interval, real_feedback_score
from operational_effect_closure import majority_certificate


def calibration_interval(level=5,alignment_count=5):
    alpha=4*sin_interval(I.rational(1,4))
    bias=I.rational(1,5)
    for _ in range(level): bias=bias*(1+alpha)/(1+alpha*bias**2)
    attenuation=bias**2/visibility_interval(alignment_count)
    flip=(1-attenuation)/2
    if flip.lo<0 or flip.hi>SCALE//2: raise ValueError("The real alignment is too weak for this calibration.")
    return bias,attenuation,flip


def calibrated_flagged_records(level=5,alignment_count=5,read_count=5,deliver_message=True):
    calibration_interval(level,alignment_count)
    alignment=majority_certificate(alignment_count)["effective_visibility_diagnostic"]
    read=majority_certificate(read_count)["effective_visibility_diagnostic"]
    flip=(1-finite_bias(level)**2/alignment)/2
    result={}
    for flag,raw in reference_flags_by_circuit(alignment_count).items():
        corrected=real_x_with_reset(raw,1,2) if deliver_message and flag<0 else raw
        calibrated=(1-flip)*corrected+flip*real_x_with_reset(corrected,1,2)
        result.update({(flag,)+key:value for key,value in records_for_reference(calibrated,read,read,read).items()})
    return result


def nonsignalling_residual(flagged):
    """Audit the full declared (h,a,b,c) outputs against outer setting changes."""
    differences=[]
    for h,a,b,x in product((-1,1),(-1,1),range(4),range(3)):
        values=[sum(flagged[h,x,z,a,b,c] for c in (-1,1)) for z in range(6)]
        differences.append(max(values)-min(values))
    for h,b,c,z in product((-1,1),range(4),(-1,1),range(6)):
        values=[sum(flagged[h,x,z,a,b,c] for a in (-1,1)) for x in range(3)]
        differences.append(max(values)-min(values))
    return max(differences)


def comparison_report():
    real=calibrated_flagged_records()
    cut=calibrated_flagged_records(deliver_message=False)
    complex_records=finite_protocol_table()
    max_error=max(abs(value-complex_records[key[1:]]/2) for key,value in real.items())
    _,_,flip=calibration_interval()
    coin_numerator=(flip.lo+flip.hi)//2
    coin=Fraction(coin_numerator,SCALE)
    coin_error=Fraction(max(coin_numerator-flip.lo,flip.hi-coin_numerator),SCALE)
    read=visibility_interval(5)
    bias,_,_=calibration_interval()
    intervention_gap=2*I.exact(2).sqrt()*read**3*bias**2
    return {"calibration_flip_interval":flip.floats(),
            "calibration_flip_lower_exact":str(Fraction(flip.lo,SCALE)),
            "calibration_flip_upper_exact":str(Fraction(flip.hi,SCALE)),
            "finite_fair_coin_bits":100,"finite_coin_probability_exact":str(coin),
            "finite_coin_record_TV_error_upper_exact":str(coin_error),
            "finite_coin_score_error_upper_exact":str(24*coin_error),
            "matched_declared_probability_entries":len(real),"max_probability_error_diagnostic":max_error,
            "exact_matching_condition":"reference_YY_correlation = round_69_private_Y_bias^2",
            "nonsignalling_residual_diagnostic":nonsignalling_residual(real),
            "real_with_message_score_diagnostic":network_score(forget_flag(real)),
            "real_message_cut_score_diagnostic":network_score(forget_flag(cut)),
            "complex_without_message_score_diagnostic":network_score(complex_records),
            "intervention_score_gap_interval":intervention_gap.floats(),
            "intervention_score_gap_lower_exact":str(Fraction(intervention_gap.lo,SCALE)),
            "matched_interface":"settings x,z; outer results a,c; Bob Bell label b and locally retained flag h",
            "matched_interface_includes_trusted_transport_logs":False,
            "matched_interface_includes_internal_pointer_and_randomness_histories":False,
            "nonsignalling_implies_no_communication":False,
            "passive_record_prediction_selects_complex_model":False}


class CommunicationCausalAuditTests(unittest.TestCase):
    def test_calibration_is_a_strictly_interior_finite_real_randomization(self):
        _,attenuation,flip=calibration_interval()
        self.assertGreater(flip.lo,0)
        self.assertLess(flip.hi,I.rational(3,100000).lo)
        self.assertGreater(attenuation.lo,0)
        self.assertLess(attenuation.hi,SCALE)

    def test_calibrated_real_model_matches_every_complex_record_and_retained_flag(self):
        real=calibrated_flagged_records()
        target=finite_protocol_table()
        self.assertEqual(len(real),576)
        for key,value in real.items(): self.assertAlmostEqual(value,target[key[1:]]/2,places=14)
        np.testing.assert_allclose(list(forget_flag(real).values()),list(target.values()),atol=1e-15)

    def test_declared_flag_is_uniform_and_statistically_independent_of_retained_network_outputs(self):
        real=calibrated_flagged_records()
        coarse=forget_flag(real)
        for key,value in real.items(): self.assertAlmostEqual(value,coarse[key[1:]]/2,places=14)
        for x,z in product(range(3),range(6)):
            for h in (-1,1):
                self.assertAlmostEqual(sum(real[h,x,z,a,b,c] for a,b,c in product((-1,1),range(4),(-1,1))),.5,places=14)

    def test_both_communicating_and_cut_models_pass_declared_nonsignalling_tests(self):
        for deliver in (False,True): self.assertLess(nonsignalling_residual(calibrated_flagged_records(deliver_message=deliver)),1e-14)

    def test_cutting_the_message_has_the_predicted_large_effect_only_on_real_feedback(self):
        report=comparison_report()
        self.assertAlmostEqual(report["real_with_message_score_diagnostic"],report["complex_without_message_score_diagnostic"],places=12)
        read=majority_certificate(5)["effective_visibility_diagnostic"]
        self.assertAlmostEqual(report["real_message_cut_score_diagnostic"],real_feedback_score(0,read,read,read),places=12)
        self.assertGreater(Fraction(report["intervention_score_gap_lower_exact"]),Fraction(2828,1000))
        self.assertAlmostEqual(report["real_with_message_score_diagnostic"]-report["real_message_cut_score_diagnostic"],report["intervention_score_gap_interval"][0],places=12)

    def test_all_positive_record_likelihood_ratios_equal_one_before_intervention(self):
        # The analytic equality in the note implies equality for any finite
        # adaptive setting policy under fresh repeated uses of this interface.
        real=calibrated_flagged_records()
        target=finite_protocol_table()
        ratios=[value/(target[key[1:]]/2) for key,value in real.items()]
        # Small-probability ratios magnify full-circuit floating accumulation;
        # exact equality is obtained algebraically from the common YY coefficient.
        np.testing.assert_allclose(ratios,np.ones(576),rtol=1e-13,atol=0)

    def test_one_hundred_fair_bits_give_a_certified_finite_precision_coin(self):
        report=comparison_report()
        _,_,interval=calibration_interval()
        coin=Fraction(report["finite_coin_probability_exact"])
        error=Fraction(report["finite_coin_record_TV_error_upper_exact"])
        self.assertLess(error,Fraction(1,10**24))
        self.assertLessEqual(abs(coin-Fraction(interval.lo,SCALE)),error)
        self.assertLessEqual(abs(coin-Fraction(interval.hi,SCALE)),error)
        self.assertEqual((coin*SCALE).denominator,1)
        self.assertEqual(Fraction(report["finite_coin_score_error_upper_exact"]),24*error)

    def test_no_network_results_are_discarded_when_the_message_is_cut(self):
        records=calibrated_flagged_records(deliver_message=False)
        self.assertGreaterEqual(min(records.values()),0)
        for x,z in product(range(3),range(6)):
            total=sum(records[h,x,z,a,b,c] for h,a,b,c in product((-1,1),(-1,1),range(4),(-1,1)))
            self.assertAlmostEqual(total,1.,places=14)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CommunicationCausalAuditTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":73,**comparison_report(),"scope":"Exact equality of a specified coarse-grained record interface; causal interventions and trusted transport records are additional data",
            "quantum_theory_derived_from_cognition":False,"automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("communication_causal_audit_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__": main()
