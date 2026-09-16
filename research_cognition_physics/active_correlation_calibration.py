"""Round 113: amplify calibration sensitivity with controls on fresh source copies."""
import argparse
import json
import math
import unittest
from fractions import Fraction as F
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I, SCALE
from correlated_environment_response import correlated_accessible
from environment_response_summary import probe_state
from nonzero_fragment_records import biased_environment
from record_correlation_frontier import biased_probe_records
from joint_relation_records import interval_fields
from unknown_correlation_learning import (calibration_visibility, exp_minus_interval,
    scheduled_radius, source_mean_interval)
from weak_relation_tradeoff import trace_distance


def one_sided_failure_bound(n, look, kappa, response_strength, target=F(1,5)):
    """Bound failure to certify the true sign at this look, for |kappa|>target."""
    kappa,a,target=map(F,(kappa,response_strength,target))
    if not target < abs(kappa) <= 1 or not 0 < a <= 1:
        raise ValueError("Positive distance beyond the activation threshold required.")
    gap=calibration_visibility(1)*I.exact(a*(abs(kappa)-target))-scheduled_radius(n,look)
    if gap.lo <= 0:
        return I.exact(1)
    # The acceptance threshold is rounded upward in the actual finite controller.
    # Reserve 4/n for integer lattice and outward interval rounding.
    conservative=F(gap.lo,SCALE)-F(4,n)
    if conservative <= 0:
        return I.exact(1)
    low_tail=exp_minus_interval(F(n,2)*conservative*conservative)
    # Very high records can yield an empty physical intersection and trigger fallback.
    # For any physical truth that upper-domain event lies beyond the radius.
    upper_conflict=exp_minus_interval(6+look)
    result=low_tail+upper_conflict
    return I(min(result.lo,SCALE),min(result.hi,SCALE))


def calibration_design_certificate(kappa=F(64,125), target=F(1,5), max_looks=12):
    output=[]
    for a in (F(19,100),F(1)):
        first=None
        for look in range(1,max_looks+1):
            n=64*2**(look-1)
            fail=one_sided_failure_bound(n,look,kappa,a,target)
            if fail.hi < SCALE//100:
                first={"look":look,"independent_samples":n,
                    "failure_to_activate_correct_sign_upper":interval_fields(fail),
                    "source_and_probe_pair_rotations":34*n,
                    "source_and_probe_extra_local_rotations":41*n,
                    "fresh_pure_slots_allocated_and_retained":7*n}
                break
        output.append({"response_strength_exact":str(a),
            "single_original_read_mean":interval_fields(source_mean_interval(kappa,a)),
            "first_scheduled_look_certifying_at_least_99_percent_correct_activation":first})
    return output


def single_source_comparison(kappa1,kappa2,count=1):
    delta=abs(F(kappa1)-F(kappa2))
    ideal=I.exact(delta)/2
    return {"source_trace_distance":interval_fields(ideal),
        "accessible_trace_distance_at_lambda1_lambda2_zero":interval_fields(ideal),
        "finite_majority_record_total_variation":interval_fields(calibration_visibility(count)*ideal),
        "finite_readout_optimal_over_all_possible_protocols":False}


class ActiveCorrelationCalibrationTests(unittest.TestCase):
    def test_full_source_distance_is_the_control_independent_ceiling(self):
        for biases in ((.6,)*3,(.2,.4,.7)):
            for k1,k2 in ((-.3,.3),(0.,.2)):
                actual=trace_distance(biased_environment(biases,k1),biased_environment(biases,k2))
                self.assertAlmostEqual(actual,abs(k1-k2)/2,places=14)

    def test_zero_unavailable_cosines_saturate_source_distinguishability(self):
        for lam0 in (0.,.5,1.):
            states=[correlated_accessible(probe_state(.8),.8,biased_environment((.6,)*3,k), (lam0,0.,0.),(0,),True)
                    for k in (-.512,.512)]
            self.assertAlmostEqual(trace_distance(*states),.512,places=13)

    def test_any_original_strength_is_below_the_fresh_control_ceiling(self):
        rng=np.random.default_rng(113)
        for _ in range(15):
            overlaps=tuple(rng.uniform(0,1,3))
            states=[correlated_accessible(probe_state(.8),.8,biased_environment((.6,)*3,k),overlaps,(0,),True)
                    for k in (-.3,.4)]
            expected=.35*math.sqrt((1-overlaps[1]**2)*(1-overlaps[2]**2))
            self.assertAlmostEqual(trace_distance(*states),expected,places=13)
            self.assertLessEqual(expected,.35)

    def test_actual_original_readouts_saturate_the_claimed_finite_probe_value(self):
        for count in (1,3):
            pair=[biased_probe_records(.8,(.6,)*3,k,(.9,0.,0.),count) for k in (-.512,.512)]
            actual=sum(abs(np.trace(pair[0][h]-pair[1][h]).real) for h in pair[0])/2
            expected=single_source_comparison(-F(64,125),F(64,125),count)
            self.assertAlmostEqual(actual,expected["finite_majority_record_total_variation"]["diagnostic"],places=13)

    def test_calibration_setting_is_separate_from_the_deployment_setting(self):
        old_signal=source_mean_interval(F(64,125),F(19,100))
        new_signal=source_mean_interval(F(64,125),F(1))
        self.assertAlmostEqual(new_signal.floats()[0]/old_signal.floats()[0],100/19,places=13)
        # Statistical sample scaling at a fixed radius is a^-2, not a^-1.
        self.assertEqual(F(1)/F(19,100)**2,F(10000,361))

    def test_finite_certified_power_improves_with_active_calibration(self):
        rows=calibration_design_certificate()
        weak,strong=[r["first_scheduled_look_certifying_at_least_99_percent_correct_activation"] for r in rows]
        self.assertIsNotNone(weak)
        self.assertIsNotNone(strong)
        self.assertLess(strong["independent_samples"],weak["independent_samples"])
        for result in (weak,strong):
            self.assertLess(F(result["failure_to_activate_correct_sign_upper"]["upper_exact"]),F(1,100))

    def test_no_extra_information_at_zero_hidden_parameter(self):
        records=biased_probe_records(.8,(.6,)*3,0.,(.9,0.,0.),1)
        for branch in records.values():
            self.assertAlmostEqual(np.trace(branch).real,.5,places=14)

    def test_margin_and_sample_requirements_are_explicit(self):
        with self.assertRaises(ValueError): one_sided_failure_bound(64,1,F(1,10),1)
        result=one_sided_failure_bound(64,1,F(64,125),1)
        self.assertEqual(result.hi,SCALE)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ActiveCorrelationCalibrationTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":113,
        "control_choice":"Fresh calibration: lambda1=lambda2=0; deployment retains lambda1=lambda2=sqrt(c)",
        "ideal_source_pair_distance_ceiling":"abs(delta_kappa)/2, attained by the accessible real probe",
        "comparison":single_source_comparison(-F(64,125),F(64,125)),
        "certified_sample_designs":calibration_design_certificate(),
        "source_assumption":"Same unknown fixed kappa, independent new preparations; provider's purification controls are not revealed to controller",
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("active_correlation_calibration_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
