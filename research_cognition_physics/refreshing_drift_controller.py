"""Round 117: finite refresh, observed drift audits, and safe fallback."""
import argparse
import json
import unittest
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I, SCALE, ceil_div
from conditional_drift_windows import predictive_interval
from drifting_source_transfer import CAPACITY, robust_action
from incompatible_relation_writes import H, calibration_state
from joint_relation_records import interval_fields
from learned_correlation_loop import learned_query_records
from unknown_correlation_learning import (source_mean_interval, confidence_set,
    matched_policy_value, uniform_risk_certificate)


def count_distribution(probabilities):
    """Certified Poisson-binomial count distribution; no identical-mean assumption."""
    lower=[SCALE];upper=[SCALE]
    for raw in probabilities:
        p=I.exact(raw);q=1-p
        if p.lo<0 or p.hi>SCALE: raise ValueError("Invalid Bernoulli interval.")
        size=len(lower);lo=[0]*(size+1);hi=[0]*(size+1)
        for k in range(size):
            lo[k]+=lower[k]*q.lo//SCALE
            hi[k]+=ceil_div(upper[k]*q.hi,SCALE)
            lo[k+1]+=lower[k]*p.lo//SCALE
            hi[k+1]+=ceil_div(upper[k]*p.hi,SCALE)
        lower,upper=lo,hi
    return [I(lo,hi) for lo,hi in zip(lower,upper)]


def source_path(index,rate=F(1,4096)):
    if not isinstance(index,int) or index<1 or F(rate)<0:
        raise ValueError("Positive source index and nonnegative rate required.")
    return max(-CAPACITY,min(CAPACITY,CAPACITY-F(rate)*(index-1)))


def audit_intervals(first,second,center_separation,rate):
    """Reject if two average-parameter confidence sets cannot satisfy the rate bound."""
    rate=F(rate);separation=F(center_separation)
    if rate<0 or separation<0: raise ValueError("Invalid audit inputs.")
    first,second=I.exact(first),I.exact(second)
    minimum_gap=max(0,second.lo-first.hi,first.lo-second.hi)
    budget=I.exact(rate*separation)
    return minimum_gap>budget.hi


def update_trusted_counts(previous_counts,previous_intervals,current_counts,current_intervals,budget):
    """An observed model conflict locks future enhanced queries to fallback."""
    if previous_counts is None:
        return [mass if not conflict else I.exact(0)
                for mass,(interval,conflict) in zip(current_counts,current_intervals)]
    output=[]
    for current_mass,(current_interval,current_conflict) in zip(current_counts,current_intervals):
        if current_conflict:
            output.append(I.exact(0));continue
        compatible_lo=compatible_hi=0
        for old_mass,(old_interval,old_conflict) in zip(previous_counts,previous_intervals):
            gap=max(0,current_interval.lo-old_interval.hi,old_interval.lo-current_interval.hi)
            if not old_conflict and gap<=budget.hi:
                compatible_lo+=old_mass.lo;compatible_hi+=old_mass.hi
        output.append(current_mass*I(compatible_lo,compatible_hi))
    return output


def complement_probability(value):
    value=1-value
    return I(max(0,value.lo),min(SCALE,value.hi))


@lru_cache(maxsize=None)
def refresh_certificate(window=512,blocks=8,rate=F(1,4096),target=F(1,10)):
    rate,target=F(rate),F(target)
    if not isinstance(window,int) or window<1 or not isinstance(blocks,int) or blocks<1:
        raise ValueError("Positive finite block sizes required.")
    if rate<0 or not 0<target<=CAPACITY: raise ValueError("Invalid drift or activation margin.")
    base,coefficient=matched_policy_value(0,0)
    rows=[];trusted=None;old_intervals=None;first_signal=None
    total_gain=I.exact(0);stale_gain=I.exact(0)
    for look in range(1,blocks+1):
        start=(look-1)*(window+1)+1
        kappas=[source_path(start+i,rate) for i in range(window)]
        current=count_distribution([(1+source_mean_interval(k,1))/2 for k in kappas])
        intervals=[confidence_set(window,k,look,CAPACITY,1) for k in range(window+1)]
        trusted=update_trusted_counts(trusted,old_intervals,current,intervals,I.exact(rate*(window+1)))
        plus=I.exact(0);minus=I.exact(0)
        for positives,mass in enumerate(trusted):
            future,conflict=predictive_interval(window,positives,look,rate)
            action=0 if conflict else robust_action(future,target)
            if action==1: plus+=mass
            if action==-1: minus+=mass
        total_trust=sum(trusted,I.exact(0))
        locked=complement_probability(total_trust)
        fallback=complement_probability(plus+minus)
        query_index=start+window
        truth=source_path(query_index,rate)
        signal=plus-minus
        if first_signal is None: first_signal=signal
        gain=coefficient*I.exact(truth)*signal
        stale=coefficient*I.exact(truth)*first_signal
        total_gain+=gain;stale_gain+=stale
        rows.append({"look":look,"calibration_first_source_index":start,"query_source_index":query_index,
            "true_query_kappa_for_evaluation_exact":str(truth),
            "positive_action_probability":interval_fields(plus),"negative_action_probability":interval_fields(minus),
            "fallback_probability_including_locked_audit":interval_fields(fallback),
            "cumulative_audit_lock_probability":interval_fields(locked),
            "refreshed_H_success":interval_fields(base+gain),
            "stale_first_policy_H_success":interval_fields(base+stale)})
        old_intervals=intervals
    return {"window":window,"blocks":blocks,"rate_per_source_exact":str(rate),"target_exact":str(target),
        "source_path":"clip(R-(source_index-1)*rate,-R,R); query preparations also advance the index",
        "controller_uses_true_path_values":False,
        "auditor":"Compare adjacent average-kappa intervals; conflict locks all later enhanced queries to fallback.",
        "rows":rows,"fallback_success":interval_fields(base),
        "mean_refreshed_H_success":interval_fields(base+total_gain/blocks),
        "mean_stale_first_policy_H_success":interval_fields(base+stale_gain/blocks),
        "total_expected_extra_correct_H_answers":interval_fields(total_gain),
        "uniform_confidence_failure_upper":interval_fields(uniform_risk_certificate()),
        "resources":{"calibration_preparations":window*blocks,"query_preparations":blocks,
            "pair_rotations_upper":34*window*blocks+38*blocks,
            "extra_local_rotations_upper":41*window*blocks+46*blocks,
            "fresh_pure_slots_allocated_and_retained":7*window*blocks+8*blocks,
            "actual_readouts":window*blocks+2*blocks,
            "calibration_result_bits":window*blocks,"policy_bits":2*blocks,"query_return_bits":blocks,
            "classical_arithmetic_and_transport_cost_modelled":False,
            "calibration_continues_on_fixed_schedule_after_audit_lock":True},
        "scope":"Independent but nonidentical source copies for this numerical path; conditional-mean theorem itself permits dependence."}


def audit_example():
    n=512
    first,conflict1=confidence_set(n,385,1,CAPACITY,1)
    second,conflict2=confidence_set(n,127,2,CAPACITY,1)
    gap=max(0,second.lo-first.hi,first.lo-second.hi)
    return {"window":n,"first_positive_count":385,"second_positive_count":127,
        "first_average_interval":interval_fields(first),"second_average_interval":interval_fields(second),
        "minimum_parameter_gap_lower_exact":str(F(gap,SCALE)),
        "permitted_gap_exact":str(F(513,4096)),
        "reject_declared_rate_or_confidence_event":conflict1 or conflict2 or audit_intervals(first,second,513,F(1,4096)),
        "passing_this_test_does_not_prove_future_stability":True}


class RefreshingDriftControllerTests(unittest.TestCase):
    def test_nonidentical_count_recursion_contains_every_exact_small_probability(self):
        ps=[F(1,3),F(3,5),F(2,7),F(4,9)]
        actual=count_distribution([I.exact(p) for p in ps])
        exact=[F(0)]*5
        for bits in product((0,1),repeat=4):
            probability=F(1)
            for bit,p in zip(bits,ps): probability*=p if bit else 1-p
            exact[sum(bits)]+=probability
        for interval,value in zip(actual,exact):
            self.assertLessEqual(F(interval.lo,SCALE),value)
            self.assertGreaterEqual(F(interval.hi,SCALE),value)

    def test_query_preparation_consumes_a_drift_step_and_average_budget_is_exact(self):
        rate=F(1,4096)
        for look in range(1,9):
            start=(look-1)*513+1
            past=[source_path(start+i) for i in range(512)]
            query=source_path(start+512)
            self.assertEqual(sum(past)/512-query,rate*F(513,2))
        self.assertEqual(source_path(4104),F(64,125)-F(4103,4096))

    def test_audit_accepts_compatible_intervals_and_flags_an_observed_jump(self):
        self.assertTrue(audit_example()["reject_declared_rate_or_confidence_event"])
        a=I.rational(1,5);b=I.rational(1,4)
        self.assertFalse(audit_intervals(a,b,513,F(1,4096)))

    def test_audit_lock_transition_matches_an_explicit_two_block_joint_tree(self):
        first=[I.rational(1,3),I.rational(2,3)]
        second=[I.rational(1,4),I.rational(3,4)]
        sets=[(I.rational(-1,4),False),(I.rational(1,4),False)]
        output=update_trusted_counts(first,sets,second,sets,I.rational(1,10))
        for interval,answer in zip(output,(F(1,12),F(1,2))):
            self.assertLessEqual(F(interval.lo,SCALE),answer)
            self.assertGreaterEqual(F(interval.hi,SCALE),answer)
        none=update_trusted_counts([I.exact(0)]*2,sets,second,sets,I.rational(1,10))
        self.assertTrue(all(x.hi==0 for x in none))

    def test_refresh_changes_sign_and_falls_back_near_the_crossing(self):
        report=refresh_certificate()
        rows=report["rows"]
        self.assertGreater(F(rows[0]["positive_action_probability"]["lower_exact"]),F(95,100))
        self.assertGreater(F(rows[-1]["negative_action_probability"]["lower_exact"]),F(6,10))
        self.assertGreater(F(rows[3]["fallback_probability_including_locked_audit"]["lower_exact"]),F(99,100))
        base=F(report["fallback_success"]["upper_exact"])
        self.assertGreater(F(report["mean_refreshed_H_success"]["lower_exact"]),base)
        self.assertLess(F(rows[-1]["stale_first_policy_H_success"]["upper_exact"]),base)
        self.assertLess(F(rows[-1]["cumulative_audit_lock_probability"]["upper_exact"]),F(1,100))

    def test_actual_future_query_tree_matches_the_declining_source_parameter(self):
        for index in (513,4104):
            kappa=source_path(index)
            action=1 if kappa>0 else -1
            success=0.
            for label in (-1,1):
                records=learned_query_records(calibration_state(H,label),float(kappa),action)
                success+=sum(np.trace(state).real for guess,state in records.values() if guess==label)/2
            expected,_=matched_policy_value(kappa,action)
            self.assertAlmostEqual(success,sum(expected.floats())/2,places=13)

    def test_all_reported_action_masses_and_cumulative_resources_are_accounted_for(self):
        report=refresh_certificate()
        for row in report["rows"]:
            fields=[row[k] for k in ("positive_action_probability","negative_action_probability",
                                    "fallback_probability_including_locked_audit")]
            self.assertLessEqual(sum(F(x["lower_exact"]) for x in fields),1)
            self.assertGreaterEqual(sum(F(x["upper_exact"]) for x in fields),1)
        self.assertEqual(report["resources"]["fresh_pure_slots_allocated_and_retained"],28736)
        self.assertEqual(report["resources"]["pair_rotations_upper"],139568)
        self.assertEqual(report["resources"]["actual_readouts"],4112)

    def test_fixed_zero_drift_instance_does_not_require_oracle_updates(self):
        report=refresh_certificate(64,2,F(0),F(1,10))
        self.assertFalse(report["controller_uses_true_path_values"])
        self.assertGreater(F(report["total_expected_extra_correct_H_answers"]["lower_exact"]),0)
        with self.assertRaises(ValueError): refresh_certificate(0,1)
        with self.assertRaises(ValueError): audit_intervals(I.exact(0),I.exact(0),1,-1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RefreshingDriftControllerTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":117,"finite_refresh_and_audit":refresh_certificate(),"observed_violation_example":audit_example(),
        "quantum_theory_derived_from_cognition":False,"automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("refreshing_drift_controller_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
