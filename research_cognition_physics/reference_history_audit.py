"""Round 123: temporal correlations, honest posteriors and paid reference repair.

The diagnostic uses fresh real Bell targets and actual old readout visibility.
It adds outcome communication; it is not free purification of a shared resource.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
import math
from pathlib import Path
import unittest
import numpy as np

from bell_network_statistics import BELL_SIGNS
from certified_intervals import Interval as I, SCALE
from independent_source_alignment import pair_reference, real_x_with_reset
from joint_relation_records import interval_fields
from one_bit_real_network import visibility_interval
from reusable_network_reference import sector_tables, retained_branch


def binomial_mass(probability,count):
    if not isinstance(count,int) or count<0: raise ValueError("Nonnegative integer count required.")
    p=I.exact(probability)
    if p.lo<0 or p.hi>SCALE: raise ValueError("Probability outside [0,1].")
    result=[I.exact(1)]
    for _ in range(count):
        updated=[I.exact(0) for _ in range(len(result)+1)]
        for k,mass in enumerate(result):
            updated[k]+=mass*(1-p)
            updated[k+1]+=mass*p
        result=updated
    return result


def diagnostic_parameters():
    g=visibility_interval(5)
    visibility=visibility_interval(5)**3/I.exact(2).sqrt()
    return g,visibility


def signed_outcome(a,b,c):
    return BELL_SIGNS[b][1]*a*c


def exact_posterior(g,visibility,count,positive):
    """Rational toy instances, independently checked against full history trees."""
    g,v=F(g),F(visibility)
    if not -1<=g<=1 or not 0<v<1 or not 0<=positive<=count: raise ValueError("Invalid parameters.")
    p=(1+v)/2
    plus=(1+g)/2*p**positive*(1-p)**(count-positive)
    minus=(1-g)/2*(1-p)**positive*p**(count-positive)
    if plus+minus==0: raise ValueError("Impossible history.")
    return (plus-minus)/(plus+minus)


def correction_certificate(count,g=None,visibility=None):
    initial_g,actual_v=diagnostic_parameters()
    g=initial_g if g is None else I.exact(g)
    v=actual_v if visibility is None else I.exact(visibility)
    if g.lo<0 or g.hi>SCALE or v.lo<=0 or v.hi>=SCALE: raise ValueError("Use 0<=g<=1, 0<v<1.")
    plus_masses=binomial_mass((1+v)/2,count)
    minus_masses=binomial_mass((1-v)/2,count)
    prior_plus,prior_minus=(1+g)/2,(1-g)/2
    risk=I.exact(0)
    optimal_lower=0
    actions=[]
    flip_probability=I.exact(0)
    for k in range(count+1):
        plus=prior_plus*plus_masses[k]
        minus=prior_minus*minus_masses[k]
        # Flip only when outward intervals prove the negative sector more likely.
        flip=minus.lo>plus.hi
        actions.append(-1 if flip else 1)
        risk+=plus if flip else minus
        optimal_lower+=min(plus.lo,minus.lo)
        if flip: flip_probability+=plus+minus
    repaired_g=1-2*risk
    return {"diagnostic_trials":count,"actions_by_positive_count":actions,
            "certified_policy_error":interval_fields(risk),
            "optimal_bayes_error_lower_exact":str(F(optimal_lower,SCALE)),
            "post_correction_mean_reference_correlation":interval_fields(repaired_g),
            "probability_of_executing_real_X":interval_fields(flip_probability),
            "resources_excluding_initial_alignment":{"fresh_target_Bell_pairs":2*count,
                "original_readouts":20*count,
                "pair_rotations_upper":23*count+(1 if any(a<0 for a in actions) else 0),
                "classical_outcome_bits_to_Charlie":2*count,
                "one_additional_pure_reset_for_possible_X":any(a<0 for a in actions)},
            "uses_actual_unknown_orientation_in_policy":False}


def variance_certificate(count):
    if not isinstance(count,int) or count<1: raise ValueError("Positive count required.")
    g,v=diagnostic_parameters()
    floor=(1-g*g)*v*v
    reused=(1-v*v)/count+floor
    fresh=(1-g*g*v*v)/count
    return {"trials":count,"reused_reference_sample_mean_variance":interval_fields(reused),
            "fresh_reference_sample_mean_variance":interval_fields(fresh),
            "persistent_covariance":interval_fields(floor),
            "persistent_standard_deviation":interval_fields(floor.sqrt())}


class ReferenceHistoryAuditTests(unittest.TestCase):
    def test_diagnostic_bit_is_sufficient_in_the_full_actual_record_table(self):
        g,v=diagnostic_parameters()
        vis=sum(v.floats())/2
        for sector in (-1,1):
            table=sector_tables()[sector]
            for a,b,c in product((-1,1),range(4),(-1,1)):
                expected=(1+sector*vis*signed_outcome(a,b,c))/16
                self.assertAlmostEqual(table[1,4,a,b,c],expected,places=14)

    def test_exact_filter_matches_rational_history_enumeration(self):
        g,v=F(3,5),F(2,3)
        for history in product((-1,1),repeat=6):
            plus,minus=(1+g)/2,(1-g)/2
            for value in history:
                plus*=(1+value*v)/2
                minus*=(1-value*v)/2
            expected=(plus-minus)/(plus+minus)
            self.assertEqual(exact_posterior(g,v,6,history.count(1)),expected)

    def test_binomial_interval_filter_contains_an_exact_bayes_risk(self):
        g,v=F(3,5),F(2,3)
        for count in (0,2,7):
            p=(1+v)/2
            expected=sum(min((1+g)/2*math.comb(count,k)*p**k*(1-p)**(count-k),
                             (1-g)/2*math.comb(count,k)*(1-p)**k*p**(count-k)) for k in range(count+1))
            report=correction_certificate(count,g,v)
            self.assertLessEqual(F(report["optimal_bayes_error_lower_exact"]),expected)
            self.assertGreaterEqual(F(report["certified_policy_error"]["upper_exact"]),expected)

    def test_actual_two_trial_reference_update_uses_public_records_only(self):
        first=(1,4,1,0,1)
        second=(1,4,-1,1,1)
        g=.6
        state=retained_branch(retained_branch(pair_reference(g),*first),*second)
        v=sum(diagnostic_parameters()[1].floats())/2
        weights=[(1+s*g)/2*np.prod([(1+s*v*signed_outcome(k[2],k[3],k[4]))/16 for k in (first,second)])
                 for s in (1,-1)]
        expected=(weights[0]-weights[1])/sum(weights)
        np.testing.assert_allclose(state/np.trace(state),pair_reference(expected),atol=3e-15,rtol=0)
        self.assertLess(expected,0)
        repaired=real_x_with_reset(state/np.trace(state),1,2)
        np.testing.assert_allclose(repaired,pair_reference(-expected),atol=3e-15,rtol=0)

    def test_persistent_covariance_and_sample_variance_follow_complete_binary_trees(self):
        g,v=F(3,5),F(2,3)
        for count in (2,4):
            mean=second=F(0)
            for history in product((-1,1),repeat=count):
                mass=sum((1+s*g)/2*math.prod((1+s*v*w)/2 for w in history) for s in (-1,1))
                value=F(sum(history),count)
                mean+=mass*value
                second+=mass*value*value
            self.assertEqual(mean,g*v)
            self.assertEqual(second-mean**2,(1-v*v)/count+(1-g*g)*v*v)

    def test_rare_adverse_history_requires_an_honest_change_of_prediction(self):
        seven=correction_certificate(7)
        eight=correction_certificate(8)
        self.assertEqual(seven["actions_by_positive_count"],[1]*8)
        self.assertEqual(eight["actions_by_positive_count"][0],-1)
        self.assertGreater(F(eight["probability_of_executing_real_X"]["lower_exact"]),0)

    def test_paid_diagnostics_improve_unconditional_reference_with_all_branches_kept(self):
        before=correction_certificate(0)
        after=correction_certificate(32)
        self.assertLess(F(after["certified_policy_error"]["upper_exact"]),
                        F(before["certified_policy_error"]["lower_exact"])/1000)
        costs=after["resources_excluding_initial_alignment"]
        self.assertEqual(costs["original_readouts"],640)
        self.assertEqual(costs["classical_outcome_bits_to_Charlie"],64)
        self.assertEqual(costs["fresh_target_Bell_pairs"],64)

    def test_reuse_variance_does_not_obey_an_iid_one_over_N_claim(self):
        report=variance_certificate(10**6)
        reused=report["reused_reference_sample_mean_variance"]
        fresh=report["fresh_reference_sample_mean_variance"]
        self.assertGreater(F(reused["lower_exact"]),5*F(fresh["upper_exact"]))
        with self.assertRaises(ValueError): variance_certificate(0)
        with self.assertRaises(ValueError): binomial_mass(F(2),1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ReferenceHistoryAuditTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    g,v=diagnostic_parameters()
    report={"round":123,"initial_reference_correlation":interval_fields(g),
            "actual_diagnostic_visibility":interval_fields(v),
            "variance_rows":[variance_certificate(n) for n in (1,1000,10**6)],
            "repair_rows":[correction_certificate(n) for n in (0,7,8,16,32,64)],
            "all_results_retained_no_postselection":True,
            "new_targets_and_communication_are_paid_resources":True,
            "quantum_theory_derived_from_cognition":False,
            "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("reference_history_audit_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
