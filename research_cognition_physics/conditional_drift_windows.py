"""Round 116: conditional-mean learning, drift budgets, and finite memory windows."""
import argparse
import json
import math
import unittest
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I, SCALE
from drifting_source_transfer import (CAPACITY, expanded_future_interval,
    robust_action, robust_gain_floor)
from joint_relation_records import interval_fields
from nonzero_fragment_records import biased_environment
from record_correlation_frontier import biased_probe_records
from unknown_correlation_learning import (confidence_set, scheduled_radius,
    calibration_visibility, uniform_risk_certificate)


def average_to_future_budget(count, rate, gap=1):
    rate=F(rate)
    if not isinstance(count,int) or count<1 or rate<0 or not isinstance(gap,int) or gap<1:
        raise ValueError("Positive window/gap and nonnegative per-source drift required.")
    return rate*(F(count-1,2)+gap)


def predictive_interval(count,positives,look,rate,gap=1):
    """The first interval covers the window's average predictable kappa, not a fixed mean."""
    average,conflict=confidence_set(count,positives,look,CAPACITY,1)
    drift=average_to_future_budget(count,rate,gap)
    return expanded_future_interval(average,drift),conflict


def window_uncertainty(count,look,rate,gap=1):
    return scheduled_radius(count,look)/calibration_visibility(1)+I.exact(
        average_to_future_budget(count,rate,gap))


@lru_cache(maxsize=None)
def optimum_window(look,rate,maximum=4096,gap=1):
    rate=F(rate)
    if not isinstance(maximum,int) or maximum<1: raise ValueError("Positive maximum required.")
    values=[window_uncertainty(n,look,rate,gap) for n in range(1,maximum+1)]
    best=min(range(maximum),key=lambda i:values[i].hi)
    possible=[i+1 for i,value in enumerate(values) if value.lo<=values[best].hi]
    return {"look":look,"rate_exact":str(rate),"gap_in_source_steps":gap,"maximum_window":maximum,
        "certified_possible_minimizers":possible,"selected_window":best+1,
        "unclipped_worst_case_radius":interval_fields(values[best]),
        "radius_if_all_maximum_samples_used":interval_fields(values[-1]),
        "optimization_scope":"This concentration-plus-drift bound, fixed look and source rate, integer window sizes only"}


def horizon_certificate(count,positives,look,rate,target=F(1,10)):
    rate,target=F(rate),F(target)
    if rate<=0 or target<0: raise ValueError("Positive rate and nonnegative target required.")
    average,conflict=confidence_set(count,positives,look,CAPACITY,1)
    if conflict: return {"action":0,"certified_future_source_steps":0}
    if average.lo>0:
        sign=1;margin=F(average.lo,SCALE)-target
    elif average.hi<0:
        sign=-1;margin=-F(average.hi,SCALE)-target
    else:
        return {"action":0,"certified_future_source_steps":0}
    horizon=margin/rate-F(count-1,2)
    steps=max(0,horizon.numerator//horizon.denominator)
    return {"action":sign if steps else 0,"certified_future_source_steps":steps,
        "rate_exact":str(rate),"target_exact":str(target),
        "steps_are_preparation_indices_not_physical_time":True}


def shared_label_predictive_parameter(history,kappa=CAPACITY):
    """Posterior conditional source parameter for the aligned shared-label family."""
    mean=calibration_visibility(1)*I.exact(F(kappa))
    weights=[]
    for sign in (1,-1):
        p=I.exact(1)
        for y in history:
            if y not in (-1,1): raise ValueError("Binary history required.")
            p*= (1+y*sign*mean)/2
        weights.append(p)
    return I.exact(F(kappa))*(weights[0]-weights[1])/(weights[0]+weights[1])


def window_record_certificate(count,positives,look,rate,target=F(1,10),gap=1):
    interval,conflict=predictive_interval(count,positives,look,rate,gap)
    action=0 if conflict else robust_action(interval,target)
    return {"window_samples":count,"positive_count":positives,"look":look,
        "rate_exact":str(rate),"average_to_future_budget_exact":str(average_to_future_budget(count,rate,gap)),
        "future_parameter_interval":interval_fields(interval),"chosen_action":action,
        "gain_floor_on_coverage_and_drift_assumption":interval_fields(robust_gain_floor(interval,action)),
        "uniform_probability_failure_bound":interval_fields(uniform_risk_certificate())}


class ConditionalDriftWindowTests(unittest.TestCase):
    def test_linear_source_path_attains_the_average_to_next_source_drift_budget(self):
        rate=F(1,4096)
        for count in (1,64,512,1024):
            path=[CAPACITY-rate*i for i in range(count)]
            future=CAPACITY-rate*count
            self.assertEqual(sum(path)/count-future,average_to_future_budget(count,rate))
        self.assertEqual(average_to_future_budget(512,rate,7)-average_to_future_budget(512,rate),6*rate)

    def test_history_dependent_source_has_conditional_not_unconditional_calibration_mean(self):
        for last in (-1,1):
            kappa=F(last,3)
            records=biased_probe_records(.8,(.6,)*3,float(kappa),(.9,0.,0.),1)
            mean=sum((1-2*h[0])*np.trace(rho).real for h,rho in records.items())
            expected=calibration_visibility(1).floats()[0]*float(kappa)
            self.assertAlmostEqual(mean,expected,places=13)

    def test_conditional_centering_preserves_exponential_bound_without_independence(self):
        alpha=calibration_visibility(1).floats()[0]
        n=10
        expectations={lam:0. for lam in (-.4,.4)}
        centered=0.;tail=0.;mass=0.
        for word in product((-1,1),repeat=n):
            probability=1.;difference=0.
            for i,y in enumerate(word):
                kappa=0. if i==0 else -.4*word[i-1]
                mean=alpha*kappa
                probability*=(1+y*mean)/2
                difference+=y-mean
            mass+=probability;centered+=probability*difference
            tail+=probability*int(abs(difference/n)>.6)
            for lam in expectations: expectations[lam]+=probability*math.exp(lam*difference)
        self.assertAlmostEqual(mass,1,places=13)
        self.assertAlmostEqual(centered,0,places=13)
        for lam,value in expectations.items(): self.assertLessEqual(value,math.exp(n*lam*lam/2)+1e-13)
        self.assertLessEqual(tail,2*math.exp(-n*.6**2/2))

    def test_shared_label_average_source_is_zero_but_prediction_moves_with_history(self):
        empty=shared_label_predictive_parameter(())
        self.assertLessEqual(empty.lo,0);self.assertGreaterEqual(empty.hi,0)
        plus=shared_label_predictive_parameter((1,1,1))
        minus=shared_label_predictive_parameter((-1,-1,-1))
        self.assertGreater(plus.lo,0);self.assertLess(minus.hi,0)
        self.assertAlmostEqual(plus.floats()[0],-minus.floats()[1],places=14)
        k=sum(plus.floats())/2
        p=(1+k/float(CAPACITY))/2
        mixture=p*biased_environment((.6,)*3,float(CAPACITY))+(1-p)*biased_environment((.6,)*3,-float(CAPACITY))
        np.testing.assert_allclose(mixture,biased_environment((.6,)*3,k),atol=1e-15)

    def test_integer_window_optimum_is_certified_and_longer_can_be_worse(self):
        report=optimum_window(1,F(1,4096))
        self.assertEqual(len(report["certified_possible_minimizers"]),1)
        n=report["selected_window"]
        self.assertLess(n,1024)
        self.assertLess(F(report["unclipped_worst_case_radius"]["upper_exact"]),
                        F(report["radius_if_all_maximum_samples_used"]["lower_exact"]))
        for neighbor in (n-1,n+1):
            self.assertLess(window_uncertainty(n,1,F(1,4096)).hi,
                            window_uncertainty(neighbor,1,F(1,4096)).lo)

    def test_zero_drift_uses_all_samples_but_nonzero_drift_has_a_finite_optimum(self):
        report=optimum_window(1,F(0),1024)
        self.assertEqual(report["selected_window"],1024)
        fast=optimum_window(1,F(1,1000),1024)
        slow=optimum_window(1,F(1,4096),1024)
        self.assertLess(fast["selected_window"],slow["selected_window"])

    def test_certified_knowledge_horizon_expires_without_new_records(self):
        report=horizon_certificate(512,385,1,F(1,4096))
        horizon=report["certified_future_source_steps"]
        self.assertGreater(horizon,1)
        before,_=predictive_interval(512,385,1,F(1,4096),horizon)
        after,_=predictive_interval(512,385,1,F(1,4096),horizon+1)
        self.assertEqual(robust_action(before,F(1,10)),1)
        self.assertEqual(robust_action(after,F(1,10)),0)

    def test_invalid_drift_and_missing_future_bridge_are_not_silently_accepted(self):
        with self.assertRaises(ValueError): average_to_future_budget(10,-1)
        with self.assertRaises(ValueError): average_to_future_budget(0,F(1,10))
        with self.assertRaises(ValueError): average_to_future_budget(10,F(1,10),0)
        with self.assertRaises(ValueError): horizon_certificate(10,5,1,0)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ConditionalDriftWindowTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    rate=F(1,4096)
    report={"round":116,
        "estimand":"window average of kappa_i(history_before_i), where E[Y_i|history]=alpha*kappa_i",
        "independence_required_for_confidence_theorem":False,
        "future_bridge_is_still_required":True,
        "per_source_rate_budget":"L*(gap+(window-1)/2)",
        "window_bound":"sqrt(2*(6+j)/window)/alpha + L*(gap+(window-1)/2)",
        "continuous_optimum_equation":"L^2*alpha^2*window^3=2*(6+j), for L>0",
        "window_examples":[optimum_window(1,l) for l in (F(0),rate,F(1,1000))],
        "record_example":window_record_certificate(512,385,1,rate),
        "expiry_example":horizon_certificate(512,385,1,rate),
        "quantum_theory_derived_from_cognition":False,"automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("conditional_drift_windows_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
