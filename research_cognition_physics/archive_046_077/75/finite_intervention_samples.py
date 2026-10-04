"""Round 75: certified finite sample risk for the fixed intervention experiment.

Perfect blocking gives likelihood ratio r^D, D=sum S. All trials, including
S=0, count. Equal-prior Bayes risk equals both conditional errors by symmetry.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

from certified_intervals import Interval as I, SCALE
from communication_causal_audit import calibration_interval
from intervention_probe_reduction import signal_interval, signal_diagnostic


def walk_distribution(count,mu):
    if type(count) is not int or count<0: raise ValueError("Use a nonnegative integer trial count.")
    states={0:mu*0+1}
    weights=((-1,(1-mu)/4),(0,Fraction(1,2)),(1,(1+mu)/4))
    for _ in range(count):
        next_states={}
        for balance,mass in states.items():
            for step,weight in weights:
                next_states[balance+step]=next_states.get(balance+step,mu*0)+mass*weight
        states=next_states
    return states


def fixed_error(count,mu):
    states=walk_distribution(count,mu)
    return sum((mass for balance,mass in states.items() if balance<0),mu*0)+states[0]/2


def binomial_distribution(count,p):
    states=[p*0+1]
    for _ in range(count):
        new=[p*0 for _ in range(len(states)+1)]
        for k,mass in enumerate(states):
            new[k]+=mass*(1-p)
            new[k+1]+=mass*p
        states=new
    return states


def ignored_flag_error(count,mu):
    first=binomial_distribution(count,(1+mu)/2)
    second=binomial_distribution(count,mu*0+Fraction(1,2))
    if isinstance(mu,I):
        return sum((I(min(a.lo,b.lo),min(a.hi,b.hi)) for a,b in zip(first,second)),I.exact(0))/2
    return sum(min(a,b) for a,b in zip(first,second))/2


def first_sample_certificate(error_function=fixed_error,target=Fraction(1,100)):
    mu=signal_interval()
    for count in range(1,201):
        if error_function(count,signal_diagnostic())<float(target): break
    else: raise RuntimeError("Diagnostic search exhausted.")
    previous,current=error_function(count-1,mu),error_function(count,mu)
    if previous.lo<I.exact(target).hi or current.hi>=I.exact(target).lo:
        raise ArithmeticError("The strict certificate did not resolve the threshold.")
    return {"minimum_fixed_trials":count,"target_error_exact":str(target),
            "previous_error_interval":previous.floats(),"error_interval":current.floats(),
            "previous_error_lower_exact":str(Fraction(previous.lo,SCALE)),
            "error_upper_exact":str(Fraction(current.hi,SCALE))}


def sample_report():
    full=first_sample_certificate()
    ignored=first_sample_certificate(ignored_flag_error)
    _,_,flip=calibration_interval()
    coin=(flip.lo+flip.hi)//2
    per_trial=Fraction(max(coin-flip.lo,flip.hi-coin),SCALE)
    finite_coin_upper=Fraction(full["error_upper_exact"])+full["minimum_fixed_trials"]*per_trial
    return {"retained_flag":full,"ignored_flag":ignored,
            "decision":"Choose complex for D>0, real feedback for D<0, and a fair coin for D=0",
            "conditional_errors_are_equal_for_exact_calibration":True,
            "finite_100_bit_calibration_each_error_upper_exact":str(finite_coin_upper),
            "finite_calibration_still_below_one_percent":finite_coin_upper<Fraction(1,100),
            "trial_count_includes_uninformative_flags":True,
            "minimum_scope":"The 18 original settings and observe/block actions, fresh calibrated trials, perfect trusted blocking, equal priors"}


class FiniteInterventionSamplesTests(unittest.TestCase):
    def test_dynamic_distribution_matches_complete_small_record_enumeration(self):
        mu=Fraction(3,5)
        weights={-1:(1-mu)/4,0:Fraction(1,2),1:(1+mu)/4}
        for n in range(6):
            expected={d:Fraction(0) for d in range(-n,n+1)}
            for record in product((-1,0,1),repeat=n):
                expected[sum(record)]+=math.prod(weights[s] for s in record)
            self.assertEqual(walk_distribution(n,mu),expected)
            self.assertEqual(sum(expected.values()),1)

    def test_bayes_risk_from_all_records_equals_the_signed_sum_rule(self):
        mu=Fraction(1,2)
        weights={-1:(1-mu)/4,0:Fraction(1,2),1:(1+mu)/4}
        for n in range(5):
            error=sum(min(math.prod(weights[s] for s in record),math.prod(weights[-s] for s in record))
                      for record in product((-1,0,1),repeat=n))/2
            self.assertEqual(error,fixed_error(n,mu))

    def test_conditioning_on_informative_count_gives_the_same_error_and_counts_every_trial(self):
        mu=Fraction(2,3)
        for n in range(7):
            total=Fraction(0)
            for m in range(n+1):
                dist=binomial_distribution(m,(1+mu)/2)
                conditional=sum(p for k,p in enumerate(dist) if 2*k<m)
                if m%2==0: conditional+=dist[m//2]/2
                total+=Fraction(math.comb(n,m),2**n)*conditional
            self.assertEqual(total,fixed_error(n,mu))

    def test_interval_dynamic_program_encloses_independent_exact_rational_values(self):
        for n in (0,1,6,13):
            value=fixed_error(n,I.rational(3,5))
            exact=fixed_error(n,Fraction(3,5))
            self.assertLessEqual(Fraction(value.lo,SCALE),exact)
            self.assertGreaterEqual(Fraction(value.hi,SCALE),exact)

    def test_minimum_sample_count_has_strict_certificates_on_both_sides(self):
        report=first_sample_certificate()
        self.assertGreater(Fraction(report["previous_error_lower_exact"]),Fraction(1,100))
        self.assertLess(Fraction(report["error_upper_exact"]),Fraction(1,100))
        self.assertGreater(report["minimum_fixed_trials"],10)
        self.assertLess(report["minimum_fixed_trials"],40)

    def test_ignoring_the_flag_strictly_increases_required_samples(self):
        report=sample_report()
        self.assertGreater(report["ignored_flag"]["minimum_fixed_trials"],report["retained_flag"]["minimum_fixed_trials"])
        self.assertEqual(ignored_flag_error(1,Fraction(3,5)),fixed_error(1,Fraction(3,5)))
        for n in (5,12): self.assertGreater(ignored_flag_error(n,Fraction(3,5)),fixed_error(n,Fraction(3,5)))

    def test_finite_fair_bit_calibration_does_not_consume_the_error_budget(self):
        self.assertTrue(sample_report()["finite_calibration_still_below_one_percent"])

    def test_zero_signal_has_no_discrimination_and_ties_are_not_dropped(self):
        for n in (0,1,7): self.assertEqual(fixed_error(n,Fraction(0)),Fraction(1,2))
        self.assertEqual(fixed_error(1,Fraction(1)),Fraction(1,4))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(FiniteInterventionSamplesTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":75,**sample_report(),"quantum_theory_derived_from_cognition":False,
            "empirical_samples_collected":0,"automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("finite_intervention_samples_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__": main()
