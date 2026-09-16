"""Round 112: finite confidence sets for a hidden source parameter."""
import argparse
import json
import math
import unittest
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I, SCALE
from joint_relation_records import interval_fields
from one_bit_real_network import visibility_interval
from record_correlation_frontier import biased_probe_records


def exp_minus_interval(x):
    """Certified exp(-x), x>=0 rational: positive Taylor series plus squaring."""
    x = F(x)
    if x < 0:
        raise ValueError("Expected a nonnegative rational.")
    depth = 0
    while x > F(1, 2):
        x /= 2
        depth += 1
    term = I.exact(1)
    total = term
    for j in range(1, 41):
        term = term * I.exact(x) / j
        total += term
    # Remainder of exp(x) <= next term / (1 - x/42).
    tail = term * I.exact(x) / 41 / (1 - I.exact(x) / 42)
    positive = I(total.lo, (total + tail).hi)
    result = 1 / positive
    for _ in range(depth):
        result *= result
    return result


@lru_cache(maxsize=None)
def calibration_visibility(count=1):
    return visibility_interval(count)


def source_mean_interval(kappa, response_strength, count=1):
    """Y is minus the observed majority sign; E[Y]=gamma*a*kappa."""
    kappa, a = F(kappa), F(response_strength)
    if not -1 <= kappa <= 1 or not 0 < a <= 1:
        raise ValueError("Nonzero response strength and physical kappa required.")
    return calibration_visibility(count) * I.exact(a * kappa)


def scheduled_radius(n, look):
    if not isinstance(n, int) or n < 1 or not isinstance(look, int) or look < 1:
        raise ValueError("Positive integer sample size and look index required.")
    # Assigned tail bound: 2 exp(-(6+j)); sum over all j>=1 is < 1/100.
    return I.rational(2 * (6 + look), n).sqrt()


def uniform_risk_certificate():
    e1 = exp_minus_interval(1)
    return 2 * exp_minus_interval(7) / (1 - e1)


def confidence_set(n, positives, look, capacity, response_strength, count=1):
    """Outward confidence hull, intersected with a known symmetric source domain."""
    capacity, a = F(capacity), F(response_strength)
    if not isinstance(positives, int) or not 0 <= positives <= n:
        raise ValueError("Invalid count.")
    if not 0 < capacity <= 1 or not 0 < a <= 1:
        raise ValueError("Positive physical capacity and response strength required.")
    radius = scheduled_radius(n, look)
    mean = I.rational(2 * positives - n, n)
    g = calibration_visibility(count) * I.exact(a)
    raw = I(mean.lo - radius.hi, mean.hi + radius.hi) / g
    domain = I.exact(capacity)
    lo, hi = max(raw.lo, -domain.hi), min(raw.hi, domain.hi)
    if lo > hi:
        # This record contradicts the confidence event; it is not a license to act.
        return I(-domain.hi, domain.hi), True
    return I(lo, hi), False


def choose_sign(n, positives, look, capacity, response_strength, target=F(1,5), count=1):
    target = F(target)
    if not 0 < target <= F(capacity):
        raise ValueError("The activation margin must fit the source family.")
    interval, conflict = confidence_set(n, positives, look, capacity, response_strength, count)
    target_i = I.exact(target)
    if conflict:
        return 0
    if interval.lo >= target_i.hi:
        return 1
    if interval.hi <= -target_i.hi:
        return -1
    return 0


def matched_policy_value(kappa, chosen_sign, c=F(4,5), d=F(4,5), m=1, q=1):
    """Exact finite expected H success for a chosen policy, including wrong signs."""
    c, d, kappa = map(F, (c, d, kappa))
    if not F(1,2) <= c < 1 or not 0 <= d < 1 or not -1 <= kappa <= 1:
        raise ValueError("Use the matched control range.")
    if chosen_sign not in (-1,0,1):
        raise ValueError("Use -1, 0 (fallback), or 1.")
    w = (1 - I.exact(d*d)).sqrt()
    baseline = (1 + calibration_visibility(q)*w*I.exact(c))/2
    coefficient = calibration_visibility(q)*w*calibration_visibility(m)*I.exact(1-c)/4
    return baseline + coefficient*I.exact(chosen_sign*kappa), coefficient


def calibration_record_certificate(n, positives, look, a, target=F(1,5)):
    cap = F(64,125)
    confidence, conflict = confidence_set(n, positives, look, cap, a)
    sign = choose_sign(n, positives, look, cap, a, target)
    baseline, coefficient = matched_policy_value(0,0)
    return {
        "fresh_independent_source_samples": n, "positive_Y_count": positives,
        "scheduled_look": look, "calibration_response_strength_exact": str(a),
        "kappa_confidence_hull": interval_fields(confidence), "confidence_domain_conflict": conflict,
        "chosen_sign_zero_means_fallback": sign,
        "activation_margin_exact": str(target),
        "simultaneous_miscoverage_upper": interval_fields(uniform_risk_certificate()),
        "matched_query_baseline": interval_fields(baseline),
        "conditional_success_gain_floor_on_coverage": interval_fields(coefficient*F(target) if sign else I.exact(0)),
        "scope": "iid unknown-kappa source; specified look schedule, matched recovery, independent future task"
    }


class UnknownCorrelationLearningTests(unittest.TestCase):
    def test_actual_biased_probe_has_the_signed_bernoulli_mean(self):
        for a, gates in ((.19,(.9,)*3), (1.,(.9,0.,0.))):
            for kappa in (-.512,0.,.3):
                for count in (1,3):
                    records = biased_probe_records(.8,(.6,)*3,kappa,gates,count)
                    mean = sum((-1 if sum(h)>count//2 else 1)*np.trace(rho).real
                               for h,rho in records.items())
                    expected = sum(calibration_visibility(count).floats())/2*a*kappa
                    self.assertAlmostEqual(mean,expected,places=13)

    def test_exponential_certificate_and_infinite_risk_sum(self):
        for value in (F(0),F(1,3),F(7),F(13)):
            lo, hi = exp_minus_interval(value).floats()
            self.assertAlmostEqual((lo+hi)/2,math.exp(-float(value)),places=14)
        self.assertLess(F(uniform_risk_certificate().hi,SCALE),F(1,100))
        self.assertLess(uniform_risk_certificate().hi-uniform_risk_certificate().lo,100000)

    def test_confidence_hull_contains_every_model_inside_the_mean_error_event(self):
        for a in (F(19,100),F(1)):
            for n,j in ((64,1),(256,3),(1024,5)):
                eps = scheduled_radius(n,j).floats()[0]
                for k in range(0,n+1,17):
                    interval, conflict = confidence_set(n,k,j,F(64,125),a)
                    for truth in (F(-64,125),F(-1,5),F(0),F(1,5),F(64,125)):
                        mu = source_mean_interval(truth,a).floats()[0]
                        if abs((2*k-n)/n-mu)<eps-1e-14:
                            self.assertFalse(conflict)
                            self.assertLessEqual(F(interval.lo,SCALE),truth)
                            self.assertGreaterEqual(F(interval.hi,SCALE),truth)

    def test_action_depends_on_evidence_and_not_an_oracle_kappa(self):
        self.assertEqual(choose_sign(256,194,3,F(64,125),F(1)),1)
        self.assertEqual(choose_sign(256,62,3,F(64,125),F(1)),-1)
        self.assertEqual(choose_sign(256,128,3,F(64,125),F(1)),0)
        self.assertEqual(choose_sign(256,194,3,F(64,125),F(19,100)),0)  # domain conflict

    def test_wrong_sign_is_costed_and_fallback_has_no_gain_or_loss(self):
        base,_ = matched_policy_value(0,0)
        right,coef = matched_policy_value(F(64,125),1)
        wrong,_ = matched_policy_value(F(64,125),-1)
        abstain,_ = matched_policy_value(F(64,125),0)
        self.assertGreater(right.lo,base.hi)
        self.assertLess(wrong.hi,base.lo)
        self.assertEqual((abstain.lo,abstain.hi),(base.lo,base.hi))
        self.assertAlmostEqual(sum((right+wrong).floats())/2,sum((2*base).floats())/2,places=14)

    def test_rounded_confidence_decisions_never_overstate_sign_margin(self):
        cap=F(64,125)
        for n,j in ((64,1),(128,2),(256,3)):
            for k in range(n+1):
                sign=choose_sign(n,k,j,cap,1)
                interval,_=confidence_set(n,k,j,cap,1)
                if sign==1: self.assertGreaterEqual(F(interval.lo,SCALE),F(1,5))
                if sign==-1: self.assertLessEqual(F(interval.hi,SCALE),-F(1,5))

    def test_fixed_probe_blocks_are_independent_only_with_new_preparations(self):
        # Two fresh sources give the product Bernoulli law, with a nonzero mixed event.
        p=(1+source_mean_interval(F(64,125),1))/2
        probabilities=[p*p,p*(1-p),(1-p)*p,(1-p)*(1-p)]
        total=sum(probabilities,I.exact(0))
        self.assertLessEqual(total.lo,SCALE)
        self.assertGreaterEqual(total.hi,SCALE)
        self.assertGreater(probabilities[1].lo,0)

    def test_degenerate_inputs_do_not_create_a_fake_confidence_certificate(self):
        for args in ((0,1),(64,0)):
            with self.assertRaises(ValueError): scheduled_radius(*args)
        with self.assertRaises(ValueError): confidence_set(64,65,1,F(1,2),1)
        with self.assertRaises(ValueError): confidence_set(64,32,1,F(1,2),0)
        with self.assertRaises(ValueError): choose_sign(64,32,1,F(1,10),1,F(1,5))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(UnknownCorrelationLearningTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":112,
        "observation_law":"P(Y=+1)=(1+gamma_l*a*kappa)/2; one fresh source per Y",
        "confidence_radius_for_mean":"sqrt(2*(6+j)/n_j), at prespecified increasing looks",
        "risk_bound_all_looks":interval_fields(uniform_risk_certificate()),
        "examples":[calibration_record_certificate(256,k,3,F(1)) for k in (194,128,62)],
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("unknown_correlation_learning_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
