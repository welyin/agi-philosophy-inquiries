"""Round 76: robust finite testing with an independently lower-bounded blocker.

Blocking success kappa is independent of current flags/results; failures deliver
the original message. Test the least favourable kappa0 and retain both signed
outcome counts. Budgets below are sufficient, not globally minimax optimal.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I, SCALE
from finite_intervention_samples import fixed_error
from intervention_probe_reduction import signal_interval,signal_diagnostic,reduced_kernel,setting_kernel,informative_settings,OUTCOMES,INCREMENTS,record_increment


def binomial_step(states,p):
    new=[p*0 for _ in range(len(states)+1)]
    for k,mass in enumerate(states):
        new[k]+=mass*(1-p)
        new[k+1]+=mass*p
    return new


def decision_boundary(m,kappa,certified=False):
    kappa=Fraction(kappa)
    if not 0<kappa<=1 or type(m) is not int or m<0: raise ValueError("Use kappa in (0,1] and a nonnegative informative count.")
    if m==0: return 1,0  # fair tie at k=0
    if kappa==1: return m//2+1,(m//2 if m%2==0 else None)
    mu=signal_diagnostic()
    v=(1-2*float(kappa))*mu
    a=(1+mu)/(1+v)
    b=(1-v)/(1-mu)
    cutoff=math.ceil(m*math.log(b)/(math.log(a)+math.log(b)))
    if certified:
        mui=signal_interval()
        vi=(1-2*kappa)*mui
        ai,bi=(1+mui)/(1+vi),(1-vi)/(1-mui)
        below=ai**(cutoff-1)-bi**(m-cutoff+1)
        above=ai**cutoff-bi**(m-cutoff)
        if below.hi>=0 or above.lo<=0: raise ArithmeticError("Unresolved likelihood ordering.")
    return cutoff,None


def conditional_error_rows(count,policy_kappa,actual_kappa=None,certified=False):
    if actual_kappa is None: actual_kappa=policy_kappa
    policy_kappa,actual_kappa=Fraction(policy_kappa),Fraction(actual_kappa)
    if not 0<policy_kappa<=actual_kappa<=1: raise ValueError("Actual success must meet the promised lower bound.")
    mu=signal_interval() if certified else signal_diagnostic()
    pc=(1+mu)/2
    pr=(1+(1-2*actual_kappa)*mu)/2
    cstates,rstates=[mu*0+1],[mu*0+1]
    rows=[]
    for m in range(count+1):
        cut,tie=decision_boundary(m,policy_kappa,certified)
        cerror=sum(cstates[:cut],mu*0)
        rerror=sum(rstates[cut:],mu*0)
        if tie is not None:
            cerror-=cstates[tie]/2
            rerror+=rstates[tie]/2
        rows.append((cerror,rerror))
        if m<count:
            cstates=binomial_step(cstates,pc)
            rstates=binomial_step(rstates,pr)
    return rows


def fixed_budget_errors(count,policy_kappa,actual_kappa=None,certified=False):
    rows=conditional_error_rows(count,policy_kappa,actual_kappa,certified)
    zero=I.exact(0) if certified else 0.
    totals=[zero,zero]
    for m,row in enumerate(rows):
        weight=Fraction(math.comb(count,m),2**count)
        if certified: weight=I.exact(weight)
        else: weight=float(weight)
        for h in range(2): totals[h]+=weight*row[h]
    return tuple(totals)


@lru_cache(maxsize=8)
def sufficient_budget(kappa):
    kappa=Fraction(kappa)
    rows=conditional_error_rows(700,kappa)
    flags=[1.]
    chosen=None
    for count in range(1,701):
        flags=binomial_step(flags,.5)
        errors=[sum(w*rows[m][h] for m,w in enumerate(flags)) for h in range(2)]
        if max(errors)<.01:
            chosen=count
            break
    if chosen is None: raise RuntimeError("Budget search exhausted.")
    errors=fixed_budget_errors(chosen,kappa,certified=True)
    if any(error.hi>=I.rational(1,100).lo for error in errors): raise ArithmeticError("Risk certificate unresolved.")
    return {"blocking_success_lower_bound_exact":str(kappa),"sufficient_fixed_trials":chosen,
            "complex_error_interval":errors[0].floats(),"real_error_at_least_favourable_blocker_interval":errors[1].floats(),
            "complex_error_upper_exact":str(Fraction(errors[0].hi,SCALE)),
            "real_error_upper_exact":str(Fraction(errors[1].hi,SCALE)),
            "both_errors_strictly_below_one_percent":True,"globally_minimum_trial_count_claimed":False}


class ImperfectInterventionTestingTests(unittest.TestCase):
    def test_mixture_of_actual_delivery_and_cut_circuits_gives_the_noisy_blocker_kernel(self):
        x,z=informative_settings()[0]
        delivered,cut=setting_kernel(x,z,False),setting_kernel(x,z,True)
        for kappa in (.25,.5,.8):
            full=kappa*cut+(1-kappa)*delivered
            reduced=np.zeros((2,3))
            for col,(h,a,b,c) in enumerate(OUTCOMES):
                reduced[:,INCREMENTS.index(record_increment(z,h,a,b,c))]+=full[:,col]
            np.testing.assert_allclose(reduced,reduced_kernel(blocking_success=kappa),atol=3e-15)

    def test_likelihood_threshold_order_is_certified_for_every_relevant_informative_count(self):
        for kappa in (Fraction(1,2),Fraction(1,4)):
            for m in range(1,120):
                cut,tie=decision_boundary(m,kappa,True)
                self.assertTrue(1<=cut<=m)
                self.assertIsNone(tie)

    def test_perfect_blocker_reduces_to_the_previous_symmetric_test(self):
        for n in (0,1,6,18):
            first,second=fixed_budget_errors(n,Fraction(1))
            self.assertAlmostEqual(first,fixed_error(n,signal_diagnostic()),places=15)
            self.assertAlmostEqual(first,second,places=15)

    def test_same_signed_sum_can_have_different_likelihoods_when_blocking_is_imperfect(self):
        kernel=reduced_kernel(blocking_success=.5)
        # Empty informative record and (+1,-1) both have balance zero.
        ratio_pair=(kernel[0,2]/kernel[1,2])*(kernel[0,0]/kernel[1,0])
        self.assertNotAlmostEqual(ratio_pair,1.)
        self.assertLess(ratio_pair,1.)

    def test_least_favourable_calibration_controls_every_larger_constant_success_probability(self):
        baseline=fixed_budget_errors(35,Fraction(1,2))
        for actual in (Fraction(3,5),Fraction(3,4),Fraction(1)):
            errors=fixed_budget_errors(35,Fraction(1,2),actual)
            self.assertAlmostEqual(errors[0],baseline[0],places=15)
            self.assertLess(errors[1],baseline[1])

    def test_each_reported_budget_certifies_both_conditional_errors(self):
        counts=[]
        for kappa in (Fraction(1),Fraction(1,2),Fraction(1,4)):
            row=sufficient_budget(kappa)
            counts.append(row["sufficient_fixed_trials"])
            self.assertLess(Fraction(row["complex_error_upper_exact"]),Fraction(1,100))
            self.assertLess(Fraction(row["real_error_upper_exact"]),Fraction(1,100))
        self.assertEqual(counts,sorted(counts))
        self.assertEqual(counts[0],18)

    def test_zero_or_unpromised_blocker_has_no_uniform_finite_separation(self):
        mu=signal_diagnostic()
        np.testing.assert_array_equal(reduced_kernel(mu,0)[0],reduced_kernel(mu,0)[1])
        for kappa in (.1,.01,.001):
            kernel=reduced_kernel(mu,kappa)
            self.assertAlmostEqual(np.abs(kernel[0]-kernel[1]).sum()/2,kappa*mu/2,places=15)
        with self.assertRaises(ValueError): conditional_error_rows(5,Fraction(1,2),Fraction(1,4))

    def test_interval_and_float_error_calculations_agree(self):
        for n,kappa in ((10,Fraction(1,2)),(25,Fraction(1,4))):
            floats=fixed_budget_errors(n,kappa)
            intervals=fixed_budget_errors(n,kappa,certified=True)
            for value,interval in zip(floats,intervals): self.assertAlmostEqual(value,sum(interval.floats())/2,places=14)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ImperfectInterventionTestingTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":76,"budgets":[sufficient_budget(k) for k in (Fraction(1),Fraction(1,2),Fraction(1,4))],
            "retained_statistics":"M=number of h=-1 trials and K=number of positive parities in those trials; all N trials count",
            "robustness_scope":"Unknown constant kappa >= certified kappa0; independent success events; failure delivers the message; complex instrument unaffected",
            "without_positive_lower_bound_no_uniform_finite_guarantee":True,
            "quantum_theory_derived_from_cognition":False,"empirical_samples_collected":0,
            "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("imperfect_intervention_testing_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__": main()
