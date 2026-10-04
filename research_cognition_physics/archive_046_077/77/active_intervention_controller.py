"""Round 77: finite predictive controller with a certified stopping rule.

Perfect trusted blocking only. Accumulate D=sum S; stop at +/-3 or a fixed cap,
then choose by sign (fair tie). Exact dynamic programming includes every early
decision and cap branch. This identifies two implementations, not a number field.
"""

import argparse
import json
import unittest
from dataclasses import dataclass, FrozenInstanceError
from fractions import Fraction
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I, SCALE
from cognitive_model_acceptance import information_gain, predictive_update
from communication_causal_audit import calibration_interval
from intervention_probe_reduction import INCREMENTS, informative_settings, reduced_kernel, signal_interval, signal_diagnostic


def stopping_statistics(threshold,cap,mu,true_hypothesis="complex"):
    if type(threshold) is not int or threshold<1 or type(cap) is not int or cap<0:
        raise ValueError("Use a positive integer threshold and nonnegative cap.")
    if true_hypothesis not in ("complex","real"): raise ValueError("Unknown hypothesis.")
    drift=mu if true_hypothesis=="complex" else -mu
    weights=((-1,(1-drift)/4),(0,Fraction(1,2)),(1,(1+drift)/4))
    states={0:mu*0+1}
    left=right=expected=mu*0
    for _ in range(cap):
        expected+=sum(states.values(),mu*0)
        next_states={}
        for balance,mass in states.items():
            for step,weight in weights:
                destination=balance+step
                if destination<=-threshold: left+=mass*weight
                elif destination>=threshold: right+=mass*weight
                else: next_states[destination]=next_states.get(destination,mu*0)+mass*weight
        states=next_states
    if true_hypothesis=="complex": error=left+sum((p for d,p in states.items() if d<0),mu*0)
    else: error=right+sum((p for d,p in states.items() if d>0),mu*0)
    error+=states.get(0,mu*0)/2
    return {"error":error,"expected_trials":expected,"unresolved_before_cap_decision":sum(states.values(),mu*0),
            "absorbed_left":left,"absorbed_right":right,"live":states}


def unbounded_statistics(threshold,mu):
    ratio=(1+mu)/(1-mu)
    error=1/(1+ratio**threshold)
    expected=2*threshold*(1-2*error)/mu
    return error,expected


def stopping_certificate(threshold=3):
    for cap in range(1,101):
        if stopping_statistics(threshold,cap,signal_diagnostic())["error"]<.01: break
    else: raise RuntimeError("No certified cap candidate found.")
    stats=stopping_statistics(threshold,cap,signal_interval())
    if stats["error"].hi>=I.rational(1,100).lo: raise ArithmeticError("Risk certificate unresolved.")
    infinite_error,infinite_expected=unbounded_statistics(threshold,signal_interval())
    _,_,flip=calibration_interval()
    coin=(flip.lo+flip.hi)//2
    coin_error=Fraction(max(coin-flip.lo,flip.hi-coin),SCALE)
    finite_coin_error=Fraction(stats["error"].hi,SCALE)+cap*coin_error
    return {"integer_evidence_threshold":threshold,"maximum_trials":cap,
            "each_conditional_error_interval":stats["error"].floats(),
            "each_error_upper_exact":str(Fraction(stats["error"].hi,SCALE)),
            "expected_trials_interval":stats["expected_trials"].floats(),
            "expected_trials_upper_exact":str(Fraction(stats["expected_trials"].hi,SCALE)),
            "cap_survival_probability_interval":stats["unresolved_before_cap_decision"].floats(),
            "unbounded_threshold_error_interval":infinite_error.floats(),
            "unbounded_threshold_expected_trials_interval":infinite_expected.floats(),
            "finite_100_bit_calibration_each_error_upper_exact":str(finite_coin_error),
            "finite_calibration_error_below_one_percent":finite_coin_error<Fraction(1,100),
            "cap_policy":"At cap, choose by sign of D and use a fair coin if D=0; all cap branches count",
            "each_record_has_99_percent_posterior_claimed":False,"globally_optimal_sequential_test_claimed":False}


@dataclass(frozen=True)
class ProbeRecord:
    trial:int
    action:str
    setting:tuple
    prediction:tuple
    increment:int
    posterior:tuple
    balance:int


class ActiveProbeController:
    def __init__(self,threshold=3,cap=20,allow_block=True):
        if type(threshold) is not int or threshold<1 or type(cap) is not int or cap<1: raise ValueError("Use positive integer limits.")
        self.threshold,self.cap,self.allow_block=threshold,cap,bool(allow_block)
        self.balance=0
        self.history=()

    @property
    def posterior(self):
        ratio=(1+signal_diagnostic())/(1-signal_diagnostic())
        complex_weight=1/(1+ratio**(-self.balance))
        return complex_weight,1-complex_weight

    def action(self):
        if abs(self.balance)>=self.threshold or len(self.history)>=self.cap: return "decide"
        if not self.allow_block: return "abstain"
        return "block"

    def observe(self,increment):
        if self.action()!="block": raise RuntimeError("No permitted pending probe.")
        if increment not in INCREMENTS: raise ValueError("Unknown reduced outcome.")
        table=reduced_kernel()
        prior=self.posterior
        prediction=tuple(np.asarray(prior)@table)
        _,posterior=predictive_update(prior,table,INCREMENTS.index(increment))
        self.balance+=increment
        record=ProbeRecord(len(self.history)+1,"block",informative_settings()[0],prediction,increment,tuple(posterior),self.balance)
        self.history=self.history+(record,)
        return record

    def decision(self):
        if self.action()!="decide": return None
        if self.balance>0: return "complex"
        if self.balance<0: return "real"
        return "fair_coin"


class ActiveInterventionControllerTests(unittest.TestCase):
    def test_capped_dynamic_program_matches_complete_small_stopping_tree(self):
        mu=Fraction(1,2)
        weights=((-1,(1-mu)/4),(0,Fraction(1,2)),(1,(1+mu)/4))
        for threshold,cap in ((1,3),(2,5)):
            leaves=[]
            def visit(balance,used,mass):
                if abs(balance)>=threshold or used==cap:
                    error=1 if balance<0 else (Fraction(1,2) if balance==0 else 0)
                    leaves.append((mass,used,error))
                    return
                for step,p in weights: visit(balance+step,used+1,mass*p)
            visit(0,0,Fraction(1))
            exact=stopping_statistics(threshold,cap,mu)
            self.assertEqual(sum(p for p,_,_ in leaves),1)
            self.assertEqual(sum(p*n for p,n,_ in leaves),exact["expected_trials"])
            self.assertEqual(sum(p*e for p,_,e in leaves),exact["error"])

    def test_both_hypotheses_have_the_same_conditional_error_and_sample_cost(self):
        first=stopping_statistics(3,20,Fraction(3,5))
        second=stopping_statistics(3,20,Fraction(3,5),"real")
        self.assertEqual(first["error"],second["error"])
        self.assertEqual(first["expected_trials"],second["expected_trials"])
        self.assertEqual(first["absorbed_left"],second["absorbed_right"])

    def test_unbounded_hitting_formulas_solve_the_absorbing_chain(self):
        mu=.6
        threshold=3
        states=list(range(-threshold+1,threshold))
        transition=np.zeros((len(states),len(states)))
        wrong=np.zeros(len(states))
        for i,d in enumerate(states):
            for step,p in ((-1,(1-mu)/4),(0,.5),(1,(1+mu)/4)):
                if d+step==-threshold: wrong[i]+=p
                elif d+step in states: transition[i,states.index(d+step)]+=p
        errors=np.linalg.solve(np.eye(len(states))-transition,wrong)
        costs=np.linalg.solve(np.eye(len(states))-transition,np.ones(len(states)))
        expected_error,expected_cost=unbounded_statistics(threshold,mu)
        self.assertAlmostEqual(errors[states.index(0)],expected_error,places=14)
        self.assertAlmostEqual(costs[states.index(0)],expected_cost,places=13)

    def test_finite_cap_certifies_risk_and_reduces_expected_trials_below_fixed_budget(self):
        report=stopping_certificate()
        self.assertLess(Fraction(report["each_error_upper_exact"]),Fraction(1,100))
        self.assertLess(Fraction(report["expected_trials_upper_exact"]),9)
        self.assertLessEqual(report["maximum_trials"],25)
        self.assertTrue(report["finite_calibration_error_below_one_percent"])

    def test_predictions_are_logged_before_updates_and_audit_records_are_immutable(self):
        controller=ActiveProbeController()
        expected=np.array([.5,.5])@reduced_kernel()
        record=controller.observe(1)
        np.testing.assert_allclose(record.prediction,expected,atol=1e-16)
        np.testing.assert_allclose(record.posterior,controller.posterior,atol=1e-15)
        with self.assertRaises(FrozenInstanceError): record.increment=-1

    def test_evidence_changes_stopping_action_and_zero_outcomes_still_use_budget(self):
        controller=ActiveProbeController(cap=20)
        for s in (1,-1,0,1,1): controller.observe(s)
        self.assertEqual(controller.action(),"block")
        self.assertEqual(len(controller.history),5)
        controller.observe(1)
        self.assertEqual(controller.action(),"decide")
        self.assertEqual(controller.decision(),"complex")
        with self.assertRaises(RuntimeError): controller.observe(1)
        tied=ActiveProbeController(cap=2)
        tied.observe(0); tied.observe(0)
        self.assertEqual(tied.decision(),"fair_coin")

    def test_hard_intervention_veto_abstains_and_cannot_be_overridden_by_prediction(self):
        controller=ActiveProbeController(allow_block=False)
        self.assertEqual(controller.action(),"abstain")
        with self.assertRaises(RuntimeError): controller.observe(1)
        self.assertIsNone(controller.decision())
        self.assertEqual(controller.posterior,(.5,.5))

    def test_chosen_action_has_positive_information_at_every_live_belief(self):
        for balance in range(-2,3):
            controller=ActiveProbeController()
            for _ in range(abs(balance)): controller.observe(1 if balance>0 else -1)
            self.assertGreater(information_gain(controller.posterior,reduced_kernel()),0)
            self.assertAlmostEqual(information_gain(controller.posterior,reduced_kernel(blocking_success=0)),0,places=14)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ActiveInterventionControllerTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":77,**stopping_certificate(),"scope":"Two specified calibrated implementations, fresh trials, equal priors, perfect trusted blocking; finite threshold controller",
            "controller_has_intervention_veto":True,"quantum_theory_derived_from_cognition":False,
            "empirical_samples_collected":0,"automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("active_intervention_controller_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__": main()
