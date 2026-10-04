"""Round 74: sufficient three-outcome record and best existing intervention setting.

The hypotheses are the two specified implementations of round 73, not all real
and complex theories. Settings are independent of hidden flags. Fresh trials
and trusted intervention semantics are explicit assumptions.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path

import numpy as np

from bell_network_statistics import BELL_SIGNS, SETTINGS
from certified_intervals import Interval as I, SCALE
from communication_causal_audit import calibration_interval, calibrated_flagged_records
from cognitive_model_acceptance import information_gain
from finite_network_separation import finite_protocol_table
from one_bit_real_network import visibility_interval


OUTCOMES=tuple(product((-1,1),(-1,1),range(4),(-1,1)))  # h,a,b,c
INCREMENTS=(-1,0,1)


def signal_interval():
    bias,_,_=calibration_interval()
    return visibility_interval(5)**3*bias**2/I.exact(2).sqrt()


def signal_diagnostic():
    return sum(signal_interval().floats())/2


def informative_settings():
    return tuple((1,z) for z,(i,j,_) in enumerate(SETTINGS) if 1 in (i,j))


def parity_label(z,a,b,c):
    i,j,sign=SETTINGS[z]
    if 1 not in (i,j): raise ValueError("Charlie setting must contain Y.")
    direction_sign=1 if i==1 else sign
    return int(a*c*BELL_SIGNS[b][1]*direction_sign)


def record_increment(z,h,a,b,c):
    if h not in (-1,1): raise ValueError("Use a binary flag.")
    return 0 if h==1 else parity_label(z,a,b,c)


def reduced_kernel(mu=None,blocking_success=1.):
    if mu is None: mu=signal_diagnostic()
    if not 0<=mu<=1 or not 0<=blocking_success<=1: raise ValueError("Use probabilities/visibility in [0,1].")
    real_bias=(1-2*blocking_success)*mu
    return np.array([[(1-mu)/4,.5,(1+mu)/4],
                     [(1-real_bias)/4,.5,(1+real_bias)/4]])


@lru_cache(maxsize=2)
def complete_hypothesis_records(blocked):
    complex_table=finite_protocol_table()
    complex_flags={(h,)+key:p/2 for h in (-1,1) for key,p in complex_table.items()}
    real_flags=calibrated_flagged_records(deliver_message=not blocked)
    return complex_flags,real_flags


def setting_kernel(x,z,blocked=True):
    tables=complete_hypothesis_records(blocked)
    return np.array([[table[h,x,z,a,b,c] for h,a,b,c in OUTCOMES] for table in tables])


def reconstruction_kernel(z,mu=None):
    """Hypothesis-independent conditional distribution of full labels given S."""
    if mu is None: mu=signal_diagnostic()
    result=np.zeros((3,len(OUTCOMES)))
    for row,s in enumerate(INCREMENTS):
        for col,(h,a,b,c) in enumerate(OUTCOMES):
            k=parity_label(z,a,b,c)
            if s==0 and h==1: result[row,col]=(1+mu*k)/16
            elif s!=0 and h==-1 and k==s: result[row,col]=1/8
    return result


def reduction_report():
    mu=signal_interval()
    rows=[]
    for x,z in product(range(3),range(6)):
        kernel=setting_kernel(x,z)
        rows.append({"alice_setting":x,"charlie_setting":z,
                     "TV_diagnostic":float(np.abs(kernel[0]-kernel[1]).sum()/2),
                     "informative":(x,z) in informative_settings()})
    return {"signal_interval":mu.floats(),"signal_lower_exact":str(Fraction(mu.lo,SCALE)),
            "signal_upper_exact":str(Fraction(mu.hi,SCALE)),"informative_settings":informative_settings(),
            "reduced_outcomes":INCREMENTS,"reduced_kernel_diagnostic":reduced_kernel().tolist(),
            "single_trial_TV_interval":(mu/2).floats(),
            "single_trial_equal_prior_error_interval":(I.rational(1,2)-mu/4).floats(),
            "setting_comparisons":rows,"zero_increment_probability":"1/2",
            "zero_increment_trials_are_counted":True,"all_existing_settings_globally_optimal_claim":True,
            "optimization_scope":"Fixed number of fresh trials; the 18 old settings and observe/block actions; no other interventions or new measurements"}


class InterventionProbeReductionTests(unittest.TestCase):
    def test_exactly_four_existing_settings_have_hypothesis_dependent_records(self):
        found=[]
        for x,z in product(range(3),range(6)):
            table=setting_kernel(x,z)
            if np.max(np.abs(table[0]-table[1]))>1e-12: found.append((x,z))
        self.assertEqual(tuple(found),informative_settings())

    def test_all_passive_settings_have_zero_information(self):
        for x,z in product(range(3),range(6)):
            table=setting_kernel(x,z,False)
            np.testing.assert_allclose(table[0],table[1],atol=1e-15)
            self.assertLess(abs(information_gain((.5,.5),table)),1e-14)

    def test_actual_full_circuits_reduce_to_the_three_outcome_likelihood(self):
        for x,z in informative_settings():
            table=setting_kernel(x,z)
            reduced=np.zeros((2,3))
            for col,(h,a,b,c) in enumerate(OUTCOMES):
                reduced[:,INCREMENTS.index(record_increment(z,h,a,b,c))]+=table[:,col]
            np.testing.assert_allclose(reduced,reduced_kernel(),atol=3e-15)

    def test_full_labels_are_reconstructed_by_a_common_stochastic_channel(self):
        for x,z in informative_settings():
            reconstruct=reconstruction_kernel(z)
            np.testing.assert_allclose(reconstruct.sum(axis=1),1,atol=1e-16)
            self.assertGreaterEqual(reconstruct.min(),0)
            np.testing.assert_allclose(reduced_kernel()@reconstruct,setting_kernel(x,z),atol=1e-15)

    def test_uninformative_experiments_are_postprocessings_of_the_informative_one(self):
        for x,z in product(range(3),range(6)):
            if (x,z) not in informative_settings():
                target=setting_kernel(x,z)
                discard=np.repeat(target[:1],3,axis=0)
                np.testing.assert_allclose(reduced_kernel()@discard,target,atol=1e-15)

    def test_likelihood_ratio_is_a_power_of_one_constant_for_perfect_blocking(self):
        mu=signal_diagnostic()
        ratio=(1+mu)/(1-mu)
        table=reduced_kernel()
        for s,c,r in zip(INCREMENTS,*table): self.assertAlmostEqual(c/r,ratio**s,places=14)
        self.assertAlmostEqual(np.abs(table[0]-table[1]).sum()/2,mu/2,places=15)

    def test_information_gain_is_positive_for_every_nontrivial_prior_and_zero_without_blocking(self):
        for prior in ((.01,.99),(.5,.5),(.95,.05)):
            self.assertGreater(information_gain(prior,reduced_kernel()),0)
            self.assertAlmostEqual(information_gain(prior,reduced_kernel(blocking_success=0)),0,places=15)

    def test_signal_certificate_retains_finite_instrument_errors(self):
        mu=signal_interval()
        ideal=I.exact(2).sqrt()/2
        self.assertGreater(mu.lo,I.rational(707,1000).hi)
        self.assertLess(mu.hi,ideal.lo)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(InterventionProbeReductionTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":74,**reduction_report(),"scope":"Two calibrated specified implementations, trusted perfect blocking, original retained interface and fresh trials",
            "all_real_models_excluded":False,"quantum_theory_derived_from_cognition":False,
            "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("intervention_probe_reduction_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in report.items() if k!='setting_comparisons'},indent=2))


if __name__ == "__main__": main()
