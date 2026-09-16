"""Round 115: identical calibration histories, opposite future response, and robust transfer."""
import argparse
import json
import math
import unittest
from fractions import Fraction as F
from itertools import product
from pathlib import Path

import numpy as np

from biased_source_compiler import mapped_pulses, phase_pulses, pulse_unitary, biased_source_preparation
from certified_intervals import Interval as I, SCALE
from independent_source_alignment import keep_systems
from joint_relation_records import interval_fields
from learned_correlation_loop import loop_distribution, learned_query_records
from incompatible_relation_writes import H, calibration_state
from nonzero_fragment_records import biased_environment, correlation_capacity
from unknown_correlation_learning import confidence_set, matched_policy_value, calibration_visibility


CAPACITY=F(64,125)


def conditional_source_pulses(biases,tau,alignment=1):
    """Label Z controls the source sign; label is NOT initialized here."""
    if not -1<=tau<=1 or alignment not in (-1,1):
        raise ValueError("Physical amplitude and sign required.")
    pulses=[]
    for i,r in enumerate(biases):
        if not 0<=r<=1: raise ValueError("Invalid bias.")
        pulses+=mapped_pulses("Y",math.acos(r),(i+1,),7)
    pulses+=phase_pulses("seed",(1,2,3),7)
    for control,other in ((4,2),(5,3)):
        pulses+=mapped_pulses("Y",math.pi/2,(control,),7)
        pulses+=phase_pulses("twirl",(control,1,other),7)
    pulses+=mapped_pulses("Y",math.pi/2,(6,),7)
    pulses+=mapped_pulses("ZY",-alignment*math.asin(tau),(0,6),7)
    pulses+=phase_pulses("cz",(6,1,2),7)
    return pulses


def shared_sign_source(calibration_copies,kappa,alignment=1):
    if not isinstance(calibration_copies,int) or not 1<=calibration_copies<=2:
        raise ValueError("Small explicit matrix check supports one or two calibration copies.")
    if alignment not in (-1,1) or not 0<=kappa<=float(CAPACITY):
        raise ValueError("Invalid parameters.")
    total=2**(3*(calibration_copies+1))
    out=np.zeros((total,total),complex)
    for z in (-1,1):
        term=np.array([[1.]],complex)
        for sign in [z]*calibration_copies+[alignment*z]:
            term=np.kron(term,biased_environment((.6,)*3,sign*kappa))
        out+=term/2
    return out


def shared_history_probability(history,kappa=CAPACITY):
    kappa=F(kappa)
    if not 0<=kappa<=CAPACITY or any(y not in (-1,1) for y in history):
        raise ValueError("Invalid source or binary history.")
    mean=calibration_visibility(1)*I.exact(kappa)
    result=I.exact(0)
    for z in (-1,1):
        term=I.exact(1)
        for y in history: term*= (1+y*z*mean)/2
        result+=term/2
    return result


def expanded_future_interval(current,radius,capacity=CAPACITY):
    current=I.exact(current);radius=F(radius);capacity=F(capacity)
    if radius<0 or not 0<capacity<=1: raise ValueError("Invalid future uncertainty.")
    r=I.exact(radius);cap=I.exact(capacity)
    lo=max(-cap.hi,current.lo-r.hi)
    hi=min(cap.hi,current.hi+r.hi)
    if lo>hi: return I(-cap.hi,cap.hi)
    return I(lo,hi)


def robust_action(future_interval,target=F(0)):
    """Minimax over mixtures of the existing +, -, and fallback policies."""
    interval=I.exact(future_interval);target=F(target)
    if target<0: raise ValueError("Nonnegative target required.")
    h=I.exact(target)
    if interval.lo>0 and interval.lo>=h.hi: return 1
    if interval.hi<0 and interval.hi<=-h.hi: return -1
    return 0


def robust_gain_floor(future_interval,action):
    _,coefficient=matched_policy_value(0,0)
    interval=I.exact(future_interval)
    if action==1: lower=F(interval.lo,SCALE)
    elif action==-1: lower=-F(interval.hi,SCALE)
    elif action==0: lower=F(0)
    else: raise ValueError("Invalid action.")
    return coefficient*I.exact(lower)


def drift_record_certificate(n,positives,look,drift,target=F(1,5)):
    current,conflict=confidence_set(n,positives,look,CAPACITY,1)
    future=expanded_future_interval(current,drift)
    action=0 if conflict else robust_action(future,target)
    return {"samples":n,"positive_count":positives,"look":look,"assumed_transfer_radius_exact":str(drift),
        "calibration_interval":interval_fields(current),"future_interval":interval_fields(future),
        "action":action,"gain_floor_on_coverage_and_bridge":interval_fields(robust_gain_floor(future,action)),
        "transfer_radius_is_an_assumption_not_inferred_from_past_data":True}


def counterexample_certificate():
    dist=loop_distribution(CAPACITY)
    base,coef=matched_policy_value(0,0)
    gain=coef*I.exact(CAPACITY)*(dist["plus"]-dist["minus"])
    return {"single_copy_unconditional_source":"sigma_0 for every calibration and future copy",
        "full_past_joint_source_and_all_past_instrument_records_identical":True,
        "aligned_future_success":interval_fields(base+gain),
        "opposite_future_success":interval_fields(base-gain),
        "fallback_success":interval_fields(base),
        "sign_rule":"hidden fair label Z: calibration kappa=Z*R, future kappa=alignment*Z*R",
        "shared_source_preparation_per_copy":{"pair_rotations":30,"local_rotations":17,"pure_slots":6},
        "shared_label_extra":{"pure_slots":1,"local_rotations":1,"label_retained_and_inaccessible":True},
        "full_1024_sample_cap_plus_query_budget":{"pair_rotations":35879,"local_rotations":44081,
                                                  "fresh_pure_slots":7177},
        "scope":"No access to source purifiers or label; future source initially independent of old AB and H label"}


class DriftingSourceTransferTests(unittest.TestCase):
    def test_correlated_wholes_have_identical_past_and_every_single_copy(self):
        pair=[shared_sign_source(2,.512,a) for a in (-1,1)]
        np.testing.assert_allclose(keep_systems(pair[0],tuple(range(6)),9),
                                   keep_systems(pair[1],tuple(range(6)),9),atol=1e-15)
        for rho in pair:
            self.assertAlmostEqual(np.trace(rho).real,1,places=14)
            for copy in range(3):
                np.testing.assert_allclose(keep_systems(rho,tuple(range(3*copy,3*copy+3)),9),
                                           biased_environment((.6,)*3,0.),atol=1e-15)
        self.assertGreater(np.linalg.norm(pair[0]-pair[1]),.01)

    def test_shared_label_is_prepared_by_real_native_gates_with_counted_cost(self):
        for alignment in (-1,1):
            pulses=conditional_source_pulses((.6,)*3,.7,alignment)
            unitary=pulse_unitary(pulses,7)
            expected=np.zeros((128,128),complex)
            expected[:64,:64]=biased_source_preparation((.6,)*3,alignment*.7)
            expected[64:,64:]=biased_source_preparation((.6,)*3,-alignment*.7)
            np.testing.assert_allclose(unitary,expected,atol=4e-15)
            counts=[sum(x!="I" for x in word) for word,angle in pulses]
            self.assertEqual((counts.count(2),counts.count(1)),(30,17))

    def test_entire_past_record_law_normalizes_and_is_not_iid_sigma_zero(self):
        total=sum((shared_history_probability(h) for h in product((-1,1),repeat=5)),I.exact(0))
        self.assertLessEqual(total.lo,SCALE);self.assertGreaterEqual(total.hi,SCALE)
        self.assertGreater(shared_history_probability((1,1)).lo,I.rational(1,4).hi)

    def test_opposite_future_reverses_actual_query_gain_not_just_a_label(self):
        for source_sign,action in ((1,1),(-1,-1)):
            values=[]
            for alignment in (-1,1):
                success=0.
                for label in (-1,1):
                    records=learned_query_records(calibration_state(H,label),alignment*source_sign*.512,action)
                    success+=sum(np.trace(rho).real for guess,rho in records.values() if guess==label)/2
                values.append(success)
            base,_=matched_policy_value(0,0)
            self.assertLess(values[0],base.floats()[0]);self.assertGreater(values[1],base.floats()[1])

    def test_full_learning_counterexample_keeps_all_error_and_timeout_mass(self):
        cert=counterexample_certificate()
        self.assertLess(F(cert["opposite_future_success"]["upper_exact"]),F(723,1000))
        self.assertGreater(F(cert["aligned_future_success"]["lower_exact"]),F(752,1000))

    def test_robust_interval_selects_or_falls_back_when_transfer_margin_expires(self):
        safe=drift_record_certificate(256,194,3,F(1,25))
        unsafe=drift_record_certificate(256,194,3,F(3,50))
        self.assertEqual(safe["action"],1);self.assertEqual(unsafe["action"],0)
        self.assertGreater(F(safe["gain_floor_on_coverage_and_bridge"]["lower_exact"]),F(6,1000))

    def test_minimax_among_randomized_existing_policies_is_attained_by_interval_rule(self):
        for lower,upper in ((F(1,10),F(2,5)),(-F(2,5),-F(1,10)),(-F(1,10),F(2,5))):
            interval=I(I.exact(lower).lo,I.exact(upper).hi)
            chosen=robust_action(interval)
            exact_best=max(F(0),lower,-upper)
            for j in range(-100,101):
                x=F(j,100)
                self.assertLessEqual(min(x*lower,x*upper),exact_best)
            self.assertEqual(min(chosen*lower,chosen*upper),exact_best)

    def test_parameter_independent_fallback_handles_unrestricted_future_sign(self):
        interval=I(-I.exact(CAPACITY).hi,I.exact(CAPACITY).hi)
        self.assertEqual(robust_action(interval),0)
        self.assertEqual(robust_gain_floor(interval,0).lo,0)
        with self.assertRaises(ValueError): expanded_future_interval(interval,-1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(DriftingSourceTransferTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":115,"counterexample":counterexample_certificate(),
        "robust_transfer_rule":"Expand the estimated conditional parameter interval by an assumed bridge radius, then choose the worst-case positive-gain policy.",
        "examples":[drift_record_certificate(256,194,3,d) for d in (F(1,25),F(3,50))],
        "quantum_theory_derived_from_cognition":False,"automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("drifting_source_transfer_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
