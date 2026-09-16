"""Round 111: exact record-correlation thresholds and a matched recovery setting."""
import argparse
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I
from biased_source_compiler import finite_biased_records
from correlated_environment_response import correlated_accessible
from environment_response_summary import probe_state,response_probe_decoder
from incompatible_relation_writes import H,calibration_state
from joint_relation_records import axis_read_memory_tree,interval_fields
from nonzero_fragment_records import (biased_environment,biased_H_factor,
    correlation_capacity,validate_biases)
from one_bit_real_network import visibility_interval
from weak_relation_tradeoff import trace_distance


def gain_margin(c,biases,kappa,lam1,lam2):
    return kappa*kappa*(1-lam1*lam1)*(1-lam2*lam2)-(1-biases[0]**2)*(c-lam1*lam2)**2


def threshold_certificate(c,biases,lam1,lam2):
    c,l1,l2=map(Fraction,(c,lam1,lam2));biases=tuple(map(Fraction,biases))
    validate_biases(biases)
    if any(not 0<=x<=1 for x in (c,l1,l2)): raise ValueError("Invalid parameters.")
    capacity_squared=math.prod(1-r*r for r in biases)
    margin=capacity_squared*(1-l1*l1)*(1-l2*l2)-(1-biases[0]**2)*(c-l1*l2)**2
    return {"c_exact":str(c),"biases_exact":list(map(str,biases)),
        "unavailable_gate_cosines_exact":[str(l1),str(l2)],
        "maximum_kappa_squared_exact":str(capacity_squared),
        "maximum_gain_test_margin_exact":str(margin),
        "some_source_in_this_family_has_strict_H_gain":margin>0,
        "scope":"H task requires d<1; comparison with same biases, same gates and kappa=0"}


def matched_certificate(c,d,r,kappa,m,q):
    c,d,r,kappa=map(Fraction,(c,d,r,kappa))
    if not Fraction(1,2)<=c<1 or not 0<=d<1 or not 0<=r<1 or kappa*kappa>(1-r*r)**3:
        raise ValueError("Use c>=1/2 and a physical uniform-bias source.")
    gamma=visibility_interval(m);eta=visibility_interval(q)
    k=I.exact(abs(kappa)*(1-c))
    factor=I.exact(c)+gamma*k/2
    w=(1-I.exact(d*d)).sqrt()
    ideal_factor=I.exact(c)+k/2
    local_G=I.exact(r)*((1-I.exact(c*c))*(1-I.exact(c))).sqrt()
    return {"c_exact":str(c),"d_exact":str(d),"uniform_bias_exact":str(r),"kappa_exact":str(kappa),
        "all_gate_cosines":"sqrt(c)","unavailable_cosine_product_exact":str(c),
        "environment_readouts":m,"final_H_readouts":q,
        "local_G_one_read_success":interval_fields((1+visibility_interval(1)*local_G)/2),
        "ideal_accessible_H_distance":interval_fields(w*ideal_factor),
        "finite_H_success":interval_fields((1+eta*w*factor)/2),
        "independent_source_same_recovery_success":interval_fields((1+eta*w*I.exact(c))/2),
        "strict_success_gain":interval_fields(eta*w*gamma*k/4),
        "ideal_X_environment_query_attains_accessible_bound":True,
        "angles_optimized_over_all_possible_protocols":False}


def biased_probe_records(c,biases,kappa,overlaps,count):
    rho=correlated_accessible(probe_state(c),c,biased_environment(biases,kappa),overlaps,(0,),True)
    decoder=response_probe_decoder(c)
    return axis_read_memory_tree(decoder@rho@decoder.conj().T,0,0.,count)


def response_error_certificate(biases,kappa1,kappa2,lam1,lam2,count):
    biases=tuple(map(Fraction,biases));k1,k2,l1,l2=map(Fraction,(kappa1,kappa2,lam1,lam2))
    validate_biases(biases)
    cap2=math.prod(1-r*r for r in biases)
    if k1*k1>cap2 or k2*k2>cap2: raise ValueError("Both sources must be physical.")
    if any(not 0<=x<=1 for x in (l1,l2)): raise ValueError("Invalid gate parameters.")
    coupling=((1-I.exact(l1*l1))*(1-I.exact(l2*l2))).sqrt()
    difference=I.exact(abs(k1-k2))*coupling/2
    capacity=I.exact(cap2).sqrt()*coupling
    return {"biases_exact":list(map(str,biases)),"kappa_1_exact":str(k1),"kappa_2_exact":str(k2),
        "sharp_channel_trace_distance":interval_fields(difference),
        "finite_probe_total_variation":interval_fields(visibility_interval(count)*difference),
        "largest_channel_difference_allowed_by_fixed_biases":interval_fields(capacity),
        "probe_total_pair_rotations":34+(count if count>1 else 0),
        "probe_extra_local_rotations":41+(3*count if count>1 else 0),
        "probe_pure_slots":7+int(count>1),"readouts":count,
        "scope":"Fixed biased-source family and known couplings; all closed fragments and purifiers remain in the whole"}


class RecordCorrelationFrontierTests(unittest.TestCase):
    def test_exact_gain_margin_matches_the_entire_closed_spectral_formula(self):
        rng=np.random.default_rng(111)
        for _ in range(90):
            c=rng.uniform(.05,.95);biases=tuple(rng.uniform(.05,.99,3))
            kappa=rng.uniform(-1,1)*correlation_capacity(biases)
            l1,l2=rng.uniform(0,.99,2)
            margin=gain_margin(c,biases,kappa,l1,l2)
            gain=biased_H_factor(c,biases[0],kappa,l1,l2)-max(c,l1*l2)
            if margin>1e-8: self.assertGreater(gain,0)
            if margin<-1e-8: self.assertAlmostEqual(gain,0,places=13)

    def test_equal_bias_transition_at_18_over_19_has_exact_rational_certificate(self):
        f=Fraction
        below=threshold_certificate(f(4,5),[f(97,100)]*3,f(9,10),f(9,10))
        above=threshold_certificate(f(4,5),[f(98,100)]*3,f(9,10),f(9,10))
        self.assertGreater(Fraction(below["maximum_gain_test_margin_exact"]),0)
        self.assertLess(Fraction(above["maximum_gain_test_margin_exact"]),0)
        # After cancelling the positive r0 factor: (1-r^2)*19 > 1.
        self.assertEqual(1-f(1,19),f(18,19))

    def test_record_strength_form_of_the_remaining_correlation_constraint(self):
        c=.8;biases=(.2,.6,.8);l1,l2=.7,.9
        v=math.sqrt(1-c*c)
        d1=v*biases[1]*math.sqrt(1-l1*l1)
        d2=v*biases[2]*math.sqrt(1-l2*l2)
        direct=correlation_capacity(biases)*math.sqrt((1-l1*l1)*(1-l2*l2))
        records=math.sqrt(1-biases[0]**2)*math.sqrt((1-l1*l1-d1*d1/(v*v))*(1-l2*l2-d2*d2/(v*v)))
        self.assertAlmostEqual(direct,records,places=14)

    def test_matched_gates_remove_nonzero_correlation_threshold_for_mixed_biases(self):
        for r in (.2,.6,.9,.99):
            c=.8;lam=math.sqrt(c)
            kappa=correlation_capacity((r,)*3)/3
            self.assertAlmostEqual(biased_H_factor(c,r,kappa,lam,lam),c+abs(kappa)*(1-c)/2,places=13)
            self.assertGreater(biased_H_factor(c,r,kappa,lam,lam),c)

    def test_matched_finite_actual_record_tree_attains_certified_formula(self):
        f=Fraction;c=.8;lam=math.sqrt(c)
        for m,q in ((1,1),(3,1)):
            cert=matched_certificate(f(4,5),f(4,5),f(3,5),f(64,125),m,q)
            success=0.
            for s in (-1,1):
                records=finite_biased_records(calibration_state(H,s),c,c,(.6,)*3,.512,(lam,)*3,m,q)
                success+=sum(np.trace(rho).real for guess,rho in records.values() if guess==s)/2
            self.assertAlmostEqual(success,cert["finite_H_success"]["diagnostic"],places=13)
            self.assertGreater(Fraction(cert["strict_success_gain"]["lower_exact"]),f(15,1000))

    def test_channel_error_sharpness_survives_nonzero_local_records_and_complex_reference(self):
        biases=(.6,)*3;overlaps=(.9,)*3;expected=.512*.19
        pair=[correlated_accessible(probe_state(.8),.8,biased_environment(biases,k),overlaps,(0,))
              for k in (-.512,.512)]
        self.assertAlmostEqual(trace_distance(*pair),expected,places=13)
        from complex_whole_real_interfaces import random_state
        rho=random_state(np.random.default_rng(211),8)
        pair=[correlated_accessible(rho,.8,biased_environment(biases,k),overlaps,(0,))
              for k in (-.512,.512)]
        self.assertLessEqual(trace_distance(*pair),expected+2e-15)

    def test_actual_signed_probe_reads_the_hidden_parameter_and_attains_variation_bound(self):
        f=Fraction
        for count in (1,3):
            pair=[biased_probe_records(.8,(.6,)*3,k,(.9,)*3,count) for k in (-.512,.512)]
            actual=sum(abs(np.trace(pair[0][h]-pair[1][h]).real) for h in pair[0])/2
            cert=response_error_certificate([f(3,5)]*3,-f(64,125),f(64,125),f(9,10),f(9,10),count)
            self.assertAlmostEqual(actual,cert["finite_probe_total_variation"]["diagnostic"],places=13)

    def test_unphysical_correlations_and_zero_hidden_capacity_are_not_counted_as_gain(self):
        f=Fraction
        with self.assertRaises(ValueError): matched_certificate(f(4,5),f(4,5),f(99,100),f(1,10),1,1)
        row=threshold_certificate(f(4,5),(f(1),f(3,5),f(3,5)),f(9,10),f(9,10))
        self.assertFalse(row["some_source_in_this_family_has_strict_H_gain"])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(RecordCorrelationFrontierTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    f=Fraction
    report={"round":111,
        "exact_gain_condition":"kappa^2*s1^2*s2^2 > (1-r0^2)*(c-lambda1*lambda2)^2",
        "equal_bias_example_threshold_squared_exact":"18/19",
        "matched_setting":"lambda1=lambda2=sqrt(c), c>=1/2; T_ideal=c+abs(kappa)*(1-c)/2",
        "threshold_examples":[threshold_certificate(f(4,5),[r]*3,f(9,10),f(9,10))
                              for r in (f(3,5),f(97,100),f(98,100))],
        "matched_examples":[matched_certificate(f(4,5),f(4,5),f(3,5),f(64,125),m,1) for m in (1,3)],
        "response_error_example":response_error_certificate([f(3,5)]*3,-f(64,125),f(64,125),f(9,10),f(9,10),1),
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("record_correlation_frontier_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
