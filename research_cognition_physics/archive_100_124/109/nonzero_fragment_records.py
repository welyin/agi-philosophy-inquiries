"""Round 109: equal nonzero local records with different hidden recovery power."""
import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I
from complex_whole_real_interfaces import random_state
from correlated_environment_response import (environment_response,correlated_joint,
    correlated_accessible,response_H_distance)
from deferred_relation_queries import memory_channel
from environment_query_recovery import interval_max
from incompatible_relation_writes import G,H,calibration_state
from independent_source_alignment import keep_systems
from joint_relation_records import axis_read_memory_tree,interval_fields
from memory_dephasing_threshold import dephase_memory
from one_bit_real_network import visibility_interval
from role_symmetry_and_swap import pauli_word
from weak_relation_tradeoff import trace_distance


def validate_biases(biases):
    if len(biases)!=3 or any(not 0<=r<=1 for r in biases):
        raise ValueError("Use three biases in [0,1].")


def correlation_capacity(biases):
    validate_biases(biases)
    return math.sqrt(math.prod(1-r*r for r in biases))


def biased_environment(biases,kappa):
    validate_biases(biases)
    if kappa*kappa>math.prod(1-r*r for r in biases)+1e-14:
        raise ValueError("The source requires kappa^2 <= product(1-r_i^2).")
    diagonal=np.ones((1,1))
    for r in biases: diagonal=np.kron(diagonal,np.eye(2)+r*pauli_word("Z"))
    return (diagonal+kappa*pauli_word("XYY"))/8


def biased_H_factor(c,r0,kappa,lam1,lam2):
    b=lam1*lam2
    k=kappa*math.sqrt((1-lam1*lam1)*(1-lam2*lam2))
    return sum(max(abs(a),math.sqrt(a*a*r0*r0+k*k)) for a in (c+b,c-b))/2


def local_record_certificate(c,r,lam,count):
    c,r,lam=map(Fraction,(c,r,lam))
    if any(not 0<=x<=1 for x in (c,r,lam)): raise ValueError("Invalid record parameter.")
    distance=I.exact(r)*((1-I.exact(c*c))*(1-I.exact(lam*lam))).sqrt()
    return {"bias_exact":str(r),"gate_cosine_exact":str(lam),"readouts":count,
        "local_G_ideal_distance":interval_fields(distance),
        "local_G_finite_success":interval_fields((1+visibility_interval(count)*distance)/2)}


def biased_certificate(c,d,biases,kappa,overlaps):
    c,d,kappa=map(Fraction,(c,d,kappa))
    biases=tuple(map(Fraction,biases));overlaps=tuple(map(Fraction,overlaps))
    validate_biases(biases)
    if len(overlaps)!=3 or any(not 0<=x<=1 for x in (c,d,*overlaps)):
        raise ValueError("Invalid write parameters.")
    capacity_squared=math.prod(1-r*r for r in biases)
    if kappa*kappa>capacity_squared: raise ValueError("Exact positivity condition fails.")
    r0=biases[0];b=overlaps[1]*overlaps[2]
    k2=kappa*kappa*(1-overlaps[1]**2)*(1-overlaps[2]**2)
    total=I.exact(0)
    for a in (c+b,c-b):
        total+=interval_max(I.exact(abs(a)),I.exact(a*a*r0*r0+k2).sqrt())/2
    w=(1-I.exact(d*d)).sqrt()
    return {"c_exact":str(c),"d_exact":str(d),"biases_exact":list(map(str,biases)),
        "kappa_exact":str(kappa),"gate_cosines_exact":list(map(str,overlaps)),
        "positive_source_capacity_squared_exact":str(capacity_squared),
        "source_positivity_margin_exact":str(capacity_squared-kappa*kappa),
        "all_environment_discarded_coherence_exact":str(math.prod(overlaps)),
        "accessible_H_ideal_distance":interval_fields(w*total),
        "independent_source_H_ideal_distance":interval_fields(w*I.exact(max(c,b))),
        "strict_H_distance_gain":interval_fields(w*(total-I.exact(max(c,b)))),
        "local_records":[local_record_certificate(c,r,lam,1) for r,lam in zip(biases,overlaps)],
        "old_records_are_readable_but_not_previously_measured_classical_histories":True}


class NonzeroFragmentRecordTests(unittest.TestCase):
    def test_complement_blocks_give_exact_positivity_condition(self):
        biases=(.6,.4,.8)
        capacity=correlation_capacity(biases)
        for kappa in (0.,capacity,-capacity):
            sigma=biased_environment(biases,kappa)
            self.assertGreaterEqual(np.linalg.eigvalsh(sigma).min(),-2e-16)
            for b in range(4):
                block=sigma[np.ix_((b,7-b),(b,7-b))]
                self.assertAlmostEqual(np.linalg.det(block).real,(capacity**2-kappa**2)/64,places=14)
        with self.assertRaises(ValueError): biased_environment(biases,capacity+.001)
        with self.assertRaises(ValueError):
            biased_certificate(Fraction(4,5),Fraction(4,5),[Fraction(3,5)]*3,Fraction(513,1000),[Fraction(9,10)]*3)

    def test_all_proper_initial_and_postcoupling_marginals_match_despite_nonzero_bias(self):
        biases=(.6,.4,.8);overlaps=(.9,.7,.8)
        state=random_state(np.random.default_rng(109),4)
        sources=[biased_environment(biases,k) for k in (0.,.2)]
        joints=[correlated_joint(state,.8,sigma,overlaps) for sigma in sources]
        for size in (1,2):
            for subset in combinations(range(3),size):
                np.testing.assert_allclose(keep_systems(sources[0],subset,3),
                                           keep_systems(sources[1],subset,3),atol=1e-15)
                np.testing.assert_allclose(keep_systems(joints[0],subset,5),
                                           keep_systems(joints[1],subset,5),atol=2e-15)

    def test_every_fragment_carries_the_same_strictly_nonzero_G_record_in_both_sources(self):
        biases=(.6,.4,.8);overlaps=(.9,.7,.8)
        for kappa in (0.,.2):
            pair=[correlated_joint(memory_channel(calibration_state(G,s),.8,.8),.8,
                    biased_environment(biases,kappa),overlaps) for s in (-1,1)]
            for i in range(3):
                local=[keep_systems(rho,(i,),5) for rho in pair]
                expected=.6*biases[i]*math.sqrt(1-overlaps[i]**2)
                self.assertAlmostEqual(trace_distance(*local),expected,places=13)
                self.assertGreater(expected,0.)

    def test_actual_noisy_single_fragment_readout_attains_its_record_certificate(self):
        f=Fraction;c=.8;biases=(.6,)*3;overlaps=(.9,)*3
        success=0.
        for s in (-1,1):
            joint=correlated_joint(memory_channel(calibration_state(G,s),c,.8),c,
                                   biased_environment(biases,.512),overlaps)
            local=keep_systems(joint,(1,),5)
            records=axis_read_memory_tree(local,0,-math.asin(.9),1)
            success+=sum(np.trace(rho).real for h,rho in records.items()
                         if (1 if sum(h)>0 else -1)==s)/2
        cert=local_record_certificate(f(4,5),f(3,5),f(9,10),1)
        self.assertAlmostEqual(success,cert["local_G_finite_success"]["diagnostic"],places=13)
        self.assertGreater(Fraction(cert["local_G_finite_success"]["lower_exact"]),f(57,100))

    def test_full_memory_channel_with_complex_reference_is_identical(self):
        state=random_state(np.random.default_rng(209),8)
        for kappa in (0.,.3,.512):
            actual=correlated_accessible(state,.8,biased_environment((.6,)*3,kappa),(.9,)*3,())
            np.testing.assert_allclose(actual,dephase_memory(state,.8,.9**3),atol=2e-15)

    def test_general_response_and_actual_H_states_match_closed_formula(self):
        for biases,kappa,overlaps in (((.6,)*3,.512,(.9,)*3),((.8,.2,.4),.3,(.7,.8,.9)),
                                     ((.2,)*3,.8,(0.,)*3),((1.,.2,.3),0.,(.4,.5,.6))):
            sigma=biased_environment(biases,kappa)
            marginal,response=environment_response(sigma,overlaps,(0,))
            expected_marginal=(np.eye(2)+biases[0]*pauli_word("Z"))/2
            k=kappa*math.sqrt((1-overlaps[1]**2)*(1-overlaps[2]**2))
            np.testing.assert_allclose(marginal,expected_marginal,atol=1e-15)
            np.testing.assert_allclose(response,overlaps[1]*overlaps[2]*marginal-k*pauli_word("X")/2,atol=1e-15)
            pair=[correlated_accessible(memory_channel(calibration_state(H,s),.8,.8),.8,
                    sigma,overlaps,(0,)) for s in (-1,1)]
            predicted=.6*biased_H_factor(.8,biases[0],kappa,*overlaps[1:])
            self.assertAlmostEqual(trace_distance(*pair),predicted,places=13)
            self.assertAlmostEqual(response_H_distance(.8,.8,marginal,response),predicted,places=13)

    def test_nonzero_record_counterexample_has_strict_rational_H_gain(self):
        f=Fraction
        cert=biased_certificate(f(4,5),f(4,5),[f(3,5)]*3,f(64,125),[f(9,10)]*3)
        self.assertGreater(Fraction(cert["strict_H_distance_gain"]["lower_exact"]),f(26,1000))
        self.assertEqual(cert["source_positivity_margin_exact"],"0")
        self.assertTrue(all(Fraction(x["local_G_ideal_distance"]["lower_exact"])>0 for x in cert["local_records"]))

    def test_any_pure_local_marginal_eliminates_this_hidden_correlation(self):
        for biases in ((1.,.6,.6),(.6,1.,.6),(.6,.6,1.)):
            self.assertEqual(correlation_capacity(biases),0)
            with self.assertRaises(ValueError): biased_environment(biases,.001)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(NonzeroFragmentRecordTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    f=Fraction
    report={"round":109,"source":"[tensor_i(I+r_i Z_i)+kappa XYY]/8",
        "exact_positivity":"kappa^2 <= product_i(1-r_i^2)",
        "local_G_distance":"sqrt(1-c^2)*r_i*sqrt(1-lambda_i^2)",
        "all_proper_environment_marginals_and_full_memory_channel_equal_across_kappa":True,
        "examples":[biased_certificate(f(4,5),f(4,5),[f(3,5)]*3,k,[f(9,10)]*3)
                    for k in (f(0),f(64,125))],
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("nonzero_fragment_records_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
