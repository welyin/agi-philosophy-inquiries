"""Round 71: exact optimum for a noisy orientation message and real correction.

Sources, alignment instrument and final outer measurements are fixed. Charlie
may apply any real qubit CPTP map to the reference conditioned on the message.
This is not an optimization over arbitrary communicating network strategies.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I, SCALE
from independent_source_alignment import pair_reference, real_x_with_reset
from one_bit_real_network import real_feedback_score, records_for_reference, reference_flags_by_circuit, visibility_interval
from operational_effect_closure import majority_certificate
from quantum_interface_audit import PAULI_Y
from bell_network_statistics import network_score


def validate_channel(channel):
    w = np.asarray(channel,dtype=float)
    if w.ndim!=2 or w.shape[0]!=2 or w.shape[1]<1 or np.min(w)<0 or not np.allclose(w.sum(axis=1),1,atol=1e-14,rtol=0):
        raise ValueError("Use two probability rows, for the + and - input flags.")
    return w


def message_distinguishability(channel):
    w = validate_channel(channel)
    return float(np.abs(w[0]-w[1]).sum()/2)


def optimal_responses(channel):
    w = validate_channel(channel)
    return np.sign(w[0]-w[1])


def response_correlation(channel,responses,alignment_visibility=1.):
    w = validate_channel(channel)
    responses = np.asarray(responses,dtype=float)
    if responses.shape!=(w.shape[1],) or np.max(np.abs(responses))>1 or not 0<=alignment_visibility<=1:
        raise ValueError("Use one response in [-1,1] per received symbol and visibility in [0,1].")
    return float(alignment_visibility*np.dot(responses,w[0]-w[1])/2)


def corrected_reference_by_channel(channel,alignment_count=3,responses=None):
    w = validate_channel(channel)
    if responses is None: responses = optimal_responses(w)
    response_correlation(w,responses)
    result = np.zeros((4,4),dtype=complex)
    for row,flag in enumerate((1,-1)):
        state = reference_flags_by_circuit(alignment_count)[flag]
        flipped = real_x_with_reset(state,1,2)
        for s,kappa in enumerate(responses):
            result += w[row,s]*((1+kappa)*state+(1-kappa)*flipped)/2
    return result


def binary_symmetric_channel(error):
    if not 0<=error<=1: raise ValueError("Use an error probability in [0,1].")
    return np.array([[1-error,error],[error,1-error]])


def message_threshold_certificate(alignment_count=5,read_count=5):
    root = I.exact(2).sqrt()
    alignment,read = visibility_interval(alignment_count),visibility_interval(read_count)
    old_bound = 6*root-1/(108+95*root)
    base = 2*root*read**2*(read*(1+alignment)+read**2)
    penalty = 4*root*read**3*alignment
    threshold = (base-old_bound)/penalty
    rows=[]
    for error in (Fraction(0),Fraction(1,2000),Fraction(1,1000),Fraction(1,2)):
        score = base-penalty*error
        rows.append({"binary_flip_probability_exact":str(error),"score_interval":score.floats(),
                     "strictly_above_round_68_bound":score.lo>old_bound.hi,
                     "strictly_below_round_68_bound":score.hi<old_bound.lo})
    return {"alignment_reads":alignment_count,"reads_per_target":read_count,
            "binary_error_threshold_interval":threshold.floats(),
            "binary_error_threshold_lower_exact":str(Fraction(threshold.lo,SCALE)),
            "binary_error_threshold_upper_exact":str(Fraction(threshold.hi,SCALE)),
            "comparisons":rows,"threshold_refers_to_conservative_round_68_bound":True}


class NoisyMessageAlignmentTests(unittest.TestCase):
    def test_all_vertices_of_the_reference_response_cube_obey_the_exact_optimum(self):
        for w in (np.array([[.7,.2,.1],[.1,.3,.6]]),np.array([[.8,.2],[.3,.7]])):
            values=[response_correlation(w,response,.8) for response in product((-1,1),repeat=w.shape[1])]
            optimum=.8*message_distinguishability(w)
            self.assertAlmostEqual(max(values),optimum,places=15)
            self.assertAlmostEqual(response_correlation(w,optimal_responses(w),.8),optimum,places=15)

    def test_noisy_channel_and_actual_real_correction_give_the_predicted_reference(self):
        for w in (binary_symmetric_channel(.07),binary_symmetric_channel(.8),np.array([[.8,.2],[.3,.7]]),np.array([[.5,.5],[.5,.5]])):
            gamma=majority_certificate(3)["effective_visibility_diagnostic"]
            np.testing.assert_allclose(corrected_reference_by_channel(w),pair_reference(gamma*message_distinguishability(w)),atol=1e-15)

    def test_real_kraus_maps_cannot_have_reference_y_gain_larger_than_one(self):
        rng=np.random.default_rng(71)
        for _ in range(12):
            isometry,_=np.linalg.qr(rng.normal(size=(6,2)))
            kraus=isometry.reshape(3,2,2)
            np.testing.assert_allclose(sum(k.T@k for k in kraus),np.eye(2),atol=8e-16)
            kappa=sum(np.linalg.det(k) for k in kraus)
            np.testing.assert_allclose(sum(k.T@PAULI_Y@k for k in kraus),kappa*PAULI_Y,atol=2e-16)
            self.assertLessEqual(abs(kappa),1.+1e-15)

    def test_general_real_reference_channels_do_not_exceed_the_response_optimum(self):
        # Nonunital maps may change the real local marginal; the measured YY
        # correlation still depends only on sum(det K), not on that translation.
        rng=np.random.default_rng(171)
        w=np.array([[.7,.3],[.2,.8]])
        states=reference_flags_by_circuit(3)
        output=np.zeros((4,4),dtype=complex)
        kappas=[]
        for s in range(2):
            isometry,_=np.linalg.qr(rng.normal(size=(6,2)))
            kraus=isometry.reshape(3,2,2)
            kappas.append(sum(np.linalg.det(k) for k in kraus))
            for row,flag in enumerate((1,-1)):
                output+=w[row,s]*sum(np.kron(np.eye(2),k)@states[flag]@np.kron(np.eye(2),k.T) for k in kraus)
        gamma=majority_certificate(3)["effective_visibility_diagnostic"]
        corr=response_correlation(w,kappas,gamma)
        self.assertAlmostEqual(np.trace(output@np.kron(PAULI_Y,PAULI_Y)).real,corr,places=14)
        self.assertAlmostEqual(network_score(records_for_reference(output)),real_feedback_score(corr),places=13)
        self.assertLessEqual(corr,gamma*message_distinguishability(w))

    def test_binary_symmetric_and_erasure_channels_have_different_distinguishability(self):
        for error in (0.,.1,.5,1.):
            self.assertAlmostEqual(message_distinguishability(binary_symmetric_channel(error)),abs(1-2*error))
            erasure=np.array([[1-error,0,error],[0,1-error,error]])
            self.assertAlmostEqual(message_distinguishability(erasure),1-error)

    def test_further_classical_postprocessing_cannot_increase_message_distinguishability(self):
        rng=np.random.default_rng(271)
        for _ in range(15):
            w=rng.dirichlet(np.ones(3),size=2)
            post=rng.dirichlet(np.ones(4),size=3)
            self.assertLessEqual(message_distinguishability(w@post),message_distinguishability(w)+1e-15)

    def test_certified_binary_noise_threshold_separates_the_two_test_errors(self):
        report=message_threshold_certificate()
        lower=Fraction(report["binary_error_threshold_lower_exact"])
        upper=Fraction(report["binary_error_threshold_upper_exact"])
        self.assertGreater(lower,Fraction(1,2000))
        self.assertLess(upper,Fraction(1,1000))
        self.assertTrue(report["comparisons"][1]["strictly_above_round_68_bound"])
        self.assertTrue(report["comparisons"][2]["strictly_below_round_68_bound"])

    def test_invalid_channels_and_responses_are_rejected(self):
        for w in ([[1,1],[0,0]],[[1.1,-.1],[0,1]],[[1,0]]):
            with self.assertRaises(ValueError): validate_channel(w)
        with self.assertRaises(ValueError): response_correlation(np.eye(2),(1,1.1))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(NoisyMessageAlignmentTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":71,"exact_optimal_reference_correlation":"alignment_visibility * TV(W(.|+),W(.|-))",
            "scope":"Fixed sources, alignment and outer instruments; any message-conditioned real CPTP map on Charlie's reference qubit",
            **message_threshold_certificate(),"quantum_theory_derived_from_cognition":False,
            "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("noisy_message_alignment_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__": main()
