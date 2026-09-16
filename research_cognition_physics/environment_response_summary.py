"""Round 108: a future-sufficient environment response and its sharp calibration error."""
import argparse
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I
from complex_whole_real_interfaces import random_state
from correlated_environment_response import (environment_response,response_channel,
    correlated_accessible,correlated_factor)
from correlated_local_recovery import correlated_environment
from independent_source_alignment import keep_systems
from joint_relation_records import memory_axis,axis_read_memory_tree,interval_fields
from one_bit_real_network import visibility_interval
from operational_effect_closure import compiled_copy_gate
from quantum_interface_audit import rotation_unitary
from role_symmetry_and_swap import pauli_word
from weak_relation_tradeoff import trace_distance


def hidden_source(kappa,epsilon):
    """XYY matters to the fixed E0 interface; ZZZ does not. Both may be present."""
    if kappa*kappa+epsilon*epsilon>1+1e-14: raise ValueError("Need kappa^2+epsilon^2 <=1.")
    return (np.eye(8)+kappa*pauli_word("XYY")+epsilon*pauli_word("ZZZ"))/8


def probe_state(c):
    vector=rotation_unitary(math.acos(c))[:,0]
    return np.outer(vector,vector.conj())


def response_probe_decoder(c):
    # A controlled real parity transfer: U^dagger Z_E U = X_E tensor L_M.
    return compiled_copy_gate()@np.kron(rotation_unitary(-math.pi/2),rotation_unitary(-math.acos(c)))


def response_probe_records(c,kappa,overlaps,count):
    state=correlated_accessible(probe_state(c),c,correlated_environment(kappa),overlaps,(0,),True)
    decode=response_probe_decoder(c)
    return axis_read_memory_tree(decode@state@decode.conj().T,0,0.,count)


def calibrated_channel_distance(kappa1,kappa2,lam1,lam2):
    return abs(kappa1-kappa2)*math.sqrt((1-lam1*lam1)*(1-lam2*lam2))/2


def calibration_certificate(kappa1,kappa2,lam1,lam2,count):
    k1,k2,l1,l2=map(Fraction,(kappa1,kappa2,lam1,lam2))
    if not all(-1<=x<=1 for x in (k1,k2)) or not all(0<=x<=1 for x in (l1,l2)):
        raise ValueError("Invalid source or coupling values.")
    q=((1-I.exact(l1*l1))*(1-I.exact(l2*l2))).sqrt()
    ideal=I.exact(abs(k1-k2))*q/2
    finite=visibility_interval(count)*ideal
    return {"kappa_1_exact":str(k1),"kappa_2_exact":str(k2),
        "unavailable_overlaps_exact":[str(l1),str(l2)],
        "maximum_accessible_channel_output_trace_distance":interval_fields(ideal),
        "compiled_probe_record_total_variation":interval_fields(finite),
        "equal_prior_source_discrimination_success":interval_fields((1+finite)/2),
        "probe_preparation":"New M in the +1 eigenspace of L, no original AB write needed",
        "readouts":count,
        "scope":"Same fixed couplings, E0 and M jointly accessible, unavailable fragments and purifiers stay closed"}


class EnvironmentResponseSummaryTests(unittest.TestCase):
    def test_distinct_positive_sources_with_same_response_give_identical_accessible_channels(self):
        rng=np.random.default_rng(108);state=random_state(rng,8)
        overlaps=(.9,.7,.8)
        sources=[hidden_source(.4,e) for e in (-.6,.6)]
        self.assertGreater(trace_distance(*sources),.5)
        for sigma in sources: self.assertGreaterEqual(np.linalg.eigvalsh(sigma).min(),-1e-15)
        summaries=[environment_response(sigma,overlaps,(0,)) for sigma in sources]
        for first,second in zip(*summaries): np.testing.assert_allclose(first,second,atol=1e-15)
        outputs=[correlated_accessible(state,.8,sigma,overlaps,(0,)) for sigma in sources]
        np.testing.assert_allclose(*outputs,atol=2e-15)

    def test_diagonal_and_coherent_probe_outputs_reconstruct_the_entire_response_pair(self):
        sigma=random_state(np.random.default_rng(208),8)
        marginal,response=environment_response(sigma,(.3,.7,.9),(0,))
        c=.8;beta=math.pi/2+math.acos(c)
        basis=rotation_unitary(beta)
        plus=np.outer(basis[:,0],basis[:,0])
        coherent=probe_state(c)  # L maps to -X in this A eigenbasis.
        first=response_channel(plus,c,marginal,response)
        np.testing.assert_allclose(keep_systems(first,(0,),2),marginal,atol=1e-15)
        u=np.kron(np.eye(2),basis)
        output=u.conj().T@response_channel(coherent,c,marginal,response)@u
        off=output.reshape(2,2,2,2)[:,0,:,1]
        np.testing.assert_allclose(-2*off,response,atol=1e-15)

    def test_accessible_channel_difference_obeys_sharp_bound_with_complex_spectators(self):
        rng=np.random.default_rng(308);overlaps=(.3,.7,.9)
        for size in (2,4,8):
            for k1,k2 in ((0.,1.),(-1.,1.),(-.4,.7)):
                rho=random_state(rng,size)
                pair=[correlated_accessible(rho,.8,correlated_environment(k),overlaps,(0,)) for k in (k1,k2)]
                self.assertLessEqual(trace_distance(*pair),
                                     calibrated_channel_distance(k1,k2,*overlaps[1:])+2e-15)

    def test_pure_L_probe_saturates_the_channel_bound_for_all_tested_couplings(self):
        for overlaps in ((.9,.9,.9),(.3,0.,0.),(.2,.7,1.)):
            pair=[correlated_accessible(probe_state(.8),.8,correlated_environment(k),overlaps,(0,))
                  for k in (-.6,.7)]
            self.assertAlmostEqual(trace_distance(*pair),calibrated_channel_distance(-.6,.7,*overlaps[1:]),places=13)

    def test_real_parity_decoder_and_actual_noisy_records_have_signed_kappa_mean(self):
        c=.8;overlaps=(.9,.7,.8);q=math.sqrt((1-.7**2)*(1-.8**2))
        l=c*pauli_word("Z")+math.sqrt(1-c*c)*pauli_word("X")
        u=response_probe_decoder(c)
        np.testing.assert_allclose(u.conj().T@pauli_word("ZI")@u,np.kron(pauli_word("X"),l),atol=1e-15)
        for kappa in (-1.,.4,1.):
            for count in (1,3):
                records=response_probe_records(c,kappa,overlaps,count)
                gamma=sum(visibility_interval(count).floats())/2
                mean=sum((1 if sum(h)>count//2 else -1)*np.trace(rho).real for h,rho in records.items())
                self.assertAlmostEqual(mean,-gamma*kappa*q,places=13)

    def test_full_noisy_probe_record_distance_attains_certified_visibility_factor(self):
        f=Fraction
        for count in (1,3):
            pair=[response_probe_records(.8,k,(.9,.9,.9),count) for k in (-1.,1.)]
            tv=sum(abs(np.trace(pair[0][h]-pair[1][h]).real) for h in pair[0])/2
            cert=calibration_certificate(-1,1,f(9,10),f(9,10),count)
            self.assertAlmostEqual(tv,cert["compiled_probe_record_total_variation"]["diagnostic"],places=13)

    def test_equal_best_H_scores_do_not_define_equal_operational_states(self):
        for kappa in (-1.,1.):
            self.assertAlmostEqual(correlated_factor(.8,kappa,0.,0.),1)
        pair=[correlated_accessible(probe_state(.8),.8,correlated_environment(k),(.9,0.,0.),(0,))
              for k in (-1.,1.)]
        self.assertAlmostEqual(trace_distance(*pair),1,places=13)
        cert=calibration_certificate(-1,1,0,0,1)
        self.assertGreater(Fraction(cert["equal_prior_source_discrimination_success"]["lower_exact"]),
                           Fraction(994,1000))

    def test_response_approximation_bound_includes_both_diagonal_and_offdiagonal_parts(self):
        rng=np.random.default_rng(408)
        summaries=[environment_response(random_state(rng,8),(.3,.7,.9),(0,)) for _ in range(2)]
        delta_sigma=summaries[0][0]-summaries[1][0]
        delta_c=summaries[0][1]-summaries[1][1]
        bound=(np.linalg.svd(delta_sigma,compute_uv=False).sum()+np.linalg.svd(delta_c,compute_uv=False).sum())/2
        for size in (2,8):
            rho=random_state(rng,size)
            pair=[response_channel(rho,.8,*summary) for summary in summaries]
            self.assertLessEqual(trace_distance(*pair),bound+1e-15)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(EnvironmentResponseSummaryTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    f=Fraction
    report={"round":108,"sufficient_and_necessary_summary":"(sigma_S,C_S) for fixed couplings and fixed permanently closed remainder",
        "general_channel_error_bound":"0.5*(trace_norm(delta_sigma_S)+trace_norm(delta_C_S))",
        "sharp_XYY_family_channel_error":"0.5*abs(delta_kappa)*sqrt((1-lambda1^2)*(1-lambda2^2))",
        "same_H_optimum_does_not_identify_future_record_statistics":True,
        "same_response_does_not_require_same_global_source":True,
        "examples":[calibration_certificate(a,b,l,l,m) for a,b,l,m in
                    ((0,1,f(9,10),1),(-1,1,f(9,10),3),(-1,1,f(0),1))],
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("environment_response_summary_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
