"""Round 106: correlated environment responses invisible to all pair marginals."""
import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np

from bilocal_record_tomography import embed_operator
from certified_intervals import Interval as I
from complex_whole_real_interfaces import random_state
from deferred_relation_queries import memory_channel
from environment_query_recovery import interval_abs,interval_max
from incompatible_relation_writes import G,H,calibration_state
from independent_source_alignment import keep_systems
from joint_relation_records import memory_axis,interval_fields
from memory_dephasing_threshold import dephasing_unitary,dephase_memory
from partial_environment_access import fragment_unitary,partition_overlaps
from quantum_interface_audit import rotation_unitary
from role_symmetry_and_swap import pauli_word
from weak_relation_tradeoff import trace_distance


def correlated_environment(kappa):
    if not -1<=kappa<=1: raise ValueError("Use kappa in [-1,1].")
    return (np.eye(8)+kappa*pauli_word("XYY"))/8


def conditional_rotation(overlaps):
    result=np.ones((1,1))
    for lam in overlaps:
        if not 0<=lam<=1: raise ValueError("Use overlaps in [0,1].")
        result=np.kron(result,rotation_unitary(2*math.acos(lam)))
    return result


def environment_response(sigma,overlaps,accessible):
    """Return sigma_S and C_S=tr_U[sigma (I_S tensor R_U^dagger)]."""
    accessible=tuple(accessible)
    partition_overlaps(overlaps,accessible)
    n=len(overlaps)
    if sigma.shape!=(2**n,2**n): raise ValueError("Environment dimension mismatch.")
    unavailable=tuple(i for i in range(n) if i not in accessible)
    ordered=keep_systems(sigma,accessible+unavailable,n)
    rest=conditional_rotation([overlaps[i] for i in unavailable])
    reduced=keep_systems(ordered,range(len(accessible)),n)
    response=keep_systems(ordered@np.kron(np.eye(2**len(accessible)),rest.conj().T),
                          range(len(accessible)),n)
    return reduced,response


def response_channel(state,c,sigma_s,response):
    axis=np.kron(memory_axis(c),np.eye(len(state)//2))
    plus=(np.eye(len(state))+axis)/2
    minus=(np.eye(len(state))-axis)/2
    return (np.kron(sigma_s,plus@state@plus+minus@state@minus)
            +np.kron(response,plus@state@minus)
            +np.kron(response.conj().T,minus@state@plus))


def correlated_joint(state,c,sigma,overlaps):
    u=fragment_unitary(c,overlaps,int(round(math.log2(len(state)))))
    return u@np.kron(sigma,state)@u.conj().T


def correlated_accessible(state,c,sigma,overlaps,accessible,undo=False):
    accessible=tuple(accessible)
    n=len(overlaps);k=len(accessible)
    rest=int(round(math.log2(len(state))))
    joint=correlated_joint(state,c,sigma,overlaps)
    reduced=keep_systems(joint,accessible+tuple(range(n,n+rest)),n+rest)
    if undo:
        u=fragment_unitary(c,[overlaps[i] for i in accessible],rest)
        reduced=u.conj().T@reduced@u
    return reduced


def response_H_distance(c,d,sigma_s,response):
    block=np.block([[c*sigma_s,response],[response.conj().T,c*sigma_s]])
    return math.sqrt(1-d*d)*np.abs(np.linalg.eigvalsh(block)).sum()/2


def correlated_factor(c,kappa,lam1,lam2):
    b=lam1*lam2
    q=kappa*math.sqrt((1-lam1*lam1)*(1-lam2*lam2))
    return (max(c,abs(b+q))+max(c,abs(b-q)))/2


def correlation_certificate(c,d,kappa,overlaps):
    c,d,kappa=map(Fraction,(c,d,kappa));overlaps=tuple(map(Fraction,overlaps))
    if len(overlaps)!=3: raise ValueError("Three fragments required.")
    if not 0<=c<=1 or not 0<=d<=1 or not -1<=kappa<=1: raise ValueError("Invalid parameters.")
    partition_overlaps(overlaps,(0,))
    _,l1,l2=overlaps
    b=I.exact(l1*l2)
    q=I.exact(kappa)*((1-I.exact(l1*l1))*(1-I.exact(l2*l2))).sqrt()
    factor=(interval_max(I.exact(c),interval_abs(b+q))
            +interval_max(I.exact(c),interval_abs(b-q)))/2
    w=(1-I.exact(d*d)).sqrt()
    return {"c_exact":str(c),"d_exact":str(d),"kappa_exact":str(kappa),
        "overlaps_exact":list(map(str,overlaps)),"available_fragment":0,
        "all_environment_discarded_coherence_exact":str(math.prod(overlaps)),
        "accessible_H_distance":interval_fields(w*factor),
        "uncorrelated_environment_accessible_H_distance":interval_fields(w*I.exact(max(c,l1*l2))),
        "all_one_and_two_fragment_environment_marginals":"maximally mixed for every input before and after coupling",
        "source_is_mixed_and_requires_preparation_or_purification":True,
        "old_AB_available_to_recovery":False}


class CorrelatedEnvironmentResponseTests(unittest.TestCase):
    def test_general_response_identity_includes_complex_environment_and_reference(self):
        rng=np.random.default_rng(106)
        sigma=random_state(rng,8);state=random_state(rng,8)
        overlaps=(.3,.7,.9)
        for subset in ((0,),(2,0),()):
            marginal,response=environment_response(sigma,overlaps,subset)
            np.testing.assert_allclose(correlated_accessible(state,.8,sigma,overlaps,subset,True),
                response_channel(state,.8,marginal,response),atol=2e-15)

    def test_three_fragment_source_is_positive_real_and_has_equal_all_proper_marginals(self):
        for kappa in (-1.,-.4,0.,1.):
            sigma=correlated_environment(kappa)
            np.testing.assert_allclose(sigma.imag,0)
            self.assertGreaterEqual(np.linalg.eigvalsh(sigma).min(),-1e-15)
            for k in (1,2):
                for subset in combinations(range(3),k):
                    np.testing.assert_allclose(keep_systems(sigma,subset,3),np.eye(2**k)/2**k,atol=1e-15)

    def test_actual_recording_leaves_all_proper_environment_marginals_identical_for_any_input(self):
        state=random_state(np.random.default_rng(206),4)
        for kappa in (0.,1.):
            joint=correlated_joint(state,.8,correlated_environment(kappa),(.9,.9,.9))
            for count in (1,2):
                for subset in combinations(range(3),count):
                    np.testing.assert_allclose(keep_systems(joint,subset,5),np.eye(2**count)/2**count,atol=2e-15)

    def test_entire_memory_channel_is_identical_despite_correlation(self):
        state=random_state(np.random.default_rng(306),8)
        for kappa in (-1.,0.,.7,1.):
            actual=correlated_accessible(state,.8,correlated_environment(kappa),(.9,.8,.7),())
            np.testing.assert_allclose(actual,dephase_memory(state,.8,.9*.8*.7),atol=2e-15)

    def test_response_contains_XYY_correlation_and_H_bound_matches_actual_states(self):
        for c,d,kappa,overlaps in ((.8,.8,0,(.9,.9,.9)),(.8,.8,1,(.9,.9,.9)),
                (.8,.8,1,(0.,0.,0.)),(.3,.7,-.6,(.2,.7,.8))):
            sigma=correlated_environment(kappa)
            marginal,response=environment_response(sigma,overlaps,(0,))
            b=overlaps[1]*overlaps[2]
            q=kappa*math.sqrt((1-overlaps[1]**2)*(1-overlaps[2]**2))
            np.testing.assert_allclose(response,(b*np.eye(2)-q*pauli_word("X"))/2,atol=1e-15)
            pair=[correlated_accessible(memory_channel(calibration_state(H,s),c,d),c,sigma,overlaps,(0,))
                  for s in (-1,1)]
            expected=math.sqrt(1-d*d)*correlated_factor(c,kappa,*overlaps[1:])
            self.assertAlmostEqual(trace_distance(*pair),expected,places=13)
            self.assertAlmostEqual(response_H_distance(c,d,marginal,response),expected,places=13)

    def test_two_active_real_qubit_rigidity_on_every_hidden_real_basis_direction(self):
        overlaps=(.6,.8)
        # Fixed single marginals leave XX,XZ,ZX,ZZ,YY. Fixed total channel fixes YY.
        for word in ("XX","XZ","ZX","ZZ"):
            delta=pauli_word(word)/4
            for subset in ((0,),(1,)):
                marginal,response=environment_response(delta,overlaps,subset)
                np.testing.assert_allclose(marginal,0,atol=1e-15)
                np.testing.assert_allclose(response,0,atol=1e-15)
        delta=pauli_word("YY")/4
        self.assertAlmostEqual(np.trace(delta@conditional_rotation(overlaps).conj().T).real,
                               -math.sqrt(1-.6**2)*math.sqrt(1-.8**2))

    def test_G_distinguishability_still_retains_its_original_value(self):
        for kappa in (0.,1.):
            pair=[correlated_accessible(memory_channel(calibration_state(G,s),.8,.8),.8,
                    correlated_environment(kappa),(.9,.9,.9),(0,)) for s in (-1,1)]
            self.assertAlmostEqual(trace_distance(*pair),.6,places=13)

    def test_counterexample_has_strict_rational_gap_and_does_not_assume_a_record_in_each_fragment(self):
        f=Fraction
        lo=correlation_certificate(f(4,5),f(4,5),0,[f(9,10)]*3)
        hi=correlation_certificate(f(4,5),f(4,5),1,[f(9,10)]*3)
        self.assertGreater(Fraction(hi["accessible_H_distance"]["lower_exact"])
                           -Fraction(lo["accessible_H_distance"]["upper_exact"]),f(53,1000))
        self.assertTrue(lo["source_is_mixed_and_requires_preparation_or_purification"])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(CorrelatedEnvironmentResponseTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    f=Fraction
    report={"round":106,"source_family":"sigma_kappa=(I+kappa*XYY)/8",
        "fixed_recovery_response":"sigma_S=I/2; C_S=(lambda1*lambda2*I-kappa*s1*s2*X)/2",
        "universal_accessible_H_distance":"w/2*(max(c,abs(B+Q))+max(c,abs(B-Q)))",
        "minimality_scope":"Three real qubit fragments, fixed controlled Ry couplings, every coupling active; fixed single marginals and full memory channel",
        "nonzero_G_record_in_each_fragment":False,
        "examples":[correlation_certificate(f(4,5),f(4,5),k,[lam]*3)
                    for lam,k in ((f(9,10),0),(f(9,10),1),(f(0),0),(f(0),1))],
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("correlated_environment_response_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
