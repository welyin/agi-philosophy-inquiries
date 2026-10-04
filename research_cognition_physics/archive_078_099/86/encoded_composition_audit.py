"""Round 86: composing independently encoded complex states is not free.

Exact real CP instrument; ideal reference-sector readout, not compiled here
into the project's noisy primitive. Both branches are retained in the audit.
"""

import argparse
import json
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from common_orientation_structure import encode_state,decode_state
from complex_control_from_reference import real_lift
from quantum_interface_audit import PAULI_X,PAULI_Y


def aligned_isometry(count):
    if count < 1: raise ValueError("At least one reference is required.")
    vector = np.ones(1,dtype=complex)
    for _ in range(count):
        vector = np.kron(vector,np.array([1,1j])/np.sqrt(2))
    return np.sqrt(2)*np.column_stack((vector.real,vector.imag))


def independent_encoding(states):
    """Order references first, then targets, without changing the state."""
    whole = np.ones((1,1))
    dimensions = []
    for state in states:
        whole = np.kron(whole,encode_state(state))
        dimensions.extend((2,len(state)))
    count = len(states)
    order = tuple(range(0,2*count,2))+tuple(range(1,2*count,2))
    full_order = order+tuple(i+2*count for i in order)
    return whole.reshape(tuple(dimensions)*2).transpose(full_order).reshape(whole.shape)


def product_state(states):
    whole = np.ones((1,1),dtype=complex)
    for state in states: whole = np.kron(whole,state)
    return whole


def two_input_instrument(first,second):
    d = len(first)*len(second)
    plus = aligned_isometry(2)
    minus = np.kron(np.eye(2),PAULI_X.real)@plus
    kraus = tuple(np.kron(v.T,np.eye(d)) for v in (plus,minus))
    input_state = independent_encoding((first,second))
    return kraus,tuple(k@input_state@k.T for k in kraus)


def aligned_success_branch(states):
    d = int(np.prod([len(state) for state in states]))
    k = np.kron(aligned_isometry(len(states)).T,np.eye(d))
    whole = independent_encoding(states)
    return k@whole@k.T


class EncodedCompositionAuditTests(unittest.TestCase):
    def test_aligned_reference_subspace_has_rank_two_for_every_tested_size(self):
        for n in range(1,6):
            v = aligned_isometry(n)
            np.testing.assert_allclose(v.T@v,np.eye(2),atol=9e-16)
            projector = v@v.T
            self.assertAlmostEqual(np.trace(projector),2)
            for i in range(1,n):
                factors = [np.eye(2)]*n
                factors[0] = factors[i] = PAULI_Y
                observable = product_state(factors)
                np.testing.assert_allclose(observable@v,v,atol=3e-16)

    def test_two_branches_form_a_complete_real_instrument(self):
        first,second = np.eye(2)/2,np.eye(3)/3
        kraus,branches = two_input_instrument(first,second)
        self.assertTrue(all(np.isrealobj(k) for k in kraus))
        np.testing.assert_allclose(sum(k.T@k for k in kraus),np.eye(24),atol=5e-16)
        self.assertAlmostEqual(sum(np.trace(b) for b in branches),1)

    def test_success_branch_is_the_desired_common_encoding_for_unknown_inputs(self):
        rng = np.random.default_rng(86)
        for da,db in ((2,2),(2,3),(3,2)):
            states=[]
            for d in (da,db):
                raw=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d))
                rho=raw@raw.conj().T
                rho=(rho+rho.conj().T)/2
                states.append(rho/np.trace(rho).real)
            _,branches=two_input_instrument(*states)
            self.assertAlmostEqual(np.trace(branches[0]),.5)
            np.testing.assert_allclose(branches[0],encode_state(np.kron(*states))/2,atol=8e-17)

    def test_other_branch_conjugates_one_logical_input_after_reference_alignment(self):
        first=np.array([[.6,.1+.2j],[.1-.2j,.4]])
        second=(np.eye(2)+.7*PAULI_Y)/2
        _,branches=two_input_instrument(first,second)
        self.assertAlmostEqual(np.trace(branches[1]),.5)
        np.testing.assert_allclose(branches[1],encode_state(np.kron(first,second.conj()))/2,atol=8e-17)

    def test_forgetting_the_branch_can_erase_an_old_logical_distinction(self):
        first=np.eye(2)/2
        outputs=[]
        for sign in (-1,1):
            second=(np.eye(2)+sign*PAULI_Y)/2
            _,branches=two_input_instrument(first,second)
            outputs.append(sum(branches))
            expected=encode_state(np.kron(first,second))
            self.assertAlmostEqual(np.abs(np.linalg.eigvalsh(sum(branches)-expected)).sum()/2,.5)
        np.testing.assert_array_equal(outputs[0],outputs[1])

    def test_retained_branch_allows_correct_interpretation_of_every_old_product_effect(self):
        first=np.array([[.6,.1+.2j],[.1-.2j,.4]])
        second=(np.eye(2)+.7*PAULI_Y)/2
        _,branches=two_input_instrument(first,second)
        rng=np.random.default_rng(186)
        for _ in range(12):
            v,w=rng.normal(size=2)+1j*rng.normal(size=2),rng.normal(size=2)+1j*rng.normal(size=2)
            v,w=v/np.linalg.norm(v),w/np.linalg.norm(w)
            a,b=np.outer(v,v.conj()),np.outer(w,w.conj())
            expected=np.trace(np.kron(first,second)@np.kron(a,b)).real
            actual=(np.trace(branches[0]@real_lift(np.kron(a,b)))+
                    np.trace(branches[1]@real_lift(np.kron(a,b.conj()))))
            self.assertAlmostEqual(actual,expected,places=14)

    def test_product_effect_reinterpretation_does_not_extend_to_all_joint_effects(self):
        vector=np.array([1.,0.,0.,1.])
        bell_effect=np.outer(vector,vector)/2
        transformed=bell_effect.reshape(2,2,2,2).transpose(0,3,2,1).reshape(4,4)
        singlet=np.array([0.,1.,-1.,0.])
        self.assertEqual(np.trace(transformed@np.outer(singlet,singlet)/2),-.5)

    def test_aligning_n_independent_encodings_has_the_exact_counted_success_probability(self):
        for n in (1,2,3,4):
            states=[(np.eye(2)+(j+1)/(n+1)*PAULI_Y)/2 for j in range(n)]
            branch=aligned_success_branch(states)
            probability=Fraction(1,2**(n-1))
            self.assertAlmostEqual(np.trace(branch),float(probability),places=14)
            np.testing.assert_allclose(branch,float(probability)*encode_state(product_state(states)),atol=1e-16)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(EncodedCompositionAuditTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={
        "round":86,
        "input":"Independent real encodings E(rho) tensor E(sigma)",
        "desired_common_encoding_branch":"E(rho tensor sigma)/2, unnormalized",
        "other_branch":"E(rho tensor conjugate(sigma))/2, unnormalized",
        "each_branch_probability_exact":"1/2",
        "all_branches_included_in_audit":True,
        "forgetting_branch_worst_example_trace_distance_exact":"1/2",
        "old_local_product_protocols_can_use_the_retained_orientation_flag":True,
        "same_reinterpretation_is_positive_for_every_joint_effect":False,
        "n_independent_inputs_this_filter_success_exact":"2^(1-n)",
        "ideal_reference_sector_readout":True,
        "compiled_into_old_noisy_primitives":False,
        "all_possible_deterministic_repair_protocols_excluded":False,
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("encoded_composition_audit_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__": main()
