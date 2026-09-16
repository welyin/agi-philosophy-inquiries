"""Round 84: preserve old protocols and certify a strict finite-task gain.

One uniform four-Bell-state input; no extra shared entangled resource.
Ideal separable/LOCC discrimination is bounded by 1/2. One old real
inverse Bell gate and two old noisy readouts strictly exceed that bound.
"""

import argparse
import json
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from bell_network_statistics import bell_basis, bell_projectors, bob_effects
from certified_intervals import Interval as I, SCALE
from one_bit_real_network import visibility_interval
from operational_effect_closure import majority_certificate, majority_effect_by_circuit
from quantum_interface_audit import ALPHA, branch_kraus


def discrimination_table(effects):
    return np.array([[np.trace(rho@effect).real for effect in effects] for rho in bell_projectors()])


def best_equal_prior_success(table):
    return np.max(table,axis=0).sum()/table.shape[0]


def finite_gain_certificate(count):
    visibility = visibility_interval(count)
    success = ((1+visibility)/2)**2
    margin = success-Fraction(1,2)
    return {
        "pointer_reads_per_bit":count,
        "old_readouts":2*count,
        "joint_inverse_bell_gates":1,
        "success_interval":success.floats(),
        "success_lower_exact":str(Fraction(success.lo,SCALE)),
        "gain_over_ideal_separable_bound_lower_exact":str(Fraction(margin.lo,SCALE)),
        "strictly_beats_one_half":margin.lo>0,
        "scope":"One promised uniform Bell input; exact old control; independent old pointer resets when count>1"}


def local_basis_povm(left,right):
    return tuple(np.kron(np.outer(left[:,i],left[:,i].conj()),
                         np.outer(right[:,j],right[:,j].conj())) for i,j in product(range(2),repeat=2))


class CognitiveBoundaryExpansionTests(unittest.TestCase):
    def test_old_protocol_probabilities_are_preserved_by_embedding_and_ignoring_partner(self):
        first = np.array([[.7,.2],[.2,.3]])
        partner = np.array([[.4,-.1],[-.1,.6]])
        for record in product((0,1),repeat=3):
            local,whole = first.copy(),np.kron(first,partner)
            for step,outcome in enumerate(record):
                kraus = branch_kraus(.2*step+.1*sum(record[:step]),outcome,.6)
                local = sum(k@local@k.conj().T for k in kraus)
                lifts = [np.kron(k,np.eye(2)) for k in kraus]
                whole = sum(k@whole@k.conj().T for k in lifts)
            np.testing.assert_allclose(whole,np.kron(local,partner),atol=2e-16)
            self.assertAlmostEqual(np.trace(whole).real,np.trace(local).real)

    def test_classical_cooperation_already_can_enlarge_a_single_observers_boundary(self):
        parity_zero = np.array([[Fraction(1,2),0],[0,Fraction(1,2)]],dtype=object)
        parity_one = np.array([[0,Fraction(1,2)],[Fraction(1,2),0]],dtype=object)
        np.testing.assert_array_equal(parity_zero.sum(axis=0),parity_one.sum(axis=0))
        np.testing.assert_array_equal(parity_zero.sum(axis=1),parity_one.sum(axis=1))
        distance = sum(abs(x-y) for x,y in zip(parity_zero.flat,parity_one.flat))/2
        self.assertEqual(distance,1)

    def test_four_bell_states_and_the_decoder_are_entirely_real(self):
        basis = bell_basis()
        np.testing.assert_array_equal(basis.imag,np.zeros((4,4)))
        np.testing.assert_allclose(basis.T@basis,np.eye(4),atol=3e-16)
        for b,rho in enumerate(bell_projectors()):
            np.testing.assert_allclose(basis.T@rho@basis,np.diag(np.eye(4)[b]),atol=3e-16)
            reduced = np.einsum("abcb->ac",rho.reshape(2,2,2,2))
            np.testing.assert_allclose(reduced,np.eye(2)/2,atol=2e-16)

    def test_product_overlap_bound_for_real_and_complex_positive_effects(self):
        rng = np.random.default_rng(84)
        for complex_entries in (False,True):
            for _ in range(20):
                x,y = rng.normal(size=(2,2)),rng.normal(size=(2,2))
                if complex_entries:
                    x = x+1j*rng.normal(size=(2,2))
                    y = y+1j*rng.normal(size=(2,2))
                effect = np.kron(x@x.conj().T,y@y.conj().T)
                for rho in bell_projectors():
                    self.assertLessEqual(np.trace(rho@effect).real,np.trace(effect).real/2+2e-14)

    def test_ideal_local_z_attains_one_half_and_other_local_bases_obey_it(self):
        self.assertAlmostEqual(best_equal_prior_success(discrimination_table(local_basis_povm(np.eye(2),np.eye(2)))),.5)
        rng = np.random.default_rng(184)
        for _ in range(10):
            left,_ = np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))
            right,_ = np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))
            self.assertLessEqual(best_equal_prior_success(discrimination_table(local_basis_povm(left,right))),.5+3e-16)

    def test_full_noisy_decoder_records_match_two_independent_bit_channels(self):
        for visibility in (0.,.4,ALPHA,1.):
            table = discrimination_table(bob_effects(visibility))
            for b,out in product(range(4),repeat=2):
                distance = (b^out).bit_count()
                expected = ((1+visibility)/2)**(2-distance)*((1-visibility)/2)**distance
                self.assertAlmostEqual(table[b,out],expected,places=15)
            np.testing.assert_allclose(table.sum(axis=1),np.ones(4),atol=4e-16)
            self.assertAlmostEqual(best_equal_prior_success(table),((1+visibility)/2)**2)

    def test_actual_finite_pointer_effects_and_strict_gain_certificates(self):
        basis = bell_basis()
        for count in (1,3,5):
            positive = majority_effect_by_circuit(count)
            effects = [basis@np.kron(positive if b//2==0 else np.eye(2)-positive,
                                     positive if b%2==0 else np.eye(2)-positive)@basis.T for b in range(4)]
            visibility = majority_certificate(count)["effective_visibility_diagnostic"]
            np.testing.assert_allclose(discrimination_table(effects),discrimination_table(bob_effects(visibility)),atol=4e-15)
            self.assertTrue(finite_gain_certificate(count)["strictly_beats_one_half"])
        self.assertGreater(Fraction(finite_gain_certificate(1)["success_lower_exact"]),Fraction(98,100))

    def test_new_distinguishable_coordinate_is_not_required_for_complex_task_gain(self):
        from role_symmetry_and_swap import pauli_word
        local_products = [pauli_word(a+b) for a,b in product("IXYZ",repeat=2)]
        pairing = np.array([[np.trace(a@b).real/4 for b in local_products] for a in local_products])
        np.testing.assert_array_equal(pairing,np.eye(16))
        # Complex local tomography is already complete, but the one-copy task
        # still has the independently proved separable upper bound of 1/2.
        self.assertGreater(best_equal_prior_success(discrimination_table(bob_effects(ALPHA))),.98)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(CognitiveBoundaryExpansionTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report = {
        "round":84,
        "candidate_principle":"Faithful embedding of old protocols, plus a strict gain on at least one declared finite task",
        "every_possible_combination_must_have_strict_gain":False,
        "old_protocol_embedding":"rho -> rho tensor sigma; K -> K tensor I; effect -> effect tensor I",
        "uniform_four_bell_input_ideal_separable_success_upper_exact":"1/2",
        "finite_real_joint_decoder":[finite_gain_certificate(m) for m in (1,3,5)],
        "strict_task_gain_requires_new_state_coordinates":False,
        "preservation_plus_strict_gain_alone_selects_complex_structure":False,
        "additional_shared_entanglement_for_separated_baseline":False,
        "postselected_trial_outputs":False,
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("cognitive_boundary_expansion_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__": main()
