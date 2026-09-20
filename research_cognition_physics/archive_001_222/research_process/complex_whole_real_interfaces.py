"""Round 87: real operational interfaces inside a declared complex whole.

Taking a real shadow is a prediction equivalence for restricted protocols,
not a physical state-erasure map. Local shadows do not compose by themselves.
"""

import argparse
import json
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from bell_network_statistics import bell_basis
from certified_intervals import SCALE
from independent_source_alignment import keep_systems
from one_bit_real_network import visibility_interval
from quantum_interface_audit import ALPHA, IDENTITY, PAULI_X, PAULI_Y, PAULI_Z
from role_symmetry_and_swap import pauli_word


def y_state(sign):
    return (IDENTITY+sign*PAULI_Y)/2


def real_shadow(state):
    """A representative of restricted-record equivalence, not a local channel."""
    return np.asarray(state).real.copy()


def pair_state(sign):
    return (np.eye(4)+sign*pauli_word("YY"))/4


def flagged_complex_whole(sign):
    """A,B carry correlated conjugate states; the existing C flag is retained."""
    return sum(np.kron(np.kron(y_state(h),y_state(sign*h)),
                       (IDENTITY+h*PAULI_Z)/2) for h in (-1,1))/2


def random_state(rng,d):
    raw=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d))
    state=raw@raw.conj().T
    state=(state+state.conj().T)/2
    return state/np.trace(state).real


def real_instrument(rng,d):
    raw=[rng.normal(size=(d,d)) for _ in range(2)]
    total=sum(k.T@k for k in raw)
    values,vectors=np.linalg.eigh(total)
    inverse=(vectors/np.sqrt(values))@vectors.T
    return [k@inverse for k in raw]


def adaptive_records(state,depth=3):
    leaves={():state}
    for step in range(depth):
        updated={}
        for history,branch in leaves.items():
            seed=8700+sum((bit+1)*3**i for i,bit in enumerate(history))
            instrument=real_instrument(np.random.default_rng(seed),len(state))
            for bit,k in enumerate(instrument):
                updated[history+(bit,)]=k@branch@k.T
        leaves=updated
    return {history:float(np.trace(branch).real) for history,branch in leaves.items()}


def finite_pair_success():
    basis=bell_basis()
    effects=[basis@np.kron(IDENTITY,(IDENTITY+s*ALPHA*PAULI_Z)/2)@basis.conj().T
             for s in (-1,1)]
    table=np.array([[np.trace(pair_state(s)@e).real for e in effects] for s in (-1,1)])
    return float(np.max(table,axis=0).sum()/2)


def finite_pair_success_lower():
    success=(1+visibility_interval(1))/2
    return Fraction(success.lo,SCALE)


class ComplexWholeRealInterfacesTests(unittest.TestCase):
    def test_real_shadow_is_positive_and_matches_every_real_effect(self):
        rng=np.random.default_rng(87)
        for d in (2,3,4):
            state=random_state(rng,d)
            shadow=real_shadow(state)
            self.assertGreater(np.linalg.eigvalsh(shadow).min(),0)
            self.assertAlmostEqual(np.trace(shadow),1)
            basis,_=np.linalg.qr(rng.normal(size=(d,d)))
            for i in range(d):
                effect=np.outer(basis[:,i],basis[:,i])
                self.assertAlmostEqual(np.trace(state@effect),np.trace(shadow@effect))

    def test_complete_finite_adaptive_real_records_match_without_state_erasure(self):
        for d in (2,3,4):
            state=random_state(np.random.default_rng(187+d),d)
            saved=state.copy()
            actual=adaptive_records(state)
            expected=adaptive_records(real_shadow(state))
            self.assertAlmostEqual(sum(actual.values()),1)
            for history in actual:
                self.assertAlmostEqual(actual[history],expected[history],places=14)
            np.testing.assert_array_equal(state,saved)
            self.assertGreater(np.linalg.norm(state.imag),0)

    def test_tensor_shadow_retains_product_of_imaginary_parts(self):
        rng=np.random.default_rng(287)
        for d,e in ((2,2),(2,3),(3,4)):
            rho,sigma=random_state(rng,d),random_state(rng,e)
            expected=np.kron(rho.real,sigma.real)-np.kron(rho.imag,sigma.imag)
            np.testing.assert_allclose(real_shadow(np.kron(rho,sigma)),expected,atol=2e-17)
            self.assertGreater(np.linalg.norm(expected-np.kron(rho.real,sigma.real)),0)

    def test_same_local_shadows_do_not_fix_the_joint_real_shadow(self):
        for sign in (-1,1):
            np.testing.assert_array_equal(real_shadow(y_state(sign)),IDENTITY/2)
            joint=real_shadow(np.kron(y_state(1),y_state(sign)))
            np.testing.assert_array_equal(joint,pair_state(sign))
            self.assertEqual(np.trace(joint@pauli_word("YY")).real,sign)
        np.testing.assert_array_equal(pair_state(1)@pair_state(-1),np.zeros((4,4)))
        self.assertGreater(np.linalg.norm(np.kron(y_state(1),y_state(1)).imag),0)

    def test_real_pair_can_be_a_marginal_of_a_nonreal_whole_in_the_fixed_basis(self):
        for sign in (-1,1):
            whole=flagged_complex_whole(sign)
            self.assertGreaterEqual(np.linalg.eigvalsh(whole).min(),-2e-16)
            self.assertEqual(np.trace(whole).real,1)
            self.assertGreater(np.linalg.norm(whole.imag),0)
            np.testing.assert_array_equal(keep_systems(whole,(0,1),3),pair_state(sign))
            for site in range(3):
                np.testing.assert_array_equal(keep_systems(whole,(site,),3),IDENTITY/2)
            self.assertEqual(np.trace(whole@pauli_word("YIZ")).real,1)
            self.assertEqual(np.trace(whole@pauli_word("IYZ")).real,sign)

    def test_access_to_retained_flag_recovers_complex_conditional_states(self):
        for sign,h in product((-1,1),repeat=2):
            whole=flagged_complex_whole(sign)
            saved=whole.copy()
            flag=(IDENTITY+h*PAULI_Z)/2
            effect=np.kron(np.eye(4),flag)
            branch=effect@whole@effect
            self.assertEqual(np.trace(branch).real,.5)
            conditional=keep_systems(branch,(0,1),3)*2
            np.testing.assert_array_equal(conditional,np.kron(y_state(h),y_state(sign*h)))
            self.assertGreater(np.linalg.norm(conditional.imag),0)
            np.testing.assert_array_equal(whole,saved)

    def test_joint_real_interface_recovers_a_relation_hidden_from_local_real_records(self):
        difference=pair_state(1)-pair_state(-1)
        for a,b in product("IXZ",repeat=2):
            self.assertEqual(np.trace(difference@pauli_word(a+b)),0)
        effects=[(np.eye(4)+s*pauli_word("YY"))/2 for s in (-1,1)]
        for sign,effect in zip((-1,1),effects):
            self.assertEqual(np.trace(pair_state(sign)@effect).real,1)
        self.assertAlmostEqual(finite_pair_success(),(1+ALPHA)/2)
        self.assertGreater(finite_pair_success_lower(),Fraction(99,100))
        basis=bell_basis()
        np.testing.assert_allclose(basis@pauli_word("IZ")@basis.conj().T,
                                   -pauli_word("YY"),atol=3e-16)

    def test_real_state_does_not_imply_a_real_theory_without_control_restrictions(self):
        real=(IDENTITY+PAULI_X)/2
        phase=np.diag([1,1j])
        np.testing.assert_array_equal(phase@real@phase.conj().T,y_state(1))
        instrument=real_instrument(np.random.default_rng(387),2)
        for k in instrument:
            output=k@real@k.T
            np.testing.assert_array_equal(output.imag,np.zeros((2,2)))
            self.assertGreaterEqual(np.linalg.eigvalsh(output).min(),-2e-16)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(ComplexWholeRealInterfacesTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={
        "round":87,
        "declared_bottom_theory":"Finite dimensional complex quantum mechanics",
        "real_interface":"Fixed basis, real Kraus operations and real symmetric effects",
        "real_shadow_preserves_all_finite_restricted_records":True,
        "real_shadow_is_applied_as_a_physical_local_channel":False,
        "local_real_shadows_determine_the_joint_shadow":False,
        "tensor_shadow_identity":"Re(rho tensor sigma) = Re(rho) tensor Re(sigma) - Im(rho) tensor Im(sigma)",
        "complex_combination_automatically_becomes_real":False,
        "real_pair_inside_complex_whole_constructed":True,
        "whole_and_flag_retained":True,
        "reference_flag_branch_probabilities_exact":["1/2","1/2"],
        "pair_alternatives_real_LOCC_success_exact":"1/2",
        "ideal_joint_real_pair_success_exact":"1",
        "old_readouts":1,
        "old_joint_inverse_bell_gates":1,
        "old_one_readout_joint_success_formula":"(1+alpha)/2",
        "old_one_readout_joint_success_diagnostic":finite_pair_success(),
        "old_one_readout_joint_success_lower_exact":str(finite_pair_success_lower()),
        "unknown_inputs_universally_conjugated":False,
        "basis_and_permission_rules_derived_from_cognition":False,
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("complex_whole_real_interfaces_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
