"""Round 89: reversible exposure of a relation and the cost of a fixed frame."""

import argparse
import json
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from bell_network_statistics import bell_basis
from complex_whole_real_interfaces import flagged_complex_whole, pair_state, random_state
from independent_source_alignment import keep_systems
from partition_interface_closure import sector_bases
from quantum_interface_audit import ALPHA, IDENTITY, PAULI_Z, effective_readout
from role_symmetry_and_swap import pauli_word


def decode(state):
    u=bell_basis()
    return u.conj().T@state@u


def local_basis(da=2,db=2):
    return tuple(np.kron(a,b) for a,b in product(sector_bases(da)[0],sector_bases(db)[0]))


def span_rank(matrices):
    # Hermitian matrices form a real vector space.
    vectors=np.array([np.concatenate((m.real.reshape(-1),m.imag.reshape(-1))) for m in matrices])
    return int(np.linalg.matrix_rank(vectors))


def input_frame(unitary,basis):
    return tuple(unitary@e@unitary.conj().T for e in basis)


def frame_report():
    original=local_basis()
    changed=input_frame(bell_basis(),original)
    merged=original+changed
    return {"original_rank":span_rank(original),"decoded_rank":span_rank(changed),
            "union_rank":span_rank(merged),
            "intersection_dimension":span_rank(original)+span_rank(changed)-span_rank(merged),
            "unrecorded_half_mixture_rank":span_rank(tuple((a+b)/2 for a,b in zip(original,changed)))}


class RedifferentiationTradeoffTests(unittest.TestCase):
    def test_old_joint_gate_moves_yy_relation_into_one_local_z_bit(self):
        for s in (-1,1):
            expected=np.kron(IDENTITY/2,(IDENTITY-s*PAULI_Z)/2)
            np.testing.assert_allclose(decode(pair_state(s)),expected,atol=2e-16)

    def test_complex_whole_and_reference_survive_the_change_of_access(self):
        extended=np.kron(bell_basis(),IDENTITY)
        for s in (-1,1):
            whole=flagged_complex_whole(s)
            changed=extended.conj().T@whole@extended
            np.testing.assert_allclose(keep_systems(changed,(0,1),3),decode(pair_state(s)),atol=1e-16)
            np.testing.assert_allclose(keep_systems(changed,(2,),3),IDENTITY/2,atol=1e-16)
            self.assertGreater(np.linalg.norm(changed.imag),0)
            np.testing.assert_allclose(extended@changed@extended.conj().T,whole,atol=2e-16)

    def test_after_splitting_again_the_old_local_readout_still_reads_the_new_relation(self):
        for s,outcome in product((-1,1),(0,1)):
            sign=2*outcome-1
            effect=np.kron(IDENTITY,effective_readout(0,outcome,1.))
            actual=np.trace(decode(pair_state(s))@effect).real
            self.assertAlmostEqual(actual,(1-sign*s*ALPHA)/2)

    def test_gate_is_reversible_for_arbitrary_complex_input_not_only_the_promise_family(self):
        rng=np.random.default_rng(89)
        for _ in range(5):
            rho=random_state(rng,4)
            u=bell_basis()
            np.testing.assert_allclose(u@decode(rho)@u.conj().T,rho,atol=2e-16)

    def test_fixed_frame_exchanges_one_of_nine_visible_directions(self):
        self.assertEqual(frame_report(),{
            "original_rank":9,"decoded_rank":9,"union_rank":10,
            "intersection_dimension":8,"unrecorded_half_mixture_rank":9})
        new=input_frame(bell_basis(),local_basis())
        old=local_basis()
        self.assertTrue(all(abs(np.trace(e@pauli_word("YY")))<1e-14 for e in old))
        self.assertTrue(all(abs(np.trace(e@pauli_word("IZ")))<1e-14 for e in new))

    def test_an_old_local_bit_becomes_invisible_under_the_same_fixed_gate(self):
        before=[(np.eye(4)+s*pauli_word("IZ"))/4 for s in (-1,1)]
        after=[decode(rho) for rho in before]
        for s,rho in zip((-1,1),before):
            self.assertEqual(np.trace(rho@pauli_word("IZ")).real,s)
        for effect in local_basis():
            self.assertAlmostEqual(np.trace((after[1]-after[0])@effect),0)
        for s,rho in zip((-1,1),after):
            np.testing.assert_allclose(rho,pair_state(s),atol=2e-16)

    def test_any_fixed_unitary_preserves_the_dimension_of_the_readable_effect_space(self):
        rng=np.random.default_rng(189)
        for da,db in ((2,2),(2,3),(3,3)):
            d=da*db
            basis=local_basis(da,db)
            raw=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d))
            u,_=np.linalg.qr(raw)
            self.assertEqual(span_rank(input_frame(u,basis)),len(basis))
            self.assertEqual(len(basis),da*(da+1)*db*(db+1)//4)

    def test_retaining_the_choice_of_frame_preserves_the_union_of_predictions(self):
        basis=tuple(pauli_word(a+b) for a,b in product("IXZ",repeat=2))
        original=basis
        changed=input_frame(bell_basis(),basis)
        real_directions=tuple(m for sector in sector_bases(4)[:1] for m in sector)
        rows=np.array([[np.trace(e@state).real for state in real_directions]
                       for e in original+changed])
        self.assertEqual(np.linalg.matrix_rank(rows),10)
        self.assertEqual(np.linalg.matrix_rank((rows[:9]+rows[9:])/2),9)
        # These are ensemble prediction coordinates, not simultaneous exact
        # readings of all noncommuting quantities on a single specimen.


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(RedifferentiationTradeoffTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={
        "round":89,
        "decoded_tau_s":"I_A/2 tensor (I_B-s Z_B)/2",
        "whole_and_existing_reference_are_retained":True,
        "old_joint_gate_is_reversible_for_arbitrary_inputs":True,
        "post_split_relation_success_formula":"(1+alpha)/2",
        "post_split_relation_success_diagnostic":(1+ALPHA)/2,
        "fixed_frame":frame_report(),
        "lost_old_direction":"I tensor Z",
        "gained_direction":"Y tensor Y",
        "one_fixed_unitary_strictly_enlarges_all_readable_effect_directions":False,
        "retained_frame_choice_covers_the_real_ten_direction_space":True,
        "ten_direction_reconstruction_is_single_specimen_exact_tomography":False,
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("redifferentiation_tradeoff_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
