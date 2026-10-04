"""Round 90: write a relation to real memory, restore the old promised whole.

The full write is reversible. Ignoring its memory dephases arbitrary inputs.
Only the declared commuting family is preserved without disturbance.
"""

import argparse
import json
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from bell_network_statistics import bell_basis
from bilocal_record_tomography import embed_operator
from certified_intervals import SCALE
from complex_whole_real_interfaces import flagged_complex_whole, random_state
from independent_source_alignment import keep_systems
from one_bit_real_network import visibility_interval
from operational_effect_closure import compiled_copy_gate
from quantum_interface_audit import ALPHA, IDENTITY, PAULI_Z, branch_kraus
from role_symmetry_and_swap import pauli_word


def sector_projector(sign):
    return (np.eye(4)+sign*pauli_word("YY"))/2


def memory_write_unitary():
    # Wire order M,A,B. M is new real memory, initialized in |0>.
    decoder=np.kron(IDENTITY,bell_basis())
    copier=embed_operator(compiled_copy_gate(),(0,2),3)
    return decoder@copier@decoder.conj().T


def memory_write_isometry():
    return memory_write_unitary()@np.kron(np.array([[1.],[0.]]),np.eye(4))


def write_whole(whole):
    # Input A,B,C -> output M,A,B,C. C is the retained original spectator.
    isometry=np.kron(memory_write_isometry(),IDENTITY)
    return isometry@whole@isometry.conj().T


def local_memory_read_branch(state,outcome):
    # Preserve M by copying its Z value into a reset N, then read only N.
    # This is an actual old noisy instrument, not an ideal measurement.
    n=int(round(np.log2(len(state))))
    p0=(IDENTITY+PAULI_Z)/2
    extended=np.kron(p0,state)
    copy=embed_operator(compiled_copy_gate(),(0,1),n+1)
    extended=copy@extended@copy.conj().T
    output=np.zeros_like(extended)
    for k in branch_kraus(0,outcome,1.):
        full=np.kron(k,np.eye(2**n))
        output+=full@extended@full.conj().T
    return keep_systems(output,tuple(range(1,n+1)),n+1)


def record_tree(state,count):
    leaves={():state}
    for _ in range(count):
        leaves={history+(r,):local_memory_read_branch(branch,r)
                for history,branch in leaves.items() for r in (0,1)}
    return leaves


def product_ray_constraint_dimension(da,db):
    def rays(d):
        eye=np.eye(d)
        return tuple(eye)+tuple((eye[i]+eye[j])/np.sqrt(2)
                               for i in range(d) for j in range(i+1,d))
    d=da*db
    rows=[]
    for a,b in product(rays(da),rays(db)):
        v=np.kron(a,b)
        rows.append(np.kron(v.reshape(1,-1),np.eye(d)-np.outer(v,v)))
    return d*d-int(np.linalg.matrix_rank(np.vstack(rows)))


def memory_read_certificate(count):
    success=(1+visibility_interval(count))/2
    return {"old_readouts":count,
            "new_memory_rebits":1,
            "reusable_read_pointer_rebits":1,
            "old_pair_interactions_to_write":3,
            "old_local_rotations_to_write":3,
            "old_pair_interactions_per_pointer_read":1,
            "old_local_rotations_per_pointer_read":3,
            "success_lower_exact":str(Fraction(success.lo,SCALE)),
            "success_diagnostic":sum(success.floats())/2}


class PersistentRelationMemoryTests(unittest.TestCase):
    def test_isometry_is_compiled_from_old_gates_with_correct_sector_labels(self):
        actual=memory_write_isometry()
        expected=(np.kron(np.array([[1.],[0.]]),sector_projector(-1))
                  +np.kron(np.array([[0.],[1.]]),sector_projector(1)))
        np.testing.assert_allclose(actual,expected,atol=4e-16)
        np.testing.assert_allclose(actual.conj().T@actual,np.eye(4),atol=5e-16)
        np.testing.assert_allclose(actual.imag,0,atol=0)

    def test_every_state_inside_one_sector_keeps_the_entire_old_correlated_whole(self):
        rng=np.random.default_rng(90)
        isometry=np.kron(memory_write_isometry(),IDENTITY)
        for s in (-1,1):
            projector=np.kron(sector_projector(s),IDENTITY)
            projected=projector@random_state(rng,8)@projector
            projected/=np.trace(projected).real
            for whole in (projected,flagged_complex_whole(s)):
                pointer=(IDENTITY-s*PAULI_Z)/2
                expected=np.kron(pointer,whole)
                np.testing.assert_allclose(isometry@whole@isometry.conj().T,expected,atol=3e-16)

    def test_actual_repeated_local_reads_keep_all_branches_and_the_original_whole(self):
        for s in (-1,1):
            whole=flagged_complex_whole(s)
            state=write_whole(whole)
            leaves=record_tree(state,3)
            self.assertEqual(len(leaves),8)
            mass=0
            success=0
            for history,branch in leaves.items():
                expected=np.prod([(1-(2*r-1)*s*ALPHA)/2 for r in history])
                probability=np.trace(branch).real
                self.assertAlmostEqual(probability,expected,places=14)
                np.testing.assert_allclose(keep_systems(branch,(1,2,3),4),
                                           expected*whole,atol=8e-16)
                mass+=probability
                predicted=1 if sum(history)<2 else -1
                if predicted==s: success+=probability
            self.assertAlmostEqual(mass,1,places=14)
            self.assertAlmostEqual(success,memory_read_certificate(3)["success_diagnostic"],places=14)
        self.assertGreater(Fraction(memory_read_certificate(3)["success_lower_exact"]),Fraction(999,1000))

    def test_forgetting_new_memory_dephases_a_general_complex_whole(self):
        rng=np.random.default_rng(190)
        yy=np.kron(pauli_word("YY"),IDENTITY)
        for _ in range(4):
            whole=random_state(rng,8)
            actual=keep_systems(write_whole(whole),(1,2,3),4)
            expected=(whole+yy@whole@yy)/2
            np.testing.assert_allclose(actual,expected,atol=2e-16)

    def test_exact_fixed_states_are_those_commuting_with_the_relation(self):
        rng=np.random.default_rng(290)
        yy=pauli_word("YY")
        for _ in range(5):
            rho=random_state(rng,4)
            fixed=(rho+yy@rho@yy)/2
            np.testing.assert_allclose(fixed@yy,yy@fixed,atol=0)
            isometry=memory_write_isometry()
            actual=keep_systems(isometry@fixed@isometry.conj().T,(1,2),3)
            np.testing.assert_allclose(actual,fixed,atol=2e-16)
            self.assertAlmostEqual(np.linalg.norm((rho+yy@rho@yy)/2-rho),
                                   np.linalg.norm(rho@yy-yy@rho)/2)

    def test_an_arbitrary_old_local_bit_can_be_disturbed_by_recording_the_relation(self):
        rho=np.diag([1.,0,0,0])
        isometry=memory_write_isometry()
        after=keep_systems(isometry@rho@isometry.conj().T,(1,2),3)
        np.testing.assert_allclose(after,np.diag([.5,0,0,.5]),atol=3e-16)
        self.assertEqual(np.trace(rho@pauli_word("IZ")).real,1)
        self.assertAlmostEqual(np.trace(after@pauli_word("IZ")).real,0)
        self.assertAlmostEqual(np.abs(np.linalg.eigvalsh(after-rho)).sum()/2,.5)

    def test_preserving_all_real_product_pure_states_forces_scalar_kraus_operators(self):
        for da,db in ((2,2),(2,3),(3,3)):
            self.assertEqual(product_ray_constraint_dimension(da,db),1)
        # The analytic proof in the note extends this finite spanning-ray
        # constraint to all finite dimensions and to complex Kraus matrices.

    def test_before_any_readout_the_enlarged_whole_can_be_exactly_uncomputed(self):
        rng=np.random.default_rng(390)
        gate=np.kron(memory_write_unitary(),IDENTITY)
        for _ in range(4):
            old=random_state(rng,8)
            initial=np.kron((IDENTITY+PAULI_Z)/2,old)
            encoded=gate@initial@gate.conj().T
            np.testing.assert_allclose(encoded,write_whole(old),atol=2e-16)
            np.testing.assert_allclose(gate.conj().T@encoded@gate,initial,atol=4e-16)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(PersistentRelationMemoryTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={
        "round":90,
        "write_isometry":"|0>_M tensor Pi_minus + |1>_M tensor Pi_plus",
        "whole_within_one_yy_sector_is_preserved_and_memory_decouples":True,
        "ignoring_memory_channel":"(rho + YY rho YY)/2, with identity on spectators",
        "fixed_point_condition":"[rho, YY tensor I_spectator] = 0",
        "noncommuting_example_old_whole_trace_distance_exact":"1/2",
        "full_write_before_readout_is_reversible":True,
        "all_finite_readout_branches_retained":True,
        "post_split_memory_read_certificates":[memory_read_certificate(m) for m in (1,3,5)],
        "universal_preservation_of_old_real_local_statistics_allows_nontrivial_record":False,
        "no_information_proof_scope":"Exact universal preservation; finite dimensional complex CPTP instrument",
        "prepared_memory_and_control_resources_derived_from_cognition":False,
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("persistent_relation_memory_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
