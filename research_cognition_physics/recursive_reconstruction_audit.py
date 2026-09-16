"""Round 83: same-class composition versus faithful local reconstruction.

A known independent shared reference restores a specified target family.
That does not make every state of the expanded real whole locally tomographic.
"""

import argparse
import json
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from bilocal_record_tomography import real_words, blocks_for_word
from independent_source_alignment import keep_systems, pair_reference
from partition_interface_closure import dimensions
from quantum_interface_audit import ALPHA
from role_symmetry_and_swap import pauli_word


LOCAL_LABELS = real_words(2)
LOCAL_PRODUCTS = tuple(a+b for a,b in product(LOCAL_LABELS,repeat=2))
GLOBAL_LABELS = real_words(4)


def trace_pairing(effects, inputs):
    first = np.array([e.reshape(-1) for e in effects])
    second = np.array([b.T.reshape(-1) for b in inputs])
    return (first @ second.T).real


def full_pairing():
    return trace_pairing([pauli_word(w) for w in LOCAL_PRODUCTS],
                         [pauli_word(w)/16 for w in GLOBAL_LABELS])


def embed_known_reference(target):
    return keep_systems(np.kron(target,pair_reference(1)),(0,2,1,3),4)


def promised_pairing(noisy=False):
    inputs = [embed_known_reference(pauli_word(w)/4) for w in real_words(2)]
    effects = [pauli_word(w) for w in LOCAL_PRODUCTS]
    pairing = trace_pairing(effects,inputs)
    if noisy:
        factors = [ALPHA**(len(blocks_for_word(w[:2]))+len(blocks_for_word(w[2:])))
                   for w in LOCAL_PRODUCTS]
        pairing = np.asarray(factors)[:,None]*pairing
    return pairing


def padded_hidden_pair(local_dimension):
    d = local_dimension
    if d < 2 or d % 2:
        raise ValueError("This explicit orthogonal family uses even local dimension.")
    j = np.kron(np.eye(d//2),np.array([[0.,1.],[-1.,0.]]))
    joint = np.kron(j,j)
    return (np.eye(d*d)+joint)/(d*d), (np.eye(d*d)-joint)/(d*d)


def transpose_first(matrix,local_dimension):
    d = local_dimension
    return matrix.reshape(d,d,d,d).transpose(2,1,0,3).reshape(d*d,d*d)


class RecursiveReconstructionAuditTests(unittest.TestCase):
    def test_unrestricted_expanded_whole_has_an_exact_36_dimensional_local_kernel(self):
        pairing = full_pairing()
        self.assertEqual(pairing.shape,(100,136))
        visible = [j for j in range(136) if np.any(pairing[:,j])]
        np.testing.assert_array_equal(pairing[:,visible],np.eye(100))
        missing = [w for j,w in enumerate(GLOBAL_LABELS) if j not in visible]
        self.assertEqual(len(missing),36)
        self.assertTrue(all(w[:2].count("Y") % 2 == w[2:].count("Y") % 2 == 1 for w in missing))

    def test_known_independent_shared_reference_makes_the_target_family_injective(self):
        pairing = promised_pairing()
        self.assertEqual(pairing.shape,(100,10))
        self.assertEqual(np.linalg.matrix_rank(pairing),10)
        # Integer Gram entries give a nonzero determinant certificate.
        gram = pairing.T@pairing
        np.testing.assert_array_equal(gram,np.eye(10))

    def test_existing_noisy_block_readings_preserve_target_reconstruction(self):
        rng = np.random.default_rng(83)
        ideal, noisy = promised_pairing(), promised_pairing(True)
        raw = rng.normal(size=(4,4))
        target = raw@raw.T
        target /= np.trace(target)
        coordinates = np.array([np.trace(target@pauli_word(w)).real for w in real_words(2)])
        record_moments = noisy@coordinates
        reconstructed = np.linalg.lstsq(noisy,record_moments,rcond=None)[0]
        np.testing.assert_allclose(reconstructed,coordinates,atol=8e-16)
        self.assertEqual(np.linalg.matrix_rank(ideal),np.linalg.matrix_rank(noisy))

    def test_any_even_finite_local_padding_still_has_orthogonal_real_hidden_states(self):
        for d in (2,4,6,8):
            first, second = padded_hidden_pair(d)
            for rho in (first,second):
                self.assertAlmostEqual(np.trace(rho),1)
                self.assertGreaterEqual(np.linalg.eigvalsh(rho).min(),-2e-17)
                np.testing.assert_array_equal(rho,rho.T)
                np.testing.assert_allclose(np.einsum("abcb->ac",rho.reshape(d,d,d,d)),np.eye(d)/d,atol=0)
            np.testing.assert_allclose(first@second,np.zeros_like(first),atol=2e-19)
            self.assertAlmostEqual(np.abs(np.linalg.eigvalsh(first-second)).sum()/2,1)
            k,l = dimensions(d)
            self.assertEqual(dimensions(d*d)[0]-k*k,l*l)

    def test_all_real_product_effects_have_zero_pairing_with_padded_difference(self):
        rng = np.random.default_rng(183)
        for d in (2,4,6):
            first,second = padded_hidden_pair(d)
            for _ in range(12):
                ka,kb = rng.normal(size=(d,d))/d,rng.normal(size=(d,d))/d
                effect = np.kron(ka.T@ka,kb.T@kb)
                self.assertAlmostEqual(np.trace((first-second)@effect),0,places=14)

    def test_partial_transpose_explains_the_whole_local_record_obstruction(self):
        for d in (2,4,6):
            first,second = padded_hidden_pair(d)
            np.testing.assert_array_equal(transpose_first(first,d),second)
            # Any real product effect is invariant under this transpose.
            vector = np.arange(1,d+1)/d
            effect = np.kron(np.outer(vector,vector),np.eye(d))
            np.testing.assert_array_equal(transpose_first(effect,d),effect)

    def test_hidden_expanded_pair_does_not_satisfy_the_same_fixed_reference_promise(self):
        first,second = padded_hidden_pair(4)
        reference_marginals = [keep_systems(rho,(1,3),4) for rho in (first,second)]
        np.testing.assert_array_equal(reference_marginals[0],pair_reference(-1))
        np.testing.assert_array_equal(reference_marginals[1],pair_reference(1))
        for rho in (first,second):
            np.testing.assert_array_equal(keep_systems(rho,(0,2),4),np.eye(4)/4)
        self.assertFalse(np.array_equal(first,embed_known_reference(np.eye(4)/4)))
        np.testing.assert_array_equal(second,embed_known_reference(np.eye(4)/4))

    def test_classical_local_tomography_needs_joint_records_and_is_not_quantum_specific(self):
        first,second = np.diag([.5,0,0,.5]),np.diag([0,.5,.5,0])
        for rho in (first,second):
            np.testing.assert_array_equal(keep_systems(rho,(0,),2),np.eye(2)/2)
            np.testing.assert_array_equal(keep_systems(rho,(1,),2),np.eye(2)/2)
        effects = [np.diag(np.eye(4)[i]) for i in range(4)]
        self.assertEqual(np.linalg.matrix_rank(trace_pairing(effects,effects)),4)
        self.assertEqual(sum(abs(np.trace((first-second)@e)) for e in effects)/2,1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(RecursiveReconstructionAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round":83,
        "user_question":"Could complex structure follow from reproducing the original cognitive structure after composition?",
        "same_class_recursive_closure_excludes_real_theory":False,
        "known_shared_reference_target_family":{"homogeneous_target_dimension":10,"local_record_rank":10,"gram_determinant_exact":1},
        "unrestricted_expanded_real_whole":{"homogeneous_dimension":136,"local_product_rank":100,"kernel_dimension":36},
        "arbitrary_finite_even_local_padding":{"dimension_deficit":"[d(d-1)/2]^2","orthogonal_hidden_pair_exists":True},
        "known_fixed_reference_recovery_contradicted_by_hidden_expanded_pair":False,
        "local_tomography_alone_forces_complex_quantum_theory":False,
        "recovering_state_statistics_equals_recovering_all_joint_actions":False,
        "scope":"Full local real instruments and classical records; no extra certified shared resource outside the whole under examination",
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}
    }
    if args.write_results:
        Path(__file__).with_name("recursive_reconstruction_audit_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    main()
