"""Round 85: a declared common real complex structure and its exact encoding.

The commuting-state cone is equivalent to complex density matrices. Merely
preserving that cone still allows conjugation; orientation-preserving Kraus
operators are a stronger requirement. None follows from task gain alone.
"""

import argparse
import json
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from ancillary_operation_equivalence import J
from complex_control_from_reference import real_lift
from partition_interface_closure import sector_bases


def orientation(d):
    return np.kron(J,np.eye(d))


def encode_state(density):
    return real_lift(density)/2


def decode_state(encoded):
    d = len(encoded)//2
    return 2*(encoded[:d,:d]+1j*encoded[:d,d:])


def decode_operator(lifted):
    d = len(lifted)//2
    return lifted[:d,:d]+1j*lifted[:d,d:]


def orientation_average(matrix):
    d = len(matrix)//2
    j = orientation(d)
    return (matrix+j@matrix@j.T)/2


def normalizer_flip(d):
    return np.diag([1.]*d+[-1.]*d)


def commuting_dimension(d,skew=False):
    basis = sector_bases(2*d)[int(skew)]
    j = orientation(d)
    constraints = np.array([(a@j-j@a).reshape(-1) for a in basis]).T
    return len(basis)-np.linalg.matrix_rank(constraints)


class CommonOrientationStructureTests(unittest.TestCase):
    def test_commuting_symmetric_states_have_exactly_d_squared_coordinates(self):
        for d in (1,2,3,4):
            self.assertEqual(commuting_dimension(d),d*d)
            np.testing.assert_array_equal(orientation(d)@orientation(d),-np.eye(2*d))

    def test_encoding_is_real_positive_and_preserves_normalization_and_decoding(self):
        rng = np.random.default_rng(85)
        for d in (2,3,4):
            for _ in range(5):
                raw = rng.normal(size=(d,d))+1j*rng.normal(size=(d,d))
                rho = raw@raw.conj().T
                rho = (rho+rho.conj().T)/2
                rho /= np.trace(rho).real
                encoded = encode_state(rho)
                self.assertTrue(np.isrealobj(encoded))
                np.testing.assert_allclose(encoded,encoded.T,atol=0)
                self.assertGreater(np.linalg.eigvalsh(encoded).min(),0)
                self.assertAlmostEqual(np.trace(encoded),1)
                np.testing.assert_array_equal(decode_state(encoded),rho)
                np.testing.assert_allclose(encoded@orientation(d),orientation(d)@encoded,atol=0)

    def test_every_averaged_real_state_decodes_to_a_valid_complex_state(self):
        rng = np.random.default_rng(185)
        for d in (2,3,4):
            raw = rng.normal(size=(2*d,2*d))
            rho = raw@raw.T
            rho /= np.trace(rho)
            averaged = orientation_average(rho)
            decoded = decode_state(averaged)
            np.testing.assert_allclose(decoded,decoded.conj().T,atol=0)
            self.assertGreaterEqual(np.linalg.eigvalsh(decoded).min(),0)
            np.testing.assert_array_equal(encode_state(decoded),averaged)
            np.testing.assert_allclose(orientation_average(averaged),averaged,atol=0)

    def test_encoded_effects_preserve_probabilities_and_complete_measurements(self):
        rng = np.random.default_rng(285)
        d = 3
        raw = rng.normal(size=(d,d))+1j*rng.normal(size=(d,d))
        rho = raw@raw.conj().T
        rho /= np.trace(rho)
        unitary,_ = np.linalg.qr(rng.normal(size=(d,d))+1j*rng.normal(size=(d,d)))
        effects = [np.outer(unitary[:,i],unitary[:,i].conj()) for i in range(d)]
        np.testing.assert_allclose(sum(real_lift(e) for e in effects),np.eye(2*d),atol=7e-16)
        for effect in effects:
            self.assertAlmostEqual(np.trace(encode_state(rho)@real_lift(effect)),
                                   np.trace(rho@effect).real)

    def test_orientation_preserving_orthogonal_control_is_the_unitary_lift(self):
        rng = np.random.default_rng(385)
        for d in (2,3,4):
            self.assertEqual(commuting_dimension(d,True),d*d)
            unitary,_ = np.linalg.qr(rng.normal(size=(d,d))+1j*rng.normal(size=(d,d)))
            lifted = real_lift(unitary)
            np.testing.assert_allclose(lifted.T@lifted,np.eye(2*d),atol=8e-16)
            np.testing.assert_allclose(lifted@orientation(d),orientation(d)@lifted,atol=0)
            np.testing.assert_array_equal(decode_operator(lifted),unitary)
            self.assertAlmostEqual(np.linalg.det(lifted),1,places=13)

    def test_commuting_real_kraus_maps_reproduce_complex_channels(self):
        d = 2
        rate = .3
        phase = np.diag([1,1j])
        kraus = (np.diag([1,np.sqrt(1-rate)])@phase,np.array([[0,np.sqrt(rate)],[0,0]])@phase)
        rho = np.array([[.6,.2+.1j],[.2-.1j,.4]])
        lifted = tuple(real_lift(k) for k in kraus)
        np.testing.assert_allclose(sum(k.T@k for k in lifted),np.eye(2*d),atol=2e-16)
        output = sum(k@encode_state(rho)@k.T for k in lifted)
        expected = encode_state(sum(k@rho@k.conj().T for k in kraus))
        np.testing.assert_allclose(output,expected,atol=1e-16)

    def test_preserving_the_state_cone_does_not_force_orientation_preservation(self):
        rho = np.array([[.5,.2j],[-.2j,.5]])
        flip = normalizer_flip(2)
        np.testing.assert_array_equal(flip@orientation(2)@flip.T,-orientation(2))
        output = flip@encode_state(rho)@flip.T
        np.testing.assert_array_equal(decode_state(output),rho.conj())
        np.testing.assert_array_equal(output@orientation(2),orientation(2)@output)
        self.assertGreater(np.linalg.norm(output-encode_state(rho)),.1)

    def test_independent_complex_partial_transpose_is_not_a_physical_channel(self):
        vector = np.array([1.,0.,0.,1.])
        bell = np.outer(vector,vector)/2
        partial = bell.reshape(2,2,2,2).transpose(0,3,2,1).reshape(4,4)
        antisymmetric = np.array([0.,1.,-1.,0.])
        witness = np.outer(antisymmetric,antisymmetric)/2
        self.assertEqual(np.trace(partial@witness),-.5)
        # Global conjugation of the entire encoded logical system is distinct
        # from treating it as a freely composable local complex channel.
        np.testing.assert_array_equal(normalizer_flip(4)@encode_state(bell)@normalizer_flip(4).T,
                                      encode_state(bell.conj()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(CommonOrientationStructureTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report = {
        "round":85,
        "declared_common_structure":"J_d = J tensor I_d, J_d^2 = -I",
        "encoded_state":"rho_real = R(rho_complex)/2",
        "state_cone":"Real positive symmetric normalized matrices commuting with J_d",
        "state_cone_homogeneous_dimension":"d^2",
        "orientation_preserving_orthogonal_group":"R(U(d))",
        "commuting_kraus_maps_match_complex_CPTP":True,
        "state_cone_preservation_alone_forces_complex_CPTP":False,
        "orientation_reversal_is_a_real_cone_preserving_operation":True,
        "logical_conjugation_is_a_universal_local_complex_channel":False,
        "partial_transpose_negative_witness_exact":"-1/2",
        "common_J_is_derived_from_boundary_expansion":False,
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("common_orientation_structure_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__": main()
