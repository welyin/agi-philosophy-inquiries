"""Round 184: fixed auxiliary type, arbitrary ensembles, and dimension scope."""

import argparse
import json
from pathlib import Path
import unittest

import numpy as np
from composition_dimension_audit import matrix_basis
from cone_duality_audit import CUBE_RAYS, DUAL_RAYS
from ensemble_extension_audit import remote_branch


def fixed_copy_data(rho, tolerance=1e-12):
    values, vectors = np.linalg.eigh(rho)
    keep = values > tolerance
    root = (vectors[:, keep] * np.sqrt(values[keep])) @ vectors[:, keep].conj().T
    inverse = (vectors[:, keep] / np.sqrt(values[keep])) @ vectors[:, keep].conj().T
    support = vectors[:, keep] @ vectors[:, keep].conj().T
    psi = root.ravel()
    return np.outer(psi, psi.conj()), root, inverse, support


def fixed_copy_effects(rho, branches):
    _, _, inverse, support = fixed_copy_data(rho)
    effects = [(inverse @ branch @ inverse).T for branch in branches]
    effects[0] = effects[0] + np.eye(len(rho)) - support.T
    return effects


def conditioning_rank(rho, complex_model):
    _, root, _, _ = fixed_copy_data(rho)
    columns = np.array([(root @ effect.T @ root).ravel()
                        for effect in matrix_basis(len(rho), complex_model)]).T
    return np.linalg.matrix_rank(np.vstack((columns.real, columns.imag)), tol=1e-10)


class UniformSteeringAuditTests(unittest.TestCase):
    def test_same_purification_same_auxiliary_for_multiple_real_and_complex_ensembles(self):
        rng = np.random.default_rng(184)
        for complex_model in (False, True):
            rho = np.array([[.6, .1], [.1, .4]], dtype=complex if complex_model else float)
            if complex_model:
                rho[0, 1], rho[1, 0] = .1j, -.1j
            joint, root, _, _ = fixed_copy_data(rho)
            for outcome_count in (2, 3, 7):
                raw = rng.normal(size=(outcome_count, 2))
                if complex_model:
                    raw = raw + 1j * rng.normal(size=raw.shape)
                frame, _ = np.linalg.qr(raw)
                source_effects = [np.outer(row.conj(), row) for row in frame]
                branches = [root @ e.T @ root for e in source_effects]
                effects = fixed_copy_effects(rho, branches)
                np.testing.assert_allclose(sum(effects), np.eye(2), atol=2e-15)
                for e, branch in zip(effects, branches):
                    self.assertGreaterEqual(np.linalg.eigvalsh(e).min(), -5e-16)
                    np.testing.assert_allclose(remote_branch(joint, e, 2), branch, atol=2e-15)

    def test_rank_deficient_states_keep_full_copy_by_padding_unused_effect_space(self):
        rho = np.diag([.7, .3, 0.])
        branches = [np.diag([.2, .1, 0.]), np.diag([.5, .2, 0.]), np.zeros((3, 3))]
        joint, _, _, _ = fixed_copy_data(rho)
        effects = fixed_copy_effects(rho, branches)
        self.assertEqual(joint.shape, (9, 9))
        np.testing.assert_allclose(sum(effects), np.eye(3), atol=6e-16)
        for e, branch in zip(effects, branches):
            self.assertGreaterEqual(np.linalg.eigvalsh(e).min(), 0)
            np.testing.assert_allclose(remote_branch(joint, e, 3), branch, atol=4e-16)

    def test_conditioning_is_invertible_only_on_full_support(self):
        for complex_model in (False, True):
            self.assertEqual(conditioning_rank(np.eye(3) / 3, complex_model), 9 if complex_model else 6)
            self.assertEqual(conditioning_rank(np.diag([.5, .5, 0]), complex_model), 4 if complex_model else 3)
            self.assertEqual(conditioning_rank(np.diag([1., 0., 0.]), complex_model), 1)

    def test_two_conditioning_isomorphisms_give_the_declared_cone_automorphism(self):
        rho = np.array([[.6, .1j], [-.1j, .4]])
        sigma = np.array([[.4, .15], [.15, .6]])
        _, root_rho, inverse_rho, _ = fixed_copy_data(rho)
        _, root_sigma, _, _ = fixed_copy_data(sigma)
        operator = root_sigma @ inverse_rho
        for target in (rho, np.array([[.5, .2j], [-.2j, .5]])):
            preimage = (inverse_rho @ target @ inverse_rho).T
            output = root_sigma @ preimage.T @ root_sigma
            np.testing.assert_allclose(output, operator @ target @ operator.conj().T, atol=6e-16)
            self.assertGreater(np.linalg.eigvalsh(preimage).min(), 0)
        np.testing.assert_allclose(operator @ rho @ operator.conj().T, sigma, atol=6e-16)
        np.testing.assert_allclose(root_rho @ inverse_rho, np.eye(2), atol=6e-16)

    def test_classical_fixed_copy_covers_zero_support_without_changing_type(self):
        p = np.array([.25, .75, 0.])
        branches = np.array([[.1, .2, 0.], [.15, .55, 0.]])
        effects = np.zeros_like(branches)
        effects[:, :2] = branches[:, :2] / p[:2]
        effects[0, 2] = 1
        np.testing.assert_allclose(effects.sum(axis=0), np.ones(3), atol=2e-16)
        np.testing.assert_allclose(effects * p, branches, atol=2e-16)
        self.assertEqual(np.linalg.matrix_rank(np.diag(p)), 2)

    def test_cube_self_copy_has_a_cone_isomorphism_obstruction(self):
        self.assertEqual(np.linalg.matrix_rank(CUBE_RAYS), 4)
        self.assertEqual(np.linalg.matrix_rank(DUAL_RAYS), 4)
        self.assertEqual((len(CUBE_RAYS), len(DUAL_RAYS)), (8, 6))
        # Equal linear dimensions do not repair unequal extreme-ray counts.

    def test_cube_center_can_be_steered_by_a_different_same_dimension_auxiliary(self):
        # A has octahedral state cone, so its effect cone is the cube cone.
        # Identity conditioning is a positive bilinear joint state in the
        # maximal composite. This concerns the center, not all B states.
        center = np.array([1., 0., 0., 0.])
        for ray in CUBE_RAYS:
            effects = [ray / 2, (2 * center - ray) / 2]
            np.testing.assert_array_equal(sum(effects), center)
            for effect in effects:
                self.assertTrue(np.all(DUAL_RAYS @ effect >= 0))
        np.testing.assert_array_equal(np.eye(4) @ center, center)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(UniformSteeringAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = {
        "round": 184,
        "quantifier_order": "for every B, one fixed auxiliary type A; for every rho, one W_rho; for every ensemble, a measurement",
        "dimension_means": "real linear operational dimension K, not Hilbert dimension or outcome count",
        "self_contained_order_isomorphism_proof_assumes_unrestricted_auxiliary_effect_interval": True,
        "binary_interval_surjectivity_with_equal_K_and_interior_marginal_implies_full_ensemble_lifting": True,
        "classical_real_and_complex_fixed_copy_constructions": True,
        "cube_self_steering_for_interior_states_possible": False,
        "cube_center_steering_by_octahedral_auxiliary_in_maximal_composite": True,
        "one_center_example_implies_uniform_steering_for_all_states": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("uniform_steering_audit_results.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
