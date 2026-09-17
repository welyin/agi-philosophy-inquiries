"""Round 179: ensemble choice on an existing extension, with its assumptions.

The matrix proof and classical conditional construction are in the note.
This is not an unknown-input purification channel or a cognition reconstruction.
"""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import unittest

import numpy as np


I = np.eye(2)
X = np.array([[0., 1.], [1., 0.]])
Y = np.array([[0., -1j], [1j, 0.]])
Z = np.diag([1., -1.])


def remote_branch(joint, effect, local_dimension):
    env_dimension = effect.shape[0]
    blocks = joint.reshape(local_dimension, env_dimension, local_dimension, env_dimension)
    return np.einsum("aibj,ji->ab", blocks, effect)


def ensemble_extension(branches, tolerance=1e-12):
    """Minimal spectral purification and POVM for positive subnormalized branches.

    Inputs here are known matrices. Numerical support uses a declared tolerance;
    the note's exact statement uses the exact support, including singular states.
    """
    rho = sum(branches)
    values, vectors = np.linalg.eigh(rho)
    support = values > tolerance
    values, vectors = values[support], vectors[:, support]
    coefficients = vectors * np.sqrt(values)
    inverse = (vectors / np.sqrt(values)).conj().T
    effects = [(inverse @ branch @ inverse.conj().T).T for branch in branches]
    psi = coefficients.ravel()
    return np.outer(psi, psi.conj()), effects


def classical_refinement(branches):
    """A fixed mixed copy p(x) delta(e,x) realizes any classical refinement."""
    p = [sum(branch[x] for branch in branches) for x in range(len(branches[0]))]
    conditional = [[branch[x] / p[x] if p[x] else Fraction(int(i == 0))
                    for x in range(len(p))] for i, branch in enumerate(branches)]
    return p, conditional


class EnsembleExtensionAuditTests(unittest.TestCase):
    def test_real_and_complex_overcomplete_ensembles_from_same_marginal(self):
        rng = np.random.default_rng(179)
        for complex_model in (False, True):
            for d, count in ((2, 5), (3, 7)):
                raw = rng.normal(size=(d, count))
                if complex_model:
                    raw = raw + 1j * rng.normal(size=raw.shape)
                raw /= np.linalg.norm(raw)
                branches = [np.outer(v, v.conj()) for v in raw.T]
                joint, effects = ensemble_extension(branches)
                np.testing.assert_allclose(sum(effects), np.eye(d), atol=4e-15)
                self.assertAlmostEqual(np.trace(joint @ joint).real, 1)
                for branch, effect in zip(branches, effects):
                    self.assertGreaterEqual(np.linalg.eigvalsh(effect).min(), -2e-15)
                    np.testing.assert_allclose(remote_branch(joint, effect, d), branch, atol=2e-15)
                if not complex_model:
                    self.assertFalse(np.iscomplexobj(joint))
                    self.assertTrue(all(not np.iscomplexobj(e) for e in effects))

    def test_singular_marginal_uses_support_and_allows_zero_branch(self):
        vectors = [np.array(v, dtype=float) for v in ((1, 0, 0), (0, 1, 0), (1, 1, 0))]
        branches = [np.outer(v, v) / 4 for v in vectors] + [np.zeros((3, 3))]
        joint, effects = ensemble_extension(branches)
        self.assertEqual(joint.shape, (6, 6))
        np.testing.assert_allclose(sum(effects), I, atol=1e-15)
        for effect, branch in zip(effects, branches):
            np.testing.assert_allclose(remote_branch(joint, effect, 3), branch, atol=1e-15)

    def test_same_real_bell_state_realizes_z_and_x_ensembles_without_signalling(self):
        psi = np.array([1., 0., 0., 1.]) / np.sqrt(2)
        joint = np.outer(psi, psi)
        for observable in (X, Z):
            effects = [(I + sign * observable) / 2 for sign in (1, -1)]
            branches = [remote_branch(joint, effect, 2) for effect in effects]
            for effect, branch in zip(effects, branches):
                np.testing.assert_allclose(branch, effect / 2, atol=2e-16)
            np.testing.assert_allclose(sum(branches), I / 2, atol=2e-16)

    def test_complex_transpose_in_steering_formula_is_required(self):
        psi = np.array([1., 0., 0., 1.]) / np.sqrt(2)
        joint = np.outer(psi, psi)
        effect = (I + Y) / 2
        np.testing.assert_allclose(remote_branch(joint, effect, 2), (I - Y) / 4, atol=2e-16)
        self.assertGreater(np.linalg.norm(effect / 2 - remote_branch(joint, effect, 2)), .7)

    def test_dephased_flag_has_same_marginal_but_cannot_prepare_x_coherence(self):
        copied_flag = np.diag([.5, 0., 0., .5])
        for observable in (X, Z, Y):
            for sign in (-1, 1):
                branch = remote_branch(copied_flag, (I + sign * observable) / 2, 2)
                self.assertEqual(branch[0, 1], 0)
        np.testing.assert_array_equal(remote_branch(copied_flag, (I + X) / 2, 2), I / 4)
        np.testing.assert_array_equal(remote_branch(copied_flag, I, 2), I / 2)

    def test_arbitrary_mixed_ensemble_branches_also_have_a_povm(self):
        branches = [np.array([[.2, .03], [.03, .1]]), np.array([[.3, -.02], [-.02, .4]])]
        joint, effects = ensemble_extension(branches)
        np.testing.assert_allclose(sum(effects), I, atol=1e-15)
        for branch, effect in zip(branches, effects):
            self.assertGreater(np.linalg.eigvalsh(effect).min(), 0)
            np.testing.assert_allclose(remote_branch(joint, effect, 2), branch, atol=1e-15)

    def test_classical_mixed_copy_realizes_every_supplied_refinement_exactly(self):
        F = Fraction
        refinements = [
            [[F(1, 3), 0, 0], [0, F(2, 3), 0]],
            [[F(1, 6), F(1, 2), 0], [F(1, 6), F(1, 6), 0]],
            [[F(1, 9), F(2, 9), 0]] * 3,
        ]
        for branches in refinements:
            p, conditional = classical_refinement(branches)
            self.assertEqual(p, [F(1, 3), F(2, 3), 0])
            for x in range(3):
                self.assertEqual(sum(c[x] for c in conditional), 1)
            for branch, q in zip(branches, conditional):
                self.assertEqual([p[x] * q[x] for x in range(3)], branch)

    def test_classical_copy_is_mixed_while_quantum_extension_is_pure(self):
        copied_flag = np.diag([.5, 0., 0., .5])
        joint, _ = ensemble_extension([np.diag([.5, 0.]), np.diag([0., .5])])
        self.assertEqual(np.trace(copied_flag @ copied_flag), .5)
        self.assertAlmostEqual(np.trace(joint @ joint), 1)
        np.testing.assert_allclose(remote_branch(joint, I, 2),
                                   remote_branch(copied_flag, I, 2), atol=2e-16)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(EnsembleExtensionAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = {
        "round": 179,
        "real_and_complex_support_formula": "E_i=(D^(-1/2) V^dagger tau_i V D^(-1/2))^T",
        "known_ensemble_matrices_are_inputs_to_construction": True,
        "classical_mixed_copy_supports_all_classical_refinements": True,
        "classical_construction_is_pure": False,
        "unconditioned_marginal_changes_with_environment_measurement": False,
        "remote_ensemble_choice_alone_selects_complex_theory": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("ensemble_extension_audit_results.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
