"""Round 173: real quantum theory has purification; purity is representation-relative.

Finite-dimensional matrix models are declared comparison theories, not deductions
from cognition. Algebraic arguments are recorded in research_note_173.md.
"""

import argparse
import json
from pathlib import Path
import unittest

import numpy as np

from common_orientation_structure import encode_state, decode_state, orientation


J = np.array([[0., -1.], [1., 0.]])
X = np.array([[0., 1.], [1., 0.]])


def purification_coefficients(rho, environment_dimension=None):
    """Rows index system A, columns environment E; C C^dagger = rho."""
    rho = np.asarray(rho)
    values, vectors = np.linalg.eigh(rho)
    if values.min() < -1e-12:
        raise ValueError("Density matrix is not positive.")
    keep = values > 1e-12
    rank = int(keep.sum())
    m = rank if environment_dimension is None else environment_dimension
    if m < rank:
        raise ValueError("Environment dimension is smaller than rank.")
    coefficients = np.zeros((len(rho), m), dtype=rho.dtype)
    coefficients[:, :rank] = vectors[:, keep] * np.sqrt(values[keep])
    return coefficients


def environment_alignment(first, second):
    """Recover U_E with second = first @ U_E.T when row Gram matrices agree."""
    if first.shape != second.shape:
        raise ValueError("Compare purifications on the same environment.")
    if not np.allclose(first @ first.conj().T, second @ second.conj().T,
                       atol=2e-12, rtol=0):
        raise ValueError("Marginals differ.")
    left, _, right_h = np.linalg.svd(first.conj().T @ second)
    return (left @ right_h).T


class RealPurificationAuditTests(unittest.TestCase):
    def test_real_and_complex_states_have_purifications_including_rank_deficient(self):
        rng = np.random.default_rng(173)
        for complex_model in (False, True):
            for d, rank in ((2, 1), (3, 2), (4, 4)):
                raw = rng.normal(size=(d, rank))
                if complex_model:
                    raw = raw + 1j * rng.normal(size=(d, rank))
                rho = raw @ raw.conj().T
                rho /= np.trace(rho)
                coefficients = purification_coefficients(rho)
                self.assertEqual(coefficients.shape[1], rank)
                np.testing.assert_allclose(coefficients @ coefficients.conj().T, rho, atol=2e-14)
                psi = coefficients.ravel()
                pure = np.outer(psi, psi.conj())
                np.testing.assert_allclose(pure @ pure, pure, atol=2e-14)
                if not complex_model:
                    self.assertFalse(np.iscomplexobj(coefficients))

    def test_essential_uniqueness_uses_environment_only_even_with_unused_dimensions(self):
        rng = np.random.default_rng(2173)
        for complex_model in (False, True):
            for rank in (1, 2, 3):
                raw = rng.normal(size=(3, rank))
                env_raw = rng.normal(size=(5, 5))
                if complex_model:
                    raw = raw + 1j * rng.normal(size=raw.shape)
                    env_raw = env_raw + 1j * rng.normal(size=env_raw.shape)
                rho = raw @ raw.conj().T
                rho /= np.trace(rho)
                first = purification_coefficients(rho, 5)
                env, _ = np.linalg.qr(env_raw)
                second = first @ env.T
                recovered = environment_alignment(first, second)
                np.testing.assert_allclose(recovered.conj().T @ recovered, np.eye(5), atol=3e-14)
                np.testing.assert_allclose(first @ recovered.T, second, atol=3e-14)

    def test_rank_is_a_necessary_environment_capacity(self):
        with self.assertRaises(ValueError):
            purification_coefficients(np.eye(3) / 3, 2)
        self.assertEqual(purification_coefficients(np.diag([.4, .6, 0.])).shape, (3, 2))

    def test_real_orthogonal_pure_seed_extension_generates_local_randomness(self):
        orthogonal = (4 * np.eye(4) + 3 * np.kron(X, J)) / 5
        np.testing.assert_allclose(orthogonal.T @ orthogonal, np.eye(4), atol=2e-16)
        psi = orthogonal[:, 0].reshape(2, 2)
        np.testing.assert_allclose(psi @ psi.T, np.diag([16 / 25, 9 / 25]), atol=1e-16)
        for rho in (np.array([[.3, .2], [.2, .7]]), np.diag([1., 0.])):
            initial = np.kron(rho, np.diag([1., 0.]))
            joint = orthogonal @ initial @ orthogonal.T
            output = np.trace(joint.reshape(2, 2, 2, 2), axis1=1, axis2=3)
            np.testing.assert_allclose(output, (16 * rho + 9 * X @ rho @ X) / 25, atol=2e-16)

    def test_real_pure_bell_and_classical_correlated_mixture_are_different_wholes(self):
        bell = np.outer([1., 0., 0., 1.], [1., 0., 0., 1.]) / 2
        classical = np.diag([.5, 0., 0., .5])
        for state in (bell, classical):
            np.testing.assert_array_equal(np.trace(state.reshape(2, 2, 2, 2), axis1=1, axis2=3), np.eye(2) / 2)
        self.assertEqual(np.trace(bell @ bell), 1)
        self.assertEqual(np.trace(classical @ classical), .5)
        self.assertEqual(np.trace(np.kron(X, X) @ bell), 1)
        self.assertEqual(np.trace(np.kron(X, X) @ classical), 0)

    def test_encoded_complex_pure_state_is_extreme_in_restricted_cone_despite_rank_two(self):
        pure = np.array([[1., -1j], [1j, 1.]]) / 2
        encoded = encode_state(pure)
        self.assertEqual(np.linalg.matrix_rank(encoded), 2)
        self.assertEqual(np.trace(encoded @ encoded), .5)
        np.testing.assert_array_equal(decode_state(encoded), pure)
        np.testing.assert_array_equal(encoded @ orientation(2), orientation(2) @ encoded)
        # The rank-one constituents of its ordinary-real spectral decomposition
        # are outside the allowed commuting cone, so do not witness mixing there.
        values, vectors = np.linalg.eigh(encoded)
        for vector in vectors[:, values > .1].T:
            projector = np.outer(vector, vector)
            self.assertGreater(np.linalg.norm(projector @ orientation(2) - orientation(2) @ projector), 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(RealPurificationAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 173,
        "declared_models": ["ordinary real quantum theory", "complex quantum theory", "equivalent restricted real encoding"],
        "real_quantum_state_purification_exists": True,
        "real_quantum_state_purification_unique_up_to_environment_orthogonal": True,
        "pure_seed_real_orthogonal_bitflip_probability": "9/25",
        "purification_selects_complex_over_ordinary_real_quantum": False,
        "encoded_complex_pure_state_ordinary_real_matrix_rank": 2,
        "operational_purity_is_extremality_in_declared_state_set": True,
        "purification_derived_from_cognition": False,
        "numerical_checks_are_examples_not_general_proofs": True,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("real_purification_audit_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
