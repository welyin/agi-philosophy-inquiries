"""Round 158: one common input for the complete three-read recovery channel.

The 36 fine branches merge into 18 signed orthogonal errors.  A 19 by 19
Gram matrix computes the full external trace distance, including the target.
Numerical eigenvalues in this module are diagnostics, not interval proofs.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import unittest

import numpy as np

from centered_history_recovery import purification
from classical_joining_history import independent_encoding, record_unitary
from coarsened_history_bound import pure_extension_metrics
from flagged_subject_joining import joining_kraus, random_state
from noisy_coherent_history_recovery import decode_label, recovery_kraus
from quantum_interface_audit import ALPHA


I2 = np.eye(2, dtype=int)
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1, -1])
BIT_ERROR = (1 - ALPHA / 2) / 2


@lru_cache(None)
def grouped_errors():
    """Exact integer W matrices and raw-code counts by Hamming distance.

    Each weight is sum_h counts[h] e**h (1-e)**(3-h) / 6.
    Signs cannot be identified: diag(I,W) and diag(I,-W) differ physically.
    """
    u = (np.kron(Z, np.kron(Z, I2)), np.kron(Z, np.eye(4, dtype=int)),
         np.kron(Z, np.kron(X, I2)), np.kron(I2, np.kron(I2, Z)),
         np.eye(8, dtype=int), np.kron(I2, np.kron(I2, X)))
    matrices, counts = [], []
    for r in range(6):
        # This order makes the grouping independent of the read probabilities.
        for g in range(6):
            w = u[g].T @ u[r]
            j = next((j for j, old in enumerate(matrices)
                      if np.array_equal(w, old)), None)
            if j is None:
                j = len(matrices)
                matrices.append(w)
                counts.append([0] * 4)
            for raw in range(8):
                if decode_label(raw) == g:
                    counts[j][(r ^ raw).bit_count()] += 1
    result = np.array(matrices), tuple(tuple(row) for row in counts)
    result[0].flags.writeable = False
    return result


def grouped_weights(error):
    return tuple(sum(c * error**h * (1-error)**(3-h)
                     for h, c in enumerate(row)) / 6
                 for row in grouped_errors()[1])


@lru_cache(None)
def gram_kernel_numerators():
    """G_ab = Tr(sigma K_ab)/8, including a final identity target.

    sigma = Re(rho_A tensor conjugate(rho_B)); sigma is real PSD, trace one.
    Products of W have only I or Z on the reference wire.  Terms involving
    the imaginary part of the encoded state therefore have zero trace.
    """
    ws = np.concatenate((grouped_errors()[0], np.eye(8, dtype=int)[None]))
    result = np.empty((len(ws), len(ws), 4, 4), dtype=int)
    for a, wa in enumerate(ws):
        for b, wb in enumerate(ws):
            partial = np.trace((wa.T @ wb).reshape(2, 4, 2, 4), axis1=0, axis2=2)
            result[a, b] = 4*np.eye(4, dtype=int) + partial + partial.T
    result.flags.writeable = False
    return result


def source_sigma(rho_a, rho_b):
    return np.kron(rho_a, rho_b.conj()).real


def bloch_state(vector):
    return (I2 + vector[0]*X + vector[1]*Y + vector[2]*Z) / 2


def bare_gram(sigma):
    return np.einsum('abij,ji->ab', gram_kernel_numerators(), sigma) / 8


def weighted_gram(sigma, error=BIT_ERROR):
    scale = np.r_[np.sqrt(np.array(grouped_weights(error), dtype=float)), 1]
    return bare_gram(sigma) * scale[:, None] * scale[None, :]


def error_from_sigma(sigma, error=BIT_ERROR):
    gram = weighted_gram(sigma, error)
    vals, vecs = np.linalg.eigh(gram)
    root = (vecs * np.sqrt(np.maximum(vals, 0))) @ vecs.T
    signs = np.r_[np.ones(len(gram)-1), -1]
    # Output minus a pure target has at most one negative eigenvalue.
    return float(max(0, -np.linalg.eigvalsh((root*signs) @ root)[0]))


def supporting_coefficients(sigma, error=BIT_ERROR):
    """Floating candidate a for B=I-sum a_j A_j; certify separately."""
    gram = weighted_gram(sigma, error)
    distance = error_from_sigma(sigma, error)
    c = np.linalg.solve(gram[:-1, :-1] + distance*np.eye(len(gram)-1), gram[:-1, -1])
    return c*np.sqrt(np.array(grouped_weights(error), dtype=float))


class JointHistoryChannelTests(unittest.TestCase):
    def test_integer_records_and_signed_groups(self):
        ws, counts = grouped_errors()
        self.assertEqual(len(ws), 18)
        self.assertEqual(tuple(np.sum(counts, axis=0)), (6, 18, 18, 6))
        for r in range(6):
            for g in range(6):
                old = record_unitary(g+1).T @ record_unitary(r+1)
                self.assertTrue(any(np.array_equal(old, w) for w in ws))
        for w in ws:
            np.testing.assert_array_equal(w.T @ w, np.eye(8))

    def test_grouping_preserves_full_channel_on_unrestricted_complex_input(self):
        rng = np.random.default_rng(158)
        c = rng.normal(size=(16, 5)) + 1j*rng.normal(size=(16, 5))
        rho = c @ c.conj().T
        rho /= np.trace(rho)
        k = joining_kraus((2, 2))
        t = np.vstack((k[(0, 0)], k[(0, 1)]))
        ws, _ = grouped_errors()
        for e in (F(0), F(1, 4), F(1, 2)):
            weights = grouped_weights(e)
            self.assertEqual(sum(weights), 1)
            merged = np.zeros_like(rho)
            for w, p in zip(ws, weights):
                a = t.T @ np.block([[np.eye(8), np.zeros((8, 8))],
                                     [np.zeros((8, 8)), w]]) @ t
                merged += float(p)*a @ rho @ a.T
            direct = sum(a @ rho @ a.T for a in recovery_kraus(e))
            np.testing.assert_allclose(merged, direct, atol=2e-16)

    def test_kernel_matches_full_purification_for_complex_independent_sources(self):
        rng = np.random.default_rng(258)
        ks = recovery_kraus(BIT_ERROR)
        for _ in range(8):
            states = (random_state(rng, 2), random_state(rng, 2))
            c = purification(independent_encoding(states))
            direct = pure_extension_metrics(ks, c)['trace_error']
            self.assertAlmostEqual(error_from_sigma(source_sigma(*states)), direct, places=13)

    def test_y_components_survive_as_a_joint_correlation(self):
        a = bloch_state((.2, .5, -.3))
        b = bloch_state((-.4, .3, .1))
        g = bare_gram(source_sigma(a, b))
        np.testing.assert_allclose(g, bare_gram(source_sigma(a.conj(), b.conj())), atol=1e-16)
        self.assertGreater(np.linalg.norm(g-bare_gram(source_sigma(a.conj(), b))), .01)

    def test_relaxation_contains_the_promised_sources_but_is_strict(self):
        rng = np.random.default_rng(358)
        for _ in range(8):
            sigma = source_sigma(random_state(rng, 2), random_state(rng, 2))
            self.assertGreaterEqual(np.linalg.eigvalsh(sigma).min(), -1e-16)
            self.assertAlmostEqual(np.trace(sigma), 1.)
        bell = np.array([1., 0., 0., 1.])/np.sqrt(2)
        sigma = np.outer(bell, bell)
        partial_transpose = sigma.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
        self.assertLess(np.linalg.eigvalsh(partial_transpose).min(), -.49)
        # Re(product) is the average of two product density operators and
        # has positive partial transpose. This real Bell state cannot be it.

    def test_zero_noise_and_complete_output_spectrum(self):
        sigma = np.eye(4)/4
        self.assertLess(error_from_sigma(sigma, 0), 1e-14)
        g = weighted_gram(sigma)
        vals, vecs = np.linalg.eigh(g)
        root = (vecs*np.sqrt(np.maximum(vals, 0))) @ vecs.T
        eig = np.linalg.eigvalsh((root*np.r_[np.ones(18), -1]) @ root)
        self.assertEqual(np.sum(eig < -1e-12), 1)
        self.assertAlmostEqual(eig.sum(), 0., places=14)
        self.assertAlmostEqual(error_from_sigma(sigma), .4361668361397797, places=13)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(JointHistoryChannelTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {'round': 158, 'signed_error_groups': 18, 'gram_dimension': 19,
              'group_hamming_counts': grouped_errors()[1],
              'weights_formula': 'sum_h counts[h] e^h (1-e)^(3-h) / 6',
              'effective_state': 'Re(rho_A tensor conjugate(rho_B))',
              'maximally_mixed_input_error_diagnostic': error_from_sigma(np.eye(4)/4),
              'all_real_4x4_density_matrices_are_only_an_upper_relaxation': True,
              'numerical_eigenvalues_are_not_global_certificates': True,
              'automated_checks': {'run': checks.testsRun, 'failures': len(checks.failures),
                                   'errors': len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name('joint_history_channel_results.json').write_text(
            json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
