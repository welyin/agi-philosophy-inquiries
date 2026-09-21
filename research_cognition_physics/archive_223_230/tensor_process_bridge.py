"""Round 224: binary tensor identification and CP endomorphism certificates.

The all-dimension proof uses positive discards, product preparations, local
tomography, and the inherited full complex state cones. Tests below exercise
the resulting normal forms; they are not a sampled proof of classification.
Only NumPy and the existing round 223 helpers are used.
"""
import argparse
import itertools
import json
from pathlib import Path
import platform
import unittest

import numpy as np

from reversible_dynamics_bridge import (
    dagger, hermitian_basis, matrix_unit, sample_generator, sample_state, unitary,
)


def partial_transpose(x, a, b, left=False, right=False):
    axes = [0, 1, 2, 3]
    if left:
        axes[0], axes[2] = axes[2], axes[0]
    if right:
        axes[1], axes[3] = axes[3], axes[1]
    return x.reshape(a, b, a, b).transpose(axes).reshape(a*b, a*b)


def partial_trace(x, a, b, keep):
    t = x.reshape(a, b, a, b)
    if keep == 'A':
        return np.trace(t, axis1=1, axis2=3)
    if keep == 'B':
        return np.trace(t, axis1=0, axis2=2)
    raise ValueError(keep)


class BinaryRepresentation:
    """Normal form pi(X,Y)=U(T^e X tensor T^f Y)U^dagger.

    encode/decode are coordinate maps, not physically implementable channels.
    In a twisted product chart decode can return a non-PSD Hermitian matrix.
    """
    def __init__(self, a, b, orientation=(0, 0)):
        self.a, self.b = a, b
        self.e, self.f = orientation
        self.u = unitary(sample_generator(a*b, seed=224), 0.317)

    def encode(self, x):
        t = partial_transpose(x, self.a, self.b, self.e, self.f)
        return self.u @ t @ dagger(self.u)

    def decode(self, x):
        t = dagger(self.u) @ x @ self.u
        return partial_transpose(t, self.a, self.b, self.e, self.f)

    def product(self, x, y):
        return self.encode(np.kron(x, y))

    def marginal(self, x, keep):
        t = dagger(self.u) @ x @ self.u
        result = partial_trace(t, self.a, self.b, keep)
        transpose = self.e if keep == 'A' else self.f
        return result.T if transpose else result

    def left(self, x):
        return self.product(x, np.eye(self.b))

    def right(self, y):
        return self.product(np.eye(self.a), y)


def frame(d, seed):
    u = unitary(sample_generator(d, seed=seed), 0.421)
    return [np.outer(u[:, j], u[:, j].conj()) for j in range(d)]


def frame_certificate(a, b, orientation):
    rep = BinaryRepresentation(a, b, orientation)
    pa, qb = frame(a, 2241), frame(b, 2242)
    products = [rep.product(p, q) for p in pa for q in qb]
    gram = np.array([[np.trace(x @ y) for y in products] for x in products])
    product_basis = [rep.product(x, y) for x in hermitian_basis(a)
                     for y in hermitian_basis(b)]
    full_gram = np.array([[np.trace(x @ y) for y in product_basis]
                          for x in product_basis])
    xa, ya = sample_generator(a, 2243), sample_generator(a, 2244)
    xb = sample_generator(b, 2245)
    ja, jb = rep.left(xa), rep.right(xb)
    marginal_errors = [np.linalg.norm(rep.marginal(rep.product(p, q), 'A')-p)
                       + np.linalg.norm(rep.marginal(rep.product(p, q), 'B')-q)
                       for p in pa for q in qb]
    # Reconstruct products using only adjoints of the marginals (left/right).
    reconstruction = max(np.linalg.norm(rep.left(p) @ rep.right(q)
                                         - rep.product(p, q)) for p in pa for q in qb)
    w = sample_state(a*b, seed=2246)
    duality = abs(np.trace(rep.left(xa) @ w)
                  - np.trace(xa @ rep.marginal(w, 'A')))
    return {
        'dimensions': [a, b], 'local_transposes': list(orientation),
        'frame_gram_error': float(np.linalg.norm(gram-np.eye(a*b))),
        'frame_sum_error': float(np.linalg.norm(sum(products)-np.eye(a*b))),
        'frame_idempotence_error': float(max(np.linalg.norm(p@p-p) for p in products)),
        'marginal_recovery_error': float(max(marginal_errors)),
        'marginal_adjoint_duality_error': float(duality),
        'product_reconstruction_error': float(reconstruction),
        'jordan_embedding_error': float(np.linalg.norm(
            rep.left((xa@ya+ya@xa)/2)-(ja@rep.left(ya)+rep.left(ya)@ja)/2)),
        'local_commutator_error': float(np.linalg.norm(ja@jb-jb@ja)),
        'full_product_basis_gram_error': float(np.linalg.norm(full_gram-np.eye(a*a*b*b))),
        'product_basis_rank': int(np.linalg.matrix_rank(full_gram)),
    }


def apply_kraus(kraus, x):
    return sum((k @ x @ dagger(k) for k in kraus),
               np.zeros((kraus[0].shape[0],)*2, dtype=complex))


def random_channel(d, count=3):
    rng = np.random.default_rng(22400+d)
    raw = [rng.normal(size=(d, d))+1j*rng.normal(size=(d, d)) for _ in range(count)]
    g = sum(dagger(k) @ k for k in raw)
    vals, vecs = np.linalg.eigh(g)
    invsqrt = (vecs*(1/np.sqrt(vals))) @ dagger(vecs)
    return [k @ invsqrt for k in raw]


def extend_left(action, x, a, b):
    """Complex-linear extension, implemented on matrix units, no CP assumption."""
    out_dim = action(matrix_unit(a, 0, 0)).shape[0]
    out = np.zeros((out_dim*b, out_dim*b), dtype=complex)
    blocks = x.reshape(a, b, a, b)
    for i in range(a):
        for j in range(a):
            out += np.kron(action(matrix_unit(a, i, j)), blocks[i, :, j, :])
    return out


def maximally_entangled(d):
    vector = np.eye(d).ravel()/np.sqrt(d)
    return np.outer(vector, vector)


def choi(action, d):
    """Unnormalized Choi matrix in output tensor input convention."""
    return d * extend_left(action, maximally_entangled(d), d, d)


def kraus_from_choi(c, input_dim, tol=1e-11):
    vals, vecs = np.linalg.eigh(c)
    if vals.min() < -tol:
        raise ValueError('Choi matrix is not positive semidefinite')
    output_dim = c.shape[0]//input_dim
    return [(np.sqrt(v)*vecs[:, i]).reshape(output_dim, input_dim)
            for i, v in enumerate(vals) if v > tol]


def cp_certificate(d, orientation):
    rep = BinaryRepresentation(d, d, orientation)
    ks = random_channel(d)
    phi = lambda x: apply_kraus(ks, x)
    oriented = lambda x: phi(x.T).T if rep.e else phi(x)
    oriented_ks = [k.conj() if rep.e else k for k in ks]
    w = sample_state(d*d, seed=2247)
    # The operational linear extension is fixed on all product inputs by F+L.
    actual = rep.encode(extend_left(phi, rep.decode(w), d, d))
    expected = rep.u @ extend_left(oriented, dagger(rep.u)@w@rep.u, d, d) @ dagger(rep.u)
    c = choi(phi, d)
    recovered = kraus_from_choi(c, d)
    basis = hermitian_basis(d)
    return {
        'dimension': d, 'local_transposes': list(orientation),
        'extension_formula_error': float(np.linalg.norm(actual-expected)),
        'extension_output_min_eigenvalue': float(np.linalg.eigvalsh(actual).min()),
        'oriented_kraus_error': float(max(np.linalg.norm(oriented(x)-apply_kraus(oriented_ks, x)) for x in basis)),
        'choi_min_eigenvalue': float(np.linalg.eigvalsh(c).min()),
        'choi_rank': len(recovered),
        'kraus_reconstruction_error': float(max(np.linalg.norm(phi(x)-apply_kraus(recovered, x)) for x in basis)),
        'trace_preservation_error': float(np.linalg.norm(sum(dagger(k)@k for k in recovered)-np.eye(d))),
    }


def noisy_transpose(d, p):
    return lambda x: p*x.T+(1-p)*np.trace(x)*np.eye(d)/d


def noise_certificate(d):
    threshold = 1/(d+1)
    rows = []
    for p in (0.0, threshold/2, threshold, (1+threshold)/2, 1.0):
        actual = float(np.linalg.eigvalsh(choi(noisy_transpose(d, p), d)/d).min())
        rows.append({'p': p, 'normalized_choi_min_eigenvalue': actual,
                     'analytic_min_eigenvalue': (1-(d+1)*p)/(d*d)})
    return {'dimension': d, 'exact_cp_threshold': f'1/{d+1}', 'samples': rows}


def orientation_witness(d=2):
    rep = BinaryRepresentation(d, d, (1, 0))
    w = rep.u @ maximally_entangled(d) @ dagger(rep.u)
    naive = rep.decode(w)
    # Negative naive eigenvalues describe a chart mismatch, not a bad global state.
    changed = rep.encode(extend_left(lambda x: x.T, naive, d, d))
    return {'dimension': d,
            'global_state_min_eigenvalue': float(np.linalg.eigvalsh(w).min()),
            'unoriented_product_chart_min_eigenvalue': float(np.linalg.eigvalsh(naive).min()),
            'roundtrip_error': float(np.linalg.norm(rep.encode(naive)-w)),
            'physical_local_transpose_output_min_eigenvalue': float(np.linalg.eigvalsh(changed).min())}


class TensorBridgeTests(unittest.TestCase):
    def test_distinguishing_frames_saturate_joint_hilbert_dimension(self):
        for a, b in ((2, 2), (2, 3), (3, 3)):
            for orientation in itertools.product((0, 1), repeat=2):
                c = frame_certificate(a, b, orientation)
                for key in ('frame_gram_error', 'frame_sum_error', 'frame_idempotence_error'):
                    self.assertLess(c[key], 2e-13, (a, b, orientation, key))

    def test_positive_marginals_and_duality_in_all_orientations(self):
        for orientation in itertools.product((0, 1), repeat=2):
            rep = BinaryRepresentation(2, 3, orientation)
            w = sample_state(6, 2248)
            for side in ('A', 'B'):
                self.assertGreaterEqual(np.linalg.eigvalsh(rep.marginal(w, side)).min(), -1e-13)
                self.assertAlmostEqual(np.trace(rep.marginal(w, side)).real, 1)
            c = frame_certificate(2, 3, orientation)
            self.assertLess(c['marginal_adjoint_duality_error'], 2e-13)
            self.assertLess(c['marginal_recovery_error'], 2e-13)

    def test_marginal_adjoints_reconstruct_products_and_jordan_embeddings(self):
        for orientation in itertools.product((0, 1), repeat=2):
            c = frame_certificate(2, 3, orientation)
            for key in ('product_reconstruction_error', 'jordan_embedding_error', 'local_commutator_error'):
                self.assertLess(c[key], 2e-13)

    def test_local_tomography_and_hilbert_schmidt_factorization(self):
        for a, b in ((2, 2), (2, 3), (3, 3)):
            c = frame_certificate(a, b, (1, 0))
            self.assertEqual(c['product_basis_rank'], a*a*b*b)
            self.assertLess(c['full_product_basis_gram_error'], 2e-13)

    def test_product_effect_probability_pairing_without_effect_permission_assumption(self):
        rep = BinaryRepresentation(2, 3, (1, 0))
        rho, sigma = sample_state(2), sample_state(3)
        e, f = frame(2, 2249)[0], frame(3, 2250)[0]
        actual = np.trace(rep.product(e, f) @ rep.product(rho, sigma))
        self.assertAlmostEqual(actual.real, (np.trace(e@rho)*np.trace(f@sigma)).real)
        self.assertGreaterEqual(np.linalg.eigvalsh(rep.product(e, f)).min(), -1e-13)

    def test_cp_endomorphism_extension_does_not_require_a_global_orientation_choice(self):
        for d in (2, 3):
            for orientation in itertools.product((0, 1), repeat=2):
                c = cp_certificate(d, orientation)
                self.assertLess(c['extension_formula_error'], 2e-13)
                self.assertLess(c['oriented_kraus_error'], 2e-13)
                self.assertGreaterEqual(c['extension_output_min_eigenvalue'], -1e-13)

    def test_choi_factorization_reconstructs_complex_kraus_and_tp(self):
        for d in (2, 3, 4):
            c = cp_certificate(d, (1, 1))
            self.assertEqual(c['choi_rank'], 3)
            self.assertLess(c['kraus_reconstruction_error'], 2e-13)
            self.assertLess(c['trace_preservation_error'], 2e-13)

    def test_success_branch_is_cp_and_trace_nonincreasing_before_normalization(self):
        k = np.diag([1, 0.5j])
        phi = lambda x: k@x@dagger(k)
        reconstructed = kraus_from_choi(choi(phi, 2), 2)
        np.testing.assert_allclose(sum(dagger(v)@v for v in reconstructed), np.diag([1, .25]), atol=2e-14)
        self.assertAlmostEqual(np.trace(phi(np.eye(2)/2)).real, .625)

    def test_noisy_transpose_threshold_is_analytic_not_an_optimization_fit(self):
        for d in (2, 3, 4):
            for row in noise_certificate(d)['samples']:
                self.assertAlmostEqual(row['normalized_choi_min_eigenvalue'], row['analytic_min_eigenvalue'])
            phi = noisy_transpose(d, 1)
            self.assertGreaterEqual(np.linalg.eigvalsh(phi(sample_state(d))).min(), -1e-13)
            with self.assertRaises(ValueError):
                kraus_from_choi(choi(phi, d), d)

    def test_active_transpose_fails_in_every_pair_orientation(self):
        for orientation in itertools.product((0, 1), repeat=2):
            rep = BinaryRepresentation(2, 2, orientation)
            w = rep.u @ maximally_entangled(2) @ dagger(rep.u)
            out = rep.encode(extend_left(lambda x: x.T, rep.decode(w), 2, 2))
            self.assertAlmostEqual(np.linalg.eigvalsh(out).min(), -.5)

    def test_coordinate_transpose_is_distinct_from_an_active_local_operation(self):
        c = orientation_witness()
        self.assertGreaterEqual(c['global_state_min_eigenvalue'], -1e-13)
        self.assertAlmostEqual(c['unoriented_product_chart_min_eigenvalue'], -.5)
        self.assertLess(c['roundtrip_error'], 2e-13)
        self.assertAlmostEqual(c['physical_local_transpose_output_min_eigenvalue'], -.5)


def report():
    return {
        'round': 224, 'date': '2026-09-20',
        'inherited': 'Full complex density-state sets; F supplies positive normalized discards, compatible product preparations and local process extensions; P supplies local tomography.',
        'analytic_result': 'Each binary product has a unitary-conjugated tensor form up to independent local transposes. Every actual A-to-A channel or unnormalized event is CP (TP or TNI respectively), independent of the binary orientation choices.',
        'additional_physical_axioms': [],
        'permission_caution': 'Marginal adjoints of arbitrary positive functionals are mathematical effects, not a grant of all measurement or gate permissions.',
        'normal_form_certificates': [frame_certificate(a, b, e) for a, b in ((2, 2), (2, 3), (3, 3)) for e in itertools.product((0, 1), repeat=2)],
        'cp_certificates': [cp_certificate(d, e) for d in (2, 3) for e in itertools.product((0, 1), repeat=2)],
        'noisy_transpose_boundary': [noise_certificate(d) for d in (2, 3, 4)],
        'orientation_witness': orientation_witness(),
        'runtime': {'python': platform.python_version(), 'numpy': np.__version__},
        'evidence_scope': 'Exact structural argument in the note plus finite numerical certificates of its normal forms and CP tests; no statistical claim of proving all-dimensional classification.',
        'not_derived': ['one coherent matrix orientation for every system and grouping', 'all A-to-B processes CP in an arbitrary independently chosen chart', 'all CP operations physically available', 'a particular Hamiltonian, physical clock or gravity'],
        'next': 'Audit orientation coherence using associativity/reordering, then cross-system CP and operational saturation separately.',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(TensorBridgeTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = report()
    data['checks'] = {'run': checks.testsRun, 'failures': 0, 'errors': 0}
    if args.write_results:
        target = Path(__file__).with_name('tensor_process_bridge_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8')) != data:
            raise RuntimeError('Existing scientific result differs; review before replacing it.')
        target.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
