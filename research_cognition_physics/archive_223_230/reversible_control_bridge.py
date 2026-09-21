"""Round 227: tensor composition excludes symplectic-only pure-state control.

The general theorem is analytic. These finite NumPy certificates check Lie
generators, pure-state orbit ranks and invariant bilinear forms independently.
Exact permission uses the explicitly stated C_path (or closed-group) branch.
"""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import platform
import unittest

import numpy as np

from reversible_dynamics_bridge import dagger, hermitian_basis, sample_state


def symplectic_form(d):
    if d < 2 or d % 2:
        raise ValueError('Even dimension at least two required')
    k = d//2
    return np.block([[np.zeros((k, k)), np.eye(k)], [-np.eye(k), np.zeros((k, k))]])


def sp_generators(d):
    k = d//2
    symplectic_form(d)
    zero = np.zeros((k, k), dtype=complex)
    result = [np.block([[1j*h, zero], [zero, (1j*h).conj()]])
              for h in hermitian_basis(k)]
    for i in range(k):
        for j in range(i, k):
            s = np.zeros((k, k), dtype=complex)
            s[i, j] = s[j, i] = 1
            for b in (s, 1j*s):
                result.append(np.block([[zero, b], [-b.conj(), zero]]))
    return result


def su_generators(d):
    return [1j*h for h in hermitian_basis(d)[1:]]


def local_product_generators(generators):
    d = len(generators[0])
    return ([np.kron(x, np.eye(d)) for x in generators]
            + [np.kron(np.eye(d), x) for x in generators])


def invariant_forms(generators):
    """Solve X^T B + B X = 0 over complex matrices using a PSD Gram operator.

    This is a bilinear form, so transpose rather than adjoint is essential.
    The row-major vectorization identity is checked separately in the tests.
    """
    n = len(generators[0])
    gram = np.zeros((n*n, n*n), dtype=complex)
    for x in generators:
        constraint = np.kron(x.T, np.eye(n))+np.kron(np.eye(n), x.T)
        gram += dagger(constraint)@constraint
    values, vectors = np.linalg.eigh(gram)
    null = values < 1e-9
    forms = [vectors[:, j].reshape(n, n) for j in np.flatnonzero(null)]
    return values, forms


def pure_orbit_rank(generators, psi):
    rho = np.outer(psi, psi.conj())
    tangents = np.column_stack([(x@rho-rho@x).ravel() for x in generators])
    return int(np.linalg.matrix_rank(np.vstack([tangents.real, tangents.imag]), tol=1e-10))


@lru_cache(None)
def certificate(family, d, copies=1):
    generators = sp_generators(d) if family == 'sp' else su_generators(d)
    if copies == 2:
        generators = local_product_generators(generators)
    n = d**copies
    values, forms = invariant_forms(generators)
    residual = max((np.linalg.norm(x.T@b+b@x) for x in generators for b in forms), default=0.)
    return {'family': family, 'local_dimension': d, 'copies': copies,
            'representation_dimension': n, 'generator_count': len(generators),
            'invariant_bilinear_form_dimension': len(forms),
            'symmetric_form_errors': [float(np.linalg.norm(b-b.T)) for b in forms],
            'antisymmetric_form_errors': [float(np.linalg.norm(b+b.T)) for b in forms],
            'max_invariance_residual': float(residual),
            'constraint_gram_smallest_positive_eigenvalue': float(min(values[values > 1e-9]))}


def orbit_certificate(d):
    generators = sp_generators(d)
    rng = np.random.default_rng(22700+d)
    vectors = [np.eye(d)[:, 0].astype(complex)]
    for _ in range(4):
        v = rng.normal(size=d)+1j*rng.normal(size=d)
        vectors.append(v/np.linalg.norm(v))
    return {'dimension': d, 'group_dimension': len(generators),
            'full_unitary_channel_group_dimension': d*d-1,
            'pure_state_orbit_tangent_ranks': [pure_orbit_rank(generators, v) for v in vectors],
            'expected_projective_dimension': 2*d-2}


def mixed_orbit_boundary():
    omega = symplectic_form(4)
    rho = np.diag([.5, 0, .5, 0])
    sigma = np.diag([.5, .5, 0, 0])
    return {'dimension': 4,
            'spectrum_error': float(np.linalg.norm(np.linalg.eigvalsh(rho)-np.linalg.eigvalsh(sigma))),
            'initial_symplectic_reality_error': float(np.linalg.norm(rho@omega-omega@rho.conj())),
            'target_symplectic_reality_error': float(np.linalg.norm(sigma@omega-omega@sigma.conj())),
            'scope': 'A single-system Sp(2) boundary example, not a model satisfying the full composite contract.'}


class ReversibleControlTests(unittest.TestCase):
    def test_sp_generators_are_skew_hermitian_traceless_and_preserve_form(self):
        for d in (2, 4, 6):
            generators = sp_generators(d)
            self.assertEqual(len(generators), (d//2)*(d+1))
            omega = symplectic_form(d)
            for x in generators:
                np.testing.assert_allclose(x+dagger(x), 0, atol=1e-14)
                self.assertAlmostEqual(abs(np.trace(x)), 0)
                np.testing.assert_allclose(x.T@omega+omega@x, 0, atol=1e-14)

    def test_bilinear_constraint_vectorization_uses_transpose(self):
        x = sp_generators(4)[-1]
        b = sample_state(4, seed=227)
        matrix = np.kron(x.T, np.eye(4))+np.kron(np.eye(4), x.T)
        np.testing.assert_allclose(matrix@b.ravel(), (x.T@b+b@x).ravel(), atol=1e-14)

    def test_sp_pure_orbits_have_full_projective_tangent_rank(self):
        for d in (2, 4, 6):
            c = orbit_certificate(d)
            self.assertEqual(c['pure_state_orbit_tangent_ranks'], [2*d-2]*5)

    def test_single_sp_has_one_antisymmetric_invariant(self):
        for d in (2, 4):
            c = certificate('sp', d)
            self.assertEqual(c['invariant_bilinear_form_dimension'], 1)
            self.assertLess(max(c['antisymmetric_form_errors']), 1e-10)
            self.assertLess(c['max_invariance_residual'], 1e-10)

    def test_double_sp_has_only_a_symmetric_invariant(self):
        for d in (2, 4):
            c = certificate('sp', d, 2)
            self.assertEqual(c['invariant_bilinear_form_dimension'], 1)
            self.assertLess(max(c['symmetric_form_errors']), 1e-10)
            self.assertAlmostEqual(c['antisymmetric_form_errors'][0], 2)
            self.assertLess(c['max_invariance_residual'], 1e-10)

    def test_tensor_invariant_is_the_product_of_local_forms(self):
        for d in (2, 4):
            omega = symplectic_form(d)
            b = np.kron(omega, omega)
            np.testing.assert_allclose(b, b.T, atol=1e-14)
            for x in local_product_generators(sp_generators(d)):
                np.testing.assert_allclose(x.T@b+b@x, 0, atol=1e-14)

    def test_su_three_has_no_invariant_bilinear_form_even_after_doubling(self):
        for copies in (1, 2):
            self.assertEqual(certificate('su', 3, copies)['invariant_bilinear_form_dimension'], 0)

    def test_pure_state_control_does_not_imply_all_single_system_gates(self):
        c = mixed_orbit_boundary()
        self.assertEqual(c['spectrum_error'], 0)
        self.assertEqual(c['initial_symplectic_reality_error'], 0)
        self.assertGreater(c['target_symplectic_reality_error'], .5)
        self.assertLess(orbit_certificate(4)['group_dimension'], 4*4-1)


def report():
    return {'round': 227, 'date': '2026-09-20',
            'baseline_result': 'Under C on identity connected groups, closure of the actual reversible group on A tensor A is PU(d_A^2). Every unitary channel on A is therefore approximable with an auxiliary copy.',
            'exact_result': 'If the identity path component is pure-state transitive on every object (C_path), or the actual reversible groups are closed, every existing object admits all unitary channels exactly.',
            'extra_condition_status': 'C_path is explicit regularity for the exact branch; not silently identified with bare connectedness or a continuous state path.',
            'orbit_certificates': [orbit_certificate(d) for d in (2, 4, 6)],
            'bilinear_certificates': [certificate(f, d, c) for f, d in [('sp', 2), ('sp', 4), ('su', 3)] for c in (1, 2)],
            'single_system_boundary': mixed_orbit_boundary(),
            'runtime': {'python': platform.python_version(), 'numpy': np.__version__},
            'not_derived': ['a rate or cost for approximate control', 'all gates from pure-state transitivity alone on one system', 'a selected physical Hamiltonian'],
            'next': 'Use exact or approximate unitary access, with the corresponding qualification, to transfer POVMs and dilate finite CP instruments.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ReversibleControlTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = report()
    data['checks'] = {'run': checks.testsRun, 'failures': 0, 'errors': 0}
    if args.write_results:
        target = Path(__file__).with_name('reversible_control_bridge_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8')) != data:
            raise RuntimeError('Existing result differs; review before replacing it.')
        target.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
