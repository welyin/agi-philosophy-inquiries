"""Round 223: finite certificates for the state-space to reversible-dynamics bridge.

Classification of affine automorphisms is a cited theorem, not inferred from
these tests. Time homogeneity and a continuous reversible clock are explicit
additional inputs. Uses the existing Python/NumPy runtime only.
"""
import argparse
import json
from pathlib import Path
import platform
import unittest

import numpy as np


def dagger(a):
    return a.conj().T


def matrix_unit(d, j, k):
    out = np.zeros((d, d), dtype=complex)
    out[j, k] = 1
    return out


def hermitian_basis(d):
    """Hilbert-Schmidt orthonormal, identity first, remaining matrices traceless."""
    basis = [np.eye(d, dtype=complex) / np.sqrt(d)]
    for j in range(d):
        for k in range(j + 1, d):
            e = matrix_unit(d, j, k)
            basis.extend([(e + dagger(e)) / np.sqrt(2),
                          -1j * (e - dagger(e)) / np.sqrt(2)])
    for k in range(1, d):
        diagonal = np.zeros(d)
        diagonal[:k], diagonal[k] = 1, -k
        basis.append(np.diag(diagonal) / np.sqrt(k * (k + 1)))
    return basis


def sample_generator(d, seed=223):
    rng = np.random.default_rng(seed + d)
    a = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    k = (a + dagger(a)) / 2
    return k - np.trace(k) * np.eye(d) / d


def sample_state(d, seed=224):
    rng = np.random.default_rng(seed + d)
    a = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    rho = a @ dagger(a)
    return rho / np.trace(rho)


def unitary(k, t):
    eigenvalues, vectors = np.linalg.eigh(k)
    return (vectors * np.exp(-1j * t * eigenvalues)) @ dagger(vectors)


def evolve(k, t, x):
    u = unitary(k, t)
    return u @ x @ dagger(u)


def derivation(k, x):
    return -1j * (k @ x - x @ k)


def real_superoperator(action, basis):
    data = np.array([[np.trace(a @ action(b)) for b in basis] for a in basis])
    if np.max(np.abs(data.imag)) > 1e-10:
        raise ValueError('Map does not preserve the Hermitian real vector space')
    return data.real


def recover_inner_generator(action, d):
    """Matrix-unit construction A=sum_j D(E_j0) E_0j, K=i antiHerm(A)."""
    a = sum((action(matrix_unit(d, j, 0)) @ matrix_unit(d, 0, j)
             for j in range(d)), np.zeros((d, d), dtype=complex))
    k = 1j * (a - dagger(a)) / 2
    return k - np.trace(k) * np.eye(d) / d


def fit_commutator(superoperator, basis):
    design = np.column_stack([
        real_superoperator(lambda x, h=h: derivation(h, x), basis).ravel()
        for h in basis[1:]])
    coefficients, _, rank, _ = np.linalg.lstsq(design, superoperator.ravel(), rcond=None)
    k = sum((a * h for a, h in zip(coefficients, basis[1:])),
            np.zeros_like(basis[0]))
    return k, int(rank), float(np.linalg.norm(design @ coefficients - superoperator.ravel()))


def leibniz_defect(action, d):
    units = [matrix_unit(d, j, k) for j in range(d) for k in range(d)]
    return float(max(np.linalg.norm(action(a @ b) - action(a) @ b - a @ action(b))
                     for a in units for b in units))


def trace_distance(a, b):
    return float(np.sum(np.abs(np.linalg.eigvalsh((a - b + dagger(a - b)) / 2))) / 2)


def clock_witness():
    k = np.diag([0.5, -0.5])
    plus = np.ones((2, 2), dtype=complex) / 2
    stationary = float(np.trace(plus @ plus).real)
    rotating = float(np.trace(plus @ evolve(k, np.pi, plus)).real)
    nonlinear = trace_distance(evolve(k, 4, plus), evolve(k, 2, plus))
    return {'same_state_space_plus_probability_at_pi':
                {'zero_generator': stationary, 'z_over_two_generator': rotating},
            'nonhomogeneous_path_t_squared_group_defect_at_s_t_1': nonlinear,
            'exact_group_defect': 'sin(1)',
            'nonhomogeneous_path_derivative_at_zero': 0}


def filter_witness():
    f = np.diag([1, 0.5])
    rho0, rho1 = np.diag([1., 0.]), np.diag([0., 1.])
    def conditional(rho):
        out = f @ rho @ f
        return out / np.trace(out)
    mixed_then_filtered = conditional((rho0 + rho1) / 2)
    filtered_then_mixed = (conditional(rho0) + conditional(rho1)) / 2
    return {'filter_diagonal': [1, 0.5],
            'success_probabilities_on_basis_states': [1, 0.25],
            'conditioned_mixture_diagonal': np.diag(mixed_then_filtered).tolist(),
            'mixture_of_conditioned_states_diagonal': np.diag(filtered_then_mixed).tolist(),
            'affinity_defect_trace_distance': trace_distance(mixed_then_filtered, filtered_then_mixed),
            'exact_affinity_defect': '3/10'}


def dimension_certificate(d):
    k, rho = sample_generator(d), sample_state(d)
    basis = hermitian_basis(d)
    action = lambda x: derivation(k, x)
    recovered = recover_inner_generator(action, d)
    eps = 1e-5
    finite_difference = lambda x: (evolve(k, eps, x) - evolve(k, -eps, x)) / (2 * eps)
    estimated, rank, residual = fit_commutator(real_superoperator(finite_difference, basis), basis)
    dissipative = lambda x: np.trace(x) * np.eye(d) / d - x
    _, _, noncommutator_residual = fit_commutator(real_superoperator(dissipative, basis), basis)
    transpose = real_superoperator(lambda x: x.T, basis)
    a, b = matrix_unit(d, 0, 1), matrix_unit(d, 1, 0)
    return {'dimension': d, 'generator_rank_mod_identity': rank,
            'inner_generator_recovery_error': float(np.linalg.norm(recovered - k)),
            'finite_difference_generator_error': float(np.linalg.norm(estimated - k)),
            'finite_difference_commutator_fit_residual': residual,
            'finite_difference_state_derivative_error': float(np.linalg.norm(finite_difference(rho) - action(rho))),
            'leibniz_defect': leibniz_defect(action, d),
            'depolarizing_generator_commutator_fit_residual': noncommutator_residual,
            'depolarizing_generator_leibniz_defect': leibniz_defect(dissipative, d),
            'transpose_real_determinant': int(round(np.linalg.det(transpose))),
            'transpose_expected_determinant': (-1) ** (d * (d - 1) // 2),
            'transpose_multiplicativity_defect': float(np.linalg.norm((a @ b).T - a.T @ b.T)),
            'transpose_antimultiplicativity_defect': float(np.linalg.norm((a @ b).T - b.T @ a.T))}


class BridgeTests(unittest.TestCase):
    def test_basis_and_real_coordinates(self):
        for d in (2, 3, 4):
            b = hermitian_basis(d)
            self.assertEqual(len(b), d*d)
            np.testing.assert_allclose([[np.trace(x @ y) for y in b] for x in b], np.eye(d*d), atol=1e-14)
            np.testing.assert_allclose([np.trace(x) for x in b[1:]], 0, atol=1e-14)

    def test_matrix_unit_reconstruction_including_scalar_and_degenerate_generators(self):
        for d in (2, 3, 4):
            for k in (sample_generator(d), 3*np.eye(d), np.diag([1.]*(d-1)+[0.])):
                recovered = recover_inner_generator(lambda x: derivation(k, x), d)
                np.testing.assert_allclose(recovered, k-np.trace(k)*np.eye(d)/d, atol=2e-14)
                self.assertLess(leibniz_defect(lambda x: derivation(k, x), d), 2e-14)

    def test_empirical_derivatives_identify_only_traceless_generator(self):
        for d in (2, 3, 4):
            c = dimension_certificate(d)
            self.assertEqual(c['generator_rank_mod_identity'], d*d-1)
            self.assertLess(c['finite_difference_generator_error'], 5e-8)
            self.assertLess(c['finite_difference_state_derivative_error'], 5e-8)

    def test_dissipative_generator_cannot_be_fitted_as_a_commutator(self):
        for d in (2, 3, 4):
            c = dimension_certificate(d)
            self.assertAlmostEqual(c['depolarizing_generator_commutator_fit_residual'], np.sqrt(d*d-1))
            self.assertGreater(c['depolarizing_generator_leibniz_defect'], 0.6)

    def test_group_inverse_trace_and_spectrum(self):
        for d in (2, 3, 4):
            k, rho = sample_generator(d), sample_state(d)
            actual = evolve(k, 0.37, evolve(k, -0.21, rho))
            np.testing.assert_allclose(actual, evolve(k, 0.16, rho), atol=3e-14)
            np.testing.assert_allclose(evolve(k, -0.37, evolve(k, 0.37, rho)), rho, atol=3e-14)
            np.testing.assert_allclose(np.linalg.eigvalsh(actual), np.linalg.eigvalsh(rho), atol=3e-14)
            self.assertAlmostEqual(np.trace(actual).real, 1)

    def test_scalar_gauge_and_clock_units(self):
        for d in (2, 3, 4):
            k, rho = sample_generator(d), sample_state(d)
            np.testing.assert_allclose(evolve(k + 7*np.eye(d), 0.63, rho), evolve(k, 0.63, rho), atol=3e-14)
            np.testing.assert_allclose(evolve(k/4, 4*0.63, rho), evolve(k, 0.63, rho), atol=3e-14)

    def test_continuous_invertible_path_does_not_supply_time_homogeneity(self):
        c = clock_witness()
        self.assertAlmostEqual(c['nonhomogeneous_path_t_squared_group_defect_at_s_t_1'], np.sin(1))
        k, rho, eps = np.diag([0.5, -0.5]), np.ones((2, 2))/2, 1e-5
        np.testing.assert_allclose(evolve(k, eps**2, rho)-evolve(k, (-eps)**2, rho), 0, atol=1e-14)
        self.assertGreater(trace_distance(evolve(k, 1, rho), rho), 0.4)

    def test_postselection_bijection_is_not_an_affine_deterministic_process(self):
        c = filter_witness()
        self.assertAlmostEqual(c['affinity_defect_trace_distance'], 0.3)
        np.testing.assert_allclose(c['conditioned_mixture_diagonal'], [0.8, 0.2])

    def test_transpose_branch_is_not_separated_by_determinant_in_all_dimensions(self):
        for d in (2, 3, 4):
            c = dimension_certificate(d)
            self.assertEqual(c['transpose_real_determinant'], c['transpose_expected_determinant'])
            self.assertAlmostEqual(c['transpose_multiplicativity_defect'], np.sqrt(2))
            self.assertEqual(c['transpose_antimultiplicativity_defect'], 0)
            rho = sample_state(d)
            np.testing.assert_allclose(np.linalg.eigvalsh(rho.T), np.linalg.eigvalsh(rho), atol=1e-14)
        self.assertEqual(dimension_certificate(4)['transpose_real_determinant'], 1)

    def test_same_state_structure_does_not_fix_the_chosen_free_evolution(self):
        p = clock_witness()['same_state_space_plus_probability_at_pi']
        self.assertAlmostEqual(p['zero_generator'], 1)
        self.assertLess(abs(p['z_over_two_generator']), 1e-14)


def report():
    return {'round': 223, 'date': '2026-09-20',
            'inherited': 'F+U+C+P and the completed complex density-state-space reconstruction.',
            'analytic_result': 'Every allowed reversible affine state map is unitary or transpose-unitary; a continuous path from identity stays unitary. A continuous time-homogeneous one-parameter group has generator -i[K,.], K Hermitian and unique modulo scalar identity.',
            'additional_clock_input': 'A declared reversible continuous one-parameter group Phi_(s+t)=Phi_s Phi_t, not supplied by pure-state transitivity alone.',
            'not_derived': ['particular K', 'physical time or energy scale', 'hbar', 'availability of all mathematical unitaries', 'all CPTP instruments', 'standard tensor identification', 'spacetime or gravity'],
            'dimensions': [dimension_certificate(d) for d in (2, 3, 4)],
            'counterexamples': {'clock': clock_witness(), 'normalized_filter': filter_witness()},
            'runtime': {'python': platform.python_version(), 'numpy': np.__version__},
            'evidence_scope': 'Finite-dimensional numerical certificates and exact analytic witnesses; the all-dimension classification is cited, not inferred from sampling.',
            'classification_source': {'url': 'https://arxiv.org/html/1911.06635', 'location': 'Section 1, Proposition 1.2 and Theorem 1.1'},
            'next': 'Audit compatibility of reconstructed composites with the standard complex tensor product, then link actual process extensions to complete positivity; separately audit process saturation/purification.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(BridgeTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = report()
    data['checks'] = {'run': checks.testsRun, 'failures': 0, 'errors': 0}
    if args.write_results:
        target = Path(__file__).with_name('reversible_dynamics_bridge_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8')) != data:
            raise RuntimeError('Existing scientific result differs; review before replacing it.')
        target.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
