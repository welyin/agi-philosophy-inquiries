"""Round 160: analytic first and second derivatives of complete recovery error.

Newton iterations discover constrained stationary points. Floating residuals
and Hessian spectra are diagnostics; the separate rational certificates give
the rigorous statements about nonconcavity and global error bounds.
"""

import argparse
from functools import lru_cache
import json
from pathlib import Path
import unittest

import numpy as np

from joint_history_channel import (I2, X, Y, Z, BIT_ERROR, weighted_gram,
    source_sigma, bloch_state, error_from_sigma, grouped_weights)
from joint_history_certificate import (promised_bloch_vectors, dual_data,
    exact_product_sigma, rayleigh_lower, LOWER_INTEGERS)
from joint_relation_records import interval_fields


@lru_cache(None)
def gram_basis():
    basis = (I2, X, Y, Z)
    return np.array([[weighted_gram(np.kron(a, b.conj()).real/4)
                      for b in basis] for a in basis])


def gram_jet(point):
    point = np.asarray(point, dtype=float)
    if point.shape != (6,):
        raise ValueError('Use the two three-coordinate Bloch vectors.')
    u, v = np.r_[1., point[:3]], np.r_[1., point[3:]]
    basis = gram_basis()
    value = np.einsum('i,j,ijab->ab', u, v, basis)
    first = np.concatenate((np.einsum('j,ijab->iab', v, basis[1:]),
                            np.einsum('i,ijab->jab', u, basis[:, 1:])))
    second = np.zeros((6, 6, 19, 19))
    second[:3, 3:] = basis[1:, 1:]
    second[3:, :3] = basis[1:, 1:].transpose(1, 0, 2, 3)
    return value, first, second


def error_jet(point, second_order=True):
    """Solve the positive Schur root, then differentiate implicitly.

    G=[[M,g],[g^T,1]], c=(M+tI)^-1 g, z=(-c,1), d=1-c^T c.
    t_i=z^T G_i z/d; t_ij=(z^T G_ij z-2 q_i^T A^-1 q_j)/d.
    """
    gram, first, second = gram_jet(point)
    m, g = gram[:-1, :-1], gram[:-1, -1]
    t = 1.
    for _ in range(30):
        a = m+t*np.eye(18)
        c = np.linalg.solve(a, g)
        d = 1-c @ c
        residual = 1-t-g @ c
        if d <= 0:
            raise ArithmeticError('The positive-root derivative denominator failed.')
        step = residual/d
        t += step
        if abs(step) < 2e-15:
            break
    else:
        raise ArithmeticError('Positive Schur root did not converge.')
    a = m+t*np.eye(18)
    c = np.linalg.solve(a, g)
    d = 1-c @ c
    z = np.r_[-c, 1.]
    gradient = np.einsum('a,iab,b->i', z, first, z)/d
    result = {'value': float(t), 'gradient': gradient, 'denominator': float(d),
              'coefficients': c*np.sqrt(np.array(grouped_weights(BIT_ERROR), dtype=float)),
              'schur_residual': float(1-t-g @ c)}
    if second_order:
        q = np.array([first[i, :-1, -1] - first[i, :-1, :-1] @ c - gradient[i]*c
                      for i in range(6)]).T
        result['hessian'] = (np.einsum('a,ijab,b->ij', z, second, z)
                             - 2*q.T @ np.linalg.solve(a, q))/d
    return result


def stationary_candidate():
    """Eight-variable Newton system for two active unit-ball constraints."""
    point = np.array(promised_bloch_vectors(), dtype=float).reshape(6)
    jet = error_jet(point)
    multipliers = np.array([jet['gradient'][:3] @ point[:3],
                            jet['gradient'][3:] @ point[3:]])
    for _ in range(8):
        jet = error_jet(point)
        normals = np.zeros((6, 2))
        normals[:3, 0], normals[3:, 1] = point[:3], point[3:]
        residual = np.r_[jet['gradient']-normals @ multipliers,
                         point[:3] @ point[:3]-1, point[3:] @ point[3:]-1]
        jacobian = np.block([[jet['hessian']-np.diag(np.repeat(multipliers, 3)), -normals],
                             [2*normals.T, np.zeros((2, 2))]])
        change = np.linalg.solve(jacobian, -residual)
        point += change[:6]
        multipliers += change[6:]
        if np.max(np.abs(residual)) < 1e-15:
            break
    jet = error_jet(point)
    tangent = np.zeros((6, 4))
    tangent[:3, :2] = np.linalg.svd(point[:3][None, :])[2][1:].T
    tangent[3:, 2:] = np.linalg.svd(point[3:][None, :])[2][1:].T
    lagrangian = jet['hessian']-np.diag(np.repeat(multipliers, 3))
    stationarity = jet['gradient']-np.repeat(multipliers, 3)*point
    return {'point': point.tolist(), 'error_diagnostic': jet['value'],
            'radial_multipliers_diagnostic': multipliers.tolist(),
            'stationarity_residual_diagnostic': float(np.max(np.abs(stationarity))),
            'constraint_residual_diagnostic': float(max(abs(point[:3] @ point[:3]-1),
                                                       abs(point[3:] @ point[3:]-1))),
            'tangent_lagrangian_eigenvalues_diagnostic': np.linalg.eigvalsh(tangent.T @ lagrangian @ tangent).tolist(),
            'full_hessian_eigenvalues_diagnostic': np.linalg.eigvalsh(jet['hessian']).tolist(),
            'stationary_point_existence_or_uniqueness_interval_certified': False}


def nonconcavity_certificate():
    """Exact Jensen counterexample, not an inference from a floating Hessian."""
    a, b = promised_bloch_vectors()
    original = exact_product_sigma(a, b)
    reflected = exact_product_sigma((a[0], -a[1], a[2]), (b[0], -b[1], b[2]))
    if not np.array_equal(original, reflected):
        raise AssertionError('Simultaneous conjugation symmetry failed.')
    midpoint = exact_product_sigma((a[0], a[1]*0, a[2]), (b[0], b[1]*0, b[2]))
    mid_coefficients = (510771285238, -158286489, 14035806177, 8713908445,
        19869529747, 7429917933, 514960878, 24005372170, 7574862422,
        19953686735, 10642695692, 27998001443, 17329930864, 6209467127,
        242304755, 3657859302, -124564716, 13794292215)
    d, numerator = dual_data(mid_coefficients)
    mid_upper = sum(midpoint[j, i]*numerator[i, j] for i in range(4) for j in range(4))/(2*d)
    endpoint_lower = rayleigh_lower(original, LOWER_INTEGERS)
    gap = endpoint_lower-mid_upper
    if gap.lo <= 0:
        raise AssertionError('The strict Jensen counterexample did not certify.')
    return {'midpoint_upper': interval_fields(mid_upper),
            'equal_endpoint_lower': interval_fields(endpoint_lower),
            'strict_Jensen_gap': interval_fields(gap),
            'midpoint_supporting_coefficient_numerators': mid_coefficients,
            'both_endpoints_and_midpoint_are_legal_independent_sources': True}


class DifferentialHistoryRecoveryTests(unittest.TestCase):
    def test_bilinear_gram_and_exact_derivative_structure(self):
        x = np.array([.2, -.3, .1, -.4, .2, .3])
        g, first, second = gram_jet(x)
        expected = weighted_gram(source_sigma(bloch_state(x[:3]), bloch_state(x[3:])))
        np.testing.assert_allclose(g, expected, atol=6e-17)
        delta = np.arange(6)*.013
        predicted = g+np.einsum('i,iab->ab', delta, first)+np.einsum('i,j,ijab->ab', delta, delta, second)/2
        np.testing.assert_allclose(gram_jet(x+delta)[0], predicted, atol=1e-16)

    def test_schur_root_matches_independent_gram_spectrum(self):
        rng = np.random.default_rng(160)
        for _ in range(12):
            x = rng.normal(size=(2, 3))
            x /= np.maximum(1, np.linalg.norm(x, axis=1))[:, None]
            jet = error_jet(x.reshape(6))
            self.assertAlmostEqual(jet['value'], error_from_sigma(source_sigma(bloch_state(x[0]), bloch_state(x[1]))), places=13)
            self.assertLess(abs(jet['schur_residual']), 5e-16)

    def test_analytic_gradient_against_centered_value_differences(self):
        x = np.array([.2, -.3, .1, -.4, .2, .3]); eps = 1e-5
        finite = [(error_jet(x+eps*v)['value']-error_jet(x-eps*v)['value'])/(2*eps) for v in np.eye(6)]
        np.testing.assert_allclose(error_jet(x)['gradient'], finite, atol=2e-10)

    def test_analytic_hessian_against_centered_gradient_differences(self):
        x = np.array([.2, -.3, .1, -.4, .2, .3]); eps = 1e-5
        finite = np.array([(error_jet(x+eps*v)['gradient']-error_jet(x-eps*v)['gradient'])/(2*eps) for v in np.eye(6)]).T
        np.testing.assert_allclose(error_jet(x)['hessian'], finite, atol=2e-10)
        np.testing.assert_allclose(error_jet(x)['hessian'], error_jet(x)['hessian'].T, atol=1e-17)

    def test_separate_concavity_and_joint_nonconcavity_diagnostics(self):
        x = np.array(promised_bloch_vectors(), dtype=float).reshape(6)
        h = error_jet(.8*x)['hessian']
        self.assertLessEqual(np.linalg.eigvalsh(h[:3, :3]).max(), 1e-15)
        self.assertLessEqual(np.linalg.eigvalsh(h[3:, 3:]).max(), 1e-15)
        self.assertGreater(np.linalg.eigvalsh(h).max(), 1e-5)

    def test_newton_stationarity_and_boundary_curvature_are_not_global_proof(self):
        report = stationary_candidate()
        self.assertLess(report['stationarity_residual_diagnostic'], 1e-14)
        self.assertLess(report['constraint_residual_diagnostic'], 1e-14)
        self.assertGreater(min(report['radial_multipliers_diagnostic']), 0)
        self.assertLess(max(report['tangent_lagrangian_eigenvalues_diagnostic']), -1e-6)
        self.assertFalse(report['stationary_point_existence_or_uniqueness_interval_certified'])

    def test_joint_nonconcavity_has_a_strict_interval_Jensen_counterexample(self):
        result = nonconcavity_certificate()
        self.assertGreater(result['strict_Jensen_gap']['diagnostic'], .000017)
        self.assertTrue(result['both_endpoints_and_midpoint_are_legal_independent_sources'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(DifferentialHistoryRecoveryTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {'round': 160, 'candidate': stationary_candidate(),
              'joint_nonconcavity_certificate': nonconcavity_certificate(),
              'analytic_gradient': 'z^T G_i z / d',
              'analytic_hessian': '(z^T G_ij z - 2 q_i^T A^-1 q_j) / d',
              'separate_concavity_follows_from_zero_same_subject_G_ij': True,
              'automated_checks': {'run': checks.testsRun, 'failures': len(checks.failures), 'errors': len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name('differential_history_recovery_results.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
