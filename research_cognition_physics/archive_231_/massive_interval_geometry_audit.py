"""Round 310: test interval geometries through states, without inverting spectral endpoints."""
from functools import lru_cache
import unittest
import numpy as np
from numpy.polynomial.legendre import leggauss
from growing_stream_audit import main
from gapped_ground_state_audit import correlation, physical_matrix, entropy_from_c, correlation_length_sites
from ground_state_modular_audit import gaussian_density, modular_matrix


@lru_cache(maxsize=32)
def elliptic_k(modulus, order=128):
    if not 0 <= modulus < 1:
        raise ValueError('Complete elliptic integral requires 0 <= modulus < 1.')
    nodes, weights = leggauss(order)
    theta = (nodes+1)*np.pi/4
    return float(np.dot(weights, 1/np.sqrt(1-modulus**2*np.sin(theta)**2))*np.pi/4)


def ctm_coefficient(delta):
    d = abs(delta)
    if not d < 1:
        raise ValueError('The singular fully dimerized modular generator is excluded.')
    ratio = (1-d)/(1+d)
    return 4*elliptic_k(np.sqrt(1-ratio**2))/(1+d)


def candidate(n, delta, shape, normalization='ctm'):
    i = np.arange(1, n)
    if shape == 'parabola':
        weight = i*(n-i)/n
    elif shape == 'triangle':
        weight = np.minimum(i, n-i)
    elif shape == 'halfline':
        weight = i
    else:
        raise ValueError('Unknown region weight.')
    coefficient = ctm_coefficient(delta) if normalization == 'ctm' else 2*np.pi
    bonds = np.diag(physical_matrix(n, delta), 1)*weight*coefficient
    return np.diag(bonds, 1)+np.diag(bonds, -1)


def covariance_of_generator(k):
    energies, modes = np.linalg.eigh(k)
    weights = np.exp(-np.logaddexp(0, energies))
    return (modes*weights)@modes.T


def distinguishability(c, k):
    predicted = covariance_of_generator(k)
    entropy = entropy_from_c(c)
    eigen_k = np.linalg.eigvalsh(k)
    relative = float(np.sum(c*k.T)+np.sum(np.logaddexp(0, -eigen_k))-entropy)
    lower = float(np.linalg.norm(c-predicted, 2))
    if relative < -1e-8:
        raise ArithmeticError('Relative entropy outside floating point tolerance.')
    return {'covariance_error_fro': float(np.linalg.norm(c-predicted)),
            'trace_distance_lower_from_number_observable': lower,
            'relative_entropy': relative,
            'pinsker_trace_distance_upper_numeric': float(np.sqrt(max(relative, 0)/2)),
            'entropy_error': entropy_from_c(predicted)-entropy,
            'stationarity_commutator_fro': float(np.linalg.norm(c@k-k@c))}


@lru_cache(maxsize=1)
def report():
    rows = []
    for d in (0., 0.02, 0.1, 0.3, -0.3):
        for n in (16, 32, 64, 128):
            c = correlation(n, d)
            rows.append({'delta': d, 'sites': n,
                         'size_over_correlation_length': 0. if d == 0 else n/correlation_length_sites(d),
                         'parabola_critical': distinguishability(c, candidate(n, d, 'parabola', 'critical')),
                         'parabola_ctm': distinguishability(c, candidate(n, d, 'parabola')),
                         'triangle_ctm': distinguishability(c, candidate(n, d, 'triangle'))})
    return {'round': 310,
            'scope': 'Finite gapped intervals: compare stated geometric ansatz states. Covariance bounds concern Fock trace distance, not full modular-generator norm.',
            'cases': rows,
            'numerical_policy': 'No logarithm of exact C at unresolved endpoints. Relative entropy is evaluated via the specified candidate generator; cancellation is checked on a small independent Fock example.',
            'normalizations': {'critical': 2*np.pi, 'ctm_delta_0_3': ctm_coefficient(0.3)},
            'conclusion': 'Mass scale, region shape, and bond cut matter. Entropy agreement alone is insufficient to identify the modular generator.'}


class Checks(unittest.TestCase):
    def test_01_critical_normalization(self):
        self.assertAlmostEqual(ctm_coefficient(0), 2*np.pi, places=13)

    def test_02_elliptic_quadrature(self):
        for k in (0.1, 0.5, 0.9):
            self.assertAlmostEqual(elliptic_k(k, 64), elliptic_k(k, 128), places=12)

    def test_03_gaussian_covariance_independent(self):
        k = candidate(4, 0.2, 'triangle')
        rho = gaussian_density(k)
        from ground_state_modular_audit import second_quantize
        c = covariance_of_generator(k)
        for i in range(4):
            for j in range(4):
                h = np.zeros((4, 4)); h[i, j] = 1
                self.assertAlmostEqual(np.trace(rho@second_quantize(h)), c[i, j], places=12)

    def test_04_relative_entropy_independent_fock(self):
        c = correlation(4, 0.2); k = candidate(4, 0.2, 'triangle')
        rho = gaussian_density(modular_matrix(c)); sigma = gaussian_density(k)
        def logm(a):
            z, u = np.linalg.eigh(a)
            return (u*np.log(z))@u.T
        direct = float(np.trace(rho@(logm(rho)-logm(sigma))))
        self.assertAlmostEqual(direct, distinguishability(c, k)['relative_entropy'], places=9)

    def test_05_trace_distance_bracket_small_fock(self):
        c = correlation(4, 0.2); k = candidate(4, 0.2, 'parabola')
        rho, sigma = gaussian_density(modular_matrix(c)), gaussian_density(k)
        distance = np.sum(abs(np.linalg.eigvalsh(rho-sigma)))/2
        b = distinguishability(c, k)
        self.assertLessEqual(b['trace_distance_lower_from_number_observable'], distance+1e-12)
        self.assertLessEqual(distance, b['pinsker_trace_distance_upper_numeric']+1e-10)

    def test_06_massless_and_massive_shapes_differ(self):
        c0, cm = correlation(32, 0), correlation(32, 0.3)
        p0 = distinguishability(c0, candidate(32, 0, 'parabola'))['covariance_error_fro']
        t0 = distinguishability(c0, candidate(32, 0, 'triangle'))['covariance_error_fro']
        pm = distinguishability(cm, candidate(32, 0.3, 'parabola'))['covariance_error_fro']
        tm = distinguishability(cm, candidate(32, 0.3, 'triangle'))['covariance_error_fro']
        self.assertLess(p0, t0/100); self.assertLess(tm, pm)

    def test_07_spectral_norm_witness(self):
        c = correlation(32, 0.02); k = candidate(32, 0.02, 'parabola')
        delta = c-covariance_of_generator(k)
        values, vectors = np.linalg.eigh(delta)
        f = vectors[:, np.argmax(abs(values))]
        self.assertAlmostEqual(abs(f@delta@f), distinguishability(c, k)['trace_distance_lower_from_number_observable'])

    def test_08_no_false_global_claim_from_entropy(self):
        b = distinguishability(correlation(64, 0.3), candidate(64, 0.3, 'triangle'))
        self.assertLess(abs(b['entropy_error']), 1e-6)
        self.assertGreater(b['covariance_error_fro'], 1e-3)

    def test_09_same_gap_cut_changes_shape_ranking(self):
        c = correlation(32, -0.3)
        p = distinguishability(c, candidate(32, -0.3, 'parabola'))['covariance_error_fro']
        t = distinguishability(c, candidate(32, -0.3, 'triangle'))['covariance_error_fro']
        self.assertLess(p, t)


if __name__ == '__main__':
    main(__name__, 'massive_interval_geometry_audit', report)
