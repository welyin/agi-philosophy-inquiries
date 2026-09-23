"""Round 305: finite quantum entropy first law and physical-generator calibration."""
import unittest
import numpy as np
from growing_stream_audit import main


def density_spectrum(rho):
    rho = np.asarray(rho, dtype=complex)
    if rho.ndim != 2 or rho.shape[0] != rho.shape[1] or not np.allclose(rho, rho.conj().T, atol=1e-12):
        raise ValueError('Hermitian square density matrix required.')
    if abs(np.trace(rho)-1) > 1e-12:
        raise ValueError('Unit trace required.')
    values, vectors = np.linalg.eigh(rho)
    if values[0] < -1e-12:
        raise ValueError('Positive semidefinite state required.')
    return np.maximum(values, 0.), vectors


def entropy(rho):
    values, _ = density_spectrum(rho)
    positive = values[values > 0]
    return float(-np.dot(positive, np.log(positive)))


def modular_hamiltonian(sigma):
    values, vectors = density_spectrum(sigma)
    if values[0] <= 0:
        raise ValueError('Full-rank reference required; support restriction is not implemented.')
    return (vectors*(-np.log(values)))@vectors.conj().T


def gibbs(hamiltonian, beta=1.):
    values, vectors = np.linalg.eigh(hamiltonian)
    weights = np.exp(-beta*(values-values[0]))
    weights /= weights.sum()
    return (vectors*weights)@vectors.conj().T


def budget(rho, sigma):
    modular = modular_hamiltonian(sigma)
    delta_s = entropy(rho)-entropy(sigma)
    delta_k = float(np.trace((rho-sigma)@modular).real)
    # Independent direct definition: Tr rho log rho - Tr rho log sigma.
    relative = float(-entropy(rho)+np.trace(rho@modular).real)
    return {'delta_S': delta_s, 'delta_K': delta_k, 'relative_entropy': relative,
            'identity_error': abs(delta_k-delta_s-relative)}


def centered(matrix):
    return matrix-np.trace(matrix)*np.eye(len(matrix))/len(matrix)


def fit_physical_generator(sigma, hamiltonian):
    kc, hc = centered(modular_hamiltonian(sigma)), centered(hamiltonian)
    denominator = float(np.trace(hc@hc).real)
    if denominator <= 1e-20:
        raise ValueError('A constant Hamiltonian cannot determine an inverse temperature.')
    beta = float(np.trace(kc@hc).real)/denominator
    residual = kc-beta*hc
    return beta, float(np.linalg.norm(residual)), residual


def tangent_hessian(sigma, delta):
    values, vectors = density_spectrum(sigma)
    transformed = vectors.conj().T@delta@vectors
    coefficient = np.empty((len(values), len(values)))
    for i, pi in enumerate(values):
        for j, pj in enumerate(values):
            coefficient[i, j] = 1/pi if abs(pi-pj) < 1e-14 else (np.log(pi)-np.log(pj))/(pi-pj)
    return float(np.sum(abs(transformed)**2*coefficient))


def example():
    hamiltonian = np.diag([0., 1., 2., 4.])
    sigma = gibbs(hamiltonian, 0.7)
    rng = np.random.default_rng(305)
    raw = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
    delta = centered((raw+raw.conj().T)/2)
    delta *= 0.4*np.linalg.eigvalsh(sigma)[0]/np.linalg.norm(delta, 2)
    return hamiltonian, sigma, delta


def unitary(generator, parameter):
    values, vectors = np.linalg.eigh(generator)
    return (vectors*np.exp(-1j*parameter*values))@vectors.conj().T


def report():
    hamiltonian, sigma, delta = example()
    coefficient = tangent_hessian(sigma, delta)/2
    rows = []
    for epsilon in (0.8, 0.4, 0.2, 0.1):
        data = budget(sigma+epsilon*delta, sigma)
        mu = np.linalg.eigvalsh(sigma)[0]-abs(epsilon)*np.linalg.norm(delta, 2)
        upper = epsilon**2*np.linalg.norm(delta)**2/(2*mu)
        rows.append({'epsilon': epsilon, **data, 'D_over_epsilon_squared': data['relative_entropy']/epsilon**2,
                     'rigorous_upper_bound': float(upper)})
    nongibbs = np.diag([0.45, 0.30, 0.20, 0.05])
    bad_beta, bad_error, _ = fit_physical_generator(nongibbs, hamiltonian)
    beta, error, _ = fit_physical_generator(sigma, hamiltonian)
    x = np.zeros((4, 4)); x[0, 1] = x[1, 0] = 1.
    restricted = gibbs(0.7*hamiltonian+0.3*x)
    mismatch = centered(modular_hamiltonian(restricted)-0.7*hamiltonian)
    return {'round': 305,
            'scope': 'Finite full-rank quantum reference and specified perturbations; entropy identity and generator calibration, without a derived physical boost, local QFT, area entropy or Unruh temperature.',
            'reference_min_eigenvalue': float(np.linalg.eigvalsh(sigma)[0]),
            'noncommuting_tangent_norm': float(np.linalg.norm(sigma@delta-delta@sigma)),
            'relative_entropy_rows': rows, 'quadratic_limit_coefficient': coefficient,
            'halving_D_ratios': [rows[i]['relative_entropy']/rows[i+1]['relative_entropy'] for i in range(3)],
            'physical_generator_fit': {'thermal_beta': beta, 'thermal_residual': error,
                                       'stationary_nongibbs_best_beta': bad_beta,
                                       'stationary_nongibbs_residual': bad_error},
            'restricted_tangent_counterexample': {'diagonal_residual_norm': float(np.linalg.norm(np.diag(mismatch))),
                                                   'full_residual_norm': float(np.linalg.norm(mismatch)),
                                                   'offdiagonal_witness': float(np.trace(x@mismatch).real)}}


class ModularTests(unittest.TestCase):
    def test_01_exact_relative_entropy_identity_and_positivity(self):
        _, sigma, delta = example()
        for epsilon in (-0.8, 0.2, 0.8):
            row = budget(sigma+epsilon*delta, sigma)
            self.assertLess(row['identity_error'], 1e-13)
            self.assertGreater(row['relative_entropy'], 0.)
        pure = np.diag([1., 0., 0., 0.])
        self.assertLess(budget(pure, sigma)['identity_error'], 1e-13)

    def test_02_first_law_and_noncommuting_quadratic_limit(self):
        _, sigma, delta = example()
        self.assertGreater(np.linalg.norm(sigma@delta-delta@sigma), 1e-4)
        expected = tangent_hessian(sigma, delta)/2
        fine = budget(sigma+0.01*delta, sigma)['relative_entropy']/0.01**2
        self.assertLess(abs(fine/expected-1), 0.002)
        for ratio in report()['halving_D_ratios']:
            self.assertTrue(3.7 < ratio < 4.3)

    def test_03_finite_remainder_bounds(self):
        _, sigma, delta = example()
        for epsilon in (-0.8, 0.4, 0.8):
            relative = budget(sigma+epsilon*delta, sigma)['relative_entropy']
            mu = np.linalg.eigvalsh(sigma)[0]-abs(epsilon)*np.linalg.norm(delta, 2)
            upper = epsilon**2*np.linalg.norm(delta)**2/(2*mu)
            lower = np.sum(abs(np.linalg.eigvalsh(epsilon*delta)))**2/2
            self.assertLessEqual(relative, upper+1e-13)
            self.assertGreaterEqual(relative+1e-13, lower)

    def test_04_thermal_generator_and_flow_agree(self):
        h, sigma, _ = example()
        beta, error, _ = fit_physical_generator(sigma, h)
        self.assertAlmostEqual(beta, 0.7)
        self.assertLess(error, 1e-13)
        observable = np.ones((4, 4))
        modular_u = unitary(modular_hamiltonian(sigma), 0.37)
        physical_u = unitary(h, beta*0.37)
        np.testing.assert_allclose(modular_u@observable@modular_u.conj().T,
                                   physical_u@observable@physical_u.conj().T, atol=1e-13)

    def test_05_stationarity_is_not_thermal_calibration(self):
        h, _, _ = example()
        sigma = np.diag([0.45, 0.30, 0.20, 0.05])
        self.assertEqual(np.linalg.norm(h@sigma-sigma@h), 0.)
        self.assertGreater(fit_physical_generator(sigma, h)[1], 0.1)
        transition_betas = [np.log(sigma[i,i]/sigma[j,j])/(h[j,j]-h[i,i]) for i,j in ((0,1),(1,2),(2,3))]
        self.assertGreater(max(transition_betas)-min(transition_betas), 0.2)

    def test_06_restricted_preparations_miss_operator_mismatch(self):
        data = report()['restricted_tangent_counterexample']
        self.assertLess(data['diagonal_residual_norm'], 1e-13)
        self.assertGreater(data['full_residual_norm'], 0.4)
        self.assertAlmostEqual(data['offdiagonal_witness'], 0.6)

    def test_07_support_and_temperature_degeneracies_are_explicit(self):
        with self.assertRaises(ValueError):
            modular_hamiltonian(np.diag([1., 0.]))
        with self.assertRaises(ValueError):
            fit_physical_generator(np.eye(2)/2, np.eye(2))
        with self.assertRaises(ValueError):
            entropy(np.diag([1.1, -0.1]))

    def test_08_unitary_covariance_and_additive_generator_constant(self):
        h, sigma, delta = example()
        rotation = unitary(np.ones((4,4)), 0.2)
        before = budget(sigma+0.4*delta, sigma)
        after = budget(rotation@(sigma+0.4*delta)@rotation.conj().T, rotation@sigma@rotation.conj().T)
        self.assertAlmostEqual(before['relative_entropy'], after['relative_entropy'])
        self.assertAlmostEqual(fit_physical_generator(sigma,h+8*np.eye(4))[0], 0.7)

    def test_09_zero_first_order_energy_can_have_finite_entropy_change(self):
        h, sigma, _ = example()
        delta = np.zeros((4,4)); delta[0,1]=delta[1,0]=0.01
        row = budget(sigma+delta, sigma)
        self.assertAlmostEqual(row['delta_K'], 0.)
        self.assertLess(row['delta_S'], 0.)
        self.assertAlmostEqual(row['delta_S'], -row['relative_entropy'])


if __name__ == '__main__':
    main(__name__, 'modular_energy_interface_audit', report)
