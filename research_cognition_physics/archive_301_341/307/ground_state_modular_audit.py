"""Round 307: derive a reduced Gaussian state from an explicit fermion ground state."""
import itertools
import unittest
import numpy as np
from growing_stream_audit import main


def chain_hamiltonian(length, antiperiodic=False):
    if length < 2 or length % 2:
        raise ValueError('An even chain is required to avoid the open-chain zero mode.')
    h = np.diag(np.full(length-1, -0.5), 1)
    h += h.T
    if antiperiodic:
        if length % 4:
            raise ValueError('Use a multiple of four for the antiperiodic half-filled ring.')
        h[0, -1] = h[-1, 0] = 0.5
    return h


def ground_projector(h):
    energies, modes = np.linalg.eigh(h)
    if np.min(abs(energies)) < 1e-12:
        raise ValueError('A zero mode needs a separate occupation prescription.')
    occupied = modes[:, energies < 0]
    return occupied@occupied.T, occupied, float(energies[energies < 0].sum())


def interval_correlation(n, ring_length=None):
    delta = np.arange(n)[:, None]-np.arange(n)
    if ring_length is None:
        # Exact infinite-chain Fourier projector: sin(pi*r/2)/(pi*r).
        return 0.5*np.sinc(delta/2)
    if ring_length % 4 or ring_length <= n:
        raise ValueError('The antiperiodic ring must be larger than the interval.')
    q = 2*np.pi*(np.arange(-ring_length//2, ring_length//2)+0.5)/ring_length
    occupied = q[abs(q) < np.pi/2]
    return np.cos(delta[:, :, None]*occupied).sum(axis=2)/ring_length


def modular_matrix(c, margin=1e-12):
    if not np.allclose(c, c.T, atol=1e-13, rtol=0):
        raise ValueError('This calibration implements real symmetric correlations.')
    z, u = np.linalg.eigh(c)
    if z[0] <= margin or z[-1] >= 1-margin:
        raise ValueError('Unresolved spectral endpoint; full logarithm refused, not clipped.')
    return (u*np.log((1-z)/z))@u.T


def second_quantize(h):
    n = len(h)
    result = np.zeros((2**n, 2**n))
    for ket in range(2**n):
        for j in range(n):
            if not (ket >> j) & 1:
                continue
            after = ket ^ (1 << j)
            sign_j = (-1)**((ket & ((1 << j)-1)).bit_count())
            for i in range(n):
                if (after >> i) & 1:
                    continue
                sign_i = (-1)**((after & ((1 << i)-1)).bit_count())
                result[after | (1 << i), ket] += h[i, j]*sign_i*sign_j
    return result


def slater_state(occupied):
    n, particles = occupied.shape
    psi = np.zeros(2**n)
    for sites in itertools.combinations(range(n), particles):
        psi[sum(1 << i for i in sites)] = np.linalg.det(occupied[list(sites), :])
    return psi


def reduced_prefix(psi, n):
    amplitudes = psi.reshape(-1, 2**n).T
    return amplitudes@amplitudes.T


def gaussian_density(k):
    values, vectors = np.linalg.eigh(second_quantize(k))
    weights = np.exp(-(values-values.min()))
    return (vectors*(weights/weights.sum()))@vectors.T


def entropy(rho):
    p = np.linalg.eigvalsh(rho)
    p = p[p > 0]
    return float(-p@np.log(p))


def small_case():
    h = chain_hamiltonian(8)
    c, occupied, energy = ground_projector(h)
    psi = slater_state(occupied)
    rho = reduced_prefix(psi, 3)
    k = modular_matrix(c[:3, :3])
    return h, c, psi, energy, rho, k


def report():
    h, c, psi, energy, rho, k = small_case()
    fock_k = second_quantize(k)
    p, u = np.linalg.eigh(rho)
    exact_k = (u*(-np.log(p)))@u.T
    difference = exact_k-fock_k
    scalar = np.trace(difference)/len(difference)
    z = np.linalg.eigvalsh(c[:3, :3])
    cinf = interval_correlation(8)
    kinf = modular_matrix(cinf)
    rows = []
    for length in (32, 64, 128, 256):
        cring = interval_correlation(8, length)
        rows.append({'ring_sites': length, 'interval_sites': 8,
                     'correlation_error_fro': float(np.linalg.norm(cring-cinf)),
                     'modular_error_op': float(np.linalg.norm(modular_matrix(cring)-kinf, 2))})
    uniform = chain_hamiltonian(8)
    fitted_beta = float(np.sum(kinf*uniform)/np.sum(uniform**2))
    return {'round': 307,
            'scope': 'Explicit 1D free-fermion branch; ground-state partial trace, not a derived spacetime or physical vacuum selection principle.',
            'inputs': ['canonical anticommutation', 'nearest-neighbour hopping -1/2',
                       'half-filled ground state', 'specified one-dimensional region and boundary'],
            'small_fock_calibration': {'chain_sites': 8, 'prefix_sites': 3,
                'ground_energy': energy, 'norm_error': float(abs(psi@psi-1)),
                'ground_equation_error': float(np.linalg.norm(second_quantize(h)@psi-energy*psi)),
                'partial_trace_vs_gaussian_fro': float(np.linalg.norm(rho-gaussian_density(k))),
                'modular_scalar': float(scalar),
                'centered_modular_error_fro': float(np.linalg.norm(difference-scalar*np.eye(len(rho)))),
                'entropy_fock': entropy(rho),
                'entropy_correlation': float(-z@np.log(z)-(1-z)@np.log(1-z))},
            'infinite_interval_counterexample': {'sites': 8,
                'long_range_k_1_4': float(kinf[0, 3]), 'uniform_best_beta': fitted_beta,
                'uniform_generator_residual_fro': float(np.linalg.norm(kinf-fitted_beta*uniform))},
            'finite_environment_convergence': rows,
            'correlation_error_halving_ratios': [rows[i]['correlation_error_fro']/rows[i+1]['correlation_error_fro'] for i in range(3)],
            'full_log_endpoint_policy': 'Reject endpoints within 1e-12; no clipping or claimed full matrix beyond resolvable sizes.'}


class Checks(unittest.TestCase):
    def test_01_ground_projector_and_gap(self):
        h = chain_hamiltonian(8)
        c, _, _ = ground_projector(h)
        np.testing.assert_allclose(c@c, c, atol=1e-14)
        self.assertAlmostEqual(np.trace(c), 4)
        with self.assertRaises(ValueError):
            ground_projector(np.zeros((2, 2)))

    def test_02_slater_ground_equation(self):
        h, _, psi, energy, _, _ = small_case()
        self.assertAlmostEqual(psi@psi, 1)
        np.testing.assert_allclose(second_quantize(h)@psi, energy*psi, atol=2e-14)

    def test_03_explicit_correlations(self):
        _, c, _, _, rho, _ = small_case()
        for i in range(3):
            for j in range(3):
                elementary = np.zeros((3, 3)); elementary[i, j] = 1
                self.assertAlmostEqual(np.trace(rho@second_quantize(elementary)), c[i, j])

    def test_04_independent_partial_trace(self):
        _, _, _, _, rho, k = small_case()
        np.testing.assert_allclose(rho, gaussian_density(k), atol=3e-14)
        self.assertAlmostEqual(np.trace(rho), 1)

    def test_05_entropy_and_modular_constant(self):
        data = report()['small_fock_calibration']
        self.assertAlmostEqual(data['entropy_fock'], data['entropy_correlation'], places=12)
        self.assertLess(data['centered_modular_error_fro'], 1e-8)

    def test_06_ring_projector_independent(self):
        c, _, _ = ground_projector(chain_hamiltonian(32, antiperiodic=True))
        np.testing.assert_allclose(c[:8, :8], interval_correlation(8, 32), atol=2e-14)

    def test_07_infinite_environment_limit(self):
        data = report()
        for ratio in data['correlation_error_halving_ratios']:
            self.assertGreater(ratio, 3.9); self.assertLess(ratio, 4.2)
        errors = [r['modular_error_op'] for r in data['finite_environment_convergence']]
        self.assertTrue(all(a > b for a, b in zip(errors, errors[1:])))

    def test_08_not_uniform_physical_generator(self):
        data = report()['infinite_interval_counterexample']
        self.assertGreater(abs(data['long_range_k_1_4']), 0.03)
        self.assertGreater(data['uniform_generator_residual_fro'], 1)

    def test_09_spectral_endpoints_rejected(self):
        with self.assertRaises(ValueError):
            modular_matrix(np.diag([0., 1.]))
        with self.assertRaises(ValueError):
            modular_matrix(interval_correlation(32))

    def test_10_half_filling_parity(self):
        k = modular_matrix(interval_correlation(8))
        parity = np.diag((-1.)**np.arange(8))
        np.testing.assert_allclose(parity@k@parity, -k, atol=1e-10)


if __name__ == '__main__':
    main(__name__, 'ground_state_modular_audit', report)
