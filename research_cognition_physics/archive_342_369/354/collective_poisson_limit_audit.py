"""Round 354: controlled collective-spin Poisson limit, not emergent geometry.

Specified symmetric coherent states and H=chi*Jz**2/N produce a spherical
classical flow for normalized collective observables on fixed time intervals.
The construction does not select the Hamiltonian, spatial locality, a metric,
or hypersurface constraints. Full-state and long-time limits are different.
"""
from math import lgamma
import unittest
import numpy as np
from growing_stream_audit import main


X = np.array([[0., 1.], [1., 0.]], dtype=complex)
Y = np.array([[0., -1j], [1j, 0.]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)


def spin_matrices(n):
    """Spin j=N/2 representation; m runs from -j to +j."""
    j = n/2
    m = np.arange(n+1)-j
    raising = np.zeros((n+1, n+1), dtype=complex)
    for k, value in enumerate(m[:-1]):
        raising[k+1, k] = np.sqrt((j-value)*(j+value+1))
    return ((raising+raising.conj().T)/2,
            (raising-raising.conj().T)/(2j), np.diag(m))


def coherent(n, z, phi=0.):
    """Symmetric-sector coefficients of identical pure qubits."""
    if n < 1 or not -1 <= z <= 1:
        raise ValueError("Need N>=1 and -1<=z<=1.")
    out = np.zeros(n+1, dtype=complex)
    if z == 1:
        out[-1] = 1.
    elif z == -1:
        out[0] = np.exp(1j*n*phi)
    else:
        k = np.arange(n+1)
        log_binomial = np.array([lgamma(n+1)-lgamma(int(a)+1) -
                                lgamma(n-int(a)+1) for a in k])
        log_probability = (log_binomial+k*np.log((1+z)/2) +
                           (n-k)*np.log((1-z)/2))
        out = np.exp(.5*log_probability+1j*(n-k)*phi)
    # Normalize lgamma floating-point residue, not the analytic formula.
    return out/np.linalg.norm(out)


def evolve_symmetric(state, time, chi=1.):
    n = len(state)-1
    m = np.arange(n+1)-n/2
    return np.exp(-1j*chi*time*m*m/n)*state


def symmetric_moments(n, z, phi, time, chi=1.):
    psi = evolve_symmetric(coherent(n, z, phi), time, chi)
    operators = [2*op/n for op in spin_matrices(n)]
    means = np.array([np.vdot(psi, op @ psi).real for op in operators])
    variances = np.array([np.vdot(op @ psi, op @ psi).real-means[k]**2
                          for k, op in enumerate(operators)])
    return means, variances


def analytic_mean(n, z, phi, time, chi=1.):
    angle = chi*time/n
    initial_plus = np.sqrt(max(0., 1-z*z))*np.exp(1j*phi)
    factor = complex(np.cos(angle), z*np.sin(angle))**(n-1)
    value = initial_plus*factor
    return np.array([value.real, value.imag, z])


def classical_mean(z, phi, time, chi=1.):
    phase = phi+chi*z*time
    radius = np.sqrt(max(0., 1-z*z))
    return np.array([radius*np.cos(phase), radius*np.sin(phase), z])


def mean_error_bound(n, z, time, chi=1.):
    radius = np.sqrt(max(0., 1-z*z))
    angle = abs(chi*time)/n
    return float(radius*(abs(z)*angle +
                         (n-1)*(1-z*z)*angle**2/2))


def total_variance(n, z, time, chi=1.):
    if n == 1:
        return 2.
    delta = (1-z*z)*np.sin(chi*time/n)**2
    if delta >= 1.:
        deficit = 1.
    else:
        deficit = -np.expm1((n-1)*np.log1p(-delta))
    return float(2/n+(1-z*z)*deficit)


def variance_bound(n, z, time, chi=1.):
    return float((2+(1-z*z)**2*(chi*time)**2)/n)


def global_fidelity(n, z, phi, time, chi=1.):
    evolved = evolve_symmetric(coherent(n, z, phi), time, chi)
    classical_product = coherent(n, z, phi+chi*z*time)
    return float(abs(np.vdot(classical_product, evolved))**2)


def tensor_qubits(n, z, phi):
    single = np.array([np.sqrt((1+z)/2),
                       np.exp(1j*phi)*np.sqrt((1-z)/2)])
    psi = np.array([1.+0j])
    for _ in range(n):
        psi = np.kron(psi, single)
    return psi


def collective_full(n):
    result = []
    for pauli in (X, Y, Z):
        collective = np.zeros((2**n, 2**n), dtype=complex)
        for site in range(n):
            term = np.array([[1.+0j]])
            for other in range(n):
                term = np.kron(term, pauli if other == site else np.eye(2))
            collective += term/2
        result.append(collective)
    return result


def full_moments(n, z, phi, time, chi=1.):
    operators = collective_full(n)
    psi = tensor_qubits(n, z, phi)
    energy = chi*np.diag(operators[2]).real**2/n
    psi = np.exp(-1j*time*energy)*psi
    return np.array([2*np.vdot(psi, op @ psi).real/n for op in operators])


def poisson_bracket(point, grad_f, grad_g):
    return float(np.dot(point, np.cross(grad_f, grad_g)))


def report():
    fixed_time = []
    for n in (16, 64, 256, 1024):
        z, phi, time, chi = .4, .3, 1.2, .8
        error = float(np.linalg.norm(analytic_mean(n, z, phi, time, chi) -
                                    classical_mean(z, phi, time, chi)))
        fixed_time.append({
            "N": n, "z": z, "phi": phi, "time": time, "chi": chi,
            "collective_mean_error": error,
            "analytic_uniform_error_bound": mean_error_bound(n, z, time, chi),
            "total_collective_variance": total_variance(n, z, time, chi),
            "variance_upper": variance_bound(n, z, time, chi)})
    collapse = []
    for n in (64, 256, 1024, 4096):
        time = np.sqrt(n)
        radius = float(np.linalg.norm(analytic_mean(n, 0., 0., time)[:2]))
        collapse.append({"N": n, "time": float(time),
                         "time_over_sqrt_N": 1.,
                         "equatorial_contrast": radius,
                         "contrast_limit": float(np.exp(-.5)),
                         "mean_error": 1-radius,
                         "total_variance": total_variance(n, 0., time),
                         "variance_limit": float(1-np.exp(-1.))})
    fidelities = []
    for n in (32, 128, 512, 2048):
        fidelity = global_fidelity(n, 0., 0., 2.)
        fidelities.append({"N": n, "fixed_time": 2., "chi": 1.,
                           "equatorial_global_fidelity": fidelity,
                           "fidelity_limit": float(1/np.sqrt(2)),
                           "pure_state_trace_distance": float(np.sqrt(1-fidelity)),
                           "collective_mean_error": float(
                               1-analytic_mean(n, 0., 0., 2.)[0])})
    return {
        "round": 354,
        "scope": {
            "kind": "specified coherent-state collective-observable classical limit",
            "proved": ["sphere Poisson brackets with effective hbar=2/N",
                       "exact one-axis-twisting first moment for every N",
                       "uniform fixed-time O(1/N) mean error and variance bound",
                       "order-one mean/concentration failure at time~sqrt(N)",
                       "local Darboux pair (phi,z), not a global plane",
                       "fixed-time collective convergence does not imply "
                       "full-state convergence to a coherent product"],
            "added_inputs": ["N labeled two-level systems and collective Pauli algebra",
                             "symmetric maximum-spin sector",
                             "identical pure coherent product initial states",
                             "specific H=chi*Jz^2/N and time normalization",
                             "restriction to fixed times and collective observations"],
            "not_proved": ["FUCP selects the state family or the Hamiltonian",
                           "all quantum states approach a classical point",
                           "channel or full-state quantum-to-classical convergence",
                           "spatial graph, three spatial dimensions, or a metric",
                           "ADM canonical metric fields or hypersurface constraints",
                           "gravitational dynamics or general-relativistic emergence"]},
        "fixed_time_rows": fixed_time,
        "sqrt_N_collapse_rows": collapse,
        "fixed_time_global_fidelity_rows": fidelities,
        "GHZ_boundary": {
            "N_at_least": 2,
            "normalized_sz_variance": 1.,
            "normalized_collective_mean": [0., 0., 0.],
            "global_X_product_expectation_pure_GHZ": 1.,
            "global_X_product_expectation_incoherent_mixture": 0.,
            "trace_distance_GHZ_to_mixture": .5,
            "interpretation": "No point concentration; classical mixtures of "
                              "sphere points are not excluded by this example."},
        "resources": {
            "physical_register": "N qubits, full Hilbert dimension 2^N",
            "symmetric_dimension": "N+1, used for exact computational compression",
            "Hamiltonian_norm": "abs(chi)*N/4",
            "interaction": "chi/4*I + chi/(2N)*sum_{a<b} Z_a Z_b",
            "no_free_geometry": "Permutation-symmetric all-to-all coupling in "
                                "labels supplies no spatial metric or local graph."},
        "sources": [
            "https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevA.47.5138/fulltext",
            "https://arxiv.org/pdf/2007.03390",
            "https://collaborate.princeton.edu/en/publications/the-classical-limit-of-quantum-spin-systems/"]}


class Checks(unittest.TestCase):
    def test_01_collective_commutator_and_effective_hbar(self):
        for n in (1, 2, 5, 16):
            spin = [2*op/n for op in spin_matrices(n)]
            for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
                np.testing.assert_allclose(spin[i] @ spin[j]-spin[j] @ spin[i],
                                           2j*spin[k]/n, atol=5e-15)

    def test_02_casimir_and_coherent_initial_covariance(self):
        for n in (1, 4, 13):
            spin = [2*op/n for op in spin_matrices(n)]
            np.testing.assert_allclose(sum(op @ op for op in spin),
                                       (1+2/n)*np.eye(n+1), atol=2e-15)
            psi = coherent(n, .3, .7)
            mean = classical_mean(.3, .7, 0.)
            for i in range(3):
                for j in range(3):
                    cov = np.vdot(psi, (spin[i] @ spin[j]+spin[j] @ spin[i]) @
                                  psi).real/2-mean[i]*mean[j]
                    self.assertAlmostEqual(cov,
                                           (float(i == j)-mean[i]*mean[j])/n,
                                           places=13)

    def test_03_exact_formula_against_symmetric_matrix_evolution(self):
        for n in (1, 2, 7, 32):
            for z in (-1., -.6, 0., .37, 1.):
                actual, _ = symmetric_moments(n, z, .43, 1.4, -.8)
                np.testing.assert_allclose(actual,
                                           analytic_mean(n, z, .43, 1.4, -.8),
                                           atol=3e-14)

    def test_04_independent_full_tensor_realization(self):
        for n in (2, 3, 5):
            np.testing.assert_allclose(full_moments(n, .27, -.8, 1.9, .7),
                                       analytic_mean(n, .27, -.8, 1.9, .7),
                                       atol=3e-14)

    def test_05_uniform_fixed_time_mean_bound(self):
        for n in (1, 3, 16, 128):
            for z in (-.9, -.4, 0., .2, .8):
                for time in np.linspace(-3., 3., 19):
                    error = np.linalg.norm(analytic_mean(n, z, .2, time) -
                                           classical_mean(z, .2, time))
                    self.assertLessEqual(error, mean_error_bound(n, z, 3.)+2e-14)

    def test_06_exact_variance_and_uniform_concentration_bound(self):
        for n in (1, 2, 9, 24):
            for z in (-.7, 0., .3, 1.):
                for time in (0., .4, 1.7):
                    _, variances = symmetric_moments(n, z, .2, time)
                    exact = total_variance(n, z, time)
                    self.assertAlmostEqual(float(sum(variances)), exact, places=12)
                    self.assertLessEqual(exact, variance_bound(n, z, time)+1e-14)

    def test_07_energy_and_z_distribution_are_preserved(self):
        n = 23
        initial = coherent(n, .4, -.3)
        final = evolve_symmetric(initial, 2.7, .8)
        np.testing.assert_allclose(abs(final)**2, abs(initial)**2, atol=1e-15)
        m = np.arange(n+1)-n/2
        energy = .8*m*m/n
        self.assertAlmostEqual(float(np.dot(abs(initial)**2, energy)),
                               float(np.dot(abs(final)**2, energy)), places=13)

    def test_08_sqrt_N_time_has_order_one_collapse(self):
        for z in (0., .4, -.6):
            radius = np.sqrt(1-z*z)
            for n in (1024, 4096):
                time = np.sqrt(n)
                mean = analytic_mean(n, z, .2, time)
                contrast = np.linalg.norm(mean[:2])/radius
                target = np.exp(-(1-z*z)/2)
                self.assertLess(abs(contrast-target), 2/n)
        self.assertGreater(1-analytic_mean(4096, 0., 0., 64.)[0], .39)

    def test_09_recurrence_prevents_irreversible_classical_claim(self):
        for n in (2, 3, 8):
            state = coherent(n, .35, .2)
            returned = evolve_symmetric(state, 2*np.pi*n)
            self.assertAlmostEqual(abs(np.vdot(state, returned)), 1., places=13)

    def test_10_sphere_poisson_flow_and_darboux_patch(self):
        for z, phi in ((-.8, .2), (0., 1.1), (.4, -2.)):
            point = classical_mean(z, phi, 0.)
            grad_phi = np.array([-point[1], point[0], 0.])/(1-z*z)
            grad_z = np.array([0., 0., 1.])
            self.assertAlmostEqual(poisson_bracket(point, grad_phi, grad_z), 1.)
            grad_h = z*grad_z
            velocity = np.array([poisson_bracket(point, np.eye(3)[i], grad_h)
                                 for i in range(3)])
            eps = 1e-6
            numerical = (classical_mean(z, phi, eps) -
                         classical_mean(z, phi, -eps))/(2*eps)
            np.testing.assert_allclose(velocity, numerical, atol=1e-10)

    def test_11_global_quantum_state_does_not_approach_product(self):
        for time in (.7, 1.3, 2.):
            target = 1/np.sqrt(1+(time/2)**2)
            for n in (128, 512, 2048):
                fidelity = global_fidelity(n, 0., 0., time)
                self.assertLess(abs(fidelity-target), 1/n)
                if time == 2.:
                    self.assertLess(fidelity, .72)
        self.assertLess(1-analytic_mean(2048, 0., 0., 2.)[0], .001)

    def test_12_ghz_is_not_a_concentrated_classical_point(self):
        for n in (2, 3, 5):
            operators = [2*op/n for op in collective_full(n)]
            ghz = np.zeros(2**n, dtype=complex)
            ghz[0] = ghz[-1] = 1/np.sqrt(2)
            self.assertAlmostEqual(np.vdot(ghz, operators[2] @
                                           operators[2] @ ghz).real, 1.)
            for op in operators:
                self.assertAlmostEqual(np.vdot(ghz, op @ ghz).real, 0.)
            parity = np.array([[1.+0j]])
            for _ in range(n):
                parity = np.kron(parity, X)
            self.assertAlmostEqual(np.vdot(ghz, parity @ ghz).real, 1.)
            self.assertEqual((parity[0, 0]+parity[-1, -1]).real/2, 0.)
            mixture = np.zeros((2**n, 2**n), dtype=complex)
            mixture[0, 0] = mixture[-1, -1] = .5
            difference = np.outer(ghz, ghz.conj())-mixture
            self.assertAlmostEqual(float(np.sum(abs(
                np.linalg.eigvalsh(difference)))/2), .5)
            energy = n*np.diag(operators[2]).real**2/4
            evolved = np.exp(-1j*.7*energy)*ghz
            self.assertAlmostEqual(abs(np.vdot(ghz, evolved)), 1.)

    def test_13_all_to_all_resource_identity(self):
        for n in (2, 3, 5):
            jz = collective_full(n)[2]
            # Diagonal product basis independently computes pair products.
            labels = np.arange(2**n)
            spins = np.array([1-2*((labels >> k) & 1) for k in range(n)])
            pair_sum = sum(spins[i]*spins[j] for i in range(n) for j in range(i))
            np.testing.assert_allclose(np.diag(jz @ jz/n),
                                       .25+pair_sum/(2*n), atol=1e-14)
            self.assertAlmostEqual(float(np.max(np.diag(jz @ jz/n).real)), n/4)

    def test_14_scaling_hamiltonian_and_finite_N_heisenberg_equation(self):
        for n in (2, 5, 13):
            jx, jy, jz = spin_matrices(n)
            plus = 2*(jx+1j*jy)/n
            sz = 2*jz/n
            h = jz @ jz/n
            actual = 1j*(h @ plus-plus @ h)
            expected = .5j*(sz @ plus+plus @ sz)
            np.testing.assert_allclose(actual, expected, atol=2e-15)


if __name__ == "__main__":
    main(__name__, "collective_poisson_limit_audit", report)
