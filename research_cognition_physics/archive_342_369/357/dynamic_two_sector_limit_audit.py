"""Round 357: changing sources and a controlled two-sector mean-field limit.

Explicit mean-square and bias estimates, not an arbitrary-input channel claim.
The microscopic interaction is all-to-all and supplies no spatial metric.
"""
from math import ceil, factorial
import unittest
import numpy as np
from growing_stream_audit import main
from collective_poisson_limit_audit import coherent, collective_full, tensor_qubits


EX = np.array([1., 0., 0.])
EZ = np.array([0., 0., 1.])
PARAMETERS = (.7, -.45, .35)  # Omega_G, Omega_M, g
INITIAL = (.3, .25, -.4, -.3)  # z_G, phi_G, z_M, phi_M


def initial_bloch(initial=INITIAL):
    result = []
    for z, phi in (initial[:2], initial[2:]):
        r = np.sqrt(1-z*z)
        result.append([r*np.cos(phi), r*np.sin(phi), z])
    return np.array(result)


def classical_rhs(vectors, parameters=PARAMETERS):
    og, om, g = parameters
    a, b = vectors
    return np.array([2*np.cross(og*EX+g*b[2]*EZ, a),
                     2*np.cross(om*EX+g*a[2]*EZ, b)])


def classical(time, parameters=PARAMETERS, initial=INITIAL, max_step=.002):
    value = initial_bloch(initial)
    steps = max(1, int(ceil(abs(time)/max_step)))
    dt = time/steps
    for _ in range(steps):
        k1 = classical_rhs(value, parameters)
        k2 = classical_rhs(value+dt*k1/2, parameters)
        k3 = classical_rhs(value+dt*k2/2, parameters)
        k4 = classical_rhs(value+dt*k3, parameters)
        value += dt*(k1+2*k2+2*k3+k4)/6
    return value


def h_action(psi, parameters=PARAMETERS):
    n = psi.shape[0]-1
    og, om, g = parameters
    z = 2*np.arange(n+1)-n
    hop = np.sqrt(np.arange(1, n+1)*np.arange(n, 0, -1))
    out = (g/n)*z[:, None]*z[None, :]*psi
    out[1:, :] += og*hop[:, None]*psi[:-1, :]
    out[:-1, :] += og*hop[:, None]*psi[1:, :]
    out[:, 1:] += om*hop[None, :]*psi[:, :-1]
    out[:, :-1] += om*hop[None, :]*psi[:, 1:]
    return out


def normalized_spin(psi, sector, axis):
    n = psi.shape[0]-1
    if sector == 1:
        return normalized_spin(psi.T, 0, axis).T
    if axis == 2:
        return ((2*np.arange(n+1)-n)/n)[:, None]*psi
    hop = np.sqrt(np.arange(1, n+1)*np.arange(n, 0, -1))/n
    out = np.zeros_like(psi)
    lower, upper = (1., 1.) if axis == 0 else (-1j, 1j)
    out[1:, :] += lower*hop[:, None]*psi[:-1, :]
    out[:-1, :] += upper*hop[:, None]*psi[1:, :]
    return out


def propagate(psi, time, parameters=PARAMETERS, order=18, max_scaled_step=.5):
    """Matrix-free Taylor action; explicit exact-arithmetic truncation bound."""
    n = psi.shape[0]-1
    norm_upper = n*sum(abs(x) for x in parameters)
    steps = max(1, int(ceil(abs(time)*norm_upper/max_scaled_step)))
    dt = time/steps
    scaled = abs(dt)*norm_upper
    remainder = np.exp(scaled)*scaled**(order+1)/factorial(order+1)
    truncation_upper = steps*remainder*np.exp(steps*remainder)
    value = psi.copy()
    for _ in range(steps):
        term, total = value.copy(), value.copy()
        for degree in range(1, order+1):
            term = (-1j*dt/degree)*h_action(term, parameters)
            total += term
        value = total
    return value, {"substeps": steps, "Taylor_order": order,
                   "scaled_step_norm_upper": scaled,
                   "exact_arithmetic_truncation_upper": float(truncation_upper),
                   "norm_error": float(abs(np.vdot(value, value).real-1))}


def initial_symmetric(n, initial=INITIAL):
    return np.outer(coherent(n, *initial[:2]), coherent(n, *initial[2:]))


def moments(psi):
    means = np.empty((2, 3))
    variances = np.empty((2, 3))
    for sector in range(2):
        for axis in range(3):
            applied = normalized_spin(psi, sector, axis)
            means[sector, axis] = np.vdot(psi, applied).real
            variances[sector, axis] = np.vdot(applied, applied).real-means[sector, axis]**2
    return means, variances


def error_budget(n, time, g):
    a = abs(g*time)
    return {"mean_square_upper": float(4*np.exp(4*a)/n),
            "bias_sum_upper": float(4*np.exp(2*a)*np.expm1(2*a)/n)}


def energy_per_N(psi, parameters=PARAMETERS):
    return float(np.vdot(psi, h_action(psi, parameters)).real/(psi.shape[0]-1))


def classical_energy(vectors, parameters=PARAMETERS):
    og, om, g = parameters
    return float(og*vectors[0, 0]+om*vectors[1, 0]+g*vectors[0, 2]*vectors[1, 2])


def fluctuation_derivative(psi, target, parameters=PARAMETERS):
    """Independent direct commutator derivative and cancellation formula."""
    g = parameters[2]
    hpsi = h_action(psi, parameters)
    target_dot = classical_rhs(target, parameters)
    direct = 0.
    for sector in range(2):
        for axis in range(3):
            delta = normalized_spin(psi, sector, axis)-target[sector, axis]*psi
            delta_h = normalized_spin(hpsi, sector, axis)-target[sector, axis]*hpsi
            # 2 Re <delta psi, delta (-i H psi)> plus explicit target derivative.
            direct += 2*np.vdot(delta, -1j*delta_h).real
            direct -= 2*target_dot[sector, axis]*np.vdot(psi, delta).real
    cancelled = 0.
    for sector in range(2):
        other = 1-sector
        dz = normalized_spin(psi, other, 2)-target[other, 2]*psi
        direction = np.cross(EZ, target[sector])
        transverse = sum(direction[axis]*(normalized_spin(psi, sector, axis)-
                            target[sector, axis]*psi) for axis in range(3))
        cancelled += 4*g*np.vdot(dz, transverse).real
    return float(direct), float(cancelled)


def full_exact(n, time, parameters=PARAMETERS, initial=INITIAL):
    single = [2*op for op in collective_full(n)]
    size = 2**n
    og, om, g = parameters
    h = og*np.kron(single[0], np.eye(size))+om*np.kron(np.eye(size), single[0])
    h += (g/n)*np.kron(single[2], single[2])
    energies, vectors = np.linalg.eigh(h)
    psi = np.kron(tensor_qubits(n, *initial[:2]), tensor_qubits(n, *initial[2:]))
    psi = vectors@(np.exp(-1j*time*energies)*(vectors.conj().T@psi))
    means = []
    variances = []
    for sector in range(2):
        row, var = [], []
        for op in single:
            operator = np.kron(op/n, np.eye(size)) if sector == 0 else np.kron(np.eye(size), op/n)
            value = operator@psi
            mean = np.vdot(psi, value).real
            row.append(mean)
            var.append(np.vdot(value, value).real-mean*mean)
        means.append(row)
        variances.append(var)
    return np.array(means), np.array(variances)


def case(n, time=.9):
    initial = initial_symmetric(n)
    final, integration = propagate(initial, time)
    target = classical(time)
    means, variances = moments(final)
    errors = np.linalg.norm(means-target, axis=1)
    mean_square = float(variances.sum()+np.sum((means-target)**2))
    budget = error_budget(n, time, PARAMETERS[2])
    return {"N_per_sector": n, "time": time, "quantum_mean_G": means[0].tolist(),
            "quantum_mean_M": means[1].tolist(), "classical_mean_G": target[0].tolist(),
            "classical_mean_M": target[1].tolist(), "bias_sum": float(errors.sum()),
            "bias_sum_upper": budget["bias_sum_upper"],
            "total_variance": float(variances.sum()), "mean_square": mean_square,
            "mean_square_upper": budget["mean_square_upper"],
            "N_times_bias_sum": float(n*errors.sum()),
            "energy_per_N_error": abs(energy_per_N(final)-energy_per_N(initial)),
            "numerical_action": integration}


def report():
    initial = initial_bloch()
    target = classical(.9)
    uncoupled = classical(.9, (PARAMETERS[0], PARAMETERS[1], 0.))
    return {"round": 357,
            "scope": {
                "proved": ["For specified pure iid sectors, total collective mean-square deviation E(t)<=4 exp(4|g|t)/N on finite time intervals.",
                           "Sum of the two mean-vector biases <=4[exp(4|g|t)-exp(2|g|t)]/N.",
                           "Both Z sources change under the noncommuting autonomous self dynamics, and each source acceleration depends on the other sector.",
                           "The same joint unitary supports finite bidirectional response and O(1/N) collective variance."],
                "inputs": ["Two N-qubit sectors, pure identical product state in each sector", "Constant Omega_G, Omega_M, g", "H=Omega_G S_Gx+Omega_M S_Mx+(g/N)S_Gz S_Mz", "Fixed times and collective observations"],
                "not_proved": ["An arbitrary-input or arbitrary-reference channel limit", "Full joint state convergence to a product", "Uniform infinite-time accuracy or an optimal validity time", "Space, dimension, local light cones, metric, HDA or Einstein equation", "FUCP selects these states or this Hamiltonian"]},
            "parameters": {"Omega_G": PARAMETERS[0], "Omega_M": PARAMETERS[1], "g": PARAMETERS[2], "initial_zG_phiG_zM_phiM": list(INITIAL)},
            "fixed_time_rows": [case(n) for n in (8, 16, 32, 64)],
            "changing_sources": {"initial_Z": initial[:, 2].tolist(), "final_classical_Z": target[:, 2].tolist(),
                                 "uncoupled_final_Z": uncoupled[:, 2].tolist(),
                                 "coupling_effect_on_Z": (target[:, 2]-uncoupled[:, 2]).tolist(),
                                 "initial_source_velocity": classical_rhs(initial)[:, 2].tolist(),
                                 "interaction_part_of_initial_Z_acceleration": [float(4*PARAMETERS[0]*PARAMETERS[2]*initial[1, 2]*initial[0, 0]), float(4*PARAMETERS[1]*PARAMETERS[2]*initial[0, 2]*initial[1, 0])]},
            "resource_account": {"physical_qubits": "2N", "full_Hilbert_dimension": "4^N", "symmetric_dimension": "(N+1)^2", "cross_pairs": "N^2", "coupling_per_pair": "g/N", "Hamiltonian_norm_upper": "N*(|Omega_G|+|Omega_M|+|g|)", "incident_cross_coupling_per_spin": "|g|", "energy_per_particle_finite": True, "spatial_locality": False},
            "primary_source": "Carollo and Lesanovsky, Phys. Rev. Lett. 133, 150401 (2024), https://arxiv.org/html/2403.17163, Eq.(1), Eq.(4), Theorem 1. N four-level cells embed this model; specific constants and bias bound are derived directly here.",
            "numerics_scope": "Matrix-free degree-18 Taylor action uses substeps with ||H dt||<=0.5 and an exact-arithmetic truncation bound. Floating-point error is checked by norm/energy and small-N full tensor diagonalization, not certified by that truncation bound. Classical ODE uses RK4 with refinement checks."}


class Checks(unittest.TestCase):
    def test_01_full_tensor_independent_evolution(self):
        for n in (1, 2, 3):
            final, _ = propagate(initial_symmetric(n), .7)
            means, variances = moments(final)
            full_mean, full_var = full_exact(n, .7)
            np.testing.assert_allclose(means, full_mean, atol=3e-14)
            np.testing.assert_allclose(variances, full_var, atol=4e-14)

    def test_02_initial_mean_square_budget(self):
        target = initial_bloch()
        for n in (1, 3, 10, 32):
            mean, var = moments(initial_symmetric(n))
            np.testing.assert_allclose(mean, target, atol=5e-14)
            self.assertAlmostEqual(var.sum(), 4/n, places=12)

    def test_03_heisenberg_source_is_not_conserved(self):
        n = 3
        psi = initial_symmetric(n)
        for sector in (0, 1):
            szpsi = normalized_spin(psi, sector, 2)
            commutator = 1j*(h_action(szpsi)-normalized_spin(h_action(psi), sector, 2))
            expected = 2*PARAMETERS[sector]*normalized_spin(psi, sector, 1)
            np.testing.assert_allclose(commutator, expected, atol=3e-15)

    def test_04_fluctuation_cancellation_identity(self):
        for n in (2, 5, 13):
            for time in (.0, .4, .9):
                final, _ = propagate(initial_symmetric(n), time)
                target = classical(time)
                direct, cancelled = fluctuation_derivative(final, target)
                self.assertAlmostEqual(direct, cancelled, places=11)
                mean, var = moments(final)
                error = var.sum()+np.sum((mean-target)**2)
                self.assertLessEqual(abs(cancelled), 4*abs(PARAMETERS[2])*error+2e-13)
        rng = np.random.default_rng(357)
        psi = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
        psi /= np.linalg.norm(psi)
        direct, cancelled = fluctuation_derivative(psi, initial_bloch())
        self.assertAlmostEqual(direct, cancelled, places=12)

    def test_05_finite_time_mean_square_and_bias_bounds(self):
        for n in (2, 8, 24):
            for time in (.2, .6, 1.):
                final, _ = propagate(initial_symmetric(n), time)
                target = classical(time)
                mean, var = moments(final)
                e = var.sum()+np.sum((mean-target)**2)
                bound = error_budget(n, time, PARAMETERS[2])
                self.assertLessEqual(e, bound["mean_square_upper"]+1e-12)
                self.assertLessEqual(np.linalg.norm(mean-target, axis=1).sum(), bound["bias_sum_upper"]+1e-12)

    def test_06_classical_lengths_and_energy(self):
        initial = initial_bloch()
        for time in (.3, .9, 1.4):
            final = classical(time)
            np.testing.assert_allclose(np.linalg.norm(final, axis=1), 1., atol=3e-12)
            self.assertAlmostEqual(classical_energy(final), classical_energy(initial), places=11)

    def test_07_time_step_refinement_and_truncation(self):
        np.testing.assert_allclose(classical(.9), classical(.9, max_step=.001), atol=3e-12)
        a, info = propagate(initial_symmetric(16), .9)
        b, _ = propagate(initial_symmetric(16), .9, order=20, max_scaled_step=.25)
        np.testing.assert_allclose(a, b, atol=4e-14)
        self.assertLess(info["exact_arithmetic_truncation_upper"], 1e-18)
        self.assertLess(info["norm_error"], 1e-13)

    def test_08_energy_and_unitary_reverse(self):
        initial = initial_symmetric(12)
        final, _ = propagate(initial, .9)
        returned, _ = propagate(final, -.9)
        np.testing.assert_allclose(returned, initial, atol=3e-14)
        self.assertAlmostEqual(energy_per_N(final), energy_per_N(initial), places=13)

    def test_09_both_dynamic_sources_respond_to_coupling(self):
        target = classical(.9)
        uncoupled = classical(.9, (PARAMETERS[0], PARAMETERS[1], 0.))
        self.assertTrue(np.all(abs(target[:, 2]-initial_bloch()[:, 2])>.03))
        self.assertTrue(np.all(abs(target[:, 2]-uncoupled[:, 2])>.03))

    def test_10_uncoupled_dynamics_exact_product_means(self):
        parameters = (PARAMETERS[0], PARAMETERS[1], 0.)
        for n in (1, 5, 17):
            final, _ = propagate(initial_symmetric(n), .9, parameters)
            means, var = moments(final)
            np.testing.assert_allclose(means, classical(.9, parameters), atol=2e-12)
            self.assertAlmostEqual(float(var.sum()), 4/n, places=12)

    def test_11_zero_transverse_fields_recover_356(self):
        n, time, g = 7, .9, .35
        final, _ = propagate(initial_symmetric(n), time, (0., 0., g))
        means, _ = moments(final)
        initial = initial_bloch()
        for sector in (0, 1):
            other = 1-sector
            multiplier = (np.cos(2*g*time/n)+1j*initial[other, 2]*np.sin(2*g*time/n))**n
            value = complex(initial[sector, 0], initial[sector, 1])*multiplier
            np.testing.assert_allclose(means[sector], [value.real, value.imag, initial[sector, 2]], atol=3e-14)

    def test_12_source_acceleration_cross_dependence(self):
        initial = initial_bloch()
        velocity = classical_rhs(initial)
        eps = 1e-5
        acceleration = (classical_rhs(initial+eps*velocity)-classical_rhs(initial-eps*velocity))/(2*eps)
        for sector in (0, 1):
            omega, g = PARAMETERS[sector], PARAMETERS[2]
            expected = 4*omega*g*initial[1-sector, 2]*initial[sector, 0]-4*omega**2*initial[sector, 2]
            self.assertAlmostEqual(acceleration[sector, 2], expected, places=10)
            for n in (2, 5):
                psi = initial_symmetric(n)
                hpsi = h_action(psi)
                h2psi = h_action(hpsi)
                exact = 2*np.vdot(hpsi, normalized_spin(hpsi, sector, 2)).real
                exact -= 2*np.vdot(psi, normalized_spin(h2psi, sector, 2)).real
                self.assertAlmostEqual(exact, expected, places=12)

    def test_13_matrix_free_generator_hermitian_and_budget(self):
        n = 3
        d = (n+1)**2
        basis = np.eye(d, dtype=complex)
        h = np.column_stack([h_action(basis[:, j].reshape(n+1, n+1)).ravel() for j in range(d)])
        np.testing.assert_allclose(h, h.conj().T, atol=1e-15)
        self.assertLessEqual(max(abs(np.linalg.eigvalsh(h))), n*sum(abs(x) for x in PARAMETERS)+1e-14)


if __name__ == "__main__":
    main(__name__, "dynamic_two_sector_limit_audit", report)
