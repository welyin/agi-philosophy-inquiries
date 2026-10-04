"""Round 373: calibrate quantum directions against an assigned Weyl flow.

The three-dimensional continuum, constant invertible E, tilt beta, position
measurement, clock and preparation access are inputs. Integrals are on R^3
with compact momentum support, not a periodic sawtooth position observable.
"""
from functools import lru_cache
import unittest
import numpy as np
from growing_stream_audit import main


I = np.eye(2, dtype=complex)
PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]],
                  [[1, 0], [0, -1]]], complex)
E = np.array([[1., .25, 0.], [0., 1.6, .1], [.15, 0., .8]])
BETA = np.array([.15, -.08, .05])
K0 = np.array([1.2, -.7, .9])
SIGMA = np.linalg.norm(E, 2)


def unit(x):
    x = np.asarray(x, float)
    return x/np.linalg.norm(x, axis=-1, keepdims=True)


def sigma(x):
    return np.einsum("...a,aij->...ij", x, PAULI)


def hamiltonian(k, e=E, beta=BETA):
    k = np.asarray(k)
    return (k@beta)[..., None, None]*I + sigma(k@e.T)


def propagator(k, time, e=E, beta=BETA):
    k = np.asarray(k)
    h = k@e.T
    radius = np.linalg.norm(h, axis=-1)
    # sin(t r)/r has a regular value at r=0.
    core = (np.cos(time*radius)[..., None, None]*I -
            1j*(time*np.sinc(time*radius/np.pi))[..., None, None]*sigma(h))
    return np.exp(-1j*time*(k@beta))[..., None, None]*core


def velocity_operators(e=E, beta=BETA):
    return beta[:, None, None]*I + sigma(e.T)


def density_and_spin(psi):
    density = np.sum(abs(psi)**2, axis=-1)
    spin = np.einsum("...i,aij,...j->...a", psi.conjugate(), PAULI, psi).real
    return density, spin


def current(psi, e=E, beta=BETA):
    density, spin = density_and_spin(psi)
    return density[..., None]*beta + spin@e


def pure_spinor(n):
    return np.linalg.eigh(sigma(unit(n)))[1][:, 1]


def upper_spinor(k):
    n = unit(np.asarray(k)@E.T)
    # The chosen compact support stays in this smooth north-pole chart.
    if np.any(1+n[..., 2] < .1):
        raise ValueError("This example requires another spinor chart here.")
    return np.stack((1+n[..., 2], n[..., 0]+1j*n[..., 1]), axis=-1)/np.sqrt(
        2*(1+n[..., 2]))[..., None]


def upper_energy(k):
    return np.asarray(k)@BETA + np.linalg.norm(np.asarray(k)@E.T, axis=-1)


def group_velocity(k):
    return BETA + unit(np.asarray(k)@E.T)@E


def envelope(k, width):
    q = (np.asarray(k)-K0)/width
    inside = np.all(abs(q) < 1, axis=-1)
    value = (315/256)**1.5/width**1.5*np.prod((1-q*q)**2, axis=-1)
    return np.where(inside, value, 0.)


@lru_cache(maxsize=None)
def quadrature(width, order=18):
    q, w = np.polynomial.legendre.leggauss(order)
    coordinates = np.stack(np.meshgrid(q, q, q, indexing="ij"), axis=-1).reshape(-1, 3)
    weights = np.prod(np.stack(np.meshgrid(w, w, w, indexing="ij"), axis=-1),
                      axis=-1).reshape(-1)*width**3
    return K0+width*coordinates, weights


def wavefunction(k, time, width, spinor=None):
    initial = upper_spinor(k) if spinor is None else np.broadcast_to(spinor, (len(k), 2))
    return envelope(k, width)[:, None]*np.einsum(
        "...ij,...j->...i", propagator(k, time), initial)


def position_moments(time, width, spinor=None, order=18, dk=1e-5):
    k, weights = quadrature(width, order)
    psi = wavefunction(k, time, width, spinor)
    centroid, second = [], []
    for axis in range(3):
        step = np.eye(3)[axis]*dk
        derivative = (wavefunction(k+step, time, width, spinor) -
                      wavefunction(k-step, time, width, spinor))/(2*dk)
        centroid.append(float(np.sum(weights*np.einsum(
            "ki,ki->k", psi.conjugate(), 1j*derivative).real)))
        second.append(float(np.sum(weights*np.sum(abs(derivative)**2, axis=1))))
    return np.array(centroid), np.array(second)


def integrated_current(time, width, spinor=None, order=18):
    k, weights = quadrature(width, order)
    return np.sum(weights[:, None]*current(wavefunction(k, time, width, spinor)), axis=0)


def initial_velocity(n):
    return BETA+E.T@np.asarray(n)


def support_bounds(width):
    epsilon = np.sqrt(3)*width
    m0 = np.linalg.norm(E@K0)
    minimum = m0-SIGMA*epsilon
    if minimum <= 0:
        raise ValueError("Support reaches the band crossing; narrow-band bound unavailable.")
    return epsilon, m0, minimum, m0+SIGMA*epsilon


@lru_cache(maxsize=None)
def calibration(step, width=.18, order=18, dk=1e-5):
    velocities, exact, raw = [], [], []
    for axis in range(3):
        pair = []
        for sign in (1., -1.):
            n = sign*np.eye(3)[axis]
            spinor = pure_spinor(n)
            centroids = np.array([position_moments(t, width, spinor, order, dk)[0]
                                  for t in (0., step, 2*step)])
            velocity = (-3*centroids[0]+4*centroids[1]-centroids[2])/(2*step)
            velocities.append(velocity)
            exact.append(initial_velocity(n))
            pair.append(centroids.tolist())
        raw.append(pair)
    velocities, exact = np.array(velocities), np.array(exact)
    learned_e = np.array([(velocities[2*a]-velocities[2*a+1])/2 for a in range(3)])
    learned_beta = velocities.mean(axis=0)
    learned_g = np.linalg.inv(learned_e.T@learned_e)
    _, _, _, maximum = support_bounds(width)
    per_velocity_bound = 4/3*step**2*maximum**2*np.linalg.norm(E, "fro")
    e_bound = np.sqrt(3)*per_velocity_bound
    return {"step": step, "width": width, "order": order, "momentum_difference": dk,
            "centroids_0_h_2h_per_axis_and_sign": raw,
            "measured_velocities": velocities.tolist(),
            "max_velocity_error": float(np.max(np.linalg.norm(velocities-exact, axis=1))),
            "finite_time_velocity_bound": float(per_velocity_bound),
            "learned_E": learned_e.tolist(), "learned_beta": learned_beta.tolist(),
            "learned_inverse_metric": learned_g.tolist(),
            "E_error_operator_norm": float(np.linalg.norm(learned_e-E, 2)),
            "E_finite_time_bound": float(e_bound),
            "beta_error": float(np.linalg.norm(learned_beta-BETA)),
            "G_error_operator_norm": float(np.linalg.norm(learned_e.T@learned_e-E.T@E, 2)),
            "G_finite_time_bound": float(2*SIGMA*e_bound+e_bound**2),
            "inverse_metric_error": float(np.linalg.norm(learned_g-np.linalg.inv(E.T@E), 2)),
            "bound_excludes": "Numerical quadrature/differences and experimental sampling/calibration errors."}


@lru_cache(maxsize=None)
def band_packet(width, time=.7, order=18, dk=1e-5):
    k, weights = quadrature(width, order)
    psi0, psit = wavefunction(k, 0., width), wavefunction(k, time, width)
    v0, energy0 = group_velocity(K0), float(upper_energy(K0))
    phase_linear = energy0+(k-K0)@v0
    rigid = np.exp(-1j*time*phase_linear)[:, None]*psi0
    state_error = np.sqrt(np.sum(weights*np.sum(abs(psit-rigid)**2, axis=1)))
    flux = integrated_current(time, width, order=order)
    x0, second0 = position_moments(0., width, order=order, dk=dk)
    xt, _ = position_moments(time, width, order=order, dk=dk)
    reduced = np.einsum("k,ki,kj->ij", weights, psi0, psi0.conjugate())
    epsilon, m0, minimum, _ = support_bounds(width)
    lipschitz = SIGMA**2/minimum
    g = np.linalg.inv(E.T@E)
    relative = flux-BETA
    return {"width": width, "time": time, "order": order,
            "support_radius_bound": float(epsilon), "band_gap_lower_bound": float(2*minimum),
            "norm": float(np.sum(weights*np.sum(abs(psit)**2, axis=1))),
            "central_group_velocity": v0.tolist(), "integrated_current": flux.tolist(),
            "centroid_initial": x0.tolist(), "centroid_final": xt.tolist(),
            "centroid_shift_over_time": ((xt-x0)/time).tolist(),
            "centroid_flux_error": float(np.linalg.norm((xt-x0)/time-flux)),
            "central_velocity_error": float(np.linalg.norm(flux-v0)),
            "central_velocity_bound": float(lipschitz*epsilon),
            "rigid_translation_state_error": float(state_error),
            "rigid_translation_bound": float(abs(time)*lipschitz*epsilon**2/2),
            "reduced_spin_purity": float(np.trace(reduced@reduced).real),
            "mean_relative_metric_speed_squared": float(relative@g@relative),
            "initial_position_variance_trace": float(second0.sum()-x0@x0),
            "position_variance_lower_bound": float(9/width**2)}


def anisotropic_angles():
    e = np.diag([1., 2., 3.])
    n, m = unit([1., 1., 0.]), unit([1., -1., 0.])
    v, w = e.T@n, e.T@m
    g = np.linalg.inv(e.T@e)
    p = float(np.trace((I+sigma(n))/2@(I+sigma(m))/2).real)
    return {"E": e.tolist(), "n": n.tolist(), "m": m.tolist(),
            "Born_probability": p, "Born_cosine": 2*p-1,
            "coordinate_ray_cosine": float(v@w/np.linalg.norm(v)/np.linalg.norm(w)),
            "metric_ray_cosine": float(v@g@w/np.sqrt((v@g@v)*(w@g@w)))}


def mixed_band_warning():
    spinor = pure_spinor([1., 0., 0.])
    width, time = .18, .7
    first = integrated_current(0., width, spinor)
    later = integrated_current(time, width, spinor)
    return {"width": width, "time": time, "initial_spin": [1., 0., 0.],
            "initial_current": first.tolist(), "later_current": later.tolist(),
            "current_change": float(np.linalg.norm(later-first)),
            "interpretation": "A fixed envelope times one pure spinor is generally a two-band superposition; the initial direction is not a persistent ray."}


def report():
    g = np.linalg.inv(E.T@E)
    packet_rows = [band_packet(width) for width in (.24, .12, .06)]
    convergence = [band_packet(.12, order=n)["centroid_shift_over_time"] for n in (12, 18, 24)]
    return {
        "round": 373,
        "scope": {
            "baseline": "Frozen through 372; only one new propagation-calibration interface.",
            "claim": "For the assigned constant Weyl generator, conserved probability flow calibrates the internal direction map and its inverse spatial metric. Upper-band narrow packets connect that flow to group velocity with bounds.",
            "inputs": ["R^3 position continuum and clock, hbar=1", "constant invertible real E and tilt beta", "natural generator H(k)=beta.k I+(E k).sigma", "repeated preparation of one momentum envelope with six labelled spin directions", "position means at 0,h,2h", "separately prepared upper-band packet for a persistent ray"],
            "not_proved": ["why space has three dimensions", "FUCP selects this generator or its coefficients", "any internal pure state is an outgoing ray", "a universal metric for all species", "events or translations emerge from a qubit cone", "Einstein dynamics or a full field theory"],
            "integration": "R^3 compact-support momentum integrals; no periodic position or assumed finite-grid canonical commutator.",
            "resources": "Calibration requires independent ensembles, known internal settings, clocks and position access. Narrow-band preparation broadens spatial extent; quadrature is not a finite-shot experiment."},
        "model": {"E": E.tolist(), "beta": BETA.tolist(), "k0": K0.tolist(),
                  "G": (E.T@E).tolist(), "inverse_metric_g": g.tolist(),
                  "tilt_squared_in_g": float(BETA@g@BETA),
                  "upper_branch_positive_for_all_nonzero_k_in_this_model": bool(BETA@g@BETA < 1)},
        "six_direction_calibrations": [calibration(step) for step in (.004, .002, .001)],
        "anisotropic_angle_witness": anisotropic_angles(),
        "upper_band_packets": packet_rows,
        "quadrature_refinement": {"orders": [12, 18, 24],
            "centroid_shift_over_time": convergence,
            "last_difference": float(np.linalg.norm(np.array(convergence[-1])-convergence[-2]))},
        "non_band_preparation": mixed_band_warning(),
        "sources": ["https://arxiv.org/html/1306.1934",
                    "https://arxiv.org/html/1507.01603",
                    "https://arxiv.org/html/1702.04624"]}


class Checks(unittest.TestCase):
    def test_01_hermitian_symbol_spectrum_and_actual_propagator(self):
        for k in [K0, np.array([-.2, .7, .3]), np.zeros(3)]:
            h = hamiltonian(k)
            energies, eigenvectors = np.linalg.eigh(h)
            np.testing.assert_allclose(energies,
                k@BETA+np.array([-1., 1.])*np.linalg.norm(E@k), atol=1e-14)
            for t in (0., .13, .7):
                expected = (eigenvectors*np.exp(-1j*t*energies))@eigenvectors.conjugate().T
                np.testing.assert_allclose(propagator(k, t), expected, atol=1e-14)
            omega = .47
            self.assertAlmostEqual(np.linalg.det(omega*I-h).real,
                (omega-BETA@k)**2-k@(E.T@E)@k, places=13)

    def test_02_local_continuity_with_independent_space_and_time_differences(self):
        modes = np.array([[.4, -.7, .2], [1.1, .3, -.5], [-.6, .1, .8]])
        amplitudes = np.array([[1., .3j], [.2+.1j, -.7], [-.4j, .25]], complex)
        def field(x, t):
            return np.sum(np.exp(1j*(modes@x))[:, None]*
                          np.einsum("kij,kj->ki", propagator(modes, t), amplitudes), axis=0)
        x, t, h = np.array([.3, -.4, .2]), .21, 2e-5
        dt_density = (density_and_spin(field(x, t+h))[0]-density_and_spin(field(x, t-h))[0])/(2*h)
        divergence = 0.
        for i in range(3):
            step = np.eye(3)[i]*h
            divergence += (current(field(x+step, t))[i]-current(field(x-step, t))[i])/(2*h)
        self.assertLess(abs(dt_density+divergence), 2e-8)

    def test_03_positive_density_maps_to_future_probability_current_cone(self):
        g = np.linalg.inv(E.T@E)
        for r in ([0., 0., 0.], [.2, -.3, .4], unit([1., 2., 3.])):
            density = (I+sigma(r))/2
            flow = np.array([np.trace(density@op).real for op in velocity_operators()])
            relative = flow-BETA
            self.assertAlmostEqual(float(relative@g@relative), float(np.asarray(r)@r), places=13)
            self.assertAlmostEqual(1-float(relative@g@relative),
                                   4*np.linalg.det(density).real, places=13)
            self.assertLessEqual(float(relative@g@relative), 1+2e-14)

    def test_04_six_initial_centroids_recover_assigned_coefficients(self):
        row = calibration(.001)
        self.assertLess(row["E_error_operator_norm"], 2e-5)
        self.assertLess(row["beta_error"], 2e-8)
        self.assertLess(row["inverse_metric_error"], 4e-5)
        self.assertLess(row["max_velocity_error"], row["finite_time_velocity_bound"])

    def test_05_forward_only_calibration_has_second_order_time_error(self):
        rows = [calibration(h) for h in (.004, .002, .001)]
        for first, second in zip(rows[:-1], rows[1:]):
            ratio = first["E_error_operator_norm"]/second["E_error_operator_norm"]
            self.assertGreater(ratio, 3.7)
            self.assertLess(ratio, 4.3)
        for row in rows:
            self.assertLess(row["E_error_operator_norm"], row["E_finite_time_bound"])
            self.assertLess(row["G_error_operator_norm"], row["G_finite_time_bound"])

    def test_06_internal_frame_changes_leave_inferred_spatial_metric_unchanged(self):
        angle = .73
        rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                             [np.sin(angle), np.cos(angle), 0], [0, 0, 1.]])
        altered = rotation@E
        n = unit([.2, -.3, .7])
        np.testing.assert_allclose(altered.T@(rotation@n), E.T@n, atol=1e-14)
        np.testing.assert_allclose(altered.T@altered, E.T@E, atol=1e-14)

    def test_07_Born_angles_equal_calibrated_metric_angles(self):
        g = np.linalg.inv(E.T@E)
        directions = [unit(v) for v in ([1, 2, 3], [-2, .3, 1], [.2, 1, -.5])]
        for n in directions:
            for m in directions:
                v, w = E.T@n, E.T@m
                probability = np.trace((I+sigma(n))/2@(I+sigma(m))/2).real
                self.assertAlmostEqual(2*probability-1, v@g@w, places=13)
                self.assertAlmostEqual(v@g@v, 1., places=13)

    def test_08_non_axis_orthogonal_directions_need_not_give_coordinate_right_angle(self):
        row = anisotropic_angles()
        self.assertAlmostEqual(row["Born_cosine"], 0., places=14)
        self.assertAlmostEqual(row["metric_ray_cosine"], 0., places=14)
        self.assertAlmostEqual(row["coordinate_ray_cosine"], -.6, places=14)

    def test_09_upper_branch_group_velocity_is_energy_gradient(self):
        for k in [K0, np.array([.4, -.3, .8])]:
            h = 2e-6
            difference = np.array([(upper_energy(k+h*a)-upper_energy(k-h*a))/(2*h)
                                   for a in np.eye(3)])
            np.testing.assert_allclose(difference, group_velocity(k), atol=3e-10)

    def test_10_projected_packet_stays_in_upper_band_and_preserves_norm(self):
        k, weights = quadrature(.12)
        for t in (0., .3, .7):
            psi = wavefunction(k, t, .12)
            np.testing.assert_allclose(np.einsum("kij,kj->ki", hamiltonian(k), psi),
                                       upper_energy(k)[:, None]*psi, atol=7e-14)
            self.assertAlmostEqual(float(np.sum(weights*np.sum(abs(psi)**2, axis=1))), 1., places=13)
        np.testing.assert_allclose(integrated_current(0., .12), integrated_current(.7, .12), atol=3e-14)

    def test_11_initial_pure_internal_state_is_not_generically_persistent_ray(self):
        row = mixed_band_warning()
        np.testing.assert_allclose(row["initial_current"], initial_velocity([1, 0, 0]), atol=2e-14)
        self.assertGreater(row["current_change"], .5)

    def test_12_narrow_band_controls_velocity_and_rigid_translation_error(self):
        rows = [band_packet(width) for width in (.24, .12, .06)]
        for row in rows:
            self.assertLess(row["central_velocity_error"], row["central_velocity_bound"])
            self.assertLess(row["rigid_translation_state_error"], row["rigid_translation_bound"])
            self.assertGreater(row["reduced_spin_purity"], .95)
            self.assertLess(row["reduced_spin_purity"], 1.)
            self.assertLess(row["mean_relative_metric_speed_squared"], 1.)
        self.assertGreater(rows[0]["central_velocity_error"], rows[-1]["central_velocity_error"])
        self.assertGreater(rows[0]["rigid_translation_state_error"],
                           rows[-1]["rigid_translation_state_error"])

    def test_13_R3_position_integral_gives_actual_centroid_motion(self):
        for width in (.24, .12, .06):
            row = band_packet(width)
            self.assertLess(row["centroid_flux_error"], 3e-7)
            displacement = np.array(row["centroid_final"])-row["centroid_initial"]
            self.assertLess(np.linalg.norm(displacement-row["time"]*np.array(
                row["central_group_velocity"])), row["time"]*row["central_velocity_bound"])

    def test_14_independent_quadrature_refinement_and_momentum_difference(self):
        rows = [band_packet(.12, order=n) for n in (12, 18, 24)]
        self.assertLess(np.linalg.norm(np.array(rows[-1]["integrated_current"])-
                                      rows[-2]["integrated_current"]), 2e-13)
        self.assertLess(np.linalg.norm(np.array(rows[-1]["centroid_shift_over_time"])-
                                      rows[-2]["centroid_shift_over_time"]), 3e-9)
        finer = band_packet(.12, dk=5e-6)
        self.assertLess(finer["centroid_flux_error"], rows[1]["centroid_flux_error"])

    def test_15_narrower_band_requires_larger_position_spread(self):
        for width in (.24, .12, .06):
            row = band_packet(width)
            self.assertGreater(row["initial_position_variance_trace"],
                               row["position_variance_lower_bound"])
        # Constant-spin envelope has exactly 3/width^2 variance per coordinate.
        width = .12
        mean, moments = position_moments(0., width, pure_spinor([1, 0, 0]))
        np.testing.assert_allclose(mean, 0., atol=1e-13)
        np.testing.assert_allclose(moments, np.full(3, 3/width**2), rtol=3e-8)


if __name__ == "__main__":
    main(__name__, "propagation_direction_calibration_audit", report)
