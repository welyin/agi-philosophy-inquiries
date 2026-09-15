"""Test a Hamiltonian quantum-classical feedback model, not Einstein gravity."""

import argparse
import json
import math
import platform
import unittest
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp

from relativistic_limit import SIGMA_X, SIGMA_Z


@dataclass(frozen=True)
class Parameters:
    sites: int = 128
    length: float = 32.0
    base_velocity: float = 0.7
    velocity_amplitude: float = 0.2
    coupling: float = 1.0
    background_speed: float = 0.45
    background_frequency: float = 0.6
    matter_mass: float = 0.6

    def __post_init__(self):
        if self.sites < 4 or self.length <= 0:
            raise ValueError("Use at least four sites and a positive interval length.")
        if not 0 <= self.velocity_amplitude < self.base_velocity:
            raise ValueError("Velocity must remain positive for every real background field.")
        if self.background_speed < 0 or self.background_frequency <= 0:
            raise ValueError("Use a nonnegative field speed and a positive restoring frequency.")

    @property
    def spacing(self):
        return self.length / self.sites

    @property
    def positions(self):
        return (np.arange(self.sites) - self.sites // 2) * self.spacing


def centered_derivative(values, spacing):
    return (np.roll(values, -1, axis=0) - np.roll(values, 1, axis=0)) / (2 * spacing)


def lattice_laplacian(values, spacing):
    return (np.roll(values, -1, axis=0) - 2 * values + np.roll(values, 1, axis=0)) / spacing**2


def velocity(field, parameters):
    return parameters.base_velocity + parameters.velocity_amplitude * np.tanh(parameters.coupling * field)


def velocity_derivative(field, parameters):
    return parameters.velocity_amplitude * parameters.coupling * (1 - np.tanh(parameters.coupling * field)**2)


def apply_matter_hamiltonian(state, field, parameters):
    speed = velocity(field, parameters)
    spatial = speed[:, None] * centered_derivative(state, parameters.spacing)
    spatial += centered_derivative(speed[:, None] * state, parameters.spacing)
    return -0.5j * (spatial @ SIGMA_Z.T) + parameters.matter_mass * (state @ SIGMA_X.T)


def kinetic_density(state, parameters):
    differentiated = centered_derivative(state, parameters.spacing) @ SIGMA_Z.T
    return np.sum(state.conj() * differentiated, axis=1).imag


def energy_components(state, field, momentum, parameters):
    gradient = (np.roll(field, -1) - field) / parameters.spacing
    field_density = (momentum**2 + parameters.background_speed**2 * gradient**2
                     + parameters.background_frequency**2 * field**2) / 2
    matter = float(np.vdot(state, apply_matter_hamiltonian(state, field, parameters)).real)
    background = float(parameters.spacing * np.sum(field_density))
    return {"matter": matter, "background": background, "total": matter + background}


def field_energy_gradient(state, field, parameters):
    return velocity_derivative(field, parameters) * kinetic_density(state, parameters)


def derivatives(state, field, momentum, parameters, feedback=True):
    state_derivative = -1j * apply_matter_hamiltonian(state, field, parameters)
    field_derivative = momentum.copy()
    momentum_derivative = (parameters.background_speed**2 * lattice_laplacian(field, parameters.spacing)
                           - parameters.background_frequency**2 * field)
    if feedback:
        momentum_derivative -= field_energy_gradient(state, field, parameters) / parameters.spacing
    return state_derivative, field_derivative, momentum_derivative


def energy_rates(state, field, momentum, parameters, feedback=True):
    state_derivative, field_derivative, momentum_derivative = derivatives(state, field, momentum, parameters, feedback)
    state_part = 2 * np.vdot(apply_matter_hamiltonian(state, field, parameters), state_derivative).real
    matter_rate = float(state_part + np.dot(field_energy_gradient(state, field, parameters), field_derivative))
    field_gradient = (-parameters.background_speed**2 * lattice_laplacian(field, parameters.spacing)
                      + parameters.background_frequency**2 * field)
    background_rate = float(parameters.spacing * (np.dot(momentum, momentum_derivative) + np.dot(field_gradient, field_derivative)))
    return {"matter_rate": matter_rate, "background_rate": background_rate, "total_rate": matter_rate + background_rate}


def initial_conditions(parameters, center=-4.0, width=1.2, wave_number=1.0):
    lattice_momentum = math.sin(wave_number * parameters.spacing) / parameters.spacing
    mode_hamiltonian = parameters.base_velocity * lattice_momentum * SIGMA_Z + parameters.matter_mass * SIGMA_X
    _, spinors = np.linalg.eigh(mode_hamiltonian)
    envelope = np.exp(-(parameters.positions - center)**2 / (4 * width**2) + 1j * wave_number * parameters.positions)
    state = envelope[:, None] * spinors[:, 1][None, :]
    state /= np.linalg.norm(state)
    return state, np.zeros(parameters.sites), np.zeros(parameters.sites)


def algebra_diagnostics():
    parameters = Parameters(sites=32, length=16)
    state, _, _ = initial_conditions(parameters, center=-2)
    field = 0.1 * np.cos(2 * math.pi * parameters.positions / parameters.length)
    momentum = 0.15 + 0.02 * np.sin(2 * math.pi * parameters.positions / parameters.length)
    return {
        "model": "Normalized one-particle spinor plus classical canonical scalar field with a shared autonomous energy functional.",
        "feedback": "The scalar source is minus the derivative of the matter energy with respect to the field, divided by the lattice spacing.",
        "energy_exchange_at_test_state": energy_rates(state, field, momentum, parameters),
        "one_way_control_at_test_state": energy_rates(state, field, momentum, parameters, feedback=False),
        "quantum_norm_rate": float(2 * np.vdot(state, derivatives(state, field, momentum, parameters)[0]).real),
        "not_established": ["Quantum dynamics of geometry", "A covariant gravity action", "Universality of gravitational coupling", "No-signalling of a fundamental hybrid theory", "A stable continuum vacuum"],
    }


def pack_state(state, field, momentum):
    return np.concatenate((state.real.ravel(), state.imag.ravel(), field, momentum))


def unpack_state(values, parameters):
    sites = parameters.sites
    state = (values[:2 * sites] + 1j * values[2 * sites:4 * sites]).reshape(sites, 2)
    return state, values[4 * sites:5 * sites], values[5 * sites:6 * sites]


def vector_field(time, values, parameters, feedback=True):
    current = unpack_state(values, parameters)
    return pack_state(*derivatives(*current, parameters, feedback))


def integrate(parameters, duration=6.0, tolerance=1e-9, method="DOP853", feedback=True, initial=None, samples=121):
    if duration <= 0 or tolerance <= 0 or samples < 2:
        raise ValueError("Use a positive duration, positive tolerance, and at least two observation times.")
    initial = initial_conditions(parameters) if initial is None else initial
    solution = solve_ivp(
        vector_field, (0.0, duration), pack_state(*initial), args=(parameters, feedback),
        method=method, rtol=tolerance, atol=tolerance * 0.01,
        t_eval=np.linspace(0, duration, samples),
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution


def trajectory_summary(solution, parameters):
    energies = [energy_components(*unpack_state(column, parameters), parameters) for column in solution.y.T]
    norms = [float(np.vdot(unpack_state(column, parameters)[0], unpack_state(column, parameters)[0]).real)
             for column in solution.y.T]
    final_state, final_field, final_momentum = unpack_state(solution.y[:, -1], parameters)
    final_density = np.sum(abs(final_state)**2, axis=1)
    return {
        "initial_energy": energies[0],
        "final_energy": energies[-1],
        "max_sampled_absolute_total_energy_drift": max(abs(energy["total"] - energies[0]["total"]) for energy in energies),
        "max_sampled_norm_squared_error": max(abs(norm - norms[0]) for norm in norms),
        "field_max_abs_at_final": float(np.max(abs(final_field))),
        "momentum_max_abs_at_final": float(np.max(abs(final_momentum))),
        "velocity_max_change_at_final": float(np.max(abs(velocity(final_field, parameters) - parameters.base_velocity))),
        "velocity_range_at_final": [float(np.min(velocity(final_field, parameters))), float(np.max(velocity(final_field, parameters)))],
        "background_response_l2_at_final": float(math.sqrt(parameters.spacing) * np.linalg.norm(final_field)),
        "matter_position_mean_at_final": float(np.dot(parameters.positions, final_density)),
        "rhs_evaluations": int(solution.nfev),
        "energy_sample_times": [float(solution.t[index]) for index in range(0, len(solution.t), 10)],
        "energy_samples": [energies[index] for index in range(0, len(energies), 10)],
        "sampling_caveat": "Errors are maxima at the specified observation times, not rigorous bounds at every internal solver step.",
    }


def free_lattice_reference(state, duration, parameters):
    wave_numbers = 2 * math.pi * np.fft.fftfreq(parameters.sites, d=parameters.spacing)
    lattice_momenta = np.sin(wave_numbers * parameters.spacing) / parameters.spacing
    transformed = np.fft.fft(state, axis=0)
    energies = np.sqrt((parameters.base_velocity * lattice_momenta)**2 + parameters.matter_mass**2)
    applied = (parameters.base_velocity * lattice_momenta[:, None] * (transformed @ SIGMA_Z.T)
               + parameters.matter_mass * (transformed @ SIGMA_X.T))
    evolved = np.cos(energies * duration)[:, None] * transformed
    evolved -= 1j * (duration * np.sinc(energies * duration / math.pi))[:, None] * applied
    return np.fft.ifft(evolved, axis=0)


def temporal_study(parameters):
    initial = initial_conditions(parameters)
    rows = []
    reference = integrate(parameters, tolerance=1e-12, initial=initial)
    for tolerance in (1e-4, 1e-7, 1e-10):
        solution = integrate(parameters, tolerance=tolerance, initial=initial)
        summary = trajectory_summary(solution, parameters)
        state, field, momentum = unpack_state(solution.y[:, -1], parameters)
        ref_state, ref_field, ref_momentum = unpack_state(reference.y[:, -1], parameters)
        rows.append({
            "relative_tolerance": tolerance,
            "absolute_tolerance": tolerance * 0.01,
            "final_state_l2_difference_from_tight_reference": float(np.linalg.norm(state - ref_state)),
            "final_field_l2_difference_from_tight_reference": float(math.sqrt(parameters.spacing) * np.linalg.norm(field - ref_field)),
            "final_momentum_l2_difference_from_tight_reference": float(math.sqrt(parameters.spacing) * np.linalg.norm(momentum - ref_momentum)),
            "max_sampled_absolute_total_energy_drift": summary["max_sampled_absolute_total_energy_drift"],
            "max_sampled_norm_squared_error": summary["max_sampled_norm_squared_error"],
            "rhs_evaluations": summary["rhs_evaluations"],
        })
    independent = integrate(parameters, tolerance=1e-10, method="RK45", initial=initial)
    independent_state, independent_field, _ = unpack_state(independent.y[:, -1], parameters)
    reference_state, reference_field, _ = unpack_state(reference.y[:, -1], parameters)
    return {
        "duration": 6.0,
        "method": "SciPy solve_ivp DOP853, real coordinates for the complex spinor and classical fields.",
        "reference_relative_tolerance": 1e-12,
        "rows": rows,
        "independent_RK45_state_l2_difference": float(np.linalg.norm(independent_state - reference_state)),
        "independent_RK45_field_l2_difference": float(math.sqrt(parameters.spacing) * np.linalg.norm(independent_field - reference_field)),
        "interpretation": "Finite-time consistency and solver comparisons, not exact discrete conservation or infinite-time stability.",
    }


def spatial_study():
    resolutions = (64, 128, 256, 512)
    final_states = []
    for sites in resolutions:
        parameters = Parameters(sites=sites)
        solution = integrate(parameters, duration=4, tolerance=1e-11, samples=41)
        final_states.append(unpack_state(solution.y[:, -1], parameters))
    rows = []
    for index, sites in enumerate(resolutions[:-1]):
        spacing = Parameters(sites=sites).spacing
        state, field, momentum = final_states[index]
        finer_state, finer_field, finer_momentum = final_states[index + 1]
        rows.append({
            "coarse_sites": sites,
            "fine_sites": resolutions[index + 1],
            "coarse_spacing": spacing,
            "state_l2_difference": float(np.linalg.norm(state - math.sqrt(2) * finer_state[::2])),
            "field_l2_difference": float(math.sqrt(spacing) * np.linalg.norm(field - finer_field[::2])),
            "momentum_l2_difference": float(math.sqrt(spacing) * np.linalg.norm(momentum - finer_momentum[::2])),
        })
    return {
        "duration": 4.0,
        "resolutions": resolutions,
        "rows": rows,
        "state_self_convergence_order": float(np.polyfit(np.log([row["coarse_spacing"] for row in rows]), np.log([row["state_l2_difference"] for row in rows]), 1)[0]),
        "field_self_convergence_order": float(np.polyfit(np.log([row["coarse_spacing"] for row in rows]), np.log([row["field_l2_difference"] for row in rows]), 1)[0]),
        "interpretation": "Smooth-packet mesh self-convergence over fixed time; no proof of a full ultraviolet-complete field theory. sqrt(2) rescales the discrete spinor when restricting the finer mesh.",
    }


def finite_time_results():
    parameters = Parameters()
    initial = initial_conditions(parameters)
    coupled = integrate(parameters, tolerance=1e-10, initial=initial)
    frozen = integrate(parameters, tolerance=1e-10, initial=initial, feedback=False)
    coupled_state = unpack_state(coupled.y[:, -1], parameters)[0]
    frozen_state = unpack_state(frozen.y[:, -1], parameters)[0]
    excited_initial = (
        initial[0],
        0.1 * np.cos(2 * math.pi * parameters.positions / parameters.length),
        0.15 + 0.02 * np.sin(2 * math.pi * parameters.positions / parameters.length),
    )
    excited_reciprocal = integrate(parameters, duration=3, tolerance=1e-10, initial=excited_initial)
    excited_one_way = integrate(parameters, duration=3, tolerance=1e-10, initial=excited_initial, feedback=False)
    return {
        "parameters": asdict(parameters),
        "initial_conditions": "q=pi=0; normalized Gaussian packet centered at -4, width 1.2, carrier k=1, positive free lattice mode spinor at that carrier.",
        "coupled": trajectory_summary(coupled, parameters),
        "zero_background_one_way_control": trajectory_summary(frozen, parameters),
        "coupled_vs_control_density_total_variation": float(0.5 * np.sum(abs(np.sum(abs(coupled_state)**2, axis=1) - np.sum(abs(frozen_state)**2, axis=1)))),
        "coupled_vs_control_state_l2_difference": float(np.linalg.norm(coupled_state - frozen_state)),
        "initially_excited_background_control": {
            "duration": 3.0,
            "field": "0.1*cos(2*pi*x/L)",
            "momentum": "0.15+0.02*sin(2*pi*x/L)",
            "reciprocal": trajectory_summary(excited_reciprocal, parameters),
            "one_way_without_force": trajectory_summary(excited_one_way, parameters),
        },
    }


def source_and_cutoff_diagnostics():
    parameters = Parameters()
    rest_state = np.ones((parameters.sites, 2), dtype=complex) / math.sqrt(2 * parameters.sites)
    zero = np.zeros(parameters.sites)
    rows = []
    for sites in (64, 128, 256, 512):
        current = replace(parameters, sites=sites)
        wave_numbers = 2 * math.pi * np.fft.fftfreq(sites, d=current.spacing)
        positive_energies = np.sqrt(
            (current.base_velocity * np.sin(wave_numbers * current.spacing) / current.spacing)**2
            + current.matter_mass**2
        )
        rows.append({
            "sites": sites,
            "spacing": current.spacing,
            "minimum_free_matter_energy": -float(positive_energies.max()),
            "energy_at_zero_momentum": current.matter_mass,
            "energy_at_brillouin_edge": float(math.hypot(
                current.base_velocity * math.sin(math.pi) / current.spacing, current.matter_mass)),
        })
    return {
        "rest_state_matter_energy": energy_components(rest_state, zero, zero, parameters)["matter"],
        "rest_state_max_field_source": float(np.max(abs(field_energy_gradient(rest_state, zero, parameters) / parameters.spacing))),
        "free_dispersion": "E_plus_minus(k)=+/-sqrt(v0^2*sin(k*a)^2/a^2+m^2), hbar=1.",
        "cutoff_rows": rows,
        "finite_grid_lower_bound": "At unit spinor norm, total energy >= -(v_max/a+abs(m)). This bound is not uniform as a tends to zero.",
        "background_principal_speed": parameters.background_speed,
        "matter_principal_speed_range": [parameters.base_velocity - parameters.velocity_amplitude,
                                         parameters.base_velocity + parameters.velocity_amplitude],
        "source_interpretation": "The source differentiates the kinetic/velocity coupling. It is not T00 and cannot simply be called a universal Newtonian potential source. Stress-like couplings can legitimately omit rest energy; a full covariant set of gravitational constraints has not been supplied.",
        "ultraviolet_interpretation": "The first-quantized Dirac branch is unbounded below as the cutoff is removed; this is not alone a proof of dynamical instability, but no second-quantized vacuum or regulator-independent completion has been built. The centered derivative also has an extra low-energy branch at the zone edge.",
    }


def ensemble_consistency_diagnostics():
    parameters = Parameters(sites=32, length=16)
    wave_number = 2 * math.pi / parameters.length
    positive = np.zeros((parameters.sites, 2), dtype=complex)
    positive[:, 0] = np.exp(1j * wave_number * parameters.positions) / math.sqrt(parameters.sites)
    negative = positive.conj()
    traveling = (positive, negative)
    standing = ((positive + negative) / math.sqrt(2), (positive - negative) / math.sqrt(2))
    zero = np.zeros(parameters.sites)
    density_matrices = []
    rows = []
    for name, states in (("traveling_wave_ensemble", traveling), ("standing_wave_ensemble", standing)):
        density_matrices.append(sum(np.outer(state.ravel(), state.ravel().conj()) for state in states) / 2)
        accelerations = [float(np.mean(derivatives(state, zero, zero, parameters)[2])) for state in states]
        duration = 0.02
        final_field_means = [float(np.mean(unpack_state(
            integrate(parameters, duration=duration, initial=(state, zero, zero), tolerance=1e-11, samples=3).y[:, -1],
            parameters)[1])) for state in states]
        rows.append({
            "ensemble": name,
            "weights": [0.5, 0.5],
            "initial_spatial_mean_field_accelerations": accelerations,
            "mean_initial_acceleration": float(np.mean(accelerations)),
            "variance_coefficient_of_mean_field_at_order_t4": float(np.var(accelerations) / 4),
            "finite_duration": duration,
            "final_spatial_mean_fields_by_branch": final_field_means,
            "finite_time_field_variance_divided_by_t4": float(np.var(final_field_means) / duration**4),
        })
    return {
        "setup": "The same mixed spinor density matrix and identical deterministic initial classical fields q=pi=0, represented by two equal-weight pure-state ensembles.",
        "parameters": asdict(parameters),
        "wave_number": wave_number,
        "initial_quantum_density_matrix_frobenius_difference": float(np.linalg.norm(density_matrices[0] - density_matrices[1])),
        "rows": rows,
        "interpretation_tested": "Each pure-state preparation branch generates its own deterministic mean-field trajectory, then classical outputs are statistically averaged after labels are discarded.",
        "conclusion_scope": "This branchwise extension does not preserve preparation equivalence: classical field variance depends on the pure-state decomposition. It disqualifies that extension as a fundamental density-matrix-based hybrid law, not the use of mean-field trajectories as a preparation-specific approximation. No universal no-go theorem for all hybrid theories or explicit signalling protocol is claimed.",
    }


class BackreactionTests(unittest.TestCase):
    def setUp(self):
        self.parameters = Parameters(sites=24, length=12)
        self.generator = np.random.default_rng(505)
        self.state = self.generator.normal(size=(24, 2)) + 1j * self.generator.normal(size=(24, 2))
        self.state /= np.linalg.norm(self.state)
        self.field = self.generator.normal(scale=0.2, size=24)
        self.momentum = self.generator.normal(scale=0.1, size=24)

    def test_real_coordinate_packing_round_trip(self):
        unpacked = unpack_state(pack_state(self.state, self.field, self.momentum), self.parameters)
        for original, restored in zip((self.state, self.field, self.momentum), unpacked):
            np.testing.assert_array_equal(restored, original)

    def test_uncoupled_integration_matches_independent_fourier_solution(self):
        parameters = Parameters(sites=48, length=24, coupling=0)
        initial = initial_conditions(parameters)
        solution = integrate(parameters, duration=2, tolerance=1e-11, initial=initial, samples=21)
        state, field, momentum = unpack_state(solution.y[:, -1], parameters)
        expected = free_lattice_reference(initial[0], 2, parameters)
        np.testing.assert_allclose(state, expected, atol=2e-9)
        np.testing.assert_allclose(field, 0, atol=1e-12)
        np.testing.assert_allclose(momentum, 0, atol=1e-12)

    def test_coupled_evolution_generates_response_and_exchanges_energy(self):
        parameters = Parameters(sites=64)
        solution = integrate(parameters, duration=4, tolerance=1e-10, samples=41)
        summary = trajectory_summary(solution, parameters)
        self.assertGreater(summary["field_max_abs_at_final"], 0.001)
        self.assertGreater(summary["final_energy"]["background"], 1e-5)
        self.assertLess(summary["max_sampled_absolute_total_energy_drift"], 1e-7)
        self.assertLess(summary["max_sampled_norm_squared_error"], 1e-7)

    def test_two_integrators_agree_over_fixed_time(self):
        parameters = Parameters(sites=32, length=16)
        high_order = integrate(parameters, duration=2, tolerance=1e-10, samples=21)
        lower_order = integrate(parameters, duration=2, tolerance=1e-10, method="RK45", samples=21)
        np.testing.assert_allclose(high_order.y[:, -1], lower_order.y[:, -1], atol=1e-8)

    def test_tighter_integration_reduces_sampled_energy_drift(self):
        parameters = Parameters(sites=64)
        loose = trajectory_summary(integrate(parameters, duration=4, tolerance=1e-4, samples=41), parameters)
        tight = trajectory_summary(integrate(parameters, duration=4, tolerance=1e-10, samples=41), parameters)
        self.assertLess(tight["max_sampled_absolute_total_energy_drift"], loose["max_sampled_absolute_total_energy_drift"])
        self.assertLess(tight["max_sampled_norm_squared_error"], loose["max_sampled_norm_squared_error"])

    def test_nonzero_rest_energy_does_not_source_this_velocity_field(self):
        parameters = Parameters(sites=24, length=12)
        state = np.ones((parameters.sites, 2), dtype=complex) / math.sqrt(2 * parameters.sites)
        field = np.zeros(parameters.sites)
        momentum = np.zeros(parameters.sites)
        self.assertAlmostEqual(energy_components(state, field, momentum, parameters)["matter"], parameters.matter_mass)
        np.testing.assert_allclose(derivatives(state, field, momentum, parameters)[2], 0, atol=1e-12)
        np.testing.assert_allclose(derivatives(state, field, momentum, parameters)[0], -1j * parameters.matter_mass * state, atol=1e-12)

    def test_source_is_not_determined_by_position_probability(self):
        parameters = Parameters(sites=32, length=16)
        wave_number = 2 * math.pi / parameters.length
        state = np.zeros((parameters.sites, 2), dtype=complex)
        state[:, 0] = np.exp(1j * wave_number * parameters.positions) / math.sqrt(parameters.sites)
        reversed_state = state.conj()
        np.testing.assert_allclose(abs(state)**2, abs(reversed_state)**2, atol=1e-12)
        positive_source = field_energy_gradient(state, np.zeros(parameters.sites), parameters)
        negative_source = field_energy_gradient(reversed_state, np.zeros(parameters.sites), parameters)
        np.testing.assert_allclose(positive_source, -negative_source, atol=1e-12)
        self.assertGreater(float(np.linalg.norm(positive_source)), 0.001)

    def test_matter_energy_has_a_fixed_lattice_bound(self):
        energy = energy_components(self.state, self.field, self.momentum, self.parameters)["matter"]
        bound = (self.parameters.base_velocity + self.parameters.velocity_amplitude) / self.parameters.spacing + abs(self.parameters.matter_mass)
        self.assertLessEqual(abs(energy), bound)

    def test_free_spectrum_matches_independent_dense_diagonalization(self):
        parameters = Parameters(sites=12, length=6)
        dimension = 2 * parameters.sites
        matrix = np.column_stack([
            apply_matter_hamiltonian(np.eye(dimension, dtype=complex)[:, index].reshape(parameters.sites, 2),
                                      np.zeros(parameters.sites), parameters).ravel()
            for index in range(dimension)
        ])
        wave_numbers = 2 * math.pi * np.fft.fftfreq(parameters.sites, d=parameters.spacing)
        energies = np.sqrt((parameters.base_velocity * np.sin(wave_numbers * parameters.spacing) / parameters.spacing)**2
                           + parameters.matter_mass**2)
        np.testing.assert_allclose(np.linalg.eigvalsh(matrix), np.sort(np.r_[-energies, energies]), atol=1e-12)

    def test_negative_energy_and_zone_edge_limitations(self):
        diagnostics = source_and_cutoff_diagnostics()
        rows = diagnostics["cutoff_rows"]
        self.assertTrue(all(later["minimum_free_matter_energy"] < earlier["minimum_free_matter_energy"]
                            for earlier, later in zip(rows, rows[1:])))
        self.assertGreater(abs(rows[-1]["minimum_free_matter_energy"]), 10)
        for row in rows:
            self.assertAlmostEqual(row["energy_at_zero_momentum"], row["energy_at_brillouin_edge"], places=12)
        self.assertAlmostEqual(diagnostics["rest_state_matter_energy"], 0.6)
        self.assertEqual(diagnostics["rest_state_max_field_source"], 0)

    def test_equivalent_ensembles_have_distinct_branchwise_field_statistics(self):
        diagnostics = ensemble_consistency_diagnostics()
        self.assertLess(diagnostics["initial_quantum_density_matrix_frobenius_difference"], 1e-12)
        traveling, standing = diagnostics["rows"]
        self.assertAlmostEqual(traveling["mean_initial_acceleration"], standing["mean_initial_acceleration"], places=12)
        self.assertGreater(traveling["variance_coefficient_of_mean_field_at_order_t4"], 1e-6)
        self.assertLess(standing["variance_coefficient_of_mean_field_at_order_t4"], 1e-24)
        self.assertLess(standing["finite_time_field_variance_divided_by_t4"], 1e-20)

    def test_ensemble_variance_matches_analytic_small_time_prediction(self):
        diagnostics = ensemble_consistency_diagnostics()
        parameters = Parameters(sites=32, length=16)
        wave_number = diagnostics["wave_number"]
        momentum = math.sin(wave_number * parameters.spacing) / parameters.spacing
        acceleration_magnitude = parameters.velocity_amplitude * parameters.coupling * momentum / parameters.length
        expected_coefficient = acceleration_magnitude**2 / 4
        traveling = diagnostics["rows"][0]
        self.assertAlmostEqual(traveling["variance_coefficient_of_mean_field_at_order_t4"], expected_coefficient, places=14)
        self.assertAlmostEqual(traveling["finite_time_field_variance_divided_by_t4"] / expected_coefficient, 1, delta=0.005)

    def test_matter_hamiltonian_is_hermitian_for_variable_field(self):
        other = self.generator.normal(size=(24, 2)) + 1j * self.generator.normal(size=(24, 2))
        first = np.vdot(self.state, apply_matter_hamiltonian(other, self.field, self.parameters))
        second = np.vdot(apply_matter_hamiltonian(self.state, self.field, self.parameters), other)
        self.assertAlmostEqual(abs(first - second), 0, places=12)

    def test_kinetic_density_matches_independent_energy_expression(self):
        mass_energy = self.parameters.matter_mass * np.vdot(self.state, self.state @ SIGMA_X.T).real
        expected = np.dot(velocity(self.field, self.parameters), kinetic_density(self.state, self.parameters)) + mass_energy
        actual = energy_components(self.state, self.field, self.momentum, self.parameters)["matter"]
        self.assertAlmostEqual(actual, expected, places=12)

    def test_feedback_force_matches_field_energy_finite_differences(self):
        step = 1e-6
        analytic = field_energy_gradient(self.state, self.field, self.parameters)
        for index in range(self.parameters.sites):
            variation = np.zeros(self.parameters.sites)
            variation[index] = step
            positive = energy_components(self.state, self.field + variation, self.momentum, self.parameters)["matter"]
            negative = energy_components(self.state, self.field - variation, self.momentum, self.parameters)["matter"]
            self.assertAlmostEqual((positive - negative) / (2 * step), analytic[index], delta=1e-9)

    def test_total_energy_and_norm_rates_vanish(self):
        rates = energy_rates(self.state, self.field, self.momentum, self.parameters)
        self.assertAlmostEqual(rates["total_rate"], 0, places=12)
        state_derivative = derivatives(self.state, self.field, self.momentum, self.parameters)[0]
        self.assertAlmostEqual(float(2 * np.vdot(self.state, state_derivative).real), 0, places=12)

    def test_energy_derivative_checked_by_full_directional_difference(self):
        step = 1e-5
        state_rate, field_rate, momentum_rate = derivatives(self.state, self.field, self.momentum, self.parameters)
        positive = energy_components(self.state + step * state_rate, self.field + step * field_rate,
                                     self.momentum + step * momentum_rate, self.parameters)["total"]
        negative = energy_components(self.state - step * state_rate, self.field - step * field_rate,
                                     self.momentum - step * momentum_rate, self.parameters)["total"]
        self.assertAlmostEqual((positive - negative) / (2 * step), 0, delta=1e-8)

    def test_one_way_control_generally_breaks_total_energy_balance(self):
        results = algebra_diagnostics()
        self.assertGreater(abs(results["one_way_control_at_test_state"]["total_rate"]), 0.001)
        self.assertLess(abs(results["energy_exchange_at_test_state"]["total_rate"]), 1e-12)

    def test_zero_coupling_removes_source_and_field_dependence(self):
        uncoupled = replace(self.parameters, coupling=0)
        np.testing.assert_allclose(field_energy_gradient(self.state, self.field, uncoupled), 0, atol=1e-12)
        np.testing.assert_allclose(apply_matter_hamiltonian(self.state, self.field, uncoupled),
                                   apply_matter_hamiltonian(self.state, np.zeros(24), uncoupled), atol=1e-12)

    def test_finite_difference_background_force_matches_energy_gradient(self):
        step = 1e-6
        _, _, acceleration = derivatives(self.state, self.field, self.momentum, self.parameters)
        for index in (0, 1, 8, 23):
            variation = np.zeros(24)
            variation[index] = step
            positive = energy_components(self.state, self.field + variation, self.momentum, self.parameters)["total"]
            negative = energy_components(self.state, self.field - variation, self.momentum, self.parameters)["total"]
            self.assertAlmostEqual((positive - negative) / (2 * step), -self.parameters.spacing * acceleration[index], delta=1e-8)

    def test_velocity_is_bounded_and_derivative_matches_difference(self):
        fields = np.linspace(-30, 30, 201)
        speeds = velocity(fields, self.parameters)
        self.assertGreaterEqual(float(speeds.min()), self.parameters.base_velocity - self.parameters.velocity_amplitude - 1e-12)
        self.assertLessEqual(float(speeds.max()), self.parameters.base_velocity + self.parameters.velocity_amplitude + 1e-12)
        step = 1e-5
        numeric = (velocity(self.field + step, self.parameters) - velocity(self.field - step, self.parameters)) / (2 * step)
        np.testing.assert_allclose(numeric, velocity_derivative(self.field, self.parameters), atol=1e-9)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(BackreactionTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    results = {
        "scope": "A proposed finite Hamiltonian feedback toy model, not a derivation of gravity or new empirical evidence.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "algebra": algebra_diagnostics(),
        "finite_time_evolution": finite_time_results(),
        "time_integration_checks": temporal_study(Parameters()),
        "spatial_self_convergence": spatial_study(),
        "source_and_cutoff_diagnostics": source_and_cutoff_diagnostics(),
        "ensemble_consistency_diagnostics": ensemble_consistency_diagnostics(),
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    encoded = json.dumps(results, indent=2)
    if arguments.write_results:
        output = Path(__file__).resolve().parent / "backreaction_results.json"
        output.write_text(encoded + "\n", encoding="utf-8")
        print(f"Results written to {output}")
    evolution = results["finite_time_evolution"]
    compact = {
        "automated_checks": results["automated_checks"],
        "algebra": results["algebra"],
        "coupled_evolution": {key: value for key, value in evolution["coupled"].items() if key not in ("energy_samples", "energy_sample_times")},
        "coupled_vs_control_density_total_variation": evolution["coupled_vs_control_density_total_variation"],
        "excited_one_way_total_energy_drift": evolution["initially_excited_background_control"]["one_way_without_force"]["max_sampled_absolute_total_energy_drift"],
        "time_integration_checks": results["time_integration_checks"],
        "spatial_self_convergence": results["spatial_self_convergence"],
        "source_and_cutoff_diagnostics": results["source_and_cutoff_diagnostics"],
        "ensemble_consistency_diagnostics": results["ensemble_consistency_diagnostics"],
        "scope": results["scope"],
    }
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()