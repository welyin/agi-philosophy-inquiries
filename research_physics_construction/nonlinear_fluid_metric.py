"""Construction 10: static spherical fluid and nonlinear metric closure.

Einstein-Hilbert dynamics and perfect-fluid matter are declared inputs.
The constant-density solution is an analytic benchmark, not a causal EOS.
"""

import argparse
import json
import math
import unittest
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from weak_metric_response import (nonlinear_einstein_from_jet, probe_unitary, probe_record_probabilities,
                                  probe_schedule, release_probe_record, RingMetric)
from massive_envelope import massive_hopping, postprocessing_table
from quantum_chain_background import initial_data
from baseline import real_lift
from common_orientation_structure import encode_state


@dataclass(frozen=True)
class UniformFluidSphere:
    mass: float = 1.0
    radius: float = 3.0
    gravitational_constant: float = 0.1

    def __post_init__(self):
        if (not all(np.isfinite(value) and value > 0 for value in (self.mass, self.radius, self.gravitational_constant))
                or self.compactness >= 8 / 9):
            raise ValueError("Use finite positive parameters and compactness below the regular-center benchmark limit.")

    @property
    def compactness(self):
        return 2 * self.gravitational_constant * self.mass / self.radius

    @property
    def density(self):
        return 3 * self.mass / (4 * math.pi * self.radius ** 3)

    def radial(self, radius):
        if not np.isfinite(radius) or radius < 0:
            raise ValueError("Use a finite nonnegative areal radius.")
        strength = self.gravitational_constant * self.mass
        if radius >= self.radius:
            lapse = math.sqrt(1 - 2 * strength / radius)
            first = strength / (radius ** 2 * lapse)
            second = -2 * strength / (radius ** 3 * lapse) - strength ** 2 / (radius ** 4 * lapse ** 3)
            factor = 1 / lapse ** 2
            factor_first = -2 * strength / (radius - 2 * strength) ** 2
            factor_second = 4 * strength / (radius - 2 * strength) ** 3
            return dict(mass=self.mass, density=0., pressure=0., pressure_first=0.,
                        lapse=lapse, lapse_first=first, lapse_second=second,
                        radial_factor=factor, radial_factor_first=factor_first, radial_factor_second=factor_second)
        coefficient = 2 * strength / self.radius ** 3
        surface = math.sqrt(1 - self.compactness)
        local = math.sqrt(1 - coefficient * radius ** 2)
        lapse = (3 * surface - local) / 2
        first = coefficient * radius / (2 * local)
        second = coefficient / (2 * local ** 3)
        pressure = self.density * (local - surface) / (3 * surface - local)
        pressure_first = -(self.density + pressure) * first / lapse
        factor = 1 / local ** 2
        factor_first = 2 * coefficient * radius / local ** 4
        factor_second = 2 * coefficient / local ** 4 + 8 * coefficient ** 2 * radius ** 2 / local ** 6
        return dict(mass=self.mass * (radius / self.radius) ** 3, density=self.density, pressure=pressure,
                    pressure_first=pressure_first, lapse=lapse, lapse_first=first, lapse_second=second,
                    radial_factor=factor, radial_factor_first=factor_first, radial_factor_second=factor_second)


def spherical_metric_jet(radius, angle, radial):
    if radius <= 0 or not 0 < angle < math.pi:
        raise ValueError("Evaluate the spherical coordinate chart away from its center and polar axes.")
    lapse, lapse_first, lapse_second = (radial[key] for key in ("lapse", "lapse_first", "lapse_second"))
    factor, factor_first, factor_second = (radial[key] for key in ("radial_factor", "radial_factor_first", "radial_factor_second"))
    sine, cosine = math.sin(angle), math.cos(angle)
    metric = np.diag([lapse ** 2, -factor, -radius ** 2, -radius ** 2 * sine ** 2])
    first = np.zeros((4, 4, 4))
    second = np.zeros((4, 4, 4, 4))
    first[1] = np.diag([2 * lapse * lapse_first, -factor_first, -2 * radius, -2 * radius * sine ** 2])
    first[2, 3, 3] = -2 * radius ** 2 * sine * cosine
    second[1, 1] = np.diag([2 * (lapse_first ** 2 + lapse * lapse_second), -factor_second, -2., -2 * sine ** 2])
    second[1, 2, 3, 3] = second[2, 1, 3, 3] = -4 * radius * sine * cosine
    second[2, 2, 3, 3] = -2 * radius ** 2 * (cosine ** 2 - sine ** 2)
    return metric, first, second


def fluid_covariant_tensor(metric, density, pressure):
    return np.diag([density * metric[0, 0], *(-pressure * np.diag(metric)[1:])])


@dataclass(frozen=True)
class LinearFluidEOS:
    surface_density: float = 0.01
    sound_speed_squared: float = 1 / 3

    def __post_init__(self):
        if (not np.isfinite(self.surface_density) or self.surface_density <= 0
                or not np.isfinite(self.sound_speed_squared) or not 0 < self.sound_speed_squared <= 1):
            raise ValueError("Use positive surface density and a causal positive barotropic sound-speed squared.")

    def density(self, pressure):
        return self.surface_density + pressure / self.sound_speed_squared

    @property
    def density_derivative(self):
        return 1 / self.sound_speed_squared

    def particle_density(self, pressure, surface_particle_mass=1.):
        if not np.isfinite(surface_particle_mass) or surface_particle_mass <= 0 or np.any(np.asarray(pressure) < 0):
            raise ValueError("Use nonnegative pressure and a positive particle-mass normalization.")
        sound = self.sound_speed_squared
        factor = 1 + (1 + sound) * np.asarray(pressure) / (sound * self.surface_density)
        return self.surface_density / surface_particle_mass * factor ** (1 / (1 + sound))

    def density_from_particles(self, number_density, surface_particle_mass=1.):
        if surface_particle_mass <= 0 or np.any(np.asarray(number_density) <= 0):
            raise ValueError("Use positive particle density and mass normalization.")
        sound = self.sound_speed_squared
        ratio = np.asarray(number_density) * surface_particle_mass / self.surface_density
        return self.surface_density / (1 + sound) * (sound + ratio ** (1 + sound))


@dataclass(frozen=True)
class ConstantDensityEOS:
    surface_density: float

    def density(self, pressure):
        return self.surface_density

    @property
    def density_derivative(self):
        return 0.


class FluidSolution:
    def __init__(self, eos=LinearFluidEOS(), central_pressure=0.0005, gravitational_constant=0.1,
                 relative_tolerance=1e-10, start_fraction=1e-6):
        if (not all(np.isfinite(value) and value > 0 for value in
                    (central_pressure, gravitational_constant, relative_tolerance, start_fraction))
                or relative_tolerance > 1e-3 or start_fraction > 1e-2):
            raise ValueError("Use positive finite physical parameters and controlled solver/center settings.")
        self.eos, self.central_pressure, self.gravitational_constant = eos, central_pressure, gravitational_constant
        self.relative_tolerance, self.start_fraction = relative_tolerance, start_fraction
        self.central_density = float(eos.density(central_pressure))
        if not np.isfinite(self.central_density) or self.central_density <= 0:
            raise ValueError("The initial density must be finite and positive.")
        scale = 1 / math.sqrt(4 * math.pi * gravitational_constant * self.central_density)
        self.start_radius = scale * start_fraction
        pressure_curvature = 2 * math.pi * gravitational_constant * (self.central_density + central_pressure) * (self.central_density / 3 + central_pressure)
        potential_curvature = 2 * math.pi * gravitational_constant * (self.central_density / 3 + central_pressure)
        self.center_pressure_coefficient, self.center_potential_coefficient = pressure_curvature, potential_curvature
        radius = self.start_radius
        initial = [4 * math.pi * self.central_density * radius ** 3 / 3,
                   central_pressure - pressure_curvature * radius ** 2, potential_curvature * radius ** 2]

        def surface_event(radius, state):
            return state[1]

        surface_event.terminal = True
        surface_event.direction = -1
        self.integration = solve_ivp(self.derivative, (radius, 30 * scale), initial, method="DOP853",
                                     rtol=relative_tolerance, atol=relative_tolerance * np.array([1e-2, 1e-5, 1e-2]),
                                     events=surface_event, dense_output=True, max_step=scale / 4)
        if not self.integration.success or len(self.integration.t_events[0]) != 1:
            raise ValueError("The declared integration did not reach a regular zero-pressure surface.")
        self.radius = float(self.integration.t_events[0][0])
        self.mass, _, surface_potential = self.integration.y_events[0][0]
        self.mass = float(self.mass)
        self.compactness = 2 * gravitational_constant * self.mass / self.radius
        if not 0 < self.compactness < 1:
            raise ValueError("A regular exterior matching surface is required.")
        self.potential_shift = math.log(1 - self.compactness) / 2 - surface_potential

    def derivative(self, radius, state):
        mass, pressure, _ = state
        density = self.eos.density(pressure)
        denominator = radius * (radius - 2 * self.gravitational_constant * mass)
        if denominator <= 0:
            raise ValueError("The static integration encountered a trapped-surface denominator.")
        potential_first = self.gravitational_constant * (mass + 4 * math.pi * radius ** 3 * pressure) / denominator
        return np.array([4 * math.pi * radius ** 2 * density,
                         -(density + pressure) * potential_first, potential_first])

    def radial(self, radius):
        if not np.isfinite(radius) or radius < 0:
            raise ValueError("Use a finite nonnegative areal radius.")
        gravity = self.gravitational_constant
        if radius >= self.radius:
            return UniformFluidSphere(self.mass, self.radius, gravity).radial(radius)
        if radius <= self.start_radius:
            pressure = self.central_pressure - self.center_pressure_coefficient * radius ** 2
            lapse = math.exp(self.potential_shift + self.center_potential_coefficient * radius ** 2)
            first = 2 * self.center_potential_coefficient * radius
            coefficient = 8 * math.pi * gravity * self.central_density / 3
            denominator = 1 - coefficient * radius ** 2
            return dict(mass=4 * math.pi * self.central_density * radius ** 3 / 3,
                        density=float(self.eos.density(pressure)), pressure=pressure,
                        pressure_first=-2 * self.center_pressure_coefficient * radius,
                        lapse=lapse, lapse_first=lapse * first,
                        lapse_second=lapse * (2 * self.center_potential_coefficient + first ** 2),
                        radial_factor=1 / denominator, radial_factor_first=2 * coefficient * radius / denominator ** 2,
                        radial_factor_second=2 * coefficient / denominator ** 2 + 8 * coefficient ** 2 * radius ** 2 / denominator ** 3)
        mass, pressure, potential = self.integration.sol(radius)
        density = self.eos.density(pressure)
        mass_first, pressure_first, potential_first = self.derivative(radius, [mass, pressure, potential])
        mass_second = 8 * math.pi * radius * density + 4 * math.pi * radius ** 2 * self.eos.density_derivative * pressure_first
        denominator = radius ** 2 - 2 * gravity * mass * radius
        denominator_first = 2 * radius - 2 * gravity * (mass + radius * mass_first)
        numerator_first = gravity * (mass_first + 12 * math.pi * radius ** 2 * pressure + 4 * math.pi * radius ** 3 * pressure_first)
        potential_second = numerator_first / denominator - potential_first * denominator_first / denominator
        lapse = math.exp(potential + self.potential_shift)
        factor_inverse = 1 - 2 * gravity * mass / radius
        inverse_first = -2 * gravity * (mass_first / radius - mass / radius ** 2)
        inverse_second = -2 * gravity * (mass_second / radius - 2 * mass_first / radius ** 2 + 2 * mass / radius ** 3)
        return dict(mass=float(mass), density=float(density), pressure=float(pressure), pressure_first=float(pressure_first),
                    lapse=lapse, lapse_first=lapse * potential_first, lapse_second=lapse * (potential_second + potential_first ** 2),
                    radial_factor=1 / factor_inverse, radial_factor_first=-inverse_first / factor_inverse ** 2,
                    radial_factor_second=2 * inverse_first ** 2 / factor_inverse ** 3 - inverse_second / factor_inverse ** 2)


@dataclass(frozen=True)
class FluidRingMetric:
    solution: FluidSolution
    length: float = 20.0
    center_offset: float = 4.0

    def __post_init__(self):
        if not np.isfinite(self.length) or self.length <= 0 or not np.isfinite(self.center_offset) or self.center_offset < 0:
            raise ValueError("Use a positive reference circle circumference and finite nonnegative offset.")
        radius = self.length / (2 * math.pi)
        if abs(self.center_offset - radius) < 1e-8:
            raise ValueError("This coordinate implementation keeps the ring away from the areal center.")

    def geometry(self, positions):
        positions = np.asarray(positions, dtype=float)
        if not np.all(np.isfinite(positions)):
            raise ValueError("Use finite reference arc coordinates.")
        circle_radius = self.length / (2 * math.pi)
        angle = 2 * math.pi * positions / self.length
        distances = np.sqrt(self.center_offset ** 2 + circle_radius ** 2 + 2 * self.center_offset * circle_radius * np.cos(angle))
        radial_tangent = -self.center_offset * np.sin(angle) / distances
        return distances, radial_tangent

    def lapse(self, positions):
        distances, _ = self.geometry(positions)
        return np.array([self.solution.radial(float(radius))["lapse"] for radius in distances.flat]).reshape(distances.shape)

    def spatial_scale(self, positions):
        distances, radial_tangent = self.geometry(positions)
        factors = np.array([self.solution.radial(float(radius))["radial_factor"] for radius in distances.flat]).reshape(distances.shape)
        return np.sqrt(1 + (factors - 1) * radial_tangent ** 2)

    def speed(self, positions):
        return self.lapse(positions) / self.spatial_scale(positions)

    @property
    def conservative_speed_floor(self):
        gravity, density, mass = self.solution.gravitational_constant, self.solution.central_density, self.solution.mass
        maximum_compactness = 2 * gravity * (4 * math.pi * density / 3) ** (1 / 3) * mass ** (2 / 3)
        if maximum_compactness >= 1:
            raise ValueError("The density/mass speed bound is insufficient for this probe schedule.")
        return self.solution.radial(0.)["lapse"] * math.sqrt(1 - maximum_compactness)


def closure_diagnostics(solution):
    residuals, balance = [], []
    for fraction in (0.05, 0.2, 0.5, 0.9, 1.1, 2.):
        radius = solution.radius * fraction
        radial = solution.radial(radius)
        jet = spherical_metric_jet(radius, 1.1, radial)
        expected = 8 * math.pi * solution.gravitational_constant * fluid_covariant_tensor(jet[0], radial["density"], radial["pressure"])
        residuals.append(float(np.linalg.norm(nonlinear_einstein_from_jet(*jet) - expected)))
        balance.append(abs(radial["pressure_first"] + (radial["density"] + radial["pressure"]) * radial["lapse_first"] / radial["lapse"]))
    derivative_residuals = []
    for fraction in (0.2, 0.5, 0.8):
        radius = solution.radius * fraction
        step = solution.radius * 1e-4
        numerical = (solution.integration.sol(radius - 2 * step) - 8 * solution.integration.sol(radius - step)
                     + 8 * solution.integration.sol(radius + step) - solution.integration.sol(radius + 2 * step)) / (12 * step)
        exact_rhs = solution.derivative(radius, solution.integration.sol(radius))
        derivative_residuals.append(float(np.max(abs(numerical - exact_rhs))))
    return {"maximum_sampled_full_Einstein_residual": max(residuals),
            "maximum_sampled_covariant_pressure_balance_residual": max(balance),
            "maximum_dense_output_five_point_ODE_residual": max(derivative_residuals),
            "scope": "curvature jets use TOV derivatives; independent dense-output derivative and refinement checks test the numerical solution, not a certified global error"}


def run_report():
    solution = FluidSolution()
    metric = FluidRingMetric(solution)
    sites, duration = 32, 2.
    prepared = np.kron(initial_data(sites), np.ones(4) / 2)
    trials = []
    for mass in (0.35, 0.55):
        gate, hamiltonian, rates = probe_unitary(sites, mass, duration, metric)
        state = gate @ prepared
        records = probe_record_probabilities(state, sites)
        flat_gate, _, _ = probe_unitary(sites, mass, duration, RingMetric(source=None))
        flat_records = probe_record_probabilities(flat_gate @ prepared, sites)
        trials.append({"proper_mass": mass, "joint_record_probabilities": records.tolist(),
                       "smooth_data_probability": float(records[0].sum()),
                       "clock_A_X_plus_probability": float(records[:, 0, :].sum()),
                       "clock_C_Y_plus_probability": float(records[:, :, 0].sum()),
                       "data_probability_change_from_flat": float(records[0].sum() - flat_records[0].sum()),
                       "full_record_TV_from_flat": float(np.abs(records - flat_records).sum() / 2),
                       "initial_data_energy": float(np.vdot(prepared, np.kron(hamiltonian, np.eye(4)) @ prepared).real),
                       "final_data_energy": float(np.vdot(state, np.kron(hamiltonian, np.eye(4)) @ state).real),
                       "clock_rates": rates.tolist()})
    refinements = []
    reference = FluidSolution(relative_tolerance=1e-12, start_fraction=5e-7)
    reference_gate, _, _ = probe_unitary(sites, 0.35, duration, FluidRingMetric(reference))
    reference_records = probe_record_probabilities(reference_gate @ prepared, sites)
    for tolerance, start in ((1e-7, 1e-6), (1e-9, 1e-6), (1e-10, 1e-6), (1e-10, 5e-7)):
        candidate = FluidSolution(relative_tolerance=tolerance, start_fraction=start)
        gate, _, _ = probe_unitary(sites, 0.35, duration, FluidRingMetric(candidate))
        records = probe_record_probabilities(gate @ prepared, sites)
        refinements.append({"rtol": tolerance, "start_fraction": start, "radius": candidate.radius, "mass": candidate.mass,
                            "radius_difference_from_refined": abs(candidate.radius - reference.radius),
                            "mass_difference_from_refined": abs(candidate.mass - reference.mass),
                            "probe_record_TV_from_refined": float(abs(records - reference_records).sum() / 2)})
    nodes, weights = np.polynomial.legendre.leggauss(128)
    radii = solution.radius * (nodes + 1) / 2
    rows = [solution.radial(float(radius)) for radius in radii]
    proper_energy = solution.radius / 2 * np.dot(weights, [4 * math.pi * radius ** 2 * row["density"] * math.sqrt(row["radial_factor"]) for radius, row in zip(radii, rows)])
    baryon_rest_mass = solution.radius / 2 * np.dot(weights, [4 * math.pi * radius ** 2 * solution.eos.particle_density(row["pressure"]) * math.sqrt(row["radial_factor"]) for radius, row in zip(radii, rows)])
    chemical_redshifts = [row["lapse"] * (row["density"] + row["pressure"]) / solution.eos.particle_density(row["pressure"]) for row in rows]
    schedule = probe_schedule(sites, duration, metric)
    clocks = metric.lapse([0., metric.length / 2])
    return {"model_version": "C10", "route_round": "10", "research_scope": "sufficiency only",
            "action_input": "Einstein-Hilbert with the C9 curvature/sign convention, plus isotropic perfect-fluid matter; G_munu=8piG T_munu is adopted nonlinear dynamics",
            "EOS_input": vars(solution.eos), "G": solution.gravitational_constant, "central_pressure": solution.central_pressure,
            "solution": {"radius": solution.radius, "mass": solution.mass, "compactness": solution.compactness,
                         "central_density": solution.central_density, "central_lapse": solution.radial(0.)["lapse"],
                         "surface_lapse": solution.radial(solution.radius)["lapse"], "proper_fluid_energy": float(proper_energy),
                         "proper_energy_minus_ADM_mass": float(proper_energy - solution.mass),
                         "particle_mass_normalization": 1.0, "baryon_rest_mass": float(baryon_rest_mass),
                         "baryon_rest_mass_minus_ADM_mass": float(baryon_rest_mass - solution.mass),
                         "redshifted_chemical_potential_range": [float(min(chemical_redshifts)), float(max(chemical_redshifts))],
                         "solver": "SciPy solve_ivp DOP853 with zero-pressure event and dense output", "rhs_evaluations": solution.integration.nfev,
                         "stored_integration_points": len(solution.integration.t), "start_radius": solution.start_radius},
            "closure_diagnostics": closure_diagnostics(solution), "refinement_diagnostics": refinements,
            "benchmark": {"uniform_density_solution": vars(UniformFluidSphere()), "causal_material_EOS": False},
            "probe_ring": {"length": metric.length, "center_offset": metric.center_offset, "chart": "areal-radius Cartesian reference coordinates",
                           "spatial_factor_rule": "b²=1+(B(r)-1)*(dr/dx)², not the C9 isotropic spatial coefficient",
                           "clock_lapses": clocks.tolist(), "near_to_far_lapse_ratio": float(clocks[1] / clocks[0]),
                           "conservative_speed_floor": metric.conservative_speed_floor},
            "probe_trials": trials, "probe_schedule": schedule,
            "fresh_probe_trials": 2, "cumulative_logical_qubit_initializations": 2 * (sites + 2),
            "resources_per_trial": {"chain_qubits": sites, "clock_qubits": 2, "raw_record_bits": sites + 2,
                                    "independent_fair_bits": 12, "threshold_table_bits": sites * 13,
                                    "bit_hops": sum(event["hops"] for event in schedule), "read_parameter": duration,
                                    "allocated_hop_interval": schedule[0]["allocated_hop_interval"],
                                    "all_records_ready_by": schedule[-1]["received_by_collector_at"], "resets": 0},
            "limits": ["nonlinear Einstein-Hilbert dynamics and EOS are inputs, not derived from cognition",
                       "causal equilibrium EOS does not prove dynamical radial stability or a microscopic matter realization",
                       "finite surface density jumps to vacuum with zero pressure; no thin shell, but not globally smooth matter",
                       "test probes and their confinement do not backreact on this classical fluid source",
                       "only metric coupling to the fluid is included for the ideal probes; material scattering and confining apparatus are omitted",
                       "areal-coordinate ring differs physically from an isotropic-coordinate circle, so this is not a same-source higher-order C9 comparison",
                       "curvature residual and tolerance comparisons are diagnostics, not interval-certified solution or probability bounds",
                       "C5 continuum bounds are not extended to the surface-crossing metric profile",
                       "no quantum source closure, lab data, sample counts, source preparation or detector energy closure"]}


class NonlinearFluidMetricTests(unittest.TestCase):
    def test_uniform_sphere_satisfies_all_nonlinear_einstein_components(self):
        sphere = UniformFluidSphere()
        for radius in (0.2, 1., 2.7, 3.4, 6.):
            for angle in (0.7, 1.2):
                radial = sphere.radial(radius)
                jet = spherical_metric_jet(radius, angle, radial)
                expected = 8 * math.pi * sphere.gravitational_constant * fluid_covariant_tensor(jet[0], radial["density"], radial["pressure"])
                np.testing.assert_allclose(nonlinear_einstein_from_jet(*jet), expected, atol=3e-14)

    def test_pressure_supplies_covariant_static_force_balance(self):
        sphere = UniformFluidSphere()
        for radius in (0., 0.5, 1.5, 2.8):
            radial = sphere.radial(radius)
            self.assertGreater(radial["pressure"], 0)
            self.assertAlmostEqual(radial["pressure_first"] + (radial["density"] + radial["pressure"]) * radial["lapse_first"] / radial["lapse"], 0)
        self.assertEqual(sphere.radial(sphere.radius)["pressure"], 0)

    def test_surface_metric_and_extrinsic_curvature_match_without_a_shell(self):
        sphere = UniformFluidSphere()
        interior = sphere.radial(sphere.radius * (1 - 1e-10))
        exterior = sphere.radial(sphere.radius)
        for key in ("lapse", "lapse_first", "radial_factor"):
            self.assertAlmostEqual(interior[key], exterior[key])

    def test_regularity_limit_and_center_pressure_are_explicit(self):
        sphere = UniformFluidSphere()
        surface = math.sqrt(1 - sphere.compactness)
        self.assertAlmostEqual(sphere.radial(0.)["pressure"] / sphere.density, (1 - surface) / (3 * surface - 1))
        with self.assertRaises(ValueError):
            UniformFluidSphere(gravitational_constant=2.)

    def test_adaptive_solver_reproduces_the_independent_uniform_density_benchmark(self):
        exact = UniformFluidSphere()
        solution = FluidSolution(ConstantDensityEOS(exact.density), exact.radial(0.)["pressure"], exact.gravitational_constant)
        self.assertAlmostEqual(solution.radius, exact.radius, places=7)
        self.assertAlmostEqual(solution.mass, exact.mass, places=7)
        for radius in (0., 0.5, 1.5, 2.9):
            for key in ("pressure", "lapse", "radial_factor"):
                self.assertAlmostEqual(solution.radial(radius)[key], exact.radial(radius)[key], places=8)

    def test_causal_EOS_solution_closes_full_curvature_and_pressure_balance(self):
        solution = FluidSolution()
        diagnostics = closure_diagnostics(solution)
        self.assertLess(diagnostics["maximum_sampled_full_Einstein_residual"], 1e-11)
        self.assertLess(diagnostics["maximum_sampled_covariant_pressure_balance_residual"], 1e-13)
        self.assertLess(diagnostics["maximum_dense_output_five_point_ODE_residual"], 1e-7)
        for radius in np.linspace(0., 0.99 * solution.radius, 15):
            state = solution.radial(float(radius))
            self.assertGreaterEqual(state["pressure"], 0.)
            self.assertLess(state["pressure"], state["density"])
            self.assertLessEqual(state["density"], solution.central_density)

    def test_causal_surface_matching_and_exterior_are_not_pressureless_interior(self):
        solution = FluidSolution()
        left = solution.radial(solution.radius * (1 - 1e-9))
        right = solution.radial(solution.radius)
        for key in ("lapse", "lapse_first", "radial_factor"):
            self.assertAlmostEqual(left[key], right[key], places=8)
        self.assertGreater(left["density"], 0.)
        self.assertEqual(right["density"], 0.)
        self.assertAlmostEqual(left["pressure"], 0., places=10)
        self.assertGreater(solution.radial(0.)["pressure"], 0.)

    def test_areal_ring_uses_tangent_contraction_of_the_full_spatial_metric(self):
        solution = FluidSolution()
        metric = FluidRingMetric(solution)
        coordinates = np.linspace(0, metric.length, 33)
        radii, radial_tangent = metric.geometry(coordinates)
        radial_factors = np.array([solution.radial(float(radius))["radial_factor"] for radius in radii])
        np.testing.assert_allclose(metric.spatial_scale(coordinates) ** 2, 1 + (radial_factors - 1) * radial_tangent ** 2)
        self.assertAlmostEqual(float(metric.spatial_scale(0.)), 1.)
        self.assertAlmostEqual(float(metric.spatial_scale(metric.length / 2)), 1.)
        self.assertGreaterEqual(metric.speed(coordinates).min(), metric.conservative_speed_floor)

    def test_same_metric_probes_preserve_unitarity_energy_and_common_J(self):
        metric = FluidRingMetric(FluidSolution())
        gate, hamiltonian, _ = probe_unitary(8, 0.35, 2., metric)
        initial = np.kron(initial_data(8), np.ones(4) / 2)
        final = gate @ initial
        np.testing.assert_allclose(gate.conj().T @ gate, np.eye(32), atol=3e-14)
        self.assertAlmostEqual(np.vdot(initial, np.kron(hamiltonian, np.eye(4)) @ initial).real,
                               np.vdot(final, np.kron(hamiltonian, np.eye(4)) @ final).real)
        lifted = real_lift(gate)
        np.testing.assert_allclose(lifted @ encode_state(np.outer(initial, initial.conj())) @ lifted.T,
                                   encode_state(np.outer(final, final.conj())), atol=3e-14)

    def test_refinement_preserves_source_size_and_probe_record(self):
        first, second = FluidSolution(relative_tolerance=1e-9), FluidSolution(relative_tolerance=1e-11, start_fraction=5e-7)
        self.assertLess(abs(first.radius - second.radius), 1e-6)
        self.assertLess(abs(first.mass - second.mass), 1e-6)
        initial = np.kron(initial_data(8), np.ones(4) / 2)
        records = [probe_record_probabilities(probe_unitary(8, 0.35, 2., FluidRingMetric(solution))[0] @ initial, 8) for solution in (first, second)]
        self.assertLess(abs(records[0] - records[1]).sum() / 2, 1e-7)

    def test_invalid_EOS_and_numerical_controls_are_rejected(self):
        for arguments in ({"sound_speed_squared": 2.}, {"surface_density": -1.}):
            with self.assertRaises(ValueError):
                LinearFluidEOS(**arguments)
        with self.assertRaises(ValueError):
            FluidSolution(central_pressure=-1.)

    def test_EOS_has_a_consistent_isentropic_energy_and_particle_number(self):
        eos = LinearFluidEOS()
        for pressure in (0., 0.0002, 0.0005):
            number = eos.particle_density(pressure)
            density = eos.density_from_particles(number)
            self.assertAlmostEqual(density, eos.density(pressure))
            derivative = (number / eos.surface_density) ** eos.sound_speed_squared
            self.assertAlmostEqual(number * derivative - density, pressure)

    def test_independent_enthalpy_lapse_relation_holds_through_the_fluid(self):
        solution = FluidSolution()
        sound, density = solution.eos.sound_speed_squared, solution.eos.surface_density
        surface_lapse = solution.radial(solution.radius)["lapse"]
        for radius in np.linspace(0., solution.radius * 0.99, 17):
            row = solution.radial(float(radius))
            enthalpy_factor = 1 + (1 + sound) * row["pressure"] / (sound * density)
            predicted_lapse = surface_lapse * enthalpy_factor ** (-sound / (1 + sound))
            self.assertLess(abs(row["lapse"] - predicted_lapse), 2e-10)

    def test_uniform_solution_has_a_controlled_Newtonian_pressure_limit(self):
        previous = None
        for gravity in (0.1, 0.05, 0.025):
            sphere = UniformFluidSphere(gravitational_constant=gravity)
            radius = sphere.radius * 0.4
            pressure_newton = 2 * math.pi / 3 * gravity * sphere.density ** 2 * (sphere.radius ** 2 - radius ** 2)
            error = abs(sphere.radial(radius)["pressure"] - pressure_newton)
            if previous is not None:
                self.assertGreater(previous / error, 3.8)
                self.assertLess(previous / error, 4.5)
            previous = error

    def test_induced_spatial_scale_equals_full_areal_metric_contraction(self):
        solution = FluidSolution()
        metric = FluidRingMetric(solution)
        circle_radius = metric.length / (2 * math.pi)
        for angle in (0.3, 1.2, 2.5):
            point = np.array([metric.center_offset + circle_radius * math.cos(angle), circle_radius * math.sin(angle), 0.])
            tangent = np.array([-math.sin(angle), math.cos(angle), 0.])
            radial_direction = point / np.linalg.norm(point)
            factor = solution.radial(float(np.linalg.norm(point)))["radial_factor"]
            spatial_metric = np.eye(3) + (factor - 1) * np.outer(radial_direction, radial_direction)
            self.assertAlmostEqual(float(metric.spatial_scale(angle * circle_radius)) ** 2, tangent @ spatial_metric @ tangent)

    def test_complete_probe_mixtures_and_both_masses_use_the_same_geometry(self):
        metric = FluidRingMetric(FluidSolution())
        gate, first_hamiltonian, _ = probe_unitary(8, 0.35, 1., metric)
        second_hamiltonian = massive_hopping(8, 0.55, metric)
        np.testing.assert_allclose(np.diag(first_hamiltonian) / 0.35, np.diag(second_hamiltonian) / 0.55)
        np.testing.assert_allclose(first_hamiltonian - np.diag(np.diag(first_hamiltonian)),
                                   second_hamiltonian - np.diag(np.diag(second_hamiltonian)))
        first, second = np.eye(8)[:, 0], np.eye(8)[:, 2]
        def output(data):
            state = gate @ np.kron(data, np.ones(4) / 2)
            return np.outer(state, state.conj())
        direct = (output(first) + output(second)) / 2
        alternate = sum(output((first + sign * second) / math.sqrt(2)) for sign in (1, -1)) / 2
        np.testing.assert_allclose(direct, alternate, atol=3e-14)

    def test_all_raw_probe_records_wait_for_conservative_metric_transport(self):
        metric = FluidRingMetric(FluidSolution())
        positions = [0] * 8
        positions[3] = 1
        schedule = probe_schedule(8, 2., metric)
        self.assertEqual(sum(event["hops"] for event in schedule), 20)
        deadline = schedule[-1]["received_by_collector_at"]
        with self.assertRaises(ValueError):
            release_probe_record(positions, [0, 1], 0, deadline - 0.01, 2., metric)
        thresholds, _ = postprocessing_table(8, center_fraction=0.5)
        for random in (0, int(thresholds[3]) - 1, int(thresholds[3]), 4095):
            self.assertEqual(release_probe_record(positions, [0, 1], random, deadline, 2., metric),
                             (int(random >= thresholds[3]), 0, 1))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    outcome = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(NonlinearFluidMetricTests))
    if not outcome.wasSuccessful():
        raise SystemExit(1)
    if arguments.write_results:
        report = run_report()
        Path(__file__).with_name("nonlinear_fluid_metric_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"tests": outcome.testsRun, "solution": report["solution"], "closure": report["closure_diagnostics"],
                          "ring": report["probe_ring"], "probes": report["probe_trials"], "resources": report["resources_per_trial"]}, indent=2))