"""Construction 11: radial perturbations of the unchanged C10 fluid sphere.

The inherited Einstein-fluid contract and adiabatic EOS are retained. Radial
linear stability is not a proof of nonlinear, nonradial or quantum stability.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np
from numpy.polynomial.chebyshev import Chebyshev
from scipy.linalg import eigh
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from nonlinear_fluid_metric import FluidSolution, LinearFluidEOS, FluidRingMetric
from massive_envelope import massive_hopping
from weak_metric_response import probe_unitary, probe_record_probabilities
from quantum_chain_background import initial_data
from clock_calibration import clock_gate
from baseline import real_lift
from common_orientation_structure import encode_state


class RadialProblem:
    def __init__(self, solution=None):
        self.solution = FluidSolution() if solution is None else solution
        if not isinstance(self.solution.eos, LinearFluidEOS):
            raise ValueError("This perturbation contract uses the C10 causal linear EOS.")

    def equilibrium(self, radius):
        if not np.isfinite(radius) or not 0 <= radius <= self.solution.radius:
            raise ValueError("Use an interior areal radius including the fluid-side surface.")
        row = self.solution.radial(radius).copy()
        if radius == self.solution.radius:
            row["density"] = self.solution.eos.surface_density
            row["pressure"] = 0.
            row["pressure_first"] = -row["density"] * row["lapse_first"] / row["lapse"]
        return row

    def coefficients(self, radius):
        if radius <= 0:
            raise ValueError("Singular Sturm-Liouville coefficients are evaluated only at positive radius.")
        row = self.equilibrium(radius)
        density, pressure = row["density"], row["pressure"]
        lapse, radial_factor = row["lapse"], row["radial_factor"]
        sound = self.solution.eos.sound_speed_squared
        common = math.sqrt(radial_factor) * lapse ** 3 / radius ** 2
        stiffness = sound * (density + pressure) * common
        potential = common * (row["pressure_first"] ** 2 / (density + pressure)
                    - 4 * row["pressure_first"] / radius
                    - 8 * math.pi * self.solution.gravitational_constant * radial_factor * pressure * (density + pressure))
        inertia = (density + pressure) * radial_factor ** 1.5 * lapse / radius ** 2
        return stiffness, potential, inertia

    def ritz(self, basis_size=10, quadrature_order=192):
        if not isinstance(basis_size, int) or basis_size < 3 or not isinstance(quadrature_order, int) or quadrature_order < 4 * basis_size:
            raise ValueError("Use at least three regular basis functions and sufficient quadrature.")
        nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
        fractions = (nodes + 1) / 2
        radii = self.solution.radius * fractions
        weights = weights * self.solution.radius / 2
        values, derivatives = regular_basis(fractions, basis_size, self.solution.radius)
        stiffness, potential, inertia = np.array([self.coefficients(float(radius)) for radius in radii]).T
        kinetic_matrix = derivatives.T @ ((weights * stiffness)[:, None] * derivatives)
        restoring_matrix = kinetic_matrix - values.T @ ((weights * potential)[:, None] * values)
        inertia_matrix = values.T @ ((weights * inertia)[:, None] * values)
        eigenvalues, eigenvectors = eigh(restoring_matrix, inertia_matrix)
        return {"eigenvalues": eigenvalues, "eigenvectors": eigenvectors,
                "restoring_matrix": restoring_matrix, "inertia_matrix": inertia_matrix,
                "basis_size": basis_size, "quadrature_order": quadrature_order}

    def sufficient_stability_bound(self):
        solution = self.solution
        gravity, central_pressure = solution.gravitational_constant, solution.central_pressure
        surface_density, sound = solution.eos.surface_density, solution.eos.sound_speed_squared
        central_density = solution.central_density
        radius_squared_upper = 3 * central_pressure / (2 * math.pi * gravity * surface_density ** 2)
        compactness_upper = 8 * math.pi * gravity * central_density * radius_squared_upper / 3
        if compactness_upper >= 1:
            raise ValueError("The simple parameter-only stability bound is inconclusive for this equilibrium.")
        enthalpy = 1 + (1 + sound) * central_pressure / (sound * surface_density)
        lapse_lower = math.sqrt(1 - compactness_upper) * enthalpy ** (-sound / (1 + sound))
        stiffness_lower = sound * surface_density * lapse_lower ** 3
        inertia_upper = (central_density + central_pressure) / (1 - compactness_upper) ** 1.5
        potential_gradient_bound = 4 * math.pi * gravity * (central_density / 3 + central_pressure) / (1 - compactness_upper)
        potential_ratio_upper = potential_gradient_bound ** 2 * radius_squared_upper + 4 * potential_gradient_bound
        kinetic_lower = math.pi ** 2 * stiffness_lower / (radius_squared_upper * inertia_upper)
        return {"radius_squared_upper": radius_squared_upper, "compactness_upper": compactness_upper,
                "lapse_lower": lapse_lower, "weighted_kinetic_lower": kinetic_lower,
                "potential_over_inertia_upper": potential_ratio_upper,
                "omega_squared_lower": kinetic_lower - potential_ratio_upper,
                "scope": "analytic conditional bound for any regular equilibrium with the declared parameters; floating evaluation is not outward rounding"}

    def shooting(self, omega_squared, relative_tolerance=1e-10, start_fraction=1e-5):
        if (not np.isfinite(omega_squared) or not 0 < relative_tolerance <= 1e-3
                or not 0 < start_fraction < 1e-2):
            raise ValueError("Use finite squared frequency and controlled integration settings.")
        radius = self.solution.radius
        center = self.equilibrium(0.)
        initial = [1., -3 * self.solution.eos.sound_speed_squared * (center["density"] + center["pressure"])]
        result = solve_ivp(lambda position, state: self.shooting_rhs(position, state, omega_squared),
                           (radius * start_fraction, radius), initial, method="DOP853", dense_output=True,
                           rtol=relative_tolerance, atol=relative_tolerance * np.array([1e-2, 1e-4]), max_step=radius / 60)
        if not result.success:
            raise ValueError("The radial shooting integration failed.")
        return result

    def shooting_rhs(self, radius, state, omega_squared):
        fractional_displacement, delta_pressure = state
        row = self.equilibrium(radius)
        density, pressure, pressure_first = row["density"], row["pressure"], row["pressure_first"]
        total = density + pressure
        sound = self.solution.eos.sound_speed_squared
        factor, lapse = row["radial_factor"], row["lapse"]
        gravity = self.solution.gravitational_constant
        first = -(3 * fractional_displacement + delta_pressure / (sound * total)) / radius - pressure_first / total * fractional_displacement
        second = fractional_displacement * (omega_squared * total * factor * radius / lapse ** 2 - 4 * pressure_first
                 + radius * pressure_first ** 2 / total - 8 * math.pi * gravity * factor * total * pressure * radius)
        second += delta_pressure * (pressure_first / total - 4 * math.pi * gravity * total * factor * radius)
        return np.array([first, second])

    def find_mode(self, index=0, relative_tolerance=1e-10, start_fraction=1e-5):
        if not isinstance(index, int) or not 0 <= index <= 2:
            raise ValueError("This report resolves only the first three radial modes.")
        estimates = self.ritz(12)["eigenvalues"]
        lower = 0. if index == 0 else float((estimates[index - 1] + estimates[index]) / 2)
        upper = float((estimates[index] + estimates[index + 1]) / 2)
        root = brentq(lambda value: self.shooting(value, relative_tolerance, start_fraction).y[1, -1],
                      lower, upper, xtol=1e-12, rtol=1e-12)
        result = self.shooting(root, relative_tolerance, start_fraction)
        grid = np.linspace(result.t[0], self.solution.radius, 1201)
        displacements = result.sol(grid)[0]
        nodes = int(np.count_nonzero(displacements[1:] * displacements[:-1] < 0))
        return {"omega_squared": root, "omega": math.sqrt(root), "shooting_solution": result,
                "node_count": nodes, "surface_pressure_residual": float(result.y[1, -1]),
                "bracket": [lower, upper], "relative_tolerance": relative_tolerance, "start_fraction": start_fraction}


class RadialMode:
    def __init__(self, problem, mode=None):
        self.problem = problem
        self.mode = problem.find_mode() if mode is None else mode
        self.solution = problem.solution
        self.radius = self.solution.radius
        self.omega = self.mode["omega"]
        self.normalization = self.mode["shooting_solution"].y[0, -1]
        if abs(self.normalization) < 1e-10:
            raise ValueError("The surface fractional displacement cannot be used to normalize this mode.")
        self.minimum_radius = self.mode["shooting_solution"].t[0]
        self.lapse_integration = solve_ivp(lambda radius, value: [self.local(radius)["delta_nu_first"]],
                                          (self.radius, self.minimum_radius), [0.], method="DOP853",
                                          rtol=1e-10, atol=1e-12, dense_output=True, max_step=self.radius / 80)
        if not self.lapse_integration.success:
            raise ValueError("The linear lapse constraint integration failed.")

    def local(self, radius):
        if not 0 < radius <= self.radius:
            raise ValueError("Use a positive fluid-side radius for the linear perturbation.")
        integration = self.mode["shooting_solution"]
        sample_radius = max(radius, integration.t[0])
        fractional, lagrangian_pressure = integration.sol(sample_radius) / self.normalization
        row = self.problem.equilibrium(radius)
        gravity, sound = self.solution.gravitational_constant, self.solution.eos.sound_speed_squared
        total = row["density"] + row["pressure"]
        displacement = radius * fractional
        eulerian_pressure = lagrangian_pressure - displacement * row["pressure_first"]
        eulerian_density = eulerian_pressure / sound
        delta_mass = -4 * math.pi * radius ** 2 * total * displacement
        delta_lambda = gravity * row["radial_factor"] * delta_mass / radius
        denominator = radius * (radius - 2 * gravity * row["mass"])
        delta_nu_first = gravity * (delta_mass + 4 * math.pi * radius ** 3 * eulerian_pressure) / denominator
        delta_nu_first += 2 * gravity * delta_mass / (radius - 2 * gravity * row["mass"]) * row["lapse_first"] / row["lapse"]
        zeta = radius ** 2 * displacement / row["lapse"]
        zeta_first = -lagrangian_pressure * radius ** 2 / (sound * total * row["lapse"])
        return {"fractional_displacement": float(fractional), "displacement": float(displacement),
                "lagrangian_pressure": float(lagrangian_pressure), "eulerian_pressure": float(eulerian_pressure),
                "eulerian_density": float(eulerian_density), "delta_mass": float(delta_mass),
                "delta_lambda": float(delta_lambda), "delta_nu_first": float(delta_nu_first),
                "zeta": float(zeta), "zeta_first": float(zeta_first)}

    def delta_nu(self, radius):
        if not np.isfinite(radius) or radius < 0:
            raise ValueError("Use a finite nonnegative fixed areal radius.")
        if radius >= self.radius:
            return 0.
        return float(self.lapse_integration.sol(max(radius, self.minimum_radius))[0])

    def conservation_report(self, quadrature_order=256):
        nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
        radii = self.radius * (nodes + 1) / 2
        weights *= self.radius / 2
        number_integrand, mass_integrand, kinetic, potential = [], [], [], []
        for radius in radii:
            row = self.problem.equilibrium(float(radius))
            perturbation = self.local(float(radius))
            number = self.solution.eos.particle_density(row["pressure"])
            delta_number = number * perturbation["eulerian_density"] / (row["density"] + row["pressure"])
            number_integrand.append(4 * math.pi * radius ** 2 * math.sqrt(row["radial_factor"]) * (delta_number + number * perturbation["delta_lambda"]))
            mass_integrand.append(4 * math.pi * radius ** 2 * perturbation["eulerian_density"])
            stiffness, source, inertia = self.problem.coefficients(float(radius))
            kinetic.append(inertia * perturbation["zeta"] ** 2)
            potential.append(stiffness * perturbation["zeta_first"] ** 2 - source * perturbation["zeta"] ** 2)
        surface = self.problem.equilibrium(self.radius)
        displacement = self.local(self.radius)["displacement"]
        number_boundary = 4 * math.pi * self.radius ** 2 * math.sqrt(surface["radial_factor"]) * self.solution.eos.particle_density(0.) * displacement
        mass_boundary = 4 * math.pi * self.radius ** 2 * surface["density"] * displacement
        inertia_integral = float(np.dot(weights, kinetic))
        stiffness_integral = float(np.dot(weights, potential))
        return {"particle_number_bulk_variation": float(np.dot(weights, number_integrand)),
                "particle_number_moving_surface_term": float(number_boundary),
                "particle_number_total_variation": float(np.dot(weights, number_integrand) + number_boundary),
                "mass_bulk_variation": float(np.dot(weights, mass_integrand)), "mass_moving_surface_term": float(mass_boundary),
                "ADM_mass_total_variation": float(np.dot(weights, mass_integrand) + mass_boundary),
                "surface_mass_constraint_total": self.local(self.radius)["delta_mass"] + mass_boundary,
                "modal_inertia": inertia_integral, "modal_restoring_energy_coefficient": stiffness_integral,
                "rayleigh_omega_squared": stiffness_integral / inertia_integral,
                "quadrature_order": quadrature_order, "normalization": "surface displacement per unit mode amplitude equals R"}

def regular_basis(fractions, size, radius):
    fractions = np.asarray(fractions)
    columns, derivatives = [], []
    argument = 2 * fractions ** 2 - 1
    for degree in range(size):
        polynomial = Chebyshev.basis(degree)
        columns.append(fractions ** 3 * polynomial(argument))
        derivatives.append((3 * fractions ** 2 * polynomial(argument)
                            + 4 * fractions ** 4 * polynomial.deriv()(argument)) / radius)
    return np.array(columns).T, np.array(derivatives).T


def rational_default_certificate():
    pi_lower, pi_upper = Fraction(157, 50), Fraction(22, 7)
    radius_squared_upper, compactness_upper = Fraction(24), Fraction(23, 100)
    lapse_lower, inertia_upper = Fraction(419, 500), Fraction(9, 500)
    gradient_upper, pi_squared_lower = Fraction(71, 10000), Fraction(197, 20)
    surface_density, central_density, central_pressure = Fraction(1, 100), Fraction(23, 2000), Fraction(1, 2000)
    gravity, sound = Fraction(1, 10), Fraction(1, 3)
    assert Fraction(75) / pi_lower < radius_squared_upper
    assert lapse_lower ** 4 < (1 - compactness_upper) ** 2 / Fraction(6, 5)
    assert inertia_upper ** 2 * (1 - compactness_upper) ** 3 > (central_density + central_pressure) ** 2
    assert 4 * pi_upper * gravity * (central_density / 3 + central_pressure) / (1 - compactness_upper) < gradient_upper
    assert pi_lower ** 2 > pi_squared_lower
    kinetic_lower = pi_squared_lower * sound * surface_density * lapse_lower ** 3 / (radius_squared_upper * inertia_upper)
    potential_upper = gradient_upper ** 2 * radius_squared_upper + 4 * gradient_upper
    difference = kinetic_lower - potential_upper
    assert difference > Fraction(3, 200)
    return {"exact_lower_fraction": str(difference), "decimal_display_only": float(difference),
            "certified_omega_squared_strictly_greater_than": "3/200",
            "method": "exact Fraction arithmetic, using 157/50 < pi < 22/7 and the derived parameter-only energy inequalities",
            "scope": "any regular C10 equilibrium with exactly G=1/10, rho_s=1/100, p_c=1/2000, sound_squared=1/3; radial linear problem only"}


def dynamic_probe_terms(mode, sites, mass, amplitude):
    if not np.isfinite(amplitude) or abs(amplitude) > 0.01:
        raise ValueError("Use the declared small linear-mode amplitude range |delta R/R|<=0.01.")
    metric = FluidRingMetric(mode.solution)
    positions = np.arange(sites) * metric.length / sites
    midpoints = positions + metric.length / (2 * sites)
    hamiltonian = massive_hopping(sites, mass, metric)

    def fractions(coordinates):
        radii, tangent = metric.geometry(coordinates)
        if np.min(abs(radii - mode.radius)) <= abs(amplitude) * mode.radius:
            raise ValueError("A moving surface can cross a sampled probe coefficient; this fixed-side linear task is not applicable.")
        lapse_changes, speed_changes = [], []
        for radius, radial_tangent in zip(radii, tangent):
            lapse_change = mode.delta_nu(float(radius))
            radial_change = mode.local(float(radius))["delta_lambda"] if radius < mode.radius else 0.
            radial_factor = mode.solution.radial(float(radius))["radial_factor"]
            spatial_squared = 1 + (radial_factor - 1) * radial_tangent ** 2
            speed_change = lapse_change - radial_factor * radial_change * radial_tangent ** 2 / spatial_squared
            lapse_changes.append(lapse_change)
            speed_changes.append(speed_change)
        return np.array(lapse_changes), np.array(speed_changes)

    lapse_changes, _ = fractions(positions)
    _, speed_changes = fractions(midpoints)
    if max(np.max(abs(amplitude * lapse_changes)), np.max(abs(amplitude * speed_changes))) >= 0.05:
        raise ValueError("The declared coefficient modulation must stay below five percent.")
    derivative = np.diag((-1.) ** np.arange(sites) * mass * metric.lapse(positions) * lapse_changes)
    for site in range(sites):
        adjacent = (site + 1) % sites
        derivative[site, adjacent] = derivative[adjacent, site] = hamiltonian[site, adjacent] * speed_changes[site]
    clock_rates = np.array([metric.lapse(0.), 1.2 * metric.lapse(metric.length / 2)])
    clock_changes = lapse_changes[[0, sites // 2]]
    speed_floor = metric.conservative_speed_floor * (1 - float(np.max(abs(amplitude * speed_changes))))
    return {"static_hamiltonian": hamiltonian, "mode_hamiltonian_derivative": derivative,
            "clock_rates": clock_rates, "clock_lapse_fractional_derivatives": clock_changes,
            "sampled_speed_fractional_derivatives": speed_changes, "metric": metric,
            "sampled_lattice_speed_floor": speed_floor}


def evolve_probe(mode, sites=32, mass=0.35, duration=2., amplitude=1e-3, relative_tolerance=1e-11):
    if not np.isfinite(duration) or duration <= 0 or not 0 < relative_tolerance <= 1e-3:
        raise ValueError("Use a positive finite interval and controlled integrator tolerance.")
    terms = dynamic_probe_terms(mode, sites, mass, amplitude)
    static = terms["static_hamiltonian"]
    derivative = terms["mode_hamiltonian_derivative"]

    def evolution(time, vector):
        gate = vector.reshape(sites, sites)
        return (-1j * (static + amplitude * math.cos(mode.omega * time) * derivative) @ gate).reshape(-1)

    integration = solve_ivp(evolution, (0., duration), np.eye(sites, dtype=complex).reshape(-1), method="DOP853",
                            rtol=relative_tolerance, atol=relative_tolerance * 1e-2, dense_output=True)
    if not integration.success:
        raise ValueError("The time-dependent probe integration failed.")
    elapsed = duration + amplitude * terms["clock_lapse_fractional_derivatives"] * math.sin(mode.omega * duration) / mode.omega
    clocks = np.kron(clock_gate(terms["clock_rates"][0], elapsed[0]), clock_gate(terms["clock_rates"][1], elapsed[1]))
    data_gate = integration.y[:, -1].reshape(sites, sites)
    full_gate = np.kron(data_gate, clocks)
    prepared_data = initial_data(sites)
    prepared = np.kron(prepared_data, np.ones(4) / 2)
    final = full_gate @ prepared
    nodes, weights = np.polynomial.legendre.leggauss(64)
    times = duration * (nodes + 1) / 2
    power = []
    for time in times:
        state = integration.sol(time).reshape(sites, sites) @ prepared_data
        power.append(-amplitude * mode.omega * math.sin(mode.omega * time) * np.vdot(state, derivative @ state).real)
    initial_energy = float(np.vdot(prepared_data, (static + amplitude * derivative) @ prepared_data).real)
    final_data = data_gate @ prepared_data
    final_energy = float(np.vdot(final_data, (static + amplitude * math.cos(mode.omega * duration) * derivative) @ final_data).real)
    return {"gate": full_gate, "integration": integration, "terms": terms,
            "records": probe_record_probabilities(final, sites),
            "unitarity_residual": float(np.linalg.norm(data_gate.conj().T @ data_gate - np.eye(sites))),
            "data_energy_initial": initial_energy, "data_energy_final": final_energy,
            "integrated_data_work": float(duration / 2 * np.dot(weights, power)),
            "clock_phases": (terms["clock_rates"] * elapsed).tolist()}


def record_schedule(sites, duration, hop_interval=1.):
    if not isinstance(sites, int) or sites < 8 or sites % 4 or duration < 0 or not np.isfinite(duration) or not np.isfinite(hop_interval) or hop_interval <= 0:
        raise ValueError("Use a valid finite protocol and positive independently supplied transport slot.")
    sources = [("position", site) for site in range(sites)] + [("clock", 0), ("clock", sites // 2)]
    receipt, schedule = duration, []
    for kind, site in sources:
        hops = min(site, sites - site)
        receipt += hops * hop_interval
        schedule.append({"kind": kind, "source_site": site, "hops": hops, "read_at": duration,
                         "received_at": receipt})
    return schedule


def release_record(positions, clocks, random_integer, receipt, duration, hop_interval=1.):
    from massive_envelope import postprocessing_table
    schedule = record_schedule(len(positions), duration, hop_interval)
    if (sum(positions) != 1 or any(bit not in (0, 1) for bit in positions) or len(clocks) != 2
            or any(bit not in (0, 1) for bit in clocks) or not isinstance(random_integer, int) or not 0 <= random_integer < 4096):
        raise ValueError("Keep all position and clock outcomes and twelve independent fair random bits.")
    if not np.isfinite(receipt) or receipt < schedule[-1]["received_at"]:
        raise ValueError("Wait for every raw record before publishing the coarse output.")
    thresholds, _ = postprocessing_table(len(positions), center_fraction=0.5)
    return int(random_integer >= thresholds[list(positions).index(1)]), *clocks


def linear_constraint_diagnostics(mode):
    euler_residuals, mass_residuals = [], []
    for fraction in (0.2, 0.5, 0.8):
        radius, step = mode.radius * fraction, mode.radius * 1e-4
        row = mode.problem.equilibrium(radius)
        perturbation = mode.local(radius)
        def derivative(key):
            return (mode.local(radius - 2 * step)[key] - 8 * mode.local(radius - step)[key]
                    + 8 * mode.local(radius + step)[key] - mode.local(radius + 2 * step)[key]) / (12 * step)
        total = row["density"] + row["pressure"]
        euler = (-mode.omega ** 2 * total * row["radial_factor"] / row["lapse"] ** 2 * perturbation["displacement"]
                 + derivative("eulerian_pressure")
                 + (perturbation["eulerian_density"] + perturbation["eulerian_pressure"]) * row["lapse_first"] / row["lapse"]
                 + total * perturbation["delta_nu_first"])
        mass = derivative("delta_mass") - 4 * math.pi * radius ** 2 * perturbation["eulerian_density"]
        euler_residuals.append(abs(euler))
        mass_residuals.append(abs(mass))
    return {"maximum_sampled_linear_Euler_residual": float(max(euler_residuals)),
            "maximum_sampled_linear_mass_constraint_residual": float(max(mass_residuals)),
            "method": "five-point differentiation of reconstructed shooting fields; numerical diagnostic, not a global error certificate"}


def run_report():
    problem = RadialProblem()
    modes = [problem.find_mode(index) for index in range(3)]
    fundamental = RadialMode(problem, modes[0])
    conservation = fundamental.conservation_report()
    refined_problem = RadialProblem(FluidSolution(relative_tolerance=1e-12, start_fraction=5e-7))
    refined_mode = refined_problem.find_mode(relative_tolerance=1e-12, start_fraction=5e-6)
    refined_fundamental = RadialMode(refined_problem, refined_mode)
    probe = evolve_probe(fundamental)
    refined_probe = evolve_probe(refined_fundamental, relative_tolerance=1e-12)
    static_gate, _, _ = probe_unitary(32, 0.35, 2., FluidRingMetric(problem.solution))
    initial = np.kron(initial_data(32), np.ones(4) / 2)
    static_records = probe_record_probabilities(static_gate @ initial, 32)
    profile = []
    for fraction in (0.05, 0.2, 0.5, 0.8, 1.):
        radius = fraction * fundamental.radius
        profile.append(dict(radius=radius, radius_fraction=fraction, delta_nu=fundamental.delta_nu(radius), **fundamental.local(radius)))
    amplitude = 1e-3
    modal_energy = amplitude ** 2 * fundamental.omega ** 2 * conservation["modal_inertia"] / 2
    return {"model_version": "C11", "route_round": "11", "research_scope": "sufficiency only",
            "equilibrium_inputs_unchanged_from_C10": {"surface_density": 0.01, "sound_speed_squared": 1 / 3, "central_pressure": 0.0005, "G": 0.1},
            "radial_equation": "(P zeta')'+(Q+omega^2 W)zeta=0, zeta=r^2 exp(-nu) xi, center zeta~r^3, fluid-side surface zeta'=0",
            "analytic_stability_bound": problem.sufficient_stability_bound(),
            "exact_parameter_stability_certificate": rational_default_certificate(),
            "modes": [{key: value for key, value in mode.items() if key != "shooting_solution"} for mode in modes],
            "ritz_refinement": [{"basis_size": size, "first_three_omega_squared": problem.ritz(size)["eigenvalues"][:3].tolist()} for size in (4, 6, 8, 12)],
            "fundamental_refinement": {"reference_omega_squared": refined_mode["omega_squared"],
                                       "absolute_omega_squared_difference": abs(modes[0]["omega_squared"] - refined_mode["omega_squared"]),
                                       "probe_record_TV_difference": float(abs(probe["records"] - refined_probe["records"]).sum() / 2)},
            "fundamental_profile": profile, "linear_conservation": conservation,
            "reconstructed_constraint_diagnostics": linear_constraint_diagnostics(fundamental),
            "mode_amplitude_surface_fraction": amplitude, "quadratic_mode_energy_in_declared_SL_normalization": modal_energy,
            "probe_task": {"duration": 2., "proper_mass": 0.35, "joint_record_probabilities": probe["records"].tolist(),
                           "static_reference_joint_record_probabilities": static_records.tolist(),
                           "record_TV_change_from_static": float(abs(probe["records"] - static_records).sum() / 2),
                           "clock_A_X_plus_probability": float(probe["records"][:, 0, :].sum()),
                           "clock_C_Y_plus_probability": float(probe["records"][:, :, 0].sum()),
                           "clock_lapse_fractional_derivatives": probe["terms"]["clock_lapse_fractional_derivatives"].tolist(),
                           "clock_phases": probe["clock_phases"],
                           "data_energy_initial": probe["data_energy_initial"], "data_energy_final": probe["data_energy_final"],
                           "integrated_data_work": probe["integrated_data_work"],
                           "unitarity_residual": probe["unitarity_residual"]},
            "resources": {"chain_qubits": 32, "clock_qubits": 2, "raw_bits": 34,
                          "independent_fair_output_bits": 12, "threshold_table_bits": 416,
                          "bit_hops": 272, "independently_allocated_hop_interval": 1., "all_records_ready_at": 274.,
                          "transport_timing_is_declared_not_derived_from_the_dynamic_metric": True, "resets": 0},
            "record_schedule": record_schedule(32, 2.),
            "limits": ["same adopted Einstein-fluid/EOS contract; perturbation is radial, adiabatic and linear only",
                       "parameter-only positive bound is a sufficient radial energy result conditional on a regular C10 equilibrium; no nonradial or nonlinear stability claim",
                       "Ritz positive eigenvalues alone are upper estimates, not the stability proof; shooting and node counts are numerical diagnostics",
                       "all moving-surface particle and ADM mass terms are retained; exterior first-order mass and lapse variation vanish",
                       "probe H(t)=H0+amplitude*cos(omega*t)*H1 defines a finite linearized metric-coefficient model; no bound on omitted amplitude-squared gravity effects",
                       "sampled sites and bonds must not be crossed by the moving fluid boundary during the protocol",
                       "probe energy changes through explicit work, but mode depletion and quantum-probe backreaction on the fluid are not included",
                       "one-unit classical transport slots are a supplied service, not a solved null-signal propagation model in the time-dependent spacetime",
                       "no interval-certified numerics, experimental samples, full quantum source, or derivation of the gravity action from cognition"]}
class RadialFluidModeTests(unittest.TestCase):
    def test_surface_uses_nonzero_interior_compressibility(self):
        current = RadialProblem()
        surface = current.equilibrium(current.solution.radius)
        self.assertEqual(surface["density"], 0.01)
        self.assertEqual(surface["pressure"], 0.)
        self.assertGreater(current.coefficients(current.solution.radius)[0], 0)
        self.assertEqual(current.solution.radial(current.solution.radius)["density"], 0.)

    def test_ritz_problem_is_symmetric_positive_inertia_and_converges(self):
        current = RadialProblem()
        previous = None
        for size in (4, 6, 10):
            result = current.ritz(size)
            np.testing.assert_allclose(result["restoring_matrix"], result["restoring_matrix"].T, atol=1e-14)
            self.assertGreater(np.linalg.eigvalsh(result["inertia_matrix"]).min(), 0)
            self.assertGreater(result["eigenvalues"][0], 0)
            if previous is not None:
                self.assertLessEqual(result["eigenvalues"][0], previous + 1e-10)
            previous = result["eigenvalues"][0]

    def test_parameter_only_energy_bound_is_positive_and_below_ritz(self):
        current = RadialProblem()
        bound = current.sufficient_stability_bound()
        self.assertGreater(bound["omega_squared_lower"], 0)
        self.assertLess(bound["omega_squared_lower"], current.ritz()["eigenvalues"][0])
        self.assertLess(current.solution.radius ** 2, bound["radius_squared_upper"])
        self.assertLess(current.solution.compactness, bound["compactness_upper"])

    def test_regular_basis_has_cubic_center_behavior(self):
        values, derivatives = regular_basis(np.array([0., 1e-6, 2e-6]), 5, 4.)
        np.testing.assert_array_equal(values[0], np.zeros(5))
        np.testing.assert_array_equal(derivatives[0], np.zeros(5))
        np.testing.assert_allclose(values[2] / values[1], 8, atol=1e-9)

    def test_first_three_shooting_modes_match_ritz_and_have_expected_nodes(self):
        problem = RadialProblem()
        estimates = problem.ritz(12)["eigenvalues"]
        for index in range(3):
            result = problem.find_mode(index)
            self.assertEqual(result["node_count"], index)
            self.assertLess(abs(result["surface_pressure_residual"]), 1e-10)
            self.assertLess(abs(result["omega_squared"] - estimates[index]), 1e-6)

    def test_shooting_and_self_adjoint_flux_equations_agree_locally(self):
        problem = RadialProblem()
        radius = problem.solution.radius * 0.6
        row = problem.equilibrium(radius)
        fractional, lagrangian_pressure = 0.4, -0.001
        omega_squared = 0.1
        first, second = problem.shooting_rhs(radius, [fractional, lagrangian_pressure], omega_squared)
        lapse_log_first = row["lapse_first"] / row["lapse"]
        radial_log_first = row["radial_factor_first"] / (2 * row["radial_factor"])
        zeta = radius ** 3 / row["lapse"] * fractional
        zeta_first = radius ** 2 / row["lapse"] * ((3 - radius * lapse_log_first) * fractional + radius * first)
        stiffness, potential, inertia = problem.coefficients(radius)
        self.assertAlmostEqual(stiffness * zeta_first, -math.sqrt(row["radial_factor"]) * row["lapse"] ** 2 * lagrangian_pressure)
        flux_first = -math.sqrt(row["radial_factor"]) * row["lapse"] ** 2 * (second + (radial_log_first + 2 * lapse_log_first) * lagrangian_pressure)
        self.assertAlmostEqual(flux_first, -(potential + omega_squared * inertia) * zeta)

    def test_moving_boundary_restores_particle_and_ADM_mass_conservation(self):
        mode = RadialMode(RadialProblem())
        report = mode.conservation_report()
        self.assertGreater(abs(report["particle_number_bulk_variation"]), 1.)
        self.assertLess(abs(report["particle_number_total_variation"]), 1e-6)
        self.assertLess(abs(report["ADM_mass_total_variation"]), 1e-6)
        self.assertAlmostEqual(report["surface_mass_constraint_total"], 0.)
        self.assertLess(abs(report["rayleigh_omega_squared"] - mode.omega ** 2), 1e-7)

    def test_lapse_perturbation_is_matched_to_zero_exterior_ADM_variation(self):
        mode = RadialMode(RadialProblem())
        self.assertEqual(mode.delta_nu(mode.radius), 0.)
        self.assertEqual(mode.delta_nu(2 * mode.radius), 0.)
        self.assertGreater(abs(mode.delta_nu(mode.radius / 3)), 1e-4)
        radius, step = mode.radius / 2, mode.radius * 1e-4
        derivative = (mode.delta_nu(radius + step) - mode.delta_nu(radius - step)) / (2 * step)
        self.assertLess(abs(derivative - mode.local(radius)["delta_nu_first"]), 1e-7)

    def test_weighted_kinetic_inequality_has_the_declared_sharp_ground_function(self):
        nodes, weights = np.polynomial.legendre.leggauss(128)
        radii, weights = (nodes + 1) / 2, weights / 2
        ground = np.sin(math.pi * radii) - math.pi * radii * np.cos(math.pi * radii)
        derivative = math.pi ** 2 * radii * np.sin(math.pi * radii)
        kinetic = np.dot(weights, derivative ** 2 / radii ** 2)
        norm = np.dot(weights, ground ** 2 / radii ** 2)
        self.assertAlmostEqual(kinetic / norm, math.pi ** 2)

    def test_time_dependent_probe_is_unitary_and_energy_change_matches_work(self):
        mode = RadialMode(RadialProblem())
        result = evolve_probe(mode, sites=8)
        self.assertLess(result["unitarity_residual"], 1e-9)
        self.assertAlmostEqual(result["data_energy_final"] - result["data_energy_initial"], result["integrated_data_work"], places=10)
        self.assertGreater(abs(result["integrated_data_work"]), 1e-10)

    def test_zero_mode_amplitude_recovers_C10_static_joint_evolution(self):
        mode = RadialMode(RadialProblem())
        result = evolve_probe(mode, sites=8, amplitude=0.)
        reference, _, _ = probe_unitary(8, 0.35, 2., FluidRingMetric(mode.solution))
        np.testing.assert_allclose(result["gate"], reference, atol=2e-10)

    def test_clock_phase_integrates_the_mode_instead_of_freezing_the_final_lapse(self):
        mode = RadialMode(RadialProblem())
        result = evolve_probe(mode, sites=8)
        terms = result["terms"]
        self.assertEqual(terms["clock_lapse_fractional_derivatives"][0], 0.)
        phase = result["clock_phases"][1]
        expected = terms["clock_rates"][1] * (2 + 1e-3 * terms["clock_lapse_fractional_derivatives"][1] * math.sin(2 * mode.omega) / mode.omega)
        self.assertAlmostEqual(phase, expected)
        self.assertAlmostEqual(result["records"][:, :, 0].sum(), (1 + math.sin(phase)) / 2)

    def test_joint_preparation_mixtures_and_common_J_remain_consistent(self):
        mode = RadialMode(RadialProblem())
        gate = evolve_probe(mode, sites=8)["gate"]
        first, second = np.eye(8)[:, 0], np.eye(8)[:, 2]
        def density(data):
            state = gate @ np.kron(data, np.ones(4) / 2)
            return np.outer(state, state.conj())
        np.testing.assert_allclose((density(first) + density(second)) / 2,
                                   sum(density((first + sign * second) / math.sqrt(2)) for sign in (1, -1)) / 2, atol=3e-14)
        initial = np.kron(initial_data(8), np.ones(4) / 2)
        lifted = real_lift(gate)
        np.testing.assert_allclose(lifted @ encode_state(np.outer(initial, initial.conj())) @ lifted.T,
                                   encode_state(density(initial_data(8))), atol=3e-14)

    def test_complete_record_delivery_and_invalid_mode_inputs(self):
        positions = [0] * 32
        positions[0] = 1
        schedule = record_schedule(32, 2.)
        self.assertEqual(schedule[-1]["received_at"], 274.)
        with self.assertRaises(ValueError):
            release_record(positions, [0, 1], 0, 273.9, 2.)
        self.assertEqual(release_record(positions, [0, 1], 0, 274., 2.)[1:], (0, 1))
        with self.assertRaises(ValueError):
            RadialProblem().find_mode(3)
        with self.assertRaises(ValueError):
            RadialProblem().ritz(2)

    def test_default_stability_sign_has_an_exact_rational_certificate(self):
        certificate = rational_default_certificate()
        self.assertGreater(Fraction(certificate["exact_lower_fraction"]), Fraction(3, 200))
        self.assertGreater(RadialProblem().sufficient_stability_bound()["omega_squared_lower"], float(Fraction(certificate["exact_lower_fraction"])))

    def test_tighter_background_and_shooting_agree_on_the_fundamental_frequency(self):
        current = RadialProblem().find_mode()
        refined = RadialProblem(FluidSolution(relative_tolerance=1e-12, start_fraction=5e-7)).find_mode(
            relative_tolerance=1e-12, start_fraction=5e-6)
        self.assertLess(abs(current["omega_squared"] - refined["omega_squared"]), 1e-8)

    def test_reconstructed_fields_satisfy_linear_Euler_and_mass_constraints(self):
        diagnostics = linear_constraint_diagnostics(RadialMode(RadialProblem()))
        self.assertLess(diagnostics["maximum_sampled_linear_Euler_residual"], 1e-7)
        self.assertLess(diagnostics["maximum_sampled_linear_mass_constraint_residual"], 1e-6)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    outcome = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RadialFluidModeTests))
    if not outcome.wasSuccessful():
        raise SystemExit(1)
    if arguments.write_results:
        report = run_report()
        Path(__file__).with_name("radial_fluid_modes_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"tests": outcome.testsRun, "bound": report["analytic_stability_bound"], "modes": report["modes"],
                          "conservation": report["linear_conservation"], "probe": report["probe_task"], "refinement": report["fundamental_refinement"]}, indent=2))