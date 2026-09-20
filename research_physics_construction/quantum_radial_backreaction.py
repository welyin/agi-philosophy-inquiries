"""Construction 12: canonical radial mode and reciprocal finite quantum model.

Gravity, the fluid EOS and quantum rules remain declared inputs. A homology
scale is a physical parameter change, not a new value of Planck's constant.
"""

import argparse
import json
import math
import unittest
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gammaln, log_ndtr
from scipy.stats import poisson

from radial_fluid_modes import (RadialMode, RadialProblem, dynamic_probe_terms,
                                evolve_probe, record_schedule, release_record)
from quantum_chain_background import oscillator_operators, initial_data, trace_distance
from massive_envelope import postprocessing_table
from clock_calibration import clock_effect
from baseline import IDENTITY, real_lift, checked_probability
from common_orientation_structure import encode_state


@lru_cache(maxsize=1)
def reference_mode():
    return RadialMode(RadialProblem())


def canonical_inertia(mode, quadrature_order=256):
    nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
    radii = mode.radius * (nodes + 1) / 2
    weights *= mode.radius / 2
    values = []
    for radius in radii:
        row = mode.problem.equilibrium(float(radius))
        displacement = mode.local(float(radius))["displacement"]
        values.append(4 * math.pi * radius ** 2 * (row["density"] + row["pressure"])
                      * row["radial_factor"] ** 1.5 / row["lapse"] * displacement ** 2)
    return float(np.dot(weights, values))


@dataclass(frozen=True)
class Parameters:
    scale: float = 1024.
    sites: int = 8
    cutoff: int = 20
    coherent_amplitude: float = 0.00025
    validity_radius: float = 0.01
    interaction_scale: float = 1.0

    def __post_init__(self):
        if (not np.isfinite(self.scale) or self.scale < 1
                or not isinstance(self.sites, int) or self.sites < 8 or self.sites % 4
                or not isinstance(self.cutoff, int) or self.cutoff < 4
                or not np.isfinite(self.coherent_amplitude)
                or not np.isfinite(self.validity_radius) or not 0 < self.validity_radius <= 0.01
                or not np.isfinite(self.interaction_scale) or not 0 <= self.interaction_scale <= 1
                or abs(self.coherent_amplitude) >= self.validity_radius):
            raise ValueError("Use scale>=1, a valid inherited ring, cutoff>=4 and a small finite mode amplitude.")


def normalization_report(parameters=Parameters()):
    mode = reference_mode()
    inertia = canonical_inertia(mode) * parameters.scale ** 3
    frequency = mode.omega / parameters.scale
    deviation = math.sqrt(1 / (2 * inertia * frequency))
    vacuum_outside = math.erfc(parameters.validity_radius / (math.sqrt(2) * deviation))
    log_vacuum_outside = math.log(2) + float(log_ndtr(-parameters.validity_radius / deviation))
    return {"physical_inertia": inertia, "frequency": frequency,
            "vacuum_fractional_radius_standard_deviation": deviation,
            "vacuum_probability_outside_declared_linear_radius": vacuum_outside,
            "log_vacuum_probability_outside_declared_linear_radius": log_vacuum_outside,
            "display_probability_underflow_is_not_zero_physical_tail": vacuum_outside == 0.,
            "coherent_alpha": parameters.coherent_amplitude / (2 * deviation),
            "scope": "canonical harmonic reduction at hbar=1; not a proof of full nonlinear quantum-fluid validity"}


def kinetic_ADM_variation(mode):
    radius = mode.minimum_radius
    row = mode.problem.equilibrium(radius)
    fractional = mode.local(radius)["fractional_displacement"]
    initial = 2 * math.pi * (row["density"] + row["pressure"]) * fractional ** 2 * radius ** 5 / (5 * row["lapse"] ** 2)

    def derivative(position, variation):
        background = mode.problem.equilibrium(position)
        displacement = mode.local(position)["displacement"]
        total = background["density"] + background["pressure"]
        coefficient = 4 * math.pi * mode.solution.gravitational_constant * position * background["radial_factor"] * total
        local_velocity_squared = background["radial_factor"] / background["lapse"] ** 2 * displacement ** 2
        return [2 * math.pi * position ** 2 * total * local_velocity_squared - coefficient * variation[0]]

    integration = solve_ivp(derivative, (radius, mode.radius), [initial], method="DOP853",
                            rtol=1e-11, atol=1e-12, max_step=mode.radius / 80)
    if not integration.success:
        raise ValueError("The fixed-particle-number kinetic mass-constraint check failed.")
    return float(integration.y[0, -1])


class QuantumRadialModel:
    def __init__(self, parameters=Parameters()):
        self.parameters = parameters
        self.mode = reference_mode()
        self.normalization = normalization_report(parameters)
        self.inertia = self.normalization["physical_inertia"]
        self.frequency = self.normalization["frequency"]
        self.sigma = self.normalization["vacuum_fractional_radius_standard_deviation"]
        self.terms = dynamic_probe_terms(self.mode, parameters.sites, 0.35, parameters.validity_radius)
        operators = oscillator_operators(parameters.cutoff)
        self.number = operators["number"]
        self.position = math.sqrt(2) * self.sigma * operators["position"]
        self.position_squared = 2 * self.sigma ** 2 * operators["position_squared_projected"]
        self.momentum = operators["momentum"] / (math.sqrt(2) * self.sigma)
        self.blocks = []
        sites, cutoff, scale = parameters.sites, parameters.cutoff, parameters.scale
        clock_rates = self.terms["clock_rates"] / scale
        clock_response = clock_rates * self.terms["clock_lapse_fractional_derivatives"]
        for first, second in ((0, 0), (0, 1), (1, 0), (1, 1)):
            signs = np.array([1 - 2 * first, 1 - 2 * second])
            matter_small = self.terms["static_hamiltonian"] / scale + (signs @ clock_rates / 2) * np.eye(sites)
            source = parameters.interaction_scale * (self.terms["mode_hamiltonian_derivative"] / scale
                                                    + (signs @ clock_response / 2) * np.eye(sites))
            matter = np.kron(matter_small, np.eye(cutoff))
            oscillator = np.kron(np.eye(sites), self.frequency * (self.number + np.eye(cutoff) / 2))
            interaction = np.kron(source, self.position)
            hamiltonian = matter + oscillator + interaction
            values, vectors = np.linalg.eigh(scale * hamiltonian)
            self.blocks.append({"matter": matter, "mode": oscillator, "interaction": interaction,
                                "hamiltonian": hamiltonian, "matter_small": matter_small, "source": source,
                                "values": values, "vectors": vectors})

    def coherent_vector(self):
        alpha = self.normalization["coherent_alpha"]
        levels = np.arange(self.parameters.cutoff)
        if alpha == 0:
            vector = np.zeros(self.parameters.cutoff)
            vector[0] = 1.
        else:
            vector = np.exp(-alpha ** 2 / 2 + levels * math.log(abs(alpha)) - gammaln(levels + 1) / 2) * np.sign(alpha) ** levels
            vector /= np.linalg.norm(vector)
        return vector

    def initial_vector(self, data=None, mode_vector=None):
        data = initial_data(self.parameters.sites) if data is None else np.asarray(data, dtype=complex)
        mode_vector = self.coherent_vector() if mode_vector is None else np.asarray(mode_vector, dtype=complex)
        if (data.shape != (self.parameters.sites,) or mode_vector.shape != (self.parameters.cutoff,)
                or not np.isclose(np.linalg.norm(data), 1, rtol=0, atol=1e-12)
                or not np.isclose(np.linalg.norm(mode_vector), 1, rtol=0, atol=1e-12)):
            raise ValueError("Use normalized data and mode vectors in the declared tensor factors.")
        return data[:, None, None] * mode_vector[None, :, None] * np.ones((1, 1, 4)) / 2

    def evolve(self, initial, reference_duration):
        shape = (self.parameters.sites, self.parameters.cutoff, 4)
        if initial.shape != shape or not np.isfinite(reference_duration):
            raise ValueError("Use the full joint vector and a finite reference time t/scale.")
        if reference_duration == 0:
            return initial.copy()
        output = np.empty_like(initial, dtype=complex)
        for index, block in enumerate(self.blocks):
            coefficients = block["vectors"].T @ initial[..., index].reshape(-1)
            output[..., index] = (block["vectors"] @ (np.exp(-1j * reference_duration * block["values"]) * coefficients)).reshape(shape[:2])
        return output

    def energies(self, state):
        return {term: float(sum(np.vdot(state[..., index].reshape(-1), block[term] @ state[..., index].reshape(-1)).real
                                for index, block in enumerate(self.blocks)))
                for term in ("matter", "mode", "interaction", "hamiltonian")}

    def reduced_states(self, state):
        data = state.reshape(self.parameters.sites, -1)
        oscillator = state.transpose(1, 0, 2).reshape(self.parameters.cutoff, -1)
        clocks = state.transpose(2, 0, 1).reshape(4, -1)
        return {name: vector @ vector.conj().T for name, vector in (("data", data), ("mode", oscillator), ("clocks", clocks))}

    def diagnostics(self, state):
        reduced = self.reduced_states(state)
        return {"energies": self.energies(state), "norm": float(np.linalg.norm(state)),
                "mode_mean_amplitude": float(np.trace(reduced["mode"] @ self.position).real),
                "mode_second_moment": float(np.trace(reduced["mode"] @ self.position_squared).real),
                "mode_number": float(np.trace(reduced["mode"] @ self.number).real),
                "mode_purity": float(np.trace(reduced["mode"] @ reduced["mode"]).real),
                "clock_A_X_plus_probability": checked_probability(reduced["clocks"], np.kron(clock_effect("X", 0), IDENTITY)),
                "clock_C_Y_plus_probability": checked_probability(reduced["clocks"], np.kron(IDENTITY, clock_effect("Y", 0)))}

    def record_probabilities(self, state):
        thresholds, denominator = postprocessing_table(self.parameters.sites, center_fraction=0.5)
        probabilities = np.zeros((2, 2, 2, 2))
        vacuum = np.zeros(self.parameters.cutoff)
        vacuum[0] = 1.
        for data_bit, mode_bit, first_clock, second_clock in np.ndindex(probabilities.shape):
            data_weight = thresholds / denominator if data_bit == 0 else 1 - thresholds / denominator
            mode_weight = vacuum if mode_bit == 0 else 1 - vacuum
            clock = np.kron(clock_effect("X", first_clock), clock_effect("Y", second_clock))
            probabilities[data_bit, mode_bit, first_clock, second_clock] = np.einsum(
                "jka,ab,jkb,j,k->", state.conj(), clock, state, data_weight, mode_weight).real
        if probabilities.min() < -1e-12 or abs(probabilities.sum() - 1) > 1e-12:
            raise ValueError("Joint records must retain all outcomes with positive normalized probabilities.")
        return probabilities

    def cutoff_budget(self, initial, reference_duration=2., samples=256):
        if reference_duration < 0 or not np.isfinite(reference_duration) or not isinstance(samples, int) or samples < 1:
            raise ValueError("Use a nonnegative finite horizon and positive quadrature count.")
        sites, cutoff = self.parameters.sites, self.parameters.cutoff
        times = reference_duration * (np.arange(samples) + 0.5) / samples
        squared, lipschitz_squared = np.zeros(samples), 0.
        for index, block in enumerate(self.blocks):
            exterior = np.zeros((sites, sites * cutoff))
            exterior[:, np.arange(sites) * cutoff + cutoff - 1] = self.parameters.scale * self.sigma * math.sqrt(cutoff) * block["source"]
            coefficients = block["vectors"].T @ initial[..., index].reshape(-1)
            exterior_modes = exterior @ block["vectors"]
            weight = float(np.vdot(coefficients, coefficients).real)
            shift = float(np.dot(abs(coefficients) ** 2, block["values"]) / weight) if weight else 0.
            frequencies = block["values"] - shift
            residual = exterior_modes @ (coefficients[:, None] * np.exp(-1j * np.outer(frequencies, times)))
            squared += np.sum(abs(residual) ** 2, axis=0)
            lipschitz = float(np.sum(abs(coefficients) * abs(frequencies) * np.linalg.norm(exterior_modes, axis=0)))
            lipschitz_squared += lipschitz ** 2
        integral = reference_duration * float(np.mean(np.sqrt(squared)))
        remainder = math.sqrt(lipschitz_squared) * reference_duration ** 2 / (4 * samples)
        return {"boundary_residual_integral": integral, "midpoint_remainder_bound": remainder,
                "same_initial_embedded_vector_error_bound": integral + remainder,
                "scope": "exact finite-spectrum Duhamel formula evaluated in floating point, not directed rounding"}

    def coherent_initial_budget(self):
        tail = float(poisson.sf(self.parameters.cutoff - 1, self.normalization["coherent_alpha"] ** 2))
        error = math.sqrt(2 * tail / (1 + math.sqrt(1 - tail)))
        return {"discarded_coherent_number_probability": tail, "initial_vector_error_bound": error}

    def validity_budget(self, data=None):
        data = initial_data(self.parameters.sites) if data is None else np.asarray(data)
        amplitude = self.parameters.coherent_amplitude
        energy = self.frequency / 2 + self.inertia * self.frequency ** 2 * amplitude ** 2 / 2
        energy += sum(np.vdot(data, (block["matter_small"] + amplitude * block["source"]) @ data).real for block in self.blocks) / 4
        matter_floor = min(float(np.linalg.eigvalsh(block["matter_small"]).min()) for block in self.blocks)
        source_norm = max(float(np.linalg.norm(block["source"], ord=2)) for block in self.blocks)
        restoring = self.inertia * self.frequency ** 2
        second_moment = 4 * (energy - matter_floor + source_norm ** 2 / restoring) / restoring
        return {"ideal_initial_total_energy": float(energy), "matter_energy_floor": matter_floor,
                "source_operator_norm": source_norm, "all_time_amplitude_second_moment_upper_bound": second_moment,
                "pointwise_time_probability_outside_linear_window_upper_bound": min(1., second_moment / self.parameters.validity_radius ** 2),
                "scope": "energy/Markov bound in the unbounded linear oscillator completion for the ideal coherent initial state; not a probability of never exiting nor a bound on omitted nonlinear terms"}


@lru_cache(maxsize=4)
def model(parameters=Parameters()):
    return QuantumRadialModel(parameters)


def embedded(state, cutoff):
    if cutoff < state.shape[1]:
        raise ValueError("The comparison embedding cannot discard levels.")
    output = np.zeros((state.shape[0], cutoff, 4), dtype=complex)
    output[:, :state.shape[1], :] = state
    return output


def total_record_budget(current, initial, reference_duration=2.):
    evolution = current.cutoff_budget(initial, reference_duration)
    preparation = current.coherent_initial_budget()
    return {"evolution": evolution, "coherent_preparation": preparation,
            "same_final_record_TV_upper_bound": min(1., evolution["same_initial_embedded_vector_error_bound"] + preparation["initial_vector_error_bound"]),
            "scope": "ideal unbounded linear-mode completion with the same final POVM; not the omitted nonlinear-fluid or other-mode error"}


def terminal_schedule(parameters=Parameters(), reference_duration=2.):
    schedule = record_schedule(parameters.sites, reference_duration * parameters.scale, parameters.scale)
    schedule.append({"kind": "collective_mode_vacuum_bit", "source_site": 0, "hops": 0,
                     "read_at": reference_duration * parameters.scale,
                     "received_at": reference_duration * parameters.scale})
    return sorted(schedule, key=lambda event: event["received_at"])


def release_joint_record(positions, clocks, mode_bit, random_integer, receipt, parameters=Parameters(), reference_duration=2.):
    if mode_bit not in (0, 1) or len(positions) != parameters.sites:
        raise ValueError("Use the complete declared ring and a retained binary collective-mode readout.")
    data_bit, first, second = release_record(positions, clocks, random_integer, receipt,
                                           reference_duration * parameters.scale, parameters.scale)
    return data_bit, mode_bit, first, second


def run_report():
    current = model()
    initial = current.initial_vector()
    final = current.evolve(initial, 2.)
    budget = total_record_budget(current, initial)
    larger = model(Parameters(cutoff=24))
    larger_initial = larger.initial_vector()
    larger_final = larger.evolve(larger_initial, 2.)
    refinements = []
    for cutoff in (12, 16, 20):
        candidate = model(Parameters(cutoff=cutoff))
        prepared = candidate.initial_vector()
        evolved = candidate.evolve(prepared, 2.)
        refinements.append({"cutoff": cutoff,
                            "vector_difference_from_D24": float(np.linalg.norm(embedded(evolved, 24) - larger_final)),
                            "budget": total_record_budget(candidate, prepared)})
    changed_data = current.evolve(current.initial_vector(data=initial_data(current.parameters.sites, -1)), 2.)
    negative_mode = current.coherent_vector() * (-1.) ** np.arange(current.parameters.cutoff)
    changed_mode = current.evolve(current.initial_vector(mode_vector=negative_mode), 2.)
    reduced = current.reduced_states(final)
    frozen = evolve_probe(current.mode, sites=current.parameters.sites, amplitude=current.parameters.coherent_amplitude)
    records = current.record_probabilities(final)
    source_norm = max(float(np.linalg.norm(block["source"], ord=2)) for block in current.blocks)
    schedule = terminal_schedule(current.parameters)
    scale = current.parameters.scale
    return {"model_version": "C12", "route_round": "12", "research_scope": "sufficiency only",
            "parameters": vars(current.parameters),
            "canonical_normalization": {"reference_SL_inertia": current.mode.conservation_report()["modal_inertia"],
                                        "reference_physical_inertia": canonical_inertia(current.mode),
                                        "reference_half_inertia": canonical_inertia(current.mode) / 2,
                                        "independent_ADM_kinetic_mass_variation": kinetic_ADM_variation(current.mode),
                                        "action": "one-half physical_inertia*(amplitude_dot^2-frequency^2*amplitude^2); physical_inertia=4pi integral W zeta_mode^2 dr",
                                        "Newtonian_limit": "physical_inertia -> 4pi integral rho r^2 xi_mode^2 dr"},
            "unscaled_linear_window_check": normalization_report(Parameters(scale=1.)),
            "scaled_normalization": current.normalization,
            "physical_homology_inputs": {"scale": scale, "G_unchanged": 0.1, "hbar_unchanged": 1.,
                                       "surface_density": 0.01 / scale ** 2, "central_pressure": 0.0005 / scale ** 2,
                                       "sound_speed_squared_unchanged": 1 / 3,
                                       "source_radius": current.mode.radius * scale, "source_ADM_mass": current.mode.solution.mass * scale,
                                       "ring_length": 20 * scale, "ring_offset": 4 * scale,
                                       "probe_proper_mass": 0.35 / scale, "bare_clock_rates": [1 / scale, 1.2 / scale],
                                       "read_parameter": 2 * scale,
                                       "note": "physical source and probe family changed at fixed G and hbar; not a mere coordinate-unit change"},
            "joint_Hamiltonian": "H_matter + P_amplitude^2/(2I) + I*omega^2*A^2/2 + A*B_matter",
            "time_rows": [dict(reference_time=time, physical_time=time * scale, **current.diagnostics(current.evolve(initial, time))) for time in (0., 0.5, 1., 2.)],
            "record_budget": budget, "D24_reference_own_budget": total_record_budget(larger, larger_initial),
            "cutoff_refinements": refinements, "linear_window_budget": current.validity_budget(),
            "source_displacement_operator_norm_bound": source_norm / (current.inertia * current.frequency ** 2),
            "reciprocal_diagnostics": {"changed_data_mode_trace_distance": trace_distance(reduced["mode"], current.reduced_states(changed_data)["mode"]),
                                       "changed_coherent_sign_data_trace_distance": trace_distance(reduced["data"], current.reduced_states(changed_mode)["data"]),
                                       "changed_coherent_sign_clock_trace_distance": trace_distance(reduced["clocks"], current.reduced_states(changed_mode)["clocks"]),
                                       "scope": "finite-cutoff diagnostics for specified preparations; default error budget not automatically reused for interventions"},
            "joint_records_data_mode_clockA_clockC": records.tolist(),
            "classical_reference_data_clock_records": frozen["records"].tolist(),
            "data_clock_record_TV_from_classical_mean_drive": float(abs(records.sum(axis=1) - frozen["records"]).sum() / 2),
            "resources": {"data_qubits": current.parameters.sites, "clock_qubits": 2,
                          "retained_collective_radial_modes": 1, "mode_cutoff": current.parameters.cutoff,
                          "occupied_complex_sector_dimension": 4 * current.parameters.sites * current.parameters.cutoff,
                          "raw_bits": current.parameters.sites + 3, "coarse_joint_outcomes": 16,
                          "independent_fair_output_bits": 12, "threshold_table_bits": current.parameters.sites * 13,
                          "bit_hops": sum(event["hops"] for event in schedule), "supplied_hop_interval": scale,
                          "all_records_ready_at": max(event["received_at"] for event in schedule), "resets": 0,
                          "full_fluid_microphysical_resource_count_closed": False,
                          "collective_mode_detector_and_transport_realization_closed": False},
            "terminal_schedule": schedule,
            "limits": ["the unscaled C10 source fails the declared weak-amplitude vacuum check and is retained as an invalid direct quantization case",
                       "the macroscopic homology, probe rescaling, initial coherent state and single-mode reduction are explicit C12 inputs",
                       "canonical kinetic normalization is checked through fixed-baryon ADM constraint variation, but the full second-order background metric is not constructed",
                       "energy is closed only within the autonomous linear single-mode joint Hamiltonian, not all Einstein constraints or all fluid modes",
                       "escape bound holds at each time separately and does not assert the state never leaves the linear window",
                       "finite probability outside the window is retained; no state projection or postselection enforces validity",
                       "Fock error budget does not cover omitted nonlinear metric terms, other radial/nonradial modes, or material interactions",
                       "one collective radial observable is supplied at the collector, not derived as a spatially local instantaneous measurement",
                       "no laboratory samples, necessity arguments, or derivation of the adopted gravitational action from cognition"]}


class QuantumRadialBackreactionTests(unittest.TestCase):
    def test_redshifted_fluid_kinetic_integral_fixes_the_angular_normalization(self):
        mode = reference_mode()
        independent = canonical_inertia(mode)
        sturm_liouville = mode.conservation_report()["modal_inertia"]
        self.assertAlmostEqual(independent, 4 * math.pi * sturm_liouville)
        self.assertGreater(independent, 90.)

    def test_unscaled_vacuum_is_not_inside_the_C11_linear_amplitude_window(self):
        result = normalization_report(Parameters(scale=1.))
        self.assertGreater(result["vacuum_fractional_radius_standard_deviation"], 0.01)
        self.assertGreater(result["vacuum_probability_outside_declared_linear_radius"], 0.9)

    def test_macroscopic_homology_reduces_zero_point_amplitude_without_changing_hbar(self):
        original, scaled = normalization_report(Parameters(scale=1.)), normalization_report()
        scale = Parameters().scale
        self.assertAlmostEqual(scaled["frequency"] * scale, original["frequency"])
        self.assertAlmostEqual(scaled["physical_inertia"] / scale ** 3, original["physical_inertia"])
        self.assertAlmostEqual(scaled["vacuum_fractional_radius_standard_deviation"] * scale,
                               original["vacuum_fractional_radius_standard_deviation"])
        self.assertLess(scaled["vacuum_fractional_radius_standard_deviation"], 0.0002)

    def test_homology_maps_the_TOV_source_and_radial_inertia_consistently(self):
        mode = reference_mode()
        scale = 7.
        radius = mode.radius * 0.4
        row = mode.problem.equilibrium(radius)
        gravity = mode.solution.gravitational_constant
        original = mode.solution.derivative(radius, [row["mass"], row["pressure"], 0.])
        scaled_radius, scaled_mass, scaled_pressure = radius * scale, row["mass"] * scale, row["pressure"] / scale ** 2
        scaled_density = row["density"] / scale ** 2
        scaled_nu_first = gravity * (scaled_mass + 4 * math.pi * scaled_radius ** 3 * scaled_pressure) / (scaled_radius * (scaled_radius - 2 * gravity * scaled_mass))
        np.testing.assert_allclose([4 * math.pi * scaled_radius ** 2 * scaled_density,
                                   -(scaled_density + scaled_pressure) * scaled_nu_first, scaled_nu_first],
                                  original / np.array([1., scale ** 3, scale]), atol=1e-15)

    def test_constrained_ADM_kinetic_variation_matches_half_the_canonical_inertia(self):
        mode = reference_mode()
        self.assertLess(abs(kinetic_ADM_variation(mode) - canonical_inertia(mode) / 2), 1e-7)

    def test_joint_evolution_conserves_total_not_separate_mode_energy(self):
        current = model()
        initial = current.initial_vector()
        final = current.evolve(initial, 2.)
        np.testing.assert_allclose(current.evolve(final, -2.), initial, atol=3e-14)
        first, last = current.energies(initial), current.energies(final)
        self.assertAlmostEqual(first["hamiltonian"], last["hamiltonian"], places=14)
        self.assertAlmostEqual(last["hamiltonian"], last["matter"] + last["mode"] + last["interaction"], places=14)
        self.assertGreater(abs(last["mode"] - first["mode"]), 1e-12)

    def test_cutoff_projection_has_only_the_correct_top_level_linear_leakage(self):
        current = model(Parameters(cutoff=8))
        extended = model(Parameters(cutoff=9))
        generator = np.random.default_rng(12)
        state = generator.normal(size=(8, 8)) + 1j * generator.normal(size=(8, 8))
        embedded = np.zeros((8, 9), dtype=complex)
        embedded[:, :8] = state
        for block, larger in zip(current.blocks, extended.blocks):
            actual = (larger["hamiltonian"] @ embedded.reshape(-1)).reshape(8, 9)
            np.testing.assert_allclose(actual[:, :8].reshape(-1), block["hamiltonian"] @ state.reshape(-1), atol=1e-15)
            np.testing.assert_allclose(actual[:, 8], current.sigma * math.sqrt(8) * block["source"] @ state[:, 7], atol=1e-15)

    def test_full_joint_records_preserve_clock_and_mode_correlations(self):
        current = model()
        state = current.evolve(current.initial_vector(), 2.)
        records = current.record_probabilities(state)
        self.assertAlmostEqual(records.sum(), 1.)
        self.assertAlmostEqual(records[:, :, :, 0].sum(), current.diagnostics(state)["clock_C_Y_plus_probability"])
        joint = records.sum(axis=(0, 2))
        self.assertGreater(np.linalg.norm(joint - np.outer(joint.sum(axis=1), joint.sum(axis=0))), 1e-8)

    def test_macro_scale_has_an_all_time_small_window_escape_allowance(self):
        current = model()
        budget = current.validity_budget()
        self.assertLess(budget["pointwise_time_probability_outside_linear_window_upper_bound"], 0.02)
        for duration in (0., 1., 2.):
            moment = current.diagnostics(current.evolve(current.initial_vector(), duration))["mode_second_moment"]
            self.assertLess(moment, budget["all_time_amplitude_second_moment_upper_bound"])

    def test_two_cutoffs_agree_with_their_independent_preparation_and_evolution_budgets(self):
        first, second = model(Parameters(cutoff=12)), model(Parameters(cutoff=20))
        initial_first, initial_second = first.initial_vector(), second.initial_vector()
        actual = np.linalg.norm(embedded(first.evolve(initial_first, 2.), 20) - second.evolve(initial_second, 2.))
        allowance = sum(total_record_budget(current, initial)["same_final_record_TV_upper_bound"]
                        for current, initial in ((first, initial_first), (second, initial_second)))
        self.assertLess(actual, allowance)

    def test_zero_interaction_keeps_the_mode_free_and_recovers_the_static_scaled_probe(self):
        current = model(Parameters(cutoff=8, interaction_scale=0.))
        initial = current.initial_vector()
        evolved = current.evolve(initial, 2.)
        for index, block in enumerate(current.blocks):
            values, vectors = np.linalg.eigh(current.parameters.scale * block["matter_small"])
            data = vectors @ (np.exp(-2j * values) * (vectors.T @ initial_data(8)))
            oscillator = np.exp(-2j * current.parameters.scale * current.frequency * (np.arange(8) + 0.5)) * current.coherent_vector()
            np.testing.assert_allclose(evolved[..., index], data[:, None] * oscillator[None, :] / 2, atol=3e-14)

    def test_both_directions_of_mode_matter_response_are_present(self):
        current = model()
        first = current.reduced_states(current.evolve(current.initial_vector(), 2.))
        changed_data = current.reduced_states(current.evolve(current.initial_vector(data=initial_data(8, -1)), 2.))
        changed_mode_vector = current.coherent_vector() * (-1.) ** np.arange(current.parameters.cutoff)
        changed_mode = current.reduced_states(current.evolve(current.initial_vector(mode_vector=changed_mode_vector), 2.))
        self.assertGreater(trace_distance(first["mode"], changed_data["mode"]), 1e-6)
        self.assertGreater(trace_distance(first["clocks"], changed_mode["clocks"]), 1e-6)

    def test_equivalent_ensembles_produce_the_same_complete_joint_density(self):
        current = model(Parameters(cutoff=6))
        first, second = np.eye(8)[:, 0], np.eye(8)[:, 2]
        def output(data):
            final = current.evolve(current.initial_vector(data), 1.).reshape(-1)
            return np.outer(final, final.conj())
        np.testing.assert_allclose((output(first) + output(second)) / 2,
                                   sum(output((first + sign * second) / math.sqrt(2)) for sign in (1, -1)) / 2, atol=3e-14)

    def test_coherent_sector_evolution_has_one_common_J_real_representation(self):
        current = model(Parameters(cutoff=4))
        initial = current.initial_vector().reshape(-1)
        gate = np.zeros((len(initial), len(initial)), dtype=complex)
        for index, block in enumerate(current.blocks):
            indices = np.arange(index, len(initial), 4)
            gate[np.ix_(indices, indices)] = (block["vectors"] * np.exp(-1j * block["values"])) @ block["vectors"].T
        final = current.evolve(current.initial_vector(), 1.).reshape(-1)
        np.testing.assert_allclose(gate @ initial, final, atol=3e-14)
        lifted = real_lift(gate)
        np.testing.assert_allclose(lifted @ encode_state(np.outer(initial, initial.conj())) @ lifted.T,
                                   encode_state(np.outer(final, final.conj())), atol=3e-14)

    def test_full_raw_records_and_random_postprocessing_have_scaled_delivery_cost(self):
        parameters = Parameters()
        schedule = terminal_schedule(parameters)
        self.assertEqual(len(schedule), 11)
        self.assertEqual(sum(event["hops"] for event in schedule), 20)
        deadline = max(event["received_at"] for event in schedule)
        self.assertEqual(deadline, 22528.)
        positions = [0] * 8
        positions[3] = 1
        with self.assertRaises(ValueError):
            release_joint_record(positions, [0, 1], 1, 0, deadline - 1)
        threshold = int(postprocessing_table(8, center_fraction=0.5)[0][3])
        for random in (0, threshold - 1, threshold, 4095):
            self.assertEqual(release_joint_record(positions, [0, 1], 1, random, deadline),
                             (int(random >= threshold), 1, 0, 1))

    def test_invalid_scales_and_amplitudes_are_rejected(self):
        for settings in ({"scale": 0.}, {"cutoff": 3}, {"coherent_amplitude": 0.02}, {"interaction_scale": -1.}):
            with self.assertRaises(ValueError):
                Parameters(**settings)

    def test_small_vacuum_tail_is_reported_in_log_form_without_claiming_exact_zero(self):
        report = normalization_report()
        self.assertTrue(np.isfinite(report["log_vacuum_probability_outside_declared_linear_radius"]))
        self.assertLess(report["log_vacuum_probability_outside_declared_linear_radius"], -100)
        if report["vacuum_probability_outside_declared_linear_radius"] == 0:
            self.assertTrue(report["display_probability_underflow_is_not_zero_physical_tail"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    outcome = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(QuantumRadialBackreactionTests))
    if not outcome.wasSuccessful():
        raise SystemExit(1)
    if arguments.write_results:
        report = run_report()
        Path(__file__).with_name("quantum_radial_backreaction_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"tests": outcome.testsRun, "normalization": report["canonical_normalization"],
                          "unscaled": report["unscaled_linear_window_check"], "scaled": report["scaled_normalization"],
                          "time_rows": report["time_rows"], "record_budget": report["record_budget"],
                          "window_budget": report["linear_window_budget"], "response": report["reciprocal_diagnostics"]}, indent=2))