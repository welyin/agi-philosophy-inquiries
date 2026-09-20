"""Construction 07: two adjacent background modes on the same massive chain.

Local mode rules, exchange, clock relocation and the optional hopping hold
are declared inputs. This finite patch is not an Einstein field theory.
"""

import argparse
import json
import math
import unittest
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

from baseline import IDENTITY, checked_probability, real_lift
from clock_calibration import clock_effect
from common_orientation_structure import encode_state
from envelope_limit import grid, hopping_matrix
from massive_envelope import BACKGROUND, DEFAULT_MASS, postprocessing_table, evolve_lattice
from quantum_chain_background import oscillator_operators, initial_data, trace_distance


@dataclass(frozen=True)
class Parameters:
    sites: int = 8
    cutoff: int = 8
    frequency: float = 1.0
    coupling: float = 0.15
    exchange: float = 0.3
    mass: float = DEFAULT_MASS
    hopping_scale: float = 1.0

    def __post_init__(self):
        grid(self.sites, BACKGROUND)
        values = (self.frequency, self.coupling, self.exchange, self.mass, self.hopping_scale)
        if (not isinstance(self.cutoff, int) or self.cutoff < 4
                or not all(np.isfinite(value) for value in values)
                or self.frequency <= 0 or not 0 <= self.coupling < 2
                or abs(self.exchange) >= self.frequency or self.mass < 0
                or not 0 <= self.hopping_scale <= 1):
            raise ValueError("Require cutoff>=4, positive frequency, |exchange|<frequency, 0<=coupling<2 and valid mass/hopping.")


def local_data_terms(parameters):
    positions, _ = grid(parameters.sites, BACKGROUND)
    hopping = parameters.hopping_scale * hopping_matrix(parameters.sites, BACKGROUND)
    diagonal = (-1.) ** np.arange(parameters.sites) * parameters.mass * BACKGROUND.lapse(positions)
    sources = []
    for center in (0, 1):
        source = np.zeros_like(hopping)
        source[center, center] = diagonal[center]
        for neighbor in ((center - 1) % parameters.sites, (center + 1) % parameters.sites):
            source[center, neighbor] = source[neighbor, center] = hopping[center, neighbor] / 2
        sources.append(source)
    return hopping + np.diag(diagonal), sources


class LocalQuantumBackground:
    def __init__(self, parameters):
        self.parameters = parameters
        sites, cutoff = parameters.sites, parameters.cutoff
        self.data_hamiltonian, self.sources = local_data_terms(parameters)
        positions, _ = grid(sites, BACKGROUND)
        self.clock_rates = np.array([BACKGROUND.lapse(positions[0]), 1.2 * BACKGROUND.lapse(positions[1])])
        self.oscillator = oscillator_operators(cutoff)
        unit = np.eye(cutoff)
        lowering = self.oscillator["lowering"]
        number = self.oscillator["number"]
        self.mode_numbers = (np.kron(number, unit), np.kron(unit, number))
        self.mode_fields = (np.kron(self.oscillator["centered_square"], unit),
                            np.kron(unit, self.oscillator["centered_square"]))
        background_free = parameters.frequency * sum(self.mode_numbers)
        background_link = parameters.exchange * (np.kron(lowering.T, lowering) + np.kron(lowering, lowering.T))
        self.background_hamiltonian = background_free + background_link
        self.blocks = []
        for first_bit, second_bit in ((0, 0), (0, 1), (1, 0), (1, 1)):
            signs = np.array([1 - 2 * first_bit, 1 - 2 * second_bit])
            clock_energy = float(np.dot(signs, self.clock_rates) / 2)
            sources = [source + sign * rate / 2 * np.eye(sites)
                       for source, sign, rate in zip(self.sources, signs, self.clock_rates)]
            terms = {"matter": np.kron(self.data_hamiltonian + clock_energy * np.eye(sites), np.eye(cutoff ** 2)),
                     "background_free": np.kron(np.eye(sites), background_free),
                     "background_link": np.kron(np.eye(sites), background_link),
                     "interaction": parameters.coupling * sum(np.kron(source, field) for source, field in zip(sources, self.mode_fields))}
            total = sum(terms.values())
            values, vectors = np.linalg.eigh(total)
            self.blocks.append(dict(terms, hamiltonian=total, sources=sources, clock_energy=clock_energy,
                                    values=values, vectors=vectors))
        source_norm = max(np.linalg.norm(source, ord=2) for block in self.blocks for source in block["sources"])
        self.stability_margin = float(parameters.frequency - abs(parameters.exchange) - 2 * parameters.coupling * source_norm)

    def initial_vector(self, data=None, levels=(0, 0)):
        sites, cutoff = self.parameters.sites, self.parameters.cutoff
        data = initial_data(sites) if data is None else np.asarray(data, dtype=complex)
        if (data.shape != (sites,) or not np.isclose(np.linalg.norm(data), 1, rtol=0, atol=1e-12)
                or len(levels) != 2 or any(not isinstance(level, int) or not 0 <= level < cutoff for level in levels)):
            raise ValueError("Use normalized single-excitation data and two valid background levels.")
        state = np.zeros((sites, cutoff, cutoff, 4), dtype=complex)
        state[:, levels[0], levels[1], :] = data[:, None] / 2
        return state

    def evolve(self, state, duration):
        expected = (self.parameters.sites, self.parameters.cutoff, self.parameters.cutoff, 4)
        if state.shape != expected or not np.isfinite(duration):
            raise ValueError("Use the full data-two-mode-two-clock vector and finite time.")
        if duration == 0:
            return state.copy()
        output = np.empty_like(state, dtype=complex)
        for index, block in enumerate(self.blocks):
            coefficients = block["vectors"].T @ state[..., index].reshape(-1)
            output[..., index] = (block["vectors"] @ (np.exp(-1j * duration * block["values"]) * coefficients)).reshape(state.shape[:-1])
        return output

    def energies(self, state):
        return {term: float(sum(np.vdot(state[..., index].reshape(-1), block[term] @ state[..., index].reshape(-1)).real
                                for index, block in enumerate(self.blocks)))
                for term in ("matter", "background_free", "background_link", "interaction", "hamiltonian")}

    def reduced_states(self, state):
        sites, cutoff = self.parameters.sites, self.parameters.cutoff
        data = state.reshape(sites, -1)
        background = state.transpose(1, 2, 0, 3).reshape(cutoff ** 2, -1)
        first_mode = state.transpose(1, 0, 2, 3).reshape(cutoff, -1)
        second_mode = state.transpose(2, 0, 1, 3).reshape(cutoff, -1)
        clocks = state.transpose(3, 0, 1, 2).reshape(4, -1)
        return {name: vector @ vector.conj().T for name, vector in
                (("data", data), ("background", background), ("mode_A", first_mode), ("mode_C", second_mode), ("clocks", clocks))}

    def remote_clock_probability(self, state, basis="Y"):
        return checked_probability(self.reduced_states(state)["clocks"], np.kron(IDENTITY, clock_effect(basis, 0)))

    def diagnostics(self, state):
        reduced = self.reduced_states(state)
        cutoff = self.parameters.cutoff
        nonvacuum = np.eye(cutoff)
        nonvacuum[0, 0] = 0
        lowering = self.oscillator["lowering"]
        coherence = np.trace(reduced["background"] @ np.kron(lowering.T, lowering))
        return {"norm": float(np.linalg.norm(state)), "energies": self.energies(state),
                "mode_numbers": [float(np.trace(reduced["background"] @ number).real) for number in self.mode_numbers],
                "mode_nonvacuum_probabilities": [checked_probability(reduced[name], nonvacuum) for name in ("mode_A", "mode_C")],
                "number_current_A_to_C": float(-2 * self.parameters.exchange * coherence.imag),
                "clock_A_X_plus_probability": checked_probability(reduced["clocks"], np.kron(clock_effect("X", 0), IDENTITY)),
                "clock_C_Y_plus_probability": self.remote_clock_probability(state),
                "data_purity": float(np.trace(reduced["data"] @ reduced["data"]).real)}

    def joint_record_probabilities(self, state):
        thresholds, denominator = postprocessing_table(self.parameters.sites, center_fraction=0.5)
        mode_weights = np.zeros((2, self.parameters.cutoff))
        mode_weights[0, 0], mode_weights[1, 1:] = 1., 1.
        probabilities = np.zeros((2, 2, 2, 2, 2))
        for data_bit, first_mode, second_mode, first_clock, second_clock in np.ndindex(probabilities.shape):
            data_weight = thresholds / denominator if data_bit == 0 else 1 - thresholds / denominator
            effect = np.kron(clock_effect("X", first_clock), clock_effect("Y", second_clock))
            probabilities[data_bit, first_mode, second_mode, first_clock, second_clock] = np.einsum(
                "jkla,ab,jklb,j,k,l->", state.conj(), effect, state, data_weight,
                mode_weights[first_mode], mode_weights[second_mode]).real
        if probabilities.min() < -1e-12 or abs(probabilities.sum() - 1) > 1e-12:
            raise ValueError("Terminal joint law must be normalized and positive.")
        return probabilities

    def cutoff_residual_matrix(self, block):
        cutoff, sites = self.parameters.cutoff, self.parameters.sites
        extended = oscillator_operators(cutoff + 2)
        inclusion = np.eye(cutoff + 2)[:, :cutoff]
        lowering = extended["lowering"][:, :cutoff]
        raising = extended["lowering"].T[:, :cutoff]
        field = extended["centered_square"][:, :cutoff]
        local = self.parameters.coupling * (
            np.kron(block["sources"][0], np.kron(field, inclusion))
            + np.kron(block["sources"][1], np.kron(inclusion, field)))
        exchange = self.parameters.exchange * np.kron(np.eye(sites),
                    np.kron(raising, lowering) + np.kron(lowering, raising))
        inside = np.zeros((cutoff + 2, cutoff + 2), dtype=bool)
        inside[:cutoff, :cutoff] = True
        return (local + exchange)[np.tile(~inside.reshape(-1), sites)]

    def cutoff_error_budget(self, initial, duration, samples=256):
        if (not np.isfinite(duration) or duration < 0 or not isinstance(samples, int) or samples < 1
                or not np.isclose(np.linalg.norm(initial), 1, rtol=0, atol=1e-12)):
            raise ValueError("Use a normalized initial vector, finite nonnegative horizon and positive quadrature count.")
        if self.stability_margin <= 0:
            raise ValueError("The infinite-oscillator comparison requires the declared positive stability margin.")
        times = duration * (np.arange(samples) + 0.5) / samples
        squared_residual = np.zeros(samples)
        squared_lipschitz = 0.
        for index, block in enumerate(self.blocks):
            residual = self.cutoff_residual_matrix(block)
            coefficients = block["vectors"].T @ initial[..., index].reshape(-1)
            exterior_modes = residual @ block["vectors"]
            weight = float(np.vdot(coefficients, coefficients).real)
            shift = float(np.dot(abs(coefficients) ** 2, block["values"]) / weight) if weight else 0.
            frequency = block["values"] - shift
            amplitudes = exterior_modes @ (coefficients[:, None] * np.exp(-1j * np.outer(frequency, times)))
            squared_residual += np.sum(abs(amplitudes) ** 2, axis=0)
            lipschitz = float(np.sum(abs(coefficients) * abs(frequency) * np.linalg.norm(exterior_modes, axis=0)))
            squared_lipschitz += lipschitz ** 2
        integral = duration * float(np.mean(np.sqrt(squared_residual)))
        remainder = math.sqrt(squared_lipschitz) * duration ** 2 / (4 * samples)
        return {"sampled_boundary_residual_integral": integral, "lipschitz_midpoint_remainder_bound": remainder,
                "norm_error_upper_bound": integral + remainder,
                "same_final_record_TV_upper_bound": min(1., integral + remainder),
                "quadrature_samples": samples, "oscillator_stability_margin": self.stability_margin,
                "scope": "two-mode rectangular Fock projection, exact spectral formula evaluated in floating point, not directed rounding"}


@lru_cache(maxsize=4)
def model(parameters=Parameters()):
    return LocalQuantumBackground(parameters)


def source_preparation(sites, occupied):
    if occupied not in (0, 2):
        raise ValueError("The local source task compares occupied site 0 with the uncoupled storage site 2.")
    return np.eye(sites, dtype=complex)[:, occupied]


def embed_cutoff_state(state, cutoff):
    previous = state.shape[1]
    if cutoff < previous:
        raise ValueError("An embedding cannot discard occupied cutoff coordinates.")
    output = np.zeros((state.shape[0], cutoff, cutoff, 4), dtype=complex)
    output[:, :previous, :previous, :] = state
    return output


def terminal_schedule(sites, duration):
    _, spacing = grid(sites, BACKGROUND)
    if not np.isfinite(duration) or duration < 0:
        raise ValueError("Use a finite nonnegative read parameter.")
    sources = [("position", site) for site in range(sites)]
    sources += [("clock", 0), ("clock", 1), ("background_vacuum_bit", 0), ("background_vacuum_bit", 1)]
    receipt = duration
    schedule = []
    for kind, source in sources:
        hops = min(source, sites - source)
        receipt += hops * spacing
        schedule.append({"kind": kind, "source_site": source, "read_at": duration,
                         "hops": hops, "received_at_site_zero": receipt})
    return schedule


def release_record(raw_positions, clock_bits, background_bits, random_integer, receipt, duration):
    sites = len(raw_positions)
    schedule = terminal_schedule(sites, duration)
    if (sum(raw_positions) != 1 or any(bit not in (0, 1) for bit in raw_positions)
            or len(clock_bits) != 2 or len(background_bits) != 2
            or any(bit not in (0, 1) for bit in (*clock_bits, *background_bits))
            or not isinstance(random_integer, int) or not 0 <= random_integer < 4096):
        raise ValueError("Retain all position, two-clock and two-background bits and twelve independent fair bits.")
    if not np.isfinite(receipt) or receipt < schedule[-1]["received_at_site_zero"]:
        raise ValueError("Do not release a decision before every scheduled raw record arrives.")
    thresholds, _ = postprocessing_table(sites, center_fraction=0.5)
    data_bit = int(random_integer >= thresholds[list(raw_positions).index(1)])
    return data_bit, *background_bits, *clock_bits


def resource_report(parameters, duration):
    schedule = terminal_schedule(parameters.sites, duration)
    _, denominator = postprocessing_table(parameters.sites, center_fraction=0.5)
    return {"chain_qubits": parameters.sites, "clock_qubits": 2, "background_modes": 2,
            "levels_per_background_mode": parameters.cutoff, "quantum_subsystems_before_detectors": parameters.sites + 4,
            "occupied_complex_sector_dimension": 4 * parameters.sites * parameters.cutoff ** 2,
            "classical_raw_outcome_bits": parameters.sites + 4, "coarse_joint_outcomes": 32,
            "classical_bit_hops": sum(event["hops"] for event in schedule),
            "final_collection_parameter": schedule[-1]["received_at_site_zero"],
            "independent_fair_output_bits": 12, "fixed_threshold_table_bits": parameters.sites * int(denominator).bit_length(),
            "resets": 0, "preparation_synchrony_detector_and_switch_energy_costs_closed": False}


@lru_cache(maxsize=2)
def source_response_report(cutoff=10, duration=2.):
    rows = []
    for exchange in (0., 0.3):
        current = model(Parameters(cutoff=cutoff, exchange=exchange, hopping_scale=0.))
        outputs, budgets = [], []
        for occupied in (0, 2):
            initial = current.initial_vector(source_preparation(8, occupied))
            outputs.append(current.evolve(initial, duration))
            budgets.append(current.cutoff_error_budget(initial, duration))
        probabilities = [current.remote_clock_probability(state) for state in outputs]
        difference = abs(probabilities[0] - probabilities[1])
        allowance = sum(budget["same_final_record_TV_upper_bound"] for budget in budgets)
        rows.append({"exchange": exchange, "hopping_scale": 0., "source_occupied_sites": [0, 2],
                     "clock_C_Y_plus_probabilities": probabilities, "absolute_remote_record_difference": difference,
                     "cutoff_budgets": budgets, "difference_lower_bound_with_cutoff_allowance": max(0., difference - allowance),
                     "mode_C_trace_distance": trace_distance(*(current.reduced_states(state)["mode_C"] for state in outputs))})
    return {"cutoff": cutoff, "duration": duration, "rows": rows,
            "control_scope": "data hopping and its bond modulation are held at zero for the entire interval; this is a declared diagnostic, not the default chain"}


def run_report():
    current = model()
    initial = current.initial_vector()
    duration = 2.
    final = current.evolve(initial, duration)
    budget = current.cutoff_error_budget(initial, duration)
    larger = model(Parameters(cutoff=10))
    larger_initial = larger.initial_vector()
    reference = larger.evolve(larger_initial, duration)
    changed_data = current.evolve(current.initial_vector(initial_data(8, -1)), duration)
    changed_background = current.evolve(current.initial_vector(levels=(1, 0)), duration)
    final_reduced = current.reduced_states(final)
    data_changed_reduced = current.reduced_states(changed_data)
    background_changed_reduced = current.reduced_states(changed_background)
    comparison = []
    for cutoff in (4, 6, 8):
        candidate = model(Parameters(cutoff=cutoff))
        prepared = candidate.initial_vector()
        evolved = candidate.evolve(prepared, duration)
        comparison.append({"cutoff": cutoff, "vector_difference_from_D10": float(np.linalg.norm(embed_cutoff_state(evolved, 10) - reference)),
                           "budget": candidate.cutoff_error_budget(prepared, duration)})
    return {"model_version": "C7", "route_round": "07", "research_scope": "sufficiency only",
            "parameters": vars(current.parameters), "background_sites_and_clock_sites": [0, 1],
            "Hamiltonian": "H_C5(hopping_scale)+H_clocks+Omega*(n_A+n_C)+exchange*(a_A^dag a_C+a_A a_C^dag)+epsilon*(S_A F_A+S_C F_C)",
            "coefficient_lower_bound": 1 - current.parameters.coupling / 2,
            "time_rows": [dict(duration=time, **current.diagnostics(current.evolve(initial, time))) for time in (0., 0.5, 1., 2.)],
            "default_cutoff_record_budget": budget,
            "cutoff_comparisons": comparison, "D10_reference_own_budget": larger.cutoff_error_budget(larger_initial, duration),
            "joint_terminal_probabilities_data_modeA_modeC_clockA_clockC": current.joint_record_probabilities(final).tolist(),
            "default_reciprocal_diagnostics": {
                "background_distance_for_changed_carrier_preparation": trace_distance(final_reduced["background"], data_changed_reduced["background"]),
                "data_distance_for_changed_mode_A_preparation": trace_distance(final_reduced["data"], background_changed_reduced["data"]),
                "clock_distance_for_changed_mode_A_preparation": trace_distance(final_reduced["clocks"], background_changed_reduced["clocks"]),
                "scope": "finite D8 state diagnostics, not individually cutoff-certified trace distances"},
            "coarser_source_response": source_response_report(8),
            "isolated_source_response": source_response_report(), "resources": resource_report(current.parameters, duration),
            "source_task_resources": resource_report(Parameters(cutoff=10, hopping_scale=0.), duration),
            "terminal_schedule": terminal_schedule(current.parameters.sites, duration),
            "limitations": ["two-site background patch, not a spatial continuum or Einstein dynamics",
                            "clock C relocated from the opposite site to adjacent site 1 as an explicit C7 input",
                            "default chain still propagates; the held-hopping task separately isolates a background path",
                            "quantum rules, oscillator law, local coupling and clock synchrony remain inputs",
                            "C5 static continuum bound not extended; no uniform stability claim under refinement",
                            "no sample counts or laboratory evidence; mode budget is not floating-point rounding certification"]}


class LocalQuantumBackgroundTests(unittest.TestCase):
    def test_data_sources_have_only_declared_local_support(self):
        parameters = Parameters(cutoff=4)
        _, sources = local_data_terms(parameters)
        for center, source in enumerate(sources):
            permitted = {(center, center)}
            for neighbor in ((center - 1) % parameters.sites, (center + 1) % parameters.sites):
                permitted.update(((center, neighbor), (neighbor, center)))
            self.assertEqual(set(zip(*np.nonzero(source))), permitted)
        current = model(parameters)
        self.assertGreater(current.stability_margin, 0)

    def test_free_background_link_transports_one_excitation_exactly(self):
        current = model(Parameters(cutoff=4, coupling=0.))
        initial = current.initial_vector(levels=(1, 0))
        duration = 1.7
        state = current.evolve(initial, duration)
        background = current.reduced_states(state)["background"]
        reference = np.zeros(16, dtype=complex)
        reference[4] = np.cos(current.parameters.exchange * duration)
        reference[1] = -1j * np.sin(current.parameters.exchange * duration)
        np.testing.assert_allclose(background, np.outer(reference, reference.conj()), atol=3e-14)

    def test_energy_and_inverse_preserve_the_complete_joint_process(self):
        current = model(Parameters(cutoff=4))
        initial = current.initial_vector()
        state = current.evolve(initial, 1.3)
        np.testing.assert_allclose(current.evolve(state, -1.3), initial, atol=3e-14)
        original, final = current.energies(initial), current.energies(state)
        self.assertAlmostEqual(original["hamiltonian"], final["hamiltonian"])
        self.assertAlmostEqual(final["hamiltonian"], sum(final[name] for name in ("matter", "background_free", "background_link", "interaction")))

    def test_held_hopping_and_cut_background_link_remove_remote_source_response(self):
        current = model(Parameters(cutoff=4, hopping_scale=0., exchange=0.))
        outputs = [current.evolve(current.initial_vector(source_preparation(8, occupied)), 2.) for occupied in (0, 2)]
        self.assertAlmostEqual(current.remote_clock_probability(outputs[0]), current.remote_clock_probability(outputs[1]))
        np.testing.assert_allclose(current.reduced_states(outputs[0])["mode_C"], current.reduced_states(outputs[1])["mode_C"], atol=3e-14)

    def test_adjacent_background_link_restores_source_to_remote_clock_response(self):
        current = model(Parameters(cutoff=4, hopping_scale=0.))
        outputs = [current.evolve(current.initial_vector(source_preparation(8, occupied)), 2.) for occupied in (0, 2)]
        difference = abs(current.remote_clock_probability(outputs[0]) - current.remote_clock_probability(outputs[1]))
        self.assertGreater(difference, 1e-6)

    def test_rectangular_residual_includes_local_pair_and_exchange_leakage(self):
        current = model(Parameters(cutoff=4))
        larger = model(Parameters(cutoff=6))
        generator = np.random.default_rng(7)
        prepared = generator.normal(size=(8, 4, 4, 4)) + 1j * generator.normal(size=(8, 4, 4, 4))
        prepared /= np.linalg.norm(prepared)
        embedded = embed_cutoff_state(prepared, 6)
        inside = np.zeros((8, 6, 6), dtype=bool)
        inside[:, :4, :4] = True
        for index, (block, extended) in enumerate(zip(current.blocks, larger.blocks)):
            actual = extended["hamiltonian"] @ embedded[..., index].reshape(-1)
            evolved = np.zeros((8, 6, 6), dtype=complex)
            evolved[:, :4, :4] = (block["hamiltonian"] @ prepared[..., index].reshape(-1)).reshape(8, 4, 4)
            residual = actual - evolved.reshape(-1)
            np.testing.assert_allclose(residual[inside.reshape(-1)], 0, atol=2e-14)
            np.testing.assert_allclose(residual[~inside.reshape(-1)], current.cutoff_residual_matrix(block) @ prepared[..., index].reshape(-1), atol=2e-14)

    def test_cutoff_difference_obeys_independent_residual_budgets(self):
        current, larger = model(Parameters(cutoff=4)), model(Parameters(cutoff=6))
        initial, reference_initial = current.initial_vector(), larger.initial_vector()
        duration = 1.
        difference = np.linalg.norm(embed_cutoff_state(current.evolve(initial, duration), 6) - larger.evolve(reference_initial, duration))
        allowance = current.cutoff_error_budget(initial, duration)["norm_error_upper_bound"] + larger.cutoff_error_budget(reference_initial, duration)["norm_error_upper_bound"]
        self.assertLess(difference, allowance)

    def test_joint_record_marginal_matches_remote_clock_without_factorization(self):
        current = model(Parameters(cutoff=4))
        state = current.evolve(current.initial_vector(), 2.)
        records = current.joint_record_probabilities(state)
        self.assertAlmostEqual(records.sum(), 1.)
        self.assertGreaterEqual(records.min(), -1e-14)
        self.assertAlmostEqual(records[..., 0].sum(), current.remote_clock_probability(state))
        background_joint = records.sum(axis=(0, 3, 4))
        self.assertGreater(np.linalg.norm(background_joint - np.outer(background_joint.sum(axis=1), background_joint.sum(axis=0))), 1e-7)

    def test_equivalent_ensembles_preserve_the_full_joint_state(self):
        current = model(Parameters(cutoff=4))
        first, second = (source_preparation(8, occupied) for occupied in (0, 2))
        direct = [current.evolve(current.initial_vector(data), 1.).reshape(-1) for data in (first, second)]
        mixed = [current.evolve(current.initial_vector((first + sign * second) / math.sqrt(2)), 1.).reshape(-1) for sign in (1, -1)]
        np.testing.assert_allclose(sum(np.outer(state, state.conj()) / 2 for state in direct),
                                   sum(np.outer(state, state.conj()) / 2 for state in mixed), atol=2e-14)

    def test_background_number_continuity_has_exchange_and_local_source_terms(self):
        current = model(Parameters(cutoff=4))
        state = current.evolve(current.initial_vector(), 0.8)
        current_value = current.diagnostics(state)["number_current_A_to_C"]
        for mode_index, orientation in ((0, -1), (1, 1)):
            observable = np.kron(np.eye(8), current.mode_numbers[mode_index])
            total_derivative, source_derivative = 0., 0.
            for index, block in enumerate(current.blocks):
                vector = state[..., index].reshape(-1)
                total_derivative += (1j * np.vdot(vector, (block["hamiltonian"] @ observable - observable @ block["hamiltonian"]) @ vector)).real
                source_derivative += (1j * np.vdot(vector, (block["interaction"] @ observable - observable @ block["interaction"]) @ vector)).real
            self.assertAlmostEqual(total_derivative, orientation * current_value + source_derivative)

    def test_complete_record_return_and_fair_bits_have_finite_cost(self):
        positions = [0] * 8
        positions[3] = 1
        schedule = terminal_schedule(8, 2.)
        self.assertEqual(len(schedule), 12)
        self.assertEqual(sum(event["hops"] for event in schedule), 18)
        self.assertEqual(schedule[-1]["received_at_site_zero"], 47.)
        with self.assertRaises(ValueError):
            release_record(positions, [0, 1], [1, 0], 0, 46.9, 2.)
        count = sum(release_record(positions, [0, 1], [1, 0], random, 47., 2.)[0] == 0 for random in range(4096))
        thresholds, _ = postprocessing_table(8, center_fraction=0.5)
        self.assertEqual(count, thresholds[3])
        self.assertEqual(release_record(positions, [0, 1], [1, 0], 0, 47., 2.)[1:], (1, 0, 0, 1))

    def test_invalid_inputs_and_unstable_completion_budget_are_rejected(self):
        for arguments in ({"cutoff": 3}, {"exchange": 1.}, {"hopping_scale": -1.}, {"coupling": 2.}):
            with self.assertRaises(ValueError):
                Parameters(**arguments)
        unstable = model(Parameters(cutoff=4, coupling=1.))
        with self.assertRaises(ValueError):
            unstable.cutoff_error_budget(unstable.initial_vector(), 1.)

    def test_zero_matter_coupling_retains_C5_propagation_and_relocated_free_clocks(self):
        current = model(Parameters(cutoff=4, coupling=0.))
        duration = 0.9
        state = current.evolve(current.initial_vector(), duration)
        data = evolve_lattice(initial_data(8), duration)
        for index, block in enumerate(current.blocks):
            expected = np.zeros((8, 4, 4), dtype=complex)
            expected[:, 0, 0] = data * np.exp(-1j * duration * block["clock_energy"]) / 2
            np.testing.assert_allclose(state[..., index], expected, atol=3e-14)

    def test_default_active_chain_retains_both_directions_of_backreaction(self):
        current = model(Parameters(cutoff=4))
        original = current.reduced_states(current.evolve(current.initial_vector(), 2.))
        changed_data = current.reduced_states(current.evolve(current.initial_vector(initial_data(8, -1)), 2.))
        changed_mode = current.reduced_states(current.evolve(current.initial_vector(levels=(1, 0)), 2.))
        self.assertGreater(trace_distance(original["background"], changed_data["background"]), 1e-5)
        self.assertGreater(trace_distance(original["data"], changed_mode["data"]), 1e-5)

    def test_coherent_clock_blocks_have_a_single_common_J_real_representation(self):
        current = model(Parameters(cutoff=4))
        initial = current.initial_vector().reshape(-1)
        duration = 0.6
        complex_gate = np.zeros((len(initial), len(initial)), dtype=complex)
        for index, block in enumerate(current.blocks):
            indices = np.arange(index, len(initial), 4)
            gate = (block["vectors"] * np.exp(-1j * duration * block["values"])) @ block["vectors"].T
            complex_gate[np.ix_(indices, indices)] = gate
        expected = current.evolve(current.initial_vector(), duration).reshape(-1)
        np.testing.assert_allclose(complex_gate @ initial, expected, atol=3e-14)
        real_gate = real_lift(complex_gate)
        actual = real_gate @ encode_state(np.outer(initial, initial.conj())) @ real_gate.T
        np.testing.assert_allclose(actual, encode_state(np.outer(expected, expected.conj())), atol=3e-14)

    def test_refined_source_record_exceeds_its_two_preparation_cutoff_budget(self):
        report = source_response_report()
        disconnected, connected = report["rows"]
        self.assertLess(disconnected["absolute_remote_record_difference"], 1e-12)
        self.assertGreater(connected["difference_lower_bound_with_cutoff_allowance"], 1e-5)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(LocalQuantumBackgroundTests)
    outcome = unittest.TextTestRunner(verbosity=2).run(suite)
    if not outcome.wasSuccessful():
        raise SystemExit(1)
    if arguments.write_results:
        report = run_report()
        Path(__file__).with_name("local_quantum_background_results.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(json.dumps({"tests": outcome.testsRun, "default_budget": report["default_cutoff_record_budget"],
                          "source_task": report["isolated_source_response"], "resources": report["resources"]}, indent=2))