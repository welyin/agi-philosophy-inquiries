"""Construction 06 / C6: one quantum background mode coupled to the C5 chain.

A finite Fock projection of q^2-1/2 modulates the same local bonds, staggered
mass and two clock rates. Clock-Z sectors allow exact block evolution while
retaining their coherences. The background is a collective mode, not a local
gravitational field. All couplings and its oscillator law are declared inputs.
"""

import argparse
import json
import math
import unittest
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

from baseline import IDENTITY, PAULI_X, PAULI_Z, real_lift, checked_probability
from common_orientation_structure import encode_state
from clock_calibration import clock_effect
from envelope_limit import grid, hopping_matrix, collection_schedule, final_record_resources
from massive_envelope import (BACKGROUND, DEFAULT_MASS, initial_modes, reconstruct,
                              massive_hopping, evolve_lattice, postprocessing_table)


@dataclass(frozen=True)
class Parameters:
    sites: int = 16
    cutoff: int = 10
    frequency: float = 1.0
    coupling: float = 0.2
    mass: float = DEFAULT_MASS

    def __post_init__(self):
        grid(self.sites, BACKGROUND)
        if (not isinstance(self.cutoff, int) or self.cutoff < 4
                or not all(np.isfinite(value) for value in (self.frequency, self.coupling, self.mass))
                or self.frequency <= 0 or not 0 <= self.coupling < 2 or self.mass < 0):
            raise ValueError("Require cutoff >=4, positive frequency, 0<=coupling<2 and nonnegative mass.")


def oscillator_operators(cutoff):
    lowering = np.diag(np.sqrt(np.arange(1, cutoff)), 1)
    number = np.diag(np.arange(cutoff, dtype=float))
    position = (lowering + lowering.T) / math.sqrt(2)
    momentum = -1j * (lowering - lowering.T) / math.sqrt(2)
    quadratic = number + np.eye(cutoff) / 2 + (lowering @ lowering + lowering.T @ lowering.T) / 2
    return {"lowering": lowering, "number": number, "position": position,
            "momentum": momentum, "position_squared_projected": quadratic,
            "centered_square": quadratic - np.eye(cutoff) / 2}


def shape(positions):
    return (1 + np.cos(BACKGROUND.frequency * (np.asarray(positions) - 0.35 * BACKGROUND.length))) / 2


def data_source(sites, mass):
    positions, spacing = grid(sites, BACKGROUND)
    source = np.zeros((sites, sites))
    for site in range(sites):
        adjacent = (site + 1) % sites
        midpoint = positions[site] + spacing / 2
        source[site, adjacent] = source[adjacent, site] = BACKGROUND.speed(midpoint) * shape(midpoint) / (2 * spacing)
    source += np.diag((-1.) ** np.arange(sites) * mass * BACKGROUND.lapse(positions) * shape(positions))
    return source


def initial_data(sites, sign=1):
    if sign not in (-1, 1):
        raise ValueError("Use a plus or minus known carrier superposition.")
    modes = (initial_modes(direction=1) + sign * initial_modes(direction=-1)) / math.sqrt(2)
    vector = reconstruct(modes, sites)
    return vector / np.linalg.norm(vector)


class QuantumChainBackground:
    def __init__(self, parameters):
        self.parameters = parameters
        sites, cutoff = parameters.sites, parameters.cutoff
        self.oscillator = oscillator_operators(cutoff)
        self.data_hamiltonian = massive_hopping(sites, parameters.mass)
        self.source = data_source(sites, parameters.mass)
        self.clock_rates = np.array([BACKGROUND.lapse(0.), 1.2 * BACKGROUND.lapse(BACKGROUND.length / 2)])
        self.clock_weights = self.clock_rates * np.array([shape(0.), shape(BACKGROUND.length / 2)])
        self.blocks = []
        for first_bit, second_bit in ((0, 0), (0, 1), (1, 0), (1, 1)):
            signs = np.array([1 - 2 * first_bit, 1 - 2 * second_bit])
            clock_energy = float(np.dot(signs, self.clock_rates) / 2)
            clock_source = float(np.dot(signs, self.clock_weights) / 2)
            source = self.source + clock_source * np.eye(sites)
            matter = np.kron(self.data_hamiltonian + clock_energy * np.eye(sites), np.eye(cutoff))
            background = np.kron(np.eye(sites), parameters.frequency * self.oscillator["number"])
            interaction = parameters.coupling * np.kron(source, self.oscillator["centered_square"])
            total = matter + background + interaction
            values, vectors = np.linalg.eigh(total)
            self.blocks.append({"source": source, "clock_energy": clock_energy,
                                "matter": matter, "background": background, "interaction": interaction,
                                "hamiltonian": total, "values": values, "vectors": vectors})
        self.source_norm = max(float(np.linalg.norm(block["source"], ord=2)) for block in self.blocks)

    def initial_vector(self, data=None, level=0):
        if not isinstance(level, int) or not 0 <= level < self.parameters.cutoff:
            raise ValueError("Initial background level is outside the finite cutoff.")
        data = initial_data(self.parameters.sites) if data is None else np.asarray(data, dtype=complex)
        if data.shape != (self.parameters.sites,) or not np.isclose(np.linalg.norm(data), 1, rtol=0, atol=1e-12):
            raise ValueError("Use a normalized vector in the declared one-excitation sector.")
        state = np.zeros((self.parameters.sites, self.parameters.cutoff, 4), dtype=complex)
        state[:, level, :] = data[:, None] / 2
        return state

    def evolve(self, initial, duration):
        if initial.shape != (self.parameters.sites, self.parameters.cutoff, 4) or not np.isfinite(duration):
            raise ValueError("Use the full data-background-two-clock vector and finite duration.")
        if duration == 0:
            return initial.copy()
        output = np.empty_like(initial, dtype=complex)
        for index, block in enumerate(self.blocks):
            vector = initial[:, :, index].reshape(-1)
            transformed = block["vectors"].T @ vector
            output[:, :, index] = (block["vectors"] @ (np.exp(-1j * duration * block["values"]) * transformed)).reshape(initial.shape[:2])
        return output

    def apply_term(self, state, term):
        output = np.empty_like(state)
        for index, block in enumerate(self.blocks):
            output[:, :, index] = (block[term] @ state[:, :, index].reshape(-1)).reshape(state.shape[:2])
        return output

    def energies(self, state):
        return {term: float(np.vdot(state, self.apply_term(state, term)).real)
                for term in ("matter", "background", "interaction", "hamiltonian")}

    def reduced_states(self, state):
        data = state.reshape(self.parameters.sites, -1)
        background = state.transpose(1, 0, 2).reshape(self.parameters.cutoff, -1)
        clocks = state.transpose(2, 0, 1).reshape(4, -1)
        return data @ data.conj().T, background @ background.conj().T, clocks @ clocks.conj().T

    def diagnostics(self, state):
        data, background, clocks = self.reduced_states(state)
        eigenvalues = np.linalg.eigvalsh(background)
        if eigenvalues.min() < -1e-12 or abs(np.trace(background) - 1) > 1e-12:
            raise ValueError("Reduced background is not a normalized positive state.")
        positive = eigenvalues[eigenvalues > 1e-14]
        entropy = max(0., float(-np.sum(positive * np.log2(positive))))
        vacuum = np.zeros_like(background)
        vacuum[0, 0] = 1
        smooth, denominator = postprocessing_table(self.parameters.sites, center_fraction=0.5)
        return {"norm": float(np.linalg.norm(state)), "energies": self.energies(state),
                "background_number": float(np.trace(background @ self.oscillator["number"]).real),
                "background_position_mean": float(np.trace(background @ self.oscillator["position"]).real),
                "background_projected_position_square": float(np.trace(background @ self.oscillator["position_squared_projected"]).real),
                "background_nonvacuum_probability": checked_probability(background, np.eye(self.parameters.cutoff) - vacuum),
                "background_entropy_bits": entropy, "data_purity": float(np.trace(data @ data).real),
                "highest_two_background_levels_probability": float(np.diag(background).real[-2:].sum()),
                "smooth_data_probability": float(np.dot(smooth / denominator, np.diag(data).real)),
                "clock_A_X_plus_probability": checked_probability(clocks, np.kron(clock_effect("X", 0), IDENTITY)),
                "clock_C_Y_plus_probability": checked_probability(clocks, np.kron(IDENTITY, clock_effect("Y", 0)))}

    def joint_record_probabilities(self, state):
        smooth, denominator = postprocessing_table(self.parameters.sites, center_fraction=0.5)
        values = np.zeros((2, 2, 2, 2))
        for data_bit in (0, 1):
            data_weight = smooth / denominator if data_bit == 0 else 1 - smooth / denominator
            for background_bit in (0, 1):
                field_weight = np.zeros(self.parameters.cutoff)
                if background_bit == 0:
                    field_weight[0] = 1
                else:
                    field_weight[1:] = 1
                for clock_A in (0, 1):
                    for clock_C in (0, 1):
                        effect = np.kron(clock_effect("X", clock_A), clock_effect("Y", clock_C))
                        values[data_bit, background_bit, clock_A, clock_C] = np.einsum(
                            "jka,ab,jkb,j,k->", state.conj(), effect, state, data_weight, field_weight).real
        if values.min() < -1e-12 or abs(values.sum() - np.vdot(state, state).real) > 1e-12:
            raise ValueError("Joint terminal record is not a normalized positive law.")
        return values

    def cutoff_error_budget(self, initial, duration, samples=256):
        if duration < 0 or not np.isfinite(duration) or not isinstance(samples, int) or samples < 1:
            raise ValueError("Use a finite nonnegative horizon and positive quadrature count.")
        stability_margin = self.parameters.frequency - 2 * self.parameters.coupling * self.source_norm
        if stability_margin <= 0:
            raise ValueError("The stated oscillator-completion budget requires the positive stability margin in C6.")
        cutoff, sites = self.parameters.cutoff, self.parameters.sites
        times = duration * (np.arange(samples) + 0.5) / samples
        residual_squared = np.zeros(samples)
        lipschitz_squared = 0.0
        for index, block in enumerate(self.blocks):
            outside = np.zeros((2 * sites, cutoff * sites))
            for target in range(sites):
                for source in range(sites):
                    entry = self.parameters.coupling * block["source"][target, source] / 2
                    outside[target, source * cutoff + cutoff - 2] = entry * math.sqrt(cutoff * (cutoff - 1))
                    outside[sites + target, source * cutoff + cutoff - 1] = entry * math.sqrt(cutoff * (cutoff + 1))
            coefficients = block["vectors"].T @ initial[:, :, index].reshape(-1)
            outside_modes = outside @ block["vectors"]
            state_norm = float(np.vdot(coefficients, coefficients).real)
            shift = float(np.dot(abs(coefficients) ** 2, block["values"]) / state_norm) if state_norm else 0.
            phases = np.exp(-1j * np.outer(block["values"] - shift, times)) * coefficients[:, None]
            residual = outside_modes @ phases
            residual_squared += np.sum(abs(residual) ** 2, axis=0)
            lipschitz = np.sum(abs(coefficients) * abs(block["values"] - shift) * np.linalg.norm(outside_modes, axis=0))
            lipschitz_squared += float(lipschitz ** 2)
        integral = duration * float(np.mean(np.sqrt(residual_squared)))
        quadrature = math.sqrt(lipschitz_squared) * duration ** 2 / (4 * samples)
        return {"sampled_boundary_residual_integral": integral,
                "lipschitz_midpoint_remainder_bound": quadrature,
                "norm_error_upper_bound": integral + quadrature,
                "same_final_record_TV_upper_bound": min(1., integral + quadrature),
                "quadrature_samples": samples,
                "sufficient_oscillator_stability_margin": stability_margin,
                "scope": "exact finite-cutoff spectral quantities; float evaluation is not directed-rounding certification"}


@lru_cache(maxsize=12)
def model(parameters=Parameters()):
    return QuantumChainBackground(parameters)


def trace_distance(first, second):
    return float(np.abs(np.linalg.eigvalsh(first - second)).sum() / 2)


def terminal_schedule(sites, duration):
    schedule = collection_schedule(sites, duration, BACKGROUND)
    schedule.append({"kind": "background_vacuum_bit", "source_site": 0,
                     "read_at": duration, "hops": 0, "received_at_site_zero": duration})
    return sorted(schedule, key=lambda event: event["received_at_site_zero"])


def release_smooth_record(raw_positions, clock_bits, background_bit, random_integer, receipt_parameter, schedule):
    if (sum(raw_positions) != 1 or any(bit not in (0, 1) for bit in raw_positions)
            or len(clock_bits) != 2 or any(bit not in (0, 1) for bit in clock_bits)
            or background_bit not in (0, 1) or not isinstance(random_integer, int) or not 0 <= random_integer < 4096):
        raise ValueError("Use complete terminal records and twelve independent fair random bits.")
    if not np.isfinite(receipt_parameter) or receipt_parameter < max(event["received_at_site_zero"] for event in schedule):
        raise ValueError("Wait until the complete terminal record has arrived.")
    thresholds, denominator = postprocessing_table(len(raw_positions), center_fraction=0.5)
    site = raw_positions.index(1)
    data_bit = 0 if random_integer < thresholds[site] else 1
    return data_bit, background_bit, *clock_bits


class QuantumChainBackgroundTests(unittest.TestCase):
    def test_projected_square_and_coefficient_positivity_include_the_cutoff_correction(self):
        for cutoff in (4, 8, 12):
            operators = oscillator_operators(cutoff)
            correction = np.zeros((cutoff, cutoff))
            correction[-1, -1] = cutoff / 2
            np.testing.assert_allclose(operators["position_squared_projected"], operators["position"] @ operators["position"] + correction, atol=2e-14)
            self.assertGreater(np.linalg.eigvalsh(operators["position_squared_projected"]).min(), 0)
            for coupling in (0.2, 1.):
                coefficient = np.eye(cutoff) + coupling * operators["centered_square"]
                self.assertGreaterEqual(np.linalg.eigvalsh(coefficient).min(), 1 - coupling / 2)

    def test_zero_coupling_is_the_unchanged_C5_data_with_free_clocks_and_mode(self):
        current = model(Parameters(sites=8, cutoff=6, coupling=0.))
        initial = current.initial_vector()
        duration = 0.7
        evolved = current.evolve(initial, duration)
        data = evolve_lattice(initial_data(8), duration)
        for index, block in enumerate(current.blocks):
            expected = np.zeros((8, 6), dtype=complex)
            expected[:, 0] = data * np.exp(-1j * duration * block["clock_energy"]) / 2
            np.testing.assert_allclose(evolved[:, :, index], expected, rtol=0, atol=2e-14)

    def test_closed_joint_dynamics_is_unitary_and_energy_is_not_just_the_data_energy(self):
        current = model()
        initial = current.initial_vector()
        original_energy = current.energies(initial)["hamiltonian"]
        evolved = current.evolve(initial, 2.)
        restored = current.evolve(evolved, -2.)
        np.testing.assert_allclose(restored, initial, rtol=0, atol=2e-14)
        self.assertAlmostEqual(np.linalg.norm(evolved), 1)
        energies = current.energies(evolved)
        self.assertAlmostEqual(energies["hamiltonian"], original_energy)
        self.assertAlmostEqual(energies["hamiltonian"], energies["matter"] + energies["background"] + energies["interaction"])
        self.assertGreater(abs(energies["background"]), 1e-5)

    def test_centered_vacuum_has_no_initial_mean_shift_but_a_nonzero_second_order_response(self):
        current = model(Parameters(sites=8, cutoff=6))
        initial = current.initial_vector()
        number = np.kron(np.eye(8), current.oscillator["number"])
        second_derivative = 0.
        expected = 0.
        for index, block in enumerate(current.blocks):
            vector = initial[:, :, index].reshape(-1)
            hamiltonian = block["hamiltonian"]
            commutator = hamiltonian @ number - number @ hamiltonian
            second_derivative -= np.vdot(vector, (hamiltonian @ commutator - commutator @ hamiltonian) @ vector).real
            expected += 2 * current.parameters.coupling ** 2 * np.vdot(initial[:, 0, index], block["source"] @ block["source"] @ initial[:, 0, index]).real
        self.assertAlmostEqual(current.oscillator["centered_square"][0, 0], 0)
        self.assertAlmostEqual(second_derivative, expected)
        self.assertGreater(second_derivative, 0)

    def test_data_changes_the_background_and_background_changes_data_and_clocks(self):
        current = model()
        positive = current.evolve(current.initial_vector(initial_data(16, 1)), 2.)
        negative = current.evolve(current.initial_vector(initial_data(16, -1)), 2.)
        first_data, first_background, first_clocks = current.reduced_states(positive)
        _, second_background, _ = current.reduced_states(negative)
        self.assertGreater(trace_distance(first_background, second_background), 1e-4)
        excited = current.evolve(current.initial_vector(initial_data(16, 1), level=1), 2.)
        other_data, _, other_clocks = current.reduced_states(excited)
        self.assertGreater(trace_distance(first_data, other_data), 1e-4)
        self.assertGreater(trace_distance(first_clocks, other_clocks), 1e-4)
        self.assertGreater(current.diagnostics(positive)["background_entropy_bits"], 1e-5)

    def test_initial_background_width_curvature_is_driven_by_the_data_source(self):
        current = model(Parameters(sites=8, cutoff=6))
        observable = np.kron(np.eye(8), current.oscillator["position_squared_projected"])
        for sign in (1, -1):
            data = initial_data(8, sign)
            initial = current.initial_vector(data)
            actual = 0.
            for index, block in enumerate(current.blocks):
                hamiltonian = block["hamiltonian"]
                commutator = hamiltonian @ observable - observable @ hamiltonian
                vector = initial[:, :, index].reshape(-1)
                actual -= np.vdot(vector, (hamiltonian @ commutator - commutator @ hamiltonian) @ vector).real
            expected = -2 * current.parameters.frequency * current.parameters.coupling * np.vdot(data, current.source @ data).real
            self.assertAlmostEqual(actual, expected)

    def test_all_joint_records_are_positive_and_keep_correlations(self):
        current = model()
        state = current.evolve(current.initial_vector(), 2.)
        probabilities = current.joint_record_probabilities(state)
        self.assertAlmostEqual(probabilities.sum(), 1)
        self.assertTrue(np.all(probabilities >= -1e-14))
        data_field = probabilities.sum(axis=(2, 3))
        independent = np.outer(data_field.sum(axis=1), data_field.sum(axis=0))
        self.assertGreater(np.max(abs(data_field - independent)), 1e-7)

    def test_equal_initial_ensembles_give_equal_complete_final_density(self):
        current = model(Parameters(sites=8, cutoff=6))
        first, second = np.eye(8)[0], np.eye(8)[1]
        alternatives = ((first, second), ((first + second) / math.sqrt(2), (first - second) / math.sqrt(2)))
        outputs = []
        for ensemble in alternatives:
            final = [current.evolve(current.initial_vector(data), 0.7).reshape(-1) for data in ensemble]
            outputs.append(sum(np.outer(vector, vector.conj()) for vector in final) / 2)
        np.testing.assert_allclose(outputs[0], outputs[1], rtol=0, atol=2e-14)

    def test_cutoff_comparison_obeys_the_boundary_residual_budget(self):
        coarse = model(Parameters(sites=8, cutoff=6))
        fine = model(Parameters(sites=8, cutoff=10))
        coarse_initial, fine_initial = coarse.initial_vector(), fine.initial_vector()
        first, second = coarse.evolve(coarse_initial, 2.), fine.evolve(fine_initial, 2.)
        padded = np.zeros_like(second)
        padded[:, :6, :] = first
        coarse_bound = coarse.cutoff_error_budget(coarse_initial, 2.)
        fine_bound = fine.cutoff_error_budget(fine_initial, 2.)
        self.assertGreater(coarse_bound["sufficient_oscillator_stability_margin"], 0)
        self.assertLessEqual(np.linalg.norm(padded - second), coarse_bound["norm_error_upper_bound"] + fine_bound["norm_error_upper_bound"] + 1e-12)
        self.assertLess(coarse_bound["norm_error_upper_bound"], 0.1)

    def test_collective_joint_vector_matches_full_matrix_and_common_J_real_form(self):
        current = model(Parameters(sites=8, cutoff=4))
        block_size = 8 * 4
        full_hamiltonian = np.zeros((4 * block_size, 4 * block_size))
        full_unitary = np.zeros_like(full_hamiltonian, dtype=complex)
        for index, block in enumerate(current.blocks):
            indices = 4 * np.arange(block_size) + index
            full_hamiltonian[np.ix_(indices, indices)] = block["hamiltonian"]
            full_unitary[np.ix_(indices, indices)] = (block["vectors"] * np.exp(-0.37j * block["values"])) @ block["vectors"].T
        initial = current.initial_vector()
        actual = current.evolve(initial, 0.37).reshape(-1)
        np.testing.assert_allclose(actual, full_unitary @ initial.reshape(-1), rtol=0, atol=2e-14)
        density = np.outer(initial.reshape(-1), initial.reshape(-1).conj())
        lifted = real_lift(full_unitary)
        np.testing.assert_allclose(lifted @ encode_state(density) @ lifted.T,
                                   encode_state(np.outer(actual, actual.conj())), rtol=0, atol=2e-14)
        np.testing.assert_array_equal(full_hamiltonian, full_hamiltonian.T)

    def test_terminal_record_return_and_finite_random_output_are_not_free_or_early(self):
        sites = 16
        schedule = terminal_schedule(sites, 2.)
        deadline = max(event["received_at_site_zero"] for event in schedule)
        self.assertEqual(sum(event["hops"] for event in schedule), 72)
        self.assertAlmostEqual(deadline, 92.)
        bits = [int(index == 6) for index in range(sites)]
        with self.assertRaises(ValueError):
            release_smooth_record(bits, (0, 1), 0, 0, deadline - 1e-6, schedule)
        thresholds, _ = postprocessing_table(sites, center_fraction=0.5)
        count = sum(release_smooth_record(bits, (0, 1), 0, draw, deadline, schedule)[0] == 0 for draw in range(4096))
        self.assertEqual(count, thresholds[6])

    def test_predicted_background_excitation_exceeds_the_declared_cutoff_error(self):
        current = model()
        initial = current.initial_vector()
        final = current.evolve(initial, 2.)
        probability = current.diagnostics(final)["background_nonvacuum_probability"]
        bound = current.cutoff_error_budget(initial, 2.)["same_final_record_TV_upper_bound"]
        self.assertGreater(probability - bound, 0.002)

    def test_invalid_mode_cutoff_and_coupling_are_rejected(self):
        for arguments in ({"cutoff": 2}, {"coupling": 2.}, {"frequency": 0}, {"sites": 10}):
            with self.assertRaises(ValueError):
                Parameters(**arguments)


def run_report():
    current = model()
    initial = current.initial_vector()
    rows = [{"duration": time, **current.diagnostics(current.evolve(initial, time))} for time in (0., 0.5, 1., 2.)]
    coupled = current.evolve(initial, 2.)
    other_data = current.evolve(current.initial_vector(initial_data(current.parameters.sites, -1)), 2.)
    other_mode = current.evolve(current.initial_vector(level=1), 2.)
    data, field, clocks = current.reduced_states(coupled)
    _, data_changed_field, _ = current.reduced_states(other_data)
    mode_changed_data, _, mode_changed_clocks = current.reduced_states(other_mode)
    uncoupled = model(Parameters(coupling=0.))
    uncoupled_final = uncoupled.evolve(uncoupled.initial_vector(), 2.)
    frozen_records = uncoupled.joint_record_probabilities(uncoupled_final)
    final_records = current.joint_record_probabilities(coupled)
    current_budget = current.cutoff_error_budget(initial, 2.)
    background_probability = current.diagnostics(coupled)["background_nonvacuum_probability"]
    comparisons = []
    reference = model(Parameters(cutoff=14))
    reference_final = reference.evolve(reference.initial_vector(), 2.)
    for cutoff in (6, 8, 10, 12):
        candidate = model(Parameters(cutoff=cutoff))
        candidate_initial = candidate.initial_vector()
        final = candidate.evolve(candidate_initial, 2.)
        padded = np.zeros_like(reference_final)
        padded[:, :cutoff, :] = final
        comparisons.append({"cutoff": cutoff, "joint_complex_sector_dimension": 4 * current.parameters.sites * cutoff,
                            "norm_difference_to_cutoff14_diagnostic": float(np.linalg.norm(padded - reference_final)),
                            "boundary_budget": candidate.cutoff_error_budget(candidate_initial, 2.)})
    resources = final_record_resources(current.parameters.sites, 2.)
    resources["qubit_slots_before_detectors"] = resources.pop("total_quantum_slots_before_detectors")
    resources.update({"additional_quantum_background_modes": 1, "mode_level_cutoff": current.parameters.cutoff,
                      "total_quantum_subsystems_including_one_background_mode": current.parameters.sites + 3,
                      "data_clock_background_complex_sector_dimension": 4 * current.parameters.sites * current.parameters.cutoff,
                      "background_binary_record_bits_at_collector": 1,
                      "classical_raw_outcome_bits": current.parameters.sites + 3,
                      "background_binary_reads": 1,
                      "fresh_independent_fair_smooth_record_bits": 12,
                      "one_collective_mode_couples_multiple_sites_not_a_local_field": True})
    schedule = terminal_schedule(current.parameters.sites, 2.)
    return {"model_version": "C6", "route_round": "06", "research_scope": "sufficiency only",
            "parameters": vars(current.parameters),
            "Hamiltonian": "H_C5 + sum_clock(omega_j Z_j/2) + Omega*n + epsilon*(B_data + sum_clock(omega_j f_j Z_j/2))*(q^2-1/2)",
            "background_deformation": "ell_j=ell0_j*(I+epsilon*f_j*(q^2-1/2)); same factor on local hopping speed; b0 fixed",
            "mode_profile": "f(x)=(1+cos(2*pi*(x-0.35L)/L))/2",
            "projected_q_squared_is_not_square_of_truncated_q": True,
            "coefficient_positivity_lower_bound": 1 - current.parameters.coupling / 2,
            "sufficient_oscillator_stability_margin": current.parameters.frequency - 2 * current.parameters.coupling * current.source_norm,
            "time_rows": rows,
            "reciprocal_effects": {"background_trace_distance_for_two_data_preparations": trace_distance(field, data_changed_field),
                                   "data_trace_distance_vacuum_vs_level_one_background": trace_distance(data, mode_changed_data),
                                   "clock_trace_distance_vacuum_vs_level_one_background": trace_distance(clocks, mode_changed_clocks)},
            "joint_terminal_probabilities_data_background_clockA_clockC": final_records.tolist(),
            "default_cutoff_record_budget": current_budget,
            "background_nonvacuum_probability_with_cutoff_allowance": [
                max(0., background_probability - current_budget["same_final_record_TV_upper_bound"]),
                min(1., background_probability + current_budget["same_final_record_TV_upper_bound"])],
            "frozen_C5_mean_background_comparison": {
                "same_initial_data_clocks_and_mode": True,
                "background_mean_initial_F": 0.,
                "full_record_TV_difference": float(np.abs(final_records - frozen_records).sum() / 2),
                "data_clock_record_TV_difference": float(np.abs(final_records.sum(axis=1) - frozen_records.sum(axis=1)).sum() / 2),
                "interpretation": "specified coupled versus uncoupled task, not a necessity or universal-optimality argument"},
            "cutoff_comparisons": comparisons,
            "cutoff14_reference_budget": reference.cutoff_error_budget(reference.initial_vector(), 2.),
            "resources": resources, "record_delivery_schedule": schedule,
            "all_terminal_records_ready_at": max(event["received_at_site_zero"] for event in schedule),
            "fine_position_records_retained_but_reported_data_bit_is_smoothed": True,
            "background_binary_read_is_vacuum_vs_nonvacuum_at_collector": True,
            "clock_interaction_is_in_joint_state_not_external_readout_only": True,
            "C5_static_continuum_error_bound_is_automatically_valid_for_C6": False,
            "finite_mode_is_a_spatially_local_gravitational_field": False,
            "Hamiltonian_or_shared_deformation_derived_from_cognition": False,
            "metric_field_equations_or_Einstein_dynamics_derived": False,
            "roundoff_errors_are_rigorously_certified": False,
            "preparation_clock_detector_and_mode_control_energy_resources_complete": False,
            "count_samples_generated": 0, "laboratory_samples": 0,
            "necessity_arguments_performed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(QuantumChainBackgroundTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    results = run_report()
    results["automated_checks"] = {"run": tests.testsRun, "failures": 0, "errors": 0}
    if args.write_results:
        Path(__file__).with_name("quantum_chain_background_results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"tests": tests.testsRun, "time_rows": results["time_rows"],
                      "reciprocal_effects": results["reciprocal_effects"],
                      "cutoffs": results["cutoff_comparisons"]}, indent=2))


if __name__ == "__main__":
    main()