"""Construction 04 / C4: a controlled two-carrier envelope limit.

The same midpoint exchange rule as C3 is put on a smooth periodic long chain.
The one-excitation sector is an invariant computational restriction of N
qubits, not log2(N) physical slots. Taylor plus unitary Duhamel estimates give
an a^2 finite-time error for fixed smooth envelope data, and hence a bound on
finite position records. Background, synchrony and preparation remain inputs.
"""

import argparse
import json
import math
import unittest
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

from baseline import embed_operator
from propagation import EXCHANGE
from clock_calibration import clock_gate, clock_effect


@dataclass(frozen=True)
class SmoothBackground:
    length: float = 20.0
    speed_scale: float = 0.6
    modulation: float = 0.25
    lapse_modulation: float = 0.15

    def __post_init__(self):
        values = (self.length, self.speed_scale, self.modulation, self.lapse_modulation)
        if (not all(math.isfinite(value) for value in values) or self.length <= 0 or self.speed_scale <= 0
                or abs(self.modulation) >= 1 or abs(self.lapse_modulation) >= 1):
            raise ValueError("Use finite positive scales and strictly positive periodic profiles.")

    @property
    def frequency(self):
        return 2 * math.pi / self.length

    @property
    def optical_length(self):
        return self.length / self.speed_scale

    def inverse_speed_derivatives(self, positions):
        angle = self.frequency * np.asarray(positions)
        return ((1 + self.modulation * np.cos(angle)) / self.speed_scale,
                -self.modulation * self.frequency * np.sin(angle) / self.speed_scale,
                -self.modulation * self.frequency ** 2 * np.cos(angle) / self.speed_scale,
                self.modulation * self.frequency ** 3 * np.sin(angle) / self.speed_scale)

    def speed(self, positions):
        return 1 / self.inverse_speed_derivatives(positions)[0]

    def speed_derivative(self, positions):
        inverse, derivative, _, _ = self.inverse_speed_derivatives(positions)
        return -derivative / inverse ** 2

    def lapse(self, positions):
        return 1 + self.lapse_modulation * np.cos(self.frequency * np.asarray(positions))

    def spatial_scale(self, positions):
        return self.lapse(positions) / self.speed(positions)

    def optical_coordinate(self, positions):
        positions = np.asarray(positions)
        return (positions + self.modulation * np.sin(self.frequency * positions) / self.frequency) / self.speed_scale

    def position_from_optical(self, coordinates):
        target = np.mod(coordinates, self.optical_length)
        position = self.speed_scale * target
        for _ in range(12):
            position -= (self.optical_coordinate(position) - target) * self.speed(position)
        if np.max(np.abs(self.optical_coordinate(position) - target)) > 1e-11:
            raise ValueError("Optical-coordinate inversion did not converge.")
        return np.mod(position, self.length)


BACKGROUND = SmoothBackground()


def grid(sites, background=BACKGROUND):
    if not isinstance(sites, int) or sites < 8 or sites % 4:
        raise ValueError("Use at least eight sites and a multiple of four for periodic carriers.")
    spacing = background.length / sites
    return np.arange(sites) * spacing, spacing


def hopping_matrix(sites, background=BACKGROUND):
    positions, spacing = grid(sites, background)
    bonds = background.speed(positions + spacing / 2) / (2 * spacing)
    hamiltonian = np.zeros((sites, sites))
    for index, bond in enumerate(bonds):
        neighbor = (index + 1) % sites
        hamiltonian[index, neighbor] = hamiltonian[neighbor, index] = bond
    return hamiltonian


@lru_cache(maxsize=12)
def spectral_data(sites, background=BACKGROUND):
    return np.linalg.eigh(hopping_matrix(sites, background))


def evolve(vector, duration, background=BACKGROUND):
    if not math.isfinite(duration):
        raise ValueError("Use finite evolution duration.")
    values, vectors = spectral_data(len(vector), background)
    return vectors @ (np.exp(-1j * duration * values) * (vectors.T @ vector))


def optical_modes(background=BACKGROUND):
    indices = np.arange(-2, 3)
    wave_numbers = 2 * math.pi * indices / background.optical_length
    coefficients = np.exp(-0.5 * (indices / 1.2) ** 2)
    coefficients /= math.sqrt(background.optical_length * float(np.sum(coefficients ** 2)))
    return wave_numbers, coefficients


def optical_wave(coordinates, derivative=0, background=BACKGROUND):
    wave_numbers, coefficients = optical_modes(background)
    center = background.optical_coordinate(background.length * 0.35)
    phases = np.exp(1j * np.outer(np.asarray(coordinates) - center, wave_numbers))
    return phases @ (coefficients * (1j * wave_numbers) ** derivative)


def envelope(positions, duration, direction, background=BACKGROUND):
    if direction not in (-1, 1):
        raise ValueError("Carrier direction must be -1 or +1.")
    inverse_speed = background.inverse_speed_derivatives(positions)[0]
    shifted = background.optical_coordinate(positions) - direction * duration
    return np.sqrt(inverse_speed) * optical_wave(shifted, background=background)


def envelope_generator(positions, duration, direction, background=BACKGROUND):
    inverse, derivative, _, _ = background.inverse_speed_derivatives(positions)
    shifted = background.optical_coordinate(positions) - direction * duration
    root = np.sqrt(inverse)
    first = optical_wave(shifted, derivative=1, background=background)
    function = optical_wave(shifted, background=background)
    spatial_derivative = derivative * function / (2 * root) + root * inverse * first
    corrected = -direction * (background.speed(positions) * spatial_derivative
                              + background.speed_derivative(positions) * root * function / 2)
    return corrected


def carrier(sites, direction):
    if direction not in (-1, 1):
        raise ValueError("Carrier direction must be -1 or +1.")
    return (-1j * direction) ** (np.arange(sites) % 4)


def sampled_solution(sites, duration, direction, background=BACKGROUND):
    positions, spacing = grid(sites, background)
    return math.sqrt(spacing) * carrier(sites, direction) * envelope(positions, duration, direction, background)


def derivative_bounds(background=BACKGROUND):
    wave_numbers, coefficients = optical_modes(background)
    moments = [float(np.sum(np.abs(coefficients) * np.abs(wave_numbers) ** order)) for order in range(4)]
    minimum = (1 - abs(background.modulation)) / background.speed_scale
    maximum = (1 + abs(background.modulation)) / background.speed_scale
    first, second, third = [abs(background.modulation) * background.frequency ** order / background.speed_scale
                             for order in (1, 2, 3)]
    roots = [math.sqrt(maximum), first / (2 * math.sqrt(minimum)),
             second / (2 * math.sqrt(minimum)) + first ** 2 / (4 * minimum ** 1.5),
             third / (2 * math.sqrt(minimum)) + 3 * first * second / (4 * minimum ** 1.5)
             + 3 * first ** 3 / (8 * minimum ** 2.5)]
    bounds = [roots[0] * moments[0],
              roots[1] * moments[0] + roots[0] * maximum * moments[1],
              roots[2] * moments[0] + (2 * roots[1] * maximum + roots[0] * first) * moments[1]
              + roots[0] * maximum ** 2 * moments[2],
              roots[3] * moments[0] + (3 * roots[2] * maximum + 3 * roots[1] * first + roots[0] * second) * moments[1]
              + (3 * roots[1] * maximum ** 2 + 3 * roots[0] * maximum * first) * moments[2]
              + roots[0] * maximum ** 3 * moments[3]]
    speed_bounds = [1 / minimum, first / minimum ** 2,
                    2 * first ** 2 / minimum ** 3 + second / minimum ** 2,
                    6 * first ** 3 / minimum ** 4 + 6 * first * second / minimum ** 3 + third / minimum ** 2]
    third_flux = (speed_bounds[0] * bounds[3] + 1.5 * speed_bounds[1] * bounds[2]
                  + 0.75 * speed_bounds[2] * bounds[1] + 0.125 * speed_bounds[3] * bounds[0])
    return {"envelope_derivative_upper_bounds": bounds, "speed_derivative_upper_bounds": speed_bounds,
            "shifted_flux_third_derivative_upper_bound": third_flux}


def vector_error_bound(sites, duration, initial_norm, amplitude_sum=1., background=BACKGROUND):
    if initial_norm <= 0 or amplitude_sum < 0 or not math.isfinite(duration):
        raise ValueError("A positive sampled initial norm and finite duration are required.")
    _, spacing = grid(sites, background)
    coefficient = derivative_bounds(background)["shifted_flux_third_derivative_upper_bound"]
    return abs(duration) * math.sqrt(background.length) * spacing ** 2 * coefficient * amplitude_sum / (6 * initial_norm)


def cell_quadrature_bound(sites, background=BACKGROUND):
    _, spacing = grid(sites, background)
    bounds = derivative_bounds(background)["envelope_derivative_upper_bounds"]
    density_second = 2 * (bounds[1] ** 2 + bounds[0] * bounds[2])
    return background.length * spacing ** 2 * density_second / 24


def continuum_cell_probabilities(sites, duration, direction, background=BACKGROUND):
    if direction not in (-1, 1):
        raise ValueError("Carrier direction must be -1 or +1.")
    positions, spacing = grid(sites, background)
    lower = background.optical_coordinate(positions - spacing / 2) - direction * duration
    upper = background.optical_coordinate(positions + spacing / 2) - direction * duration
    width = upper - lower
    center = background.optical_coordinate(background.length * 0.35)
    midpoint = (upper + lower) / 2 - center
    wave_numbers, coefficients = optical_modes(background)
    differences = wave_numbers[:, None] - wave_numbers[None, :]
    weights = coefficients[:, None] * coefficients.conj()[None, :]
    oscillation = np.exp(1j * midpoint[:, None, None] * differences)
    integrals = width[:, None, None] * oscillation * np.sinc(width[:, None, None] * differences / (2 * math.pi))
    probabilities = np.sum(integrals * weights, axis=(1, 2)).real
    if probabilities.min() < -1e-12 or abs(probabilities.sum() - 1) > 1e-12:
        raise ValueError("Integrated continuum probabilities failed positivity or normalization.")
    return probabilities


def normalized_trace_distance(first, second):
    first = first / np.linalg.norm(first)
    second = second / np.linalg.norm(second)
    overlap = min(1., float(abs(np.vdot(first, second)) ** 2))
    return math.sqrt(max(0., 1 - overlap))


def compare(sites, duration=2., direction=1, background=BACKGROUND):
    initial = sampled_solution(sites, 0., direction, background)
    initial_norm = float(np.linalg.norm(initial))
    exact = evolve(initial / initial_norm, duration, background)
    reference = sampled_solution(sites, duration, direction, background) / initial_norm
    vector_error = float(np.linalg.norm(exact - reference))
    bound = vector_error_bound(sites, duration, initial_norm, background=background)
    normalized_reference = reference / np.linalg.norm(reference)
    positions, spacing = grid(sites, background)
    window = (positions >= 0.4 * background.length) & (positions < 0.65 * background.length)
    exact_probability = float(np.sum(abs(exact[window]) ** 2))
    reference_probability = float(np.sum(abs(normalized_reference[window]) ** 2))
    cell_probabilities = continuum_cell_probabilities(sites, duration, direction, background)
    quadrature_bound = cell_quadrature_bound(sites, background)
    return {"sites": sites, "spacing": spacing, "direction": direction, "duration": duration,
            "sampled_initial_norm": initial_norm, "sampled_reference_norm": float(np.linalg.norm(reference)),
            "exact_norm": float(np.linalg.norm(exact)), "vector_error_diagnostic": vector_error,
            "analytic_vector_error_upper_bound": bound,
            "trace_distance_diagnostic": normalized_trace_distance(exact, normalized_reference),
            "analytic_all_final_record_TV_upper_bound": min(1., 2 * bound),
            "position_record_TV_diagnostic": float(np.sum(np.abs(abs(exact) ** 2 - abs(normalized_reference) ** 2)) / 2),
            "window_probability_exact": exact_probability, "window_probability_envelope": reference_probability,
            "window_probability_error": abs(exact_probability - reference_probability),
            "continuum_cell_record_TV_diagnostic": float(np.abs(abs(exact) ** 2 - cell_probabilities).sum() / 2),
            "sample_to_integrated_cell_TV_upper_bound": quadrature_bound,
            "analytic_integrated_cell_record_TV_upper_bound": min(1., 2 * bound + quadrature_bound),
            "unitary_energy_drift": float(abs(np.vdot(exact, hopping_matrix(sites, background) @ exact)
                                             - np.vdot(initial / initial_norm, hopping_matrix(sites, background) @ (initial / initial_norm))))}


def final_record_resources(sites, duration, background=BACKGROUND):
    _, spacing = grid(sites, background)
    bit_hops = sum(min(index, sites - index) for index in range(sites)) + sites // 2
    edge_delay = spacing
    return {"chain_qubit_slots": sites, "occupied_excitation_sector_dimension": sites,
            "retained_clock_slots": 2, "total_quantum_slots_before_detectors": sites + 2,
            "local_destructive_position_readouts": sites, "clock_reads": 2,
            "classical_raw_outcome_bits": sites + 2,
            "minimum_fixed_label_bits_for_this_single_excitation_support": math.ceil(math.log2(sites)),
            "classical_bit_hops_for_position_bits": sites * sites // 4,
            "classical_bit_hops_including_clock_bits": bit_hops,
            "declared_serial_hop_parameter": edge_delay,
            "latest_record_collection_parameter": duration + bit_hops * edge_delay,
            "record_transport_is_serial_no_queue_optimization": True, "resets": 0,
            "initial_smooth_carrier_preparation_is_a_declared_input": True,
            "full_detector_memory_clock_energy_and_preparation_gate_costs_closed": False}


def collection_schedule(sites, duration, background=BACKGROUND):
    _, spacing = grid(sites, background)
    current = duration
    events = []
    sources = [("position", index) for index in range(sites)] + [("clock", 0), ("clock", sites // 2)]
    for kind, source in sources:
        hops = min(source, sites - source)
        current += hops * spacing
        events.append({"kind": kind, "source_site": source, "read_at": duration,
                       "hops": hops, "received_at_site_zero": current})
    return events


def classify_position_record(bits, received_at, schedule, decision_table):
    if len(bits) != len(decision_table) or any(bit not in (0, 1) for bit in bits) or sum(bits) != 1:
        raise ValueError("Use the complete single-excitation position record.")
    if received_at < max(event["received_at_site_zero"] for event in schedule):
        raise ValueError("The complete record and clock bits have not arrived.")
    return int(decision_table[bits.index(1)])


def final_clock_probabilities(duration, background=BACKGROUND):
    probabilities = []
    for location, bare, basis in ((0., 1., "X"), (background.length / 2, 1.2, "Y")):
        rate = bare * background.lapse(location)
        gate = clock_gate(rate, duration)
        state = gate @ (np.ones((2, 2)) / 2) @ gate.conj().T
        probabilities.append([float(np.trace(clock_effect(basis, outcome) @ state).real) for outcome in (0, 1)])
    return np.outer(probabilities[0], probabilities[1])


def finite_record_report(sites=128, duration=2., background=BACKGROUND):
    records, approximations, bounds = [], [], []
    clocks = final_clock_probabilities(duration, background)
    for direction in (1, -1):
        initial = sampled_solution(sites, 0., direction, background)
        normalization = np.linalg.norm(initial)
        evolved = evolve(initial / normalization, duration, background)
        records.append(abs(evolved) ** 2)
        approximations.append(continuum_cell_probabilities(sites, duration, direction, background))
        bounds.append(min(1., 2 * vector_error_bound(sites, duration, normalization, background=background)
                  + cell_quadrature_bound(sites, background)))
    decisions = np.where(approximations[0] >= approximations[1], 1, -1)
    exact_success = float((records[0][decisions == 1].sum() + records[1][decisions == -1].sum()) / 2)
    predicted_success = float((approximations[0][decisions == 1].sum() + approximations[1][decisions == -1].sum()) / 2)
    schedule = collection_schedule(sites, duration, background)
    joint_record_errors = [float(np.sum(np.abs(np.einsum("n,ab->nab", actual - approximate, clocks))) / 2)
                           for actual, approximate in zip(records, approximations)]
    return {"sites": sites, "duration": duration, "task": "equal-prior classification of two declared directional preparations after terminal position readout",
            "decision_table": decisions.tolist(), "policy_source": "integrated continuum cell probabilities only",
            "predicted_success": predicted_success, "exact_lattice_success": exact_success,
            "success_error": abs(exact_success - predicted_success),
            "success_error_upper_bound": sum(bounds) / 2,
            "position_and_two_clock_joint_record_TV": joint_record_errors,
            "clock_joint_probabilities": clocks.tolist(), "resources": final_record_resources(sites, duration, background),
            "record_collection_schedule": schedule,
            "decision_not_allowed_before": max(event["received_at_site_zero"] for event in schedule),
            "local_clock_preparation_and_synchronous_readout_are_inputs": True,
            "approximation_used_only_before_terminal_measurement": True,
            "sample_counts_generated": 0}


class EnvelopeLimitTests(unittest.TestCase):
    def test_midpoint_hopping_is_hermitian_local_and_extends_the_same_exchange_term(self):
        sites = 8
        positions, spacing = grid(sites)
        hamiltonian = hopping_matrix(sites)
        np.testing.assert_array_equal(hamiltonian, hamiltonian.T)
        self.assertEqual(np.count_nonzero(hamiltonian), 2 * sites)
        full = np.zeros((2 ** sites, 2 ** sites), dtype=complex)
        for index in range(sites):
            full += BACKGROUND.speed(positions[index] + spacing / 2) / (2 * spacing) * embed_operator(EXCHANGE, (index, (index + 1) % sites), sites)
        sector = [1 << (sites - 1 - index) for index in range(sites)]
        np.testing.assert_allclose(full[np.ix_(sector, sector)], hamiltonian, rtol=0, atol=1e-15)
        outside = sorted(set(range(2 ** sites)) - set(sector))
        np.testing.assert_array_equal(full[np.ix_(outside, sector)], np.zeros((len(outside), sites)))

    def test_carrier_demodulation_matches_exact_difference_generator_for_both_directions(self):
        sites = 32
        positions, spacing = grid(sites)
        function = envelope(positions, 0.3, 1)
        forward = BACKGROUND.speed(positions + spacing / 2) * np.roll(function, -1)
        backward = BACKGROUND.speed(positions - spacing / 2) * np.roll(function, 1)
        for direction in (-1, 1):
            modulation = carrier(sites, direction)
            actual = (-1j * (hopping_matrix(sites) @ (modulation * function))) / modulation
            np.testing.assert_allclose(actual, -direction * (forward - backward) / (2 * spacing), rtol=0, atol=2e-15)

    def test_optical_solution_satisfies_variable_speed_equation_including_half_derivative(self):
        positions, _ = grid(64)
        inverse = BACKGROUND.inverse_speed_derivatives(positions)[0]
        for direction in (-1, 1):
            exact_time_derivative = -direction * np.sqrt(inverse) * optical_wave(BACKGROUND.optical_coordinate(positions) - direction * 0.7, derivative=1)
            np.testing.assert_allclose(envelope_generator(positions, 0.7, direction), exact_time_derivative, rtol=0, atol=2e-15)
        self.assertGreater(np.max(np.abs(BACKGROUND.speed_derivative(positions) * envelope(positions, 0.7, 1))), 0.01)

    def test_taylor_residual_bound_covers_both_carriers_without_fitting_a_rate(self):
        coefficient = derivative_bounds()["shifted_flux_third_derivative_upper_bound"]
        for sites in (32, 64, 128):
            positions, spacing = grid(sites)
            for direction in (-1, 1):
                reference = sampled_solution(sites, 0.71, direction)
                derivative = math.sqrt(spacing) * carrier(sites, direction) * envelope_generator(positions, 0.71, direction)
                residual = derivative + 1j * hopping_matrix(sites) @ reference
                self.assertLessEqual(np.linalg.norm(residual), math.sqrt(BACKGROUND.length) * spacing ** 2 * coefficient / 6 + 1e-12)

    def test_finite_time_convergence_and_record_bounds(self):
        for direction in (-1, 1):
            rows = [compare(sites, direction=direction) for sites in (32, 64, 128)]
            for row in rows:
                self.assertAlmostEqual(row["exact_norm"], 1)
                self.assertLessEqual(row["vector_error_diagnostic"], row["analytic_vector_error_upper_bound"])
                self.assertLessEqual(row["position_record_TV_diagnostic"], row["analytic_all_final_record_TV_upper_bound"])
                self.assertLessEqual(row["window_probability_error"], row["analytic_all_final_record_TV_upper_bound"])
                self.assertLessEqual(row["continuum_cell_record_TV_diagnostic"], row["analytic_integrated_cell_record_TV_upper_bound"])
                self.assertLess(row["unitary_energy_drift"], 1e-12)
            ratios = [rows[index]["vector_error_diagnostic"] / rows[index + 1]["vector_error_diagnostic"] for index in (0, 1)]
            self.assertTrue(all(3.5 < ratio < 4.5 for ratio in ratios))

    def test_coherent_two_carrier_preparation_obeys_the_sum_residual_bound(self):
        sites, duration = 64, 2.
        weights = (math.sqrt(0.6), 1j * math.sqrt(0.4))
        initial = sum(weight * sampled_solution(sites, 0., direction) for weight, direction in zip(weights, (1, -1)))
        norm = np.linalg.norm(initial)
        reference = sum(weight * sampled_solution(sites, duration, direction) for weight, direction in zip(weights, (1, -1))) / norm
        evolved = evolve(initial / norm, duration)
        bound = vector_error_bound(sites, duration, norm, amplitude_sum=sum(abs(weight) for weight in weights))
        self.assertLessEqual(np.linalg.norm(evolved - reference), bound)
        self.assertAlmostEqual(np.linalg.norm(evolved), 1)

    def test_cross_carrier_overlap_has_a_smooth_alternating_sum_bound(self):
        bounds = derivative_bounds()["envelope_derivative_upper_bounds"]
        product_second = 2 * bounds[0] * bounds[2] + 2 * bounds[1] ** 2
        for sites in (32, 64, 128):
            _, spacing = grid(sites)
            overlap = abs(np.vdot(sampled_solution(sites, 0.7, 1), sampled_solution(sites, 0.7, -1)))
            self.assertLessEqual(overlap, BACKGROUND.length * spacing ** 2 * product_second / 4 + 1e-12)

    def test_integrated_cell_probabilities_are_normalized_and_obey_midpoint_error_bound(self):
        for sites in (32, 64, 128):
            for direction in (-1, 1):
                integrated = continuum_cell_probabilities(sites, 2., direction)
                sampled = abs(sampled_solution(sites, 2., direction)) ** 2
                budget = cell_quadrature_bound(sites)
                self.assertAlmostEqual(integrated.sum(), 1)
                self.assertGreaterEqual(integrated.min(), 0)
                self.assertLessEqual(np.abs(integrated - sampled).sum(), budget + 1e-12)
                self.assertLessEqual(np.abs(integrated - sampled / sampled.sum()).sum() / 2, budget + 1e-12)

    def test_optical_fourier_integral_matches_resolved_numerical_cell_quadrature(self):
        sites = 32
        positions, spacing = grid(sites)
        integrated = continuum_cell_probabilities(sites, 0.7, 1)
        for index in (0, 11, 31):
            fine = np.linspace(positions[index] - spacing / 2, positions[index] + spacing / 2, 4097)
            reference = np.trapezoid(abs(envelope(fine, 0.7, 1)) ** 2, fine)
            self.assertAlmostEqual(integrated[index], reference, delta=2e-10)

    def test_characteristics_and_clock_rates_use_one_positive_metric_profile(self):
        positions, _ = grid(64)
        np.testing.assert_allclose(BACKGROUND.lapse(positions) / BACKGROUND.spatial_scale(positions), BACKGROUND.speed(positions), rtol=0, atol=2e-15)
        position = 7.
        duration = 1e-4
        optical = BACKGROUND.optical_coordinate(position)
        for direction in (-1, 1):
            future = BACKGROUND.position_from_optical(optical + direction * duration)
            past = BACKGROUND.position_from_optical(optical - direction * duration)
            self.assertAlmostEqual((future - past) / (2 * duration), direction * BACKGROUND.speed(position), delta=1e-8)
        for bare in (1., 1.2):
            rate = bare * BACKGROUND.lapse(position)
            plus = np.ones((2, 2)) / 2
            gate = clock_gate(rate, 0.4)
            self.assertAlmostEqual(np.trace(clock_effect("Y", 0) @ gate @ plus @ gate.conj().T).real,
                                   (1 + math.sin(rate * 0.4)) / 2)

    def test_finite_position_measurements_and_retained_records_fit_probability_error(self):
        sites = 64
        row = compare(sites)
        duration = row["duration"]
        initial = sampled_solution(sites, 0, 1)
        evolved = evolve(initial / np.linalg.norm(initial), duration)
        reference = sampled_solution(sites, duration, 1)
        reference /= np.linalg.norm(reference)
        true_records = abs(evolved) ** 2
        model_records = abs(reference) ** 2
        self.assertAlmostEqual(true_records.sum(), 1)
        self.assertLessEqual(np.abs(true_records - model_records).sum() / 2, row["analytic_all_final_record_TV_upper_bound"])
        resources = final_record_resources(sites, duration)
        self.assertEqual(resources["classical_bit_hops_including_clock_bits"], sites * sites // 4 + sites // 2)
        self.assertGreater(resources["latest_record_collection_parameter"], duration)
        self.assertEqual(resources["chain_qubit_slots"], sites)
        self.assertEqual(resources["resets"], 0)

    def test_terminal_clock_and_position_law_keeps_the_same_record_error_bound(self):
        result = finite_record_report(128)
        self.assertAlmostEqual(sum(sum(row) for row in result["clock_joint_probabilities"]), 1)
        self.assertLessEqual(result["success_error"], result["success_error_upper_bound"])
        self.assertTrue(all(value <= result["success_error_upper_bound"] for value in result["position_and_two_clock_joint_record_TV"]))
        self.assertGreater(result["exact_lattice_success"], 0.5)

    def test_record_routing_does_not_allow_an_early_or_incomplete_decision(self):
        sites = 32
        result = finite_record_report(sites)
        schedule = result["record_collection_schedule"]
        deadline = result["decision_not_allowed_before"]
        self.assertAlmostEqual(deadline, result["resources"]["latest_record_collection_parameter"])
        self.assertEqual(sum(event["hops"] for event in schedule), sites * sites // 4 + sites // 2)
        for index in (0, sites // 3, sites - 1):
            bits = [int(site == index) for site in range(sites)]
            with self.assertRaises(ValueError):
                classify_position_record(bits, deadline - 1e-8, schedule, result["decision_table"])
            self.assertEqual(classify_position_record(bits, deadline, schedule, result["decision_table"]), result["decision_table"][index])
        with self.assertRaises(ValueError):
            classify_position_record([0] * sites, deadline, schedule, result["decision_table"])

    def test_grid_and_background_contracts_reject_invalid_inputs(self):
        for sites in (0, 7, 18):
            with self.assertRaises(ValueError):
                grid(sites)
        with self.assertRaises(ValueError):
            SmoothBackground(modulation=1)
        with self.assertRaises(ValueError):
            sampled_solution(32, 0, 0)


def report():
    rows = [compare(sites, direction=direction) for direction in (1, -1) for sites in (32, 64, 128, 256)]
    return {"model_version": "C4", "route_round": "04", "research_scope": "sufficiency only",
            "parent_rule": "C3 nearest-neighbor midpoint exchange g=v/(2a), now on a periodic long-chain family",
            "background_inputs": vars(BACKGROUND),
            "profile": "v=v0/(1+eta*cos(2*pi*x/L)); ell=1+0.15*cos(2*pi*x/L); b=ell/v",
            "carrier_convention": "c_s(j)=exp(-i*s*pi*j/2), s=+1 right-moving, s=-1 left-moving",
            "derived_envelope_equation": "partial_t chi_s=-s*(v partial_x+v_prime/2)chi_s",
            "optical_solution": "chi_s(x,t)=f(y(x)-s*t)/sqrt(v(x)), y_prime=1/v",
            "prepared_optical_fourier_modes": [-2, -1, 0, 1, 2],
            "derivative_bounds": derivative_bounds(),
            "vector_error_theorem": "||U_a(t)psi_a(0)-psi_a_cont(t)|| <= |t|sqrt(L)a^2*M3/(6*initial_sample_norm)",
            "final_record_bound": "TV <= min(1,2*vector_bound) after normalizing sampled continuum reference",
            "integrated_cell_record_bound": "TV <= min(1,2*vector_bound + L*a^2*sup_abs_density_second_derivative/24)",
            "rows": rows,
            "resource_examples": [final_record_resources(sites, 2.) for sites in (32, 64, 128, 256)],
            "finite_position_clock_and_decision_task": finite_record_report(),
            "physical_slots_are_not_log2_sector_dimension": True,
            "mass_coupling_added": False,
            "all_lattice_states_or_times_uniformly_approximated": False,
            "two_envelopes_are_two_disjoint_copies_of_the_full_lattice": False,
            "fixed_finite_lattice_has_an_exact_relativistic_causal_cone": False,
            "scope_of_geometry": "massless two-direction transport on a specified static 1+1 metric; no dynamical gravity",
            "background_source_or_Einstein_equations_derived": False,
            "detector_synchrony_and_clock_dynamics_are_working_inputs": True,
            "no_new_count_samples_collected": True,
            "necessity_arguments_performed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(EnvelopeLimitTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    results = report()
    results["automated_checks"] = {"run": checks.testsRun, "failures": 0, "errors": 0}
    if args.write_results:
        Path(__file__).with_name("envelope_limit_results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"tests": checks.testsRun, "bounds": results["derivative_bounds"], "rows": results["rows"]}, indent=2))


if __name__ == "__main__":
    main()