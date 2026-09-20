"""Construction 05 / C5: staggered mass on the same midpoint hopping chain.

A declared optical-coordinate lapse makes the mass coefficient a finite
Fourier function. A Hermitian Galerkin reference has a separately bounded
boundary defect. Sobolev derivative bounds and the C4 Taylor argument control
lattice propagation and specified smooth terminal records. This is a static
single-particle sufficiency construction, not a gravity or necessity proof.
"""

import argparse
import json
import math
import unittest
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

from envelope_limit import (SmoothBackground, grid, hopping_matrix, carrier, optical_modes,
                            sampled_solution, collection_schedule, final_record_resources,
                            final_clock_probabilities)
from baseline import IDENTITY, PAULI_X, PAULI_Z, embed_operator, binary_effect, real_lift
from common_orientation_structure import encode_state
from propagation import EXCHANGE
from clock_calibration import clock_gate, clock_effect


@dataclass(frozen=True)
class MassiveBackground(SmoothBackground):
    def lapse(self, positions):
        return 1 + self.lapse_modulation * np.cos(2 * math.pi * self.optical_coordinate(positions) / self.optical_length)


BACKGROUND = MassiveBackground()
DEFAULT_MASS = 0.35
DEFAULT_CUTOFF = 12


def validate_mass(mass):
    if not np.isfinite(mass) or mass < 0:
        raise ValueError("Use a finite nonnegative declared proper-mass coefficient.")


def massive_hopping(sites, mass=DEFAULT_MASS, background=BACKGROUND):
    validate_mass(mass)
    positions, _ = grid(sites, background)
    return hopping_matrix(sites, background) + np.diag((-1.) ** np.arange(sites) * mass * background.lapse(positions))


@lru_cache(maxsize=12)
def lattice_spectrum(sites, mass=DEFAULT_MASS, background=BACKGROUND):
    return np.linalg.eigh(massive_hopping(sites, mass, background))


def evolve_lattice(initial, duration, mass=DEFAULT_MASS, background=BACKGROUND):
    if not np.isfinite(duration):
        raise ValueError("Use finite evolution duration.")
    values, vectors = lattice_spectrum(len(initial), mass, background)
    return vectors @ (np.exp(-1j * duration * values) * (vectors.T @ initial))


def galerkin_hamiltonian(cutoff=DEFAULT_CUTOFF, mass=DEFAULT_MASS, background=BACKGROUND):
    validate_mass(mass)
    if not isinstance(cutoff, int) or cutoff < 3:
        raise ValueError("Use cutoff >=3 around the initial modes -2..2.")
    modes = np.arange(-cutoff, cutoff + 1)
    wave_numbers = 2 * math.pi * modes / background.optical_length
    hamiltonian = np.zeros((2 * len(modes), 2 * len(modes)), dtype=complex)
    neighbor = mass * background.lapse_modulation / 2
    for index, wave_number in enumerate(wave_numbers):
        block = slice(2 * index, 2 * index + 2)
        hamiltonian[block, block] = wave_number * PAULI_Z + mass * PAULI_X
        if index + 1 < len(modes):
            next_block = slice(2 * index + 2, 2 * index + 4)
            hamiltonian[block, next_block] = neighbor * PAULI_X
            hamiltonian[next_block, block] = neighbor * PAULI_X
    return hamiltonian


def initial_modes(cutoff=DEFAULT_CUTOFF, direction=1, background=BACKGROUND):
    if direction not in (-1, 1) or cutoff < 3:
        raise ValueError("Use a declared carrier direction and cutoff >=3.")
    result = np.zeros((2 * cutoff + 1, 2), dtype=complex)
    wave_numbers, coefficients = optical_modes(background)
    center = background.optical_coordinate(background.length * 0.35)
    result[cutoff - 2:cutoff + 3, 0 if direction == 1 else 1] = (
        np.sqrt(background.optical_length) * coefficients * np.exp(-1j * wave_numbers * center))
    return result


@lru_cache(maxsize=12)
def galerkin_spectrum(cutoff=DEFAULT_CUTOFF, mass=DEFAULT_MASS, background=BACKGROUND):
    return np.linalg.eigh(galerkin_hamiltonian(cutoff, mass, background))


def reference_modes(duration, cutoff=DEFAULT_CUTOFF, mass=DEFAULT_MASS, direction=1, background=BACKGROUND):
    values, vectors = galerkin_spectrum(cutoff, mass, background)
    initial = initial_modes(cutoff, direction, background).reshape(-1)
    evolved = vectors @ (np.exp(-1j * duration * values) * (vectors.conj().T @ initial))
    return evolved.reshape(-1, 2)


def optical_spinor(positions, modes, derivative=0, background=BACKGROUND):
    cutoff = (len(modes) - 1) // 2
    wave_numbers = 2 * math.pi * np.arange(-cutoff, cutoff + 1) / background.optical_length
    phases = np.exp(1j * np.outer(background.optical_coordinate(positions), wave_numbers))
    return phases @ ((1j * wave_numbers)[:, None] ** derivative * modes) / np.sqrt(background.optical_length)


def reconstruct(modes, sites, background=BACKGROUND):
    positions, spacing = grid(sites, background)
    spinor = optical_spinor(positions, modes, background=background) / np.sqrt(background.speed(positions))[:, None]
    return math.sqrt(spacing) * (carrier(sites, 1) * spinor[:, 0] + carrier(sites, -1) * spinor[:, 1])


def regularity_bounds(duration, mass=DEFAULT_MASS, background=BACKGROUND):
    validate_mass(mass)
    horizon = abs(duration)
    wave_numbers, coefficients = optical_modes(background)
    optical_length = background.optical_length
    initial_norms = [math.sqrt(optical_length * float(np.sum(abs(coefficients) ** 2 * abs(wave_numbers) ** (2 * order))))
                     for order in range(5)]
    fundamental = 2 * math.pi / optical_length
    modulation = abs(mass * background.lapse_modulation)
    polynomials = []
    for order in range(5):
        polynomial = np.zeros(order + 1)
        polynomial[0] = initial_norms[order]
        for derivative in range(1, order + 1):
            lower = polynomials[order - derivative]
            forcing = math.comb(order, derivative) * modulation * fundamental ** derivative
            polynomial[1:len(lower) + 1] += forcing * lower / np.arange(1, len(lower) + 1)
        polynomials.append(polynomial)
    sobolev = [float(np.polynomial.polynomial.polyval(horizon, polynomial)) for polynomial in polynomials]
    embedding = math.sqrt(0.5 / math.tanh(optical_length / 2))
    optical_bounds = [embedding * math.hypot(sobolev[order], sobolev[order + 1]) for order in range(4)]
    lower = (1 - abs(background.modulation)) / background.speed_scale
    upper = (1 + abs(background.modulation)) / background.speed_scale
    first, second, third = [abs(background.modulation) * background.frequency ** order / background.speed_scale
                            for order in (1, 2, 3)]
    roots = [math.sqrt(upper), first / (2 * math.sqrt(lower)),
             second / (2 * math.sqrt(lower)) + first ** 2 / (4 * lower ** 1.5),
             third / (2 * math.sqrt(lower)) + 3 * first * second / (4 * lower ** 1.5) + 3 * first ** 3 / (8 * lower ** 2.5)]
    spatial = [roots[0] * optical_bounds[0],
               roots[1] * optical_bounds[0] + roots[0] * upper * optical_bounds[1],
               roots[2] * optical_bounds[0] + (2 * roots[1] * upper + roots[0] * first) * optical_bounds[1]
               + roots[0] * upper ** 2 * optical_bounds[2],
               roots[3] * optical_bounds[0] + (3 * roots[2] * upper + 3 * roots[1] * first + roots[0] * second) * optical_bounds[1]
               + (3 * roots[1] * upper ** 2 + 3 * roots[0] * upper * first) * optical_bounds[2]
               + roots[0] * upper ** 3 * optical_bounds[3]]
    speed = [1 / lower, first / lower ** 2,
             2 * first ** 2 / lower ** 3 + second / lower ** 2,
             6 * first ** 3 / lower ** 4 + 6 * first * second / lower ** 3 + third / lower ** 2]
    flux = speed[0] * spatial[3] + 1.5 * speed[1] * spatial[2] + 0.75 * speed[2] * spatial[1] + 0.125 * speed[3] * spatial[0]
    return {"optical_L2_derivative_bounds_0_to_4": sobolev,
            "optical_sup_derivative_bounds_0_to_3": optical_bounds,
            "spatial_spinor_sup_bounds_0_to_3": spatial,
            "flux_third_derivative_bound": flux,
            "maximum_inverse_speed": upper}


def reference_tail_bound(duration, cutoff=DEFAULT_CUTOFF, mass=DEFAULT_MASS, background=BACKGROUND):
    distance = cutoff - 2
    if distance < 1:
        raise ValueError("Leave at least one mode beyond the initial support.")
    strength = abs(mass * background.lapse_modulation)
    argument = strength * abs(duration)
    tail = math.exp(argument) * argument ** distance / math.factorial(distance)
    return {"boundary_mode_norm_upper_bound": tail,
            "continuum_L2_reference_error_upper_bound": abs(duration) * strength * tail / 2}


def evolution_error_budget(sites, duration, initial_norm, cutoff=DEFAULT_CUTOFF, mass=DEFAULT_MASS, background=BACKGROUND):
    _, spacing = grid(sites, background)
    bounds = regularity_bounds(duration, mass, background)
    tail = reference_tail_bound(duration, cutoff, mass, background)
    strength = abs(mass * background.lapse_modulation)
    spatial_residual = math.sqrt(2 * background.length) * spacing ** 2 * bounds["flux_third_derivative_bound"] / 6
    reference_residual = strength * math.sqrt(background.length * bounds["maximum_inverse_speed"] / background.optical_length) * tail["boundary_mode_norm_upper_bound"]
    vector_bound = abs(duration) * (spatial_residual + reference_residual) / initial_norm
    return {"lattice_vs_sampled_reference_vector_bound": vector_bound,
            "all_same_final_POVM_TV_bound": min(1., 2 * vector_bound),
            "reference_truncation": tail, "regularity": bounds}


def smooth_effect(positions, center_fraction=0.35, background=BACKGROUND):
    return (1 + np.cos(background.frequency * (np.asarray(positions) - center_fraction * background.length))) / 2


def smooth_record_budget(sites, duration, mass=DEFAULT_MASS, background=BACKGROUND):
    _, spacing = grid(sites, background)
    amplitude, first, second, _ = regularity_bounds(duration, mass, background)["spatial_spinor_sup_bounds_0_to_3"]
    def defect(effect_0, effect_1, effect_2):
        density_second = (effect_2 * amplitude ** 2 + 4 * effect_1 * amplitude * first
                          + 2 * effect_0 * (first ** 2 + amplitude * second))
        interference_second = (effect_2 * amplitude ** 2 / 2 + 2 * effect_1 * amplitude * first
                               + effect_0 * (amplitude * second + first ** 2))
        return background.length * spacing ** 2 * (density_second / 24 + interference_second / 2)
    normalization = defect(1., 0., 0.)
    measurement = defect(1., background.frequency / 2, background.frequency ** 2 / 2)
    return {"unnormalized_sample_mass_error_bound": normalization,
            "unnormalized_smooth_record_error_bound": measurement,
            "normalized_smooth_record_error_bound": min(1., (measurement + normalization) / (1 - normalization))
            if normalization < 1 else 1.}


def smooth_continuum_probability(modes, center_fraction=0.35, points=4096, background=BACKGROUND):
    positions = np.arange(points) * background.length / points
    spinor = optical_spinor(positions, modes, background=background) / np.sqrt(background.speed(positions))[:, None]
    density = np.sum(abs(spinor) ** 2, axis=1)
    return float(background.length / points * np.sum(smooth_effect(positions, center_fraction, background) * density))


def comparison(sites, duration=2., mass=DEFAULT_MASS, direction=1, cutoff=DEFAULT_CUTOFF, background=BACKGROUND, center_fraction=0.35):
    start_modes = initial_modes(cutoff, direction, background)
    initial = reconstruct(start_modes, sites, background)
    normalization = float(np.linalg.norm(initial))
    exact = evolve_lattice(initial / normalization, duration, mass, background)
    modes = reference_modes(duration, cutoff, mass, direction, background)
    sampled = reconstruct(modes, sites, background)
    reference = sampled / normalization
    budget = evolution_error_budget(sites, duration, normalization, cutoff, mass, background)
    record_budget = smooth_record_budget(sites, duration, mass, background)
    positions, _ = grid(sites, background)
    discrete_probability = float(np.sum(smooth_effect(positions, center_fraction, background) * abs(exact) ** 2))
    continuum_probability = smooth_continuum_probability(modes, center_fraction, background=background)
    continuum_coarse = smooth_continuum_probability(modes, center_fraction, points=2048, background=background)
    amplitude, first, second, _ = budget["regularity"]["spatial_spinor_sup_bounds_0_to_3"]
    weighted_density_second = (background.frequency ** 2 * amplitude ** 2 / 2
                               + 2 * background.frequency * amplitude * first
                               + 2 * (first ** 2 + amplitude * second))
    integration_bound = background.length * (background.length / 4096) ** 2 * weighted_density_second / 24
    final_bound = min(1., budget["all_same_final_POVM_TV_bound"]
                      + record_budget["normalized_smooth_record_error_bound"]
                      + budget["reference_truncation"]["continuum_L2_reference_error_upper_bound"] + integration_bound)
    return {"sites": sites, "spacing": background.length / sites, "duration": duration, "proper_mass": mass,
            "initial_direction": direction, "cutoff": cutoff, "smooth_window_center_fraction": center_fraction,
            "exact_norm": float(np.linalg.norm(exact)), "initial_sample_norm": normalization,
            "sampled_reference_norm": float(np.linalg.norm(reference)),
            "vector_error_diagnostic": float(np.linalg.norm(exact - reference)),
            "vector_error_upper_bound": budget["lattice_vs_sampled_reference_vector_bound"],
            "continuum_reference_error_upper_bound": budget["reference_truncation"]["continuum_L2_reference_error_upper_bound"],
            "opposite_direction_probability_Galerkin": float(np.sum(abs(modes[:, 1 if direction == 1 else 0]) ** 2)),
            "smooth_record_probability_lattice": discrete_probability,
            "smooth_record_probability_continuum_reference": continuum_probability,
            "smooth_record_error_diagnostic": abs(discrete_probability - continuum_probability),
            "smooth_record_error_upper_bound": final_bound,
            "record_sampling_budget": record_budget,
            "continuum_probability_quadrature_upper_bound": integration_bound,
            "quadrature_halving_difference_diagnostic": abs(continuum_probability - continuum_coarse),
            "energy_drift": float(abs(np.vdot(exact, massive_hopping(sites, mass, background) @ exact)
                                      - np.vdot(initial / normalization, massive_hopping(sites, mass, background) @ (initial / normalization))))}


def postprocessing_table(sites, random_bits=12, center_fraction=0.5, background=BACKGROUND):
    if not isinstance(random_bits, int) or not 1 <= random_bits <= 24:
        raise ValueError("Use a finite random-bit budget between 1 and 24.")
    positions, _ = grid(sites, background)
    denominator = 2 ** random_bits
    thresholds = np.rint(denominator * smooth_effect(positions, center_fraction, background)).astype(np.int64)
    return thresholds, denominator


def classify_collected_record(bits, receipt_time, schedule, random_integer, thresholds, denominator):
    if (len(bits) != len(thresholds) or any(bit not in (0, 1) for bit in bits) or sum(bits) != 1
            or not isinstance(random_integer, int) or not 0 <= random_integer < denominator):
        raise ValueError("Use a complete one-excitation record and a valid independent random draw.")
    if not np.isfinite(receipt_time) or receipt_time < max(event["received_at_site_zero"] for event in schedule):
        raise ValueError("Wait until all position and clock records have arrived.")
    return 1 if random_integer < int(thresholds[bits.index(1)]) else -1


def terminal_task(sites=512, duration=2., mass=DEFAULT_MASS, background=BACKGROUND):
    thresholds, denominator = postprocessing_table(sites, background=background)
    rows, probabilities = [], []
    for direction in (1, -1):
        row = comparison(sites, duration, mass, direction, background=background, center_fraction=0.5)
        initial = reconstruct(initial_modes(direction=direction, background=background), sites, background)
        evolved = evolve_lattice(initial / np.linalg.norm(initial), duration, mass, background)
        probabilities.append(float(np.sum(thresholds / denominator * abs(evolved) ** 2)))
        rows.append(row)
    quantization = 1 / (2 * denominator)
    success = (probabilities[0] + 1 - probabilities[1]) / 2
    predicted = (rows[0]["smooth_record_probability_continuum_reference"]
                 + 1 - rows[1]["smooth_record_probability_continuum_reference"]) / 2
    bound = min(1., sum(row["smooth_record_error_upper_bound"] for row in rows) / 2 + quantization)
    resources = final_record_resources(sites, duration, background)
    resources.update({"staggered_onsite_controls": sites, "fresh_independent_fair_postprocessing_bits": 12,
                      "threshold_table_storage_bits_excluding_metadata": 13 * sites,
                      "smooth_output_bits": 1,
                      "all_fine_record_TV_against_continuum_spinor_density_is_bounded": False})
    schedule = collection_schedule(sites, duration, background)
    return {"sites": sites, "duration": duration, "rows": rows,
            "coarse_record_plus_probabilities_lattice": probabilities,
            "coarse_record_plus_probabilities_reference": [row["smooth_record_probability_continuum_reference"] for row in rows],
            "fair_bit_denominator": denominator, "quantization_error_bound": quantization,
            "decision_rule": "plus output guesses initial + carrier, minus output guesses initial - carrier; priors 1/2",
            "exact_expected_classification_success": success, "predicted_classification_success": predicted,
            "classification_success_difference": abs(success - predicted),
            "classification_success_error_upper_bound": bound,
            "conservative_success_lower_bound": max(0., predicted - bound),
            "clock_joint_probabilities": final_clock_probabilities(duration, background).tolist(),
            "resources": resources, "record_collection_schedule": schedule,
            "decision_not_allowed_before": max(event["received_at_site_zero"] for event in schedule),
            "raw_position_record_kept_but_bound_applies_to_released_smooth_bit_and_clock_bits_only": True,
            "sample_counts_generated": 0}


class MassiveEnvelopeTests(unittest.TestCase):
    def test_staggered_onsite_term_exchanges_carriers_exactly(self):
        sites = 32
        positions, _ = grid(sites)
        mass_profile = DEFAULT_MASS * BACKGROUND.lapse(positions)
        for direction in (-1, 1):
            np.testing.assert_array_equal((-1.) ** np.arange(sites) * carrier(sites, direction), carrier(sites, -direction))
            np.testing.assert_array_equal((-1.) ** np.arange(sites) * mass_profile * carrier(sites, direction), mass_profile * carrier(sites, -direction))

    def test_full_spin_chain_restriction_includes_the_declared_local_mass(self):
        sites = 8
        positions, spacing = grid(sites)
        full = np.zeros((2 ** sites, 2 ** sites), dtype=complex)
        for index in range(sites):
            full += BACKGROUND.speed(positions[index] + spacing / 2) / (2 * spacing) * embed_operator(EXCHANGE, (index, (index + 1) % sites), sites)
            full += (-1) ** index * DEFAULT_MASS * BACKGROUND.lapse(positions[index]) * embed_operator(binary_effect(1), (index,), sites)
        sector = [1 << (sites - index - 1) for index in range(sites)]
        np.testing.assert_allclose(full[np.ix_(sector, sector)], massive_hopping(sites), rtol=0, atol=1e-14)
        outside = sorted(set(range(2 ** sites)) - set(sector))
        np.testing.assert_array_equal(full[np.ix_(outside, sector)], np.zeros((len(outside), sites)))

    def test_constant_background_plane_wave_has_exact_massive_two_band_dispersion(self):
        background = MassiveBackground(modulation=0., lapse_modulation=0.)
        sites = 32
        positions, spacing = grid(sites, background)
        wave_number = 2 * math.pi / background.length
        vectors = np.column_stack([carrier(sites, direction) * np.exp(1j * wave_number * positions) / math.sqrt(sites)
                                   for direction in (1, -1)])
        kinetic = background.speed_scale * math.sin(wave_number * spacing) / spacing
        reduced = kinetic * PAULI_Z + DEFAULT_MASS * PAULI_X
        np.testing.assert_allclose(massive_hopping(sites, DEFAULT_MASS, background) @ vectors, vectors @ reduced, rtol=0, atol=2e-14)
        np.testing.assert_allclose(np.linalg.eigvalsh(reduced), [-math.hypot(kinetic, DEFAULT_MASS), math.hypot(kinetic, DEFAULT_MASS)], atol=1e-14)

    def test_zero_mass_reference_reduces_to_C4_optical_solution(self):
        sites = 64
        for direction in (-1, 1):
            modes = reference_modes(0.7, mass=0., direction=direction)
            np.testing.assert_allclose(reconstruct(modes, sites), sampled_solution(sites, 0.7, direction, BACKGROUND), rtol=0, atol=2e-14)

    def test_reference_is_hermitian_norm_preserving_and_changes_carrier_content(self):
        matrix = galerkin_hamiltonian()
        np.testing.assert_array_equal(matrix, matrix.conj().T)
        modes = reference_modes(2.)
        self.assertAlmostEqual(np.linalg.norm(modes), 1)
        self.assertGreater(np.sum(abs(modes[:, 1]) ** 2), 0.1)

    def test_regularities_cover_reference_derivatives_at_all_checked_times(self):
        horizon = 2.
        bounds = regularity_bounds(horizon)
        waves = 2 * math.pi * np.arange(-DEFAULT_CUTOFF, DEFAULT_CUTOFF + 1) / BACKGROUND.optical_length
        for time in (0., 0.5, 1., 2.):
            modes = reference_modes(time)
            for order in range(5):
                norm = np.linalg.norm((1j * waves)[:, None] ** order * modes)
                self.assertLessEqual(norm, bounds["optical_L2_derivative_bounds_0_to_4"][order] + 1e-12)

    def test_galerkin_boundary_defect_is_separately_bounded(self):
        cutoff = 6
        bounds = reference_tail_bound(2., cutoff)
        for time in (0.5, 1., 2.):
            modes = reference_modes(time, cutoff)
            self.assertLessEqual(np.linalg.norm(modes[[0, -1]]), bounds["boundary_mode_norm_upper_bound"] + 1e-14)
        coarse = reference_modes(2., cutoff)
        fine = reference_modes(2., cutoff + 4)
        embedded = np.zeros_like(fine)
        embedded[4:-4] = coarse
        self.assertLessEqual(np.linalg.norm(embedded - fine),
                             bounds["continuum_L2_reference_error_upper_bound"] + reference_tail_bound(2., cutoff + 4)["continuum_L2_reference_error_upper_bound"] + 1e-13)

    def test_combined_spatial_and_reference_residual_contains_the_mass_term(self):
        horizon, cutoff = 2., 6
        matrix = galerkin_hamiltonian(cutoff)
        for sites in (32, 64):
            residual_bound = evolution_error_budget(sites, horizon, 1., cutoff)["lattice_vs_sampled_reference_vector_bound"] / horizon
            for duration in (0.3, 1., 2.):
                modes = reference_modes(duration, cutoff)
                time_derivative = (-1j * matrix @ modes.reshape(-1)).reshape(-1, 2)
                defect = reconstruct(time_derivative, sites) + 1j * massive_hopping(sites) @ reconstruct(modes, sites)
                self.assertLessEqual(np.linalg.norm(defect), residual_bound + 1e-12)

    def test_variable_mass_chain_converges_with_finite_time_vector_and_record_bounds(self):
        rows = [comparison(sites) for sites in (32, 64, 128)]
        for row in rows:
            self.assertAlmostEqual(row["exact_norm"], 1)
            self.assertLessEqual(row["vector_error_diagnostic"], row["vector_error_upper_bound"])
            self.assertLessEqual(row["smooth_record_error_diagnostic"], row["smooth_record_error_upper_bound"])
            self.assertLess(row["energy_drift"], 1e-12)
            self.assertLessEqual(row["quadrature_halving_difference_diagnostic"], 5 * row["continuum_probability_quadrature_upper_bound"])
        self.assertTrue(all(3.2 < rows[index]["vector_error_diagnostic"] / rows[index + 1]["vector_error_diagnostic"] < 4.8 for index in (0, 1)))

    def test_mass_and_clock_background_share_the_same_declared_lapse(self):
        positions, _ = grid(32)
        profile = DEFAULT_MASS * BACKGROUND.lapse(positions)
        np.testing.assert_allclose(profile / DEFAULT_MASS, BACKGROUND.lapse(positions), atol=1e-15)
        np.testing.assert_allclose(BACKGROUND.lapse(positions) / BACKGROUND.spatial_scale(positions), BACKGROUND.speed(positions), atol=1e-15)
        for bare in (1., 1.2):
            rate = bare * BACKGROUND.lapse(0.)
            gate = clock_gate(rate, 0.4)
            plus = np.ones((2, 2)) / 2
            self.assertAlmostEqual(np.trace(clock_effect("Y", 0) @ gate @ plus @ gate.conj().T).real, (1 + math.sin(rate * 0.4)) / 2)

    def test_finite_random_bit_postprocessing_and_delayed_record_rule(self):
        sites = 32
        thresholds, denominator = postprocessing_table(sites)
        positions, _ = grid(sites)
        self.assertLessEqual(np.max(abs(thresholds / denominator - smooth_effect(positions, 0.5))), 1 / (2 * denominator) + 1e-15)
        schedule = collection_schedule(sites, 2.)
        deadline = max(event["received_at_site_zero"] for event in schedule)
        for site in (0, 7, 16, 31):
            bits = [int(index == site) for index in range(sites)]
            with self.assertRaises(ValueError):
                classify_collected_record(bits, deadline - 1e-6, schedule, 0, thresholds, denominator)
            plus_count = sum(classify_collected_record(bits, deadline, schedule, draw, thresholds, denominator) == 1
                             for draw in range(denominator))
            self.assertEqual(plus_count, thresholds[site])

    def test_smooth_terminal_task_has_a_nontrivial_sufficient_error_and_resource_certificate(self):
        task = terminal_task()
        self.assertLessEqual(task["classification_success_difference"], task["classification_success_error_upper_bound"])
        self.assertGreater(task["conservative_success_lower_bound"], 0.5)
        self.assertAlmostEqual(np.sum(task["clock_joint_probabilities"]), 1)
        self.assertEqual(task["resources"]["total_quantum_slots_before_detectors"], 514)
        self.assertEqual(task["resources"]["classical_bit_hops_including_clock_bits"], 65792)
        self.assertAlmostEqual(task["decision_not_allowed_before"], 2572.)

    def test_massive_lattice_unitary_matches_the_inherited_common_J_real_representation(self):
        sites = 16
        values, vectors = lattice_spectrum(sites)
        unitary = (vectors * np.exp(-0.7j * values)) @ vectors.T
        initial = reconstruct(initial_modes(), sites)
        initial /= np.linalg.norm(initial)
        density = np.outer(initial, initial.conj())
        lifted = real_lift(unitary)
        np.testing.assert_allclose(lifted @ encode_state(density) @ lifted.T,
                                   encode_state(unitary @ density @ unitary.conj().T), rtol=0, atol=2e-14)

    def test_invalid_parameters_fail_without_changing_old_models(self):
        with self.assertRaises(ValueError):
            massive_hopping(32, -1)
        with self.assertRaises(ValueError):
            galerkin_hamiltonian(2)
        with self.assertRaises(ValueError):
            initial_modes(direction=0)


def report():
    rows = [comparison(sites, direction=direction) for direction in (1, -1) for sites in (32, 64, 128, 256)]
    return {"model_version": "C5", "route_round": "05", "scope": "sufficiency only",
            "background_inputs": vars(BACKGROUND),
            "lapse_change_from_C4": "ell(x)=1+0.15*cos(2*pi*y(x)/Y); same C4 speed, b=ell/v",
            "declared_proper_mass": DEFAULT_MASS,
            "new_local_term": "sum_j (-1)^j*m*ell(x_j)*n_j",
            "derived_effective_equation": "i partial_t chi = -i sigma_z(v partial_x+v_prime/2)chi + m ell(x) sigma_x chi",
            "optical_reference": "i partial_t F=(-i sigma_z partial_y + m(1+0.15cos(2pi y/Y))*sigma_x)F",
            "galerkin_cutoff": DEFAULT_CUTOFF,
            "initial_modes": [-2, -1, 0, 1, 2],
            "regularity_bounds_T2": regularity_bounds(2.),
            "reference_tail_T2": reference_tail_bound(2.),
            "rows": rows,
            "terminal_record_task": terminal_task(),
            "smooth_record_interpretation": "terminal position read then declared bounded cosine-window random postprocessing; not a bound for raw microscopic position records versus spinor density",
            "geometry_interpretation": "massive single-particle density-weighted Dirac transport on a declared static 1+1 metric",
            "proper_mass_lapse_coupling_is_derived_from_cognition": False,
            "C3_count_uncertainty_propagated_into_C5": False,
            "full_quantum_field_or_background_Einstein_dynamics_derived": False,
            "float_values_are_directed_rounding_certificates": False,
            "all_grid_states_and_unbounded_times_covered": False,
            "necessity_arguments_performed": False,
            "sample_counts_generated": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(MassiveEnvelopeTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    results = report()
    results["automated_checks"] = {"run": tests.testsRun, "failures": 0, "errors": 0}
    if args.write_results:
        Path(__file__).with_name("massive_envelope_results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"tests": tests.testsRun, "reference_tail": results["reference_tail_T2"], "rows": results["rows"]}, indent=2))


if __name__ == "__main__":
    main()