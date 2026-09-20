"""Construction 08: a declared gradient oscillator field and controlled limits.

The C7 number-conserving background exchange is explicitly replaced. Quantum
rules, the new generator, its mass and speed remain inputs, not derived GR.
"""

import argparse
import json
import math
import unittest
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from clock_calibration import clock_effect
from massive_envelope import BACKGROUND, DEFAULT_MASS
from quantum_chain_background import oscillator_operators


@dataclass(frozen=True)
class FieldParameters:
    length: float = 20.0
    speed: float = 0.6
    gap: float = 0.8
    coupling: float = 0.15

    def __post_init__(self):
        if (not all(np.isfinite(value) for value in (self.length, self.speed, self.gap, self.coupling))
                or self.length <= 0 or self.speed <= 0 or self.gap <= 0 or not 0 <= self.coupling < 2):
            raise ValueError("Use finite positive length, speed and gap, and 0<=coupling<2.")


PARAMETERS = FieldParameters()


def field_grid(sites, parameters=PARAMETERS):
    if not isinstance(sites, int) or sites < 8 or sites % 4:
        raise ValueError("Use at least eight sites, with the inherited four-site carrier periodicity.")
    spacing = parameters.length / sites
    return np.arange(sites) * spacing, spacing


def stiffness(sites, parameters=PARAMETERS):
    _, spacing = field_grid(sites, parameters)
    identity = np.eye(sites)
    laplacian = 2 * identity - np.roll(identity, 1, axis=0) - np.roll(identity, -1, axis=0)
    return parameters.gap ** 2 * identity + parameters.speed ** 2 / spacing ** 2 * laplacian


def frequencies(sites, modes, parameters=PARAMETERS):
    _, spacing = field_grid(sites, parameters)
    wave_numbers = 2 * np.pi * np.asarray(modes) / parameters.length
    lattice = np.sqrt(parameters.gap ** 2 + 4 * parameters.speed ** 2 / spacing ** 2 * np.sin(wave_numbers * spacing / 2) ** 2)
    continuum = np.sqrt(parameters.gap ** 2 + parameters.speed ** 2 * wave_numbers ** 2)
    error = parameters.speed ** 2 * spacing ** 2 * wave_numbers ** 4 / (24 * parameters.gap)
    return lattice, continuum, error


def harmonic_flow(matrix, duration):
    if not np.isfinite(duration) or matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Use a square stiffness and finite duration.")
    if not np.allclose(matrix, matrix.T, rtol=0, atol=1e-12):
        raise ValueError("The stiffness must be symmetric.")
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    if eigenvalues.min() <= 0:
        raise ValueError("The declared Gaussian ground state requires positive stiffness.")
    frequency = np.sqrt(eigenvalues)
    cosine = (eigenvectors * np.cos(frequency * duration)) @ eigenvectors.T
    sine_over = (eigenvectors * (np.sin(frequency * duration) / frequency)) @ eigenvectors.T
    minus_sine = (eigenvectors * (-frequency * np.sin(frequency * duration))) @ eigenvectors.T
    return np.block([[cosine, sine_over], [minus_sine, cosine]])


def positive_root(matrix):
    values, vectors = np.linalg.eigh(matrix)
    if values.min() <= 0:
        raise ValueError("Positive stiffness is required.")
    return (vectors * np.sqrt(values)) @ vectors.T


def ground_covariance(matrix):
    root = positive_root(matrix)
    zeros = np.zeros_like(root)
    return np.block([[np.linalg.inv(root) / 2, zeros], [zeros, root / 2]])


def gaussian_branch(matrix, initial_width, duration, constant=0.):
    values, vectors = np.linalg.eigh(matrix)
    if values.min() <= 0 or not np.isfinite(duration):
        raise ValueError("Use positive stiffness and finite time.")
    frequency = np.sqrt(values)
    root = (vectors * frequency) @ vectors.T
    relative = float(np.linalg.norm(initial_width - root, ord=2) / frequency.min())
    if relative >= 1:
        raise ValueError("The branch-continuous determinant formula requires the declared relative-width bound below one.")
    cosine = (vectors * np.cos(frequency * duration)) @ vectors.T
    sine_over = (vectors * (np.sin(frequency * duration) / frequency)) @ vectors.T
    denominator = cosine + 1j * sine_over @ initial_width
    numerator = cosine @ initial_width + 1j * matrix @ sine_over
    width = np.linalg.solve(denominator.T, numerator.T).T
    correction = (vectors * (1j * np.sin(frequency * duration) * np.exp(-1j * frequency * duration) / frequency)) @ vectors.T @ (initial_width - root)
    log_determinant = 1j * duration * frequency.sum() + np.log(1 + np.linalg.eigvals(correction)).sum()
    log_factor = -0.5 * log_determinant - 1j * constant * duration
    return {"width": width, "log_factor": log_factor, "relative_width_bound": relative}


def gaussian_overlap(bra, ket, initial_width):
    combined = bra["width"].conj() + ket["width"]
    eigenvalues = np.linalg.eigvals(combined)
    if eigenvalues.real.min() <= 0:
        raise ValueError("The Gaussian overlap must have an integrable positive real quadratic form.")
    log_value = (bra["log_factor"].conj() + ket["log_factor"]
                 + np.linalg.slogdet(initial_width)[1] / 2 + len(initial_width) * math.log(2) / 2
                 - np.log(eigenvalues).sum() / 2)
    return np.exp(log_value)


class HeldSourceField:
    def __init__(self, sites=8, parameters=PARAMETERS, link_scale=1.):
        positions, spacing = field_grid(sites, parameters)
        if parameters.length != BACKGROUND.length or not 0 <= link_scale <= 1:
            raise ValueError("Use the inherited physical ring length and a declared link scale between zero and one.")
        self.sites, self.parameters = sites, parameters
        self.base = stiffness(sites, parameters)
        self.width = positive_root(self.base)
        self.covariance = ground_covariance(self.base)
        self.link_scale = link_scale
        self.background_matrix = parameters.gap ** 2 * np.eye(sites) + link_scale * (self.base - parameters.gap ** 2 * np.eye(sites))
        self.mass_values = (-1.) ** np.arange(sites) * DEFAULT_MASS * BACKGROUND.lapse(positions)
        self.clock_rates = np.array([BACKGROUND.lapse(0.), 1.2 * BACKGROUND.lapse(spacing)])

    def sectors(self, duration):
        sectors = []
        for occupied in (0, 2):
            for first, second in ((0, 0), (0, 1), (1, 0), (1, 1)):
                signs = np.array([1 - 2 * first, 1 - 2 * second])
                source = np.zeros(self.sites)
                source[occupied] = self.mass_values[occupied]
                source[:2] += signs * self.clock_rates / 2
                matrix = self.background_matrix + 2 * self.parameters.coupling * np.diag(source)
                constant = float(self.mass_values[occupied] + signs @ self.clock_rates / 2
                                 - self.parameters.coupling * source.sum() / 2)
                branch = gaussian_branch(matrix, self.width, duration, constant)
                branch.update({"matrix": matrix, "source": source, "constant": constant,
                               "occupied": occupied, "clock_bits": (first, second)})
                sectors.append(branch)
        return sectors

    def reduced_matter(self, data_density, duration):
        data_density = np.asarray(data_density, dtype=complex)
        if (data_density.shape != (2, 2) or not np.allclose(data_density, data_density.conj().T, atol=1e-12)
                or abs(np.trace(data_density) - 1) > 1e-12 or np.linalg.eigvalsh(data_density).min() < -1e-12):
            raise ValueError("Use a normalized positive preparation on the two declared data locations.")
        initial = np.kron(data_density, np.ones((4, 4)) / 4)
        branches = self.sectors(duration)
        overlap = np.array([[gaussian_overlap(bra, ket, self.width) for bra in branches] for ket in branches])
        return initial * overlap

    def clock_record(self, occupied, duration):
        if occupied not in (0, 2):
            raise ValueError("Source records use the declared fresh preparations at 0 or 2.")
        data = np.diag([int(occupied == 0), int(occupied == 2)])
        density = self.reduced_matter(data, duration)
        probabilities = np.zeros((2, 2, 2))
        for label, first, second in np.ndindex(probabilities.shape):
            effect = np.kron(np.diag([int(label == 0), int(label == 1)]),
                             np.kron(clock_effect("X", first), clock_effect("Y", second)))
            probabilities[label, first, second] = np.trace(density @ effect).real
        if probabilities.min() < -1e-12 or abs(probabilities.sum() - 1) > 1e-12:
            raise ValueError("The complete data and clock record law is not normalized and positive.")
        return probabilities

    def joint_density_on_span(self, data_density, duration):
        branches = self.sectors(duration)
        gram = np.array([[gaussian_overlap(bra, ket, self.width) for ket in branches] for bra in branches])
        values, vectors = np.linalg.eigh(gram)
        if values.min() < -1e-12:
            raise ValueError("Gaussian branch Gram matrix must be positive.")
        coordinates = np.sqrt(np.maximum(values, 0))[:, None] * vectors.conj().T
        inclusion = np.zeros((len(branches) ** 2, len(branches)), dtype=complex)
        for index in range(len(branches)):
            inclusion[index * len(branches):(index + 1) * len(branches), index] = coordinates[:, index]
        initial = np.kron(data_density, np.ones((4, 4)) / 4)
        return inclusion @ initial @ inclusion.conj().T

    def energy_report(self, occupied, duration):
        rows = []
        for branch in self.sectors(duration):
            if branch["occupied"] != occupied:
                continue
            flow = harmonic_flow(branch["matrix"], duration)
            covariance = flow @ self.covariance @ flow.T
            position, momentum = covariance[:self.sites, :self.sites], covariance[self.sites:, self.sites:]
            matter = self.mass_values[occupied] + np.dot(1 - 2 * np.array(branch["clock_bits"]), self.clock_rates) / 2
            background = np.trace(momentum + self.background_matrix @ position) / 2
            interaction = self.parameters.coupling * np.dot(branch["source"], np.diag(position) - 0.5)
            rows.append([matter, background, interaction, matter + background + interaction])
        return dict(zip(("matter", "background", "interaction", "total"), np.mean(rows, axis=0).tolist()))


def normal_cdf(value):
    return 0.5 * math.erfc(-value / math.sqrt(2))


def smooth_modes(center_fraction, parameters=PARAMETERS):
    modes = np.arange(-2, 3)
    coefficients = np.exp(-0.5 * (modes / 1.2) ** 2) * np.exp(-2j * np.pi * modes * center_fraction)
    return modes, coefficients / np.linalg.norm(coefficients)


def smooth_record(sites, duration=1., bits=16, parameters=PARAMETERS, weight_bits=24, saturation_target=None):
    positions, spacing = field_grid(sites, parameters)
    if (not np.isfinite(duration) or duration < 0 or not isinstance(bits, int) or not 2 <= bits <= 32
            or not isinstance(weight_bits, int) or not 2 <= weight_bits <= 48
            or (saturation_target is not None and (not np.isfinite(saturation_target) or not 0 < saturation_target < 1))):
        raise ValueError("Use finite nonnegative time, 2..32 detector bits, 2..48 weight bits and a valid optional tail target.")
    modes, preparation = smooth_modes(0.35, parameters)
    _, weight = smooth_modes(0.5, parameters)
    preparation = 1.5 * preparation
    lattice, continuum, frequency_bound = frequencies(sites, modes, parameters)
    mean_lattice = float(np.vdot(weight, preparation * np.cos(lattice * duration)).real)
    mean_continuum = float(np.vdot(weight, preparation * np.cos(continuum * duration)).real)
    variance_lattice = float(np.sum(abs(weight) ** 2 / (2 * lattice)))
    variance_continuum = float(np.sum(abs(weight) ** 2 / (2 * continuum)))
    mean_bound = float(duration * np.sum(abs(weight.conj() * preparation) * frequency_bound))
    variance_bound = float(np.sum(abs(weight) ** 2 * frequency_bound / (2 * parameters.gap ** 2)))
    variance_floor = float(np.sum(abs(weight) ** 2) / (2 * continuum.max()))
    probability_bound = (mean_bound / math.sqrt(variance_floor)
                         + abs(mean_continuum) * variance_bound / (2 * variance_floor ** 1.5)) / math.sqrt(2 * math.pi)
    lattice_probability = normal_cdf(mean_lattice / math.sqrt(variance_lattice))
    continuum_probability = normal_cdf(mean_continuum / math.sqrt(variance_continuum))
    wave = np.exp(2j * np.pi * np.outer(positions / parameters.length, modes)) / math.sqrt(sites)
    local_mean = (wave @ (preparation * np.cos(lattice * duration))).real
    local_weights = (wave @ weight).real
    all_frequencies = frequencies(sites, np.arange(sites), parameters)[0]
    local_variance = float(np.mean(1 / (2 * all_frequencies)))
    tail_standard_deviations = 6. if saturation_target is None else math.sqrt(2 * math.log(2 * sites / saturation_target))
    radius = float(np.max(abs(local_mean)) + tail_standard_deviations * math.sqrt(local_variance))
    interior_bins = 2 ** bits - 2
    step = 2 * radius / interior_bins
    weight_denominator = 2 ** weight_bits
    weight_integers = np.rint(local_weights * weight_denominator).astype(int)
    weight_rounding = float(radius * np.sum(abs(local_weights - weight_integers / weight_denominator)))
    rounding_radius = float(step * np.sum(abs(local_weights)) / 2 + weight_rounding)
    saturation = float(sum(normal_cdf((-radius - mean) / math.sqrt(local_variance))
                           + normal_cdf((mean - radius) / math.sqrt(local_variance)) for mean in local_mean))
    rounding_probability = min(1., 2 * rounding_radius / math.sqrt(2 * math.pi * variance_lattice) + saturation)
    full_budget = min(1., probability_bound + rounding_probability)
    return {"sites": sites, "duration": duration, "bits_per_local_record": bits,
            "mean_lattice": mean_lattice, "mean_continuum": mean_continuum, "mean_error_bound": mean_bound,
            "variance_lattice": variance_lattice, "variance_continuum": variance_continuum, "variance_error_bound": variance_bound,
            "ideal_lattice_sign_probability": lattice_probability, "continuum_sign_probability": continuum_probability,
            "spatial_sign_probability_error_bound": min(1., probability_bound),
            "quadrature_radius": radius, "interior_bin_width": step, "rounding_radius": rounding_radius,
            "weight_rounding_radius": weight_rounding,
            "weight_fractional_bits": weight_bits, "declared_saturation_target": saturation_target,
            "saturation_union_bound": saturation, "finite_bin_probability_error_bound": rounding_probability,
            "complete_sign_probability_error_bound": full_budget,
            "finite_record_probability_interval": [max(0., continuum_probability - full_budget), min(1., continuum_probability + full_budget)],
            "raw_bits": bits * sites + 2, "serial_bit_hops": bits * sites ** 2 // 4 + 1,
            "all_records_ready_at": duration + spacing * (bits * sites ** 2 // 4 + 1),
            "local_weight_integer_table": weight_integers.tolist(), "weight_denominator": weight_denominator,
            "weight_table_bits": (weight_bits + 2) * sites, "integer_accumulator_bits": weight_bits + bits + math.ceil(math.log2(sites)) + 3,
            "chain_qubits_retained_unmeasured": sites, "background_modes": sites, "clock_qubits": 2,
            "ground_zero_point_energy": float(all_frequencies.sum() / 2),
            "smooth_displacement_energy": float(np.sum(lattice ** 2 * abs(preparation) ** 2) / 2),
            "scope": "free field, five fixed Fourier modes, own lattice/continuum ground states, local commuting finite-bin quadratures; float evaluation not directed rounding"}


def quadrature_bins(values, radius, bits):
    values = np.asarray(values)
    if not np.all(np.isfinite(values)) or not np.isfinite(radius) or radius <= 0 or not isinstance(bits, int) or not 2 <= bits <= 32:
        raise ValueError("Use finite quadratures, positive range and 2..32 bits.")
    count = 2 ** bits - 2
    interior = np.floor((values + radius) * count / (2 * radius)).astype(int)
    return np.where(values < -radius, 0, np.where(values >= radius, count + 1, 1 + np.clip(interior, 0, count - 1)))


def release_field_record(codes, clock_bits, receipt, report):
    codes = np.asarray(codes)
    count = 2 ** report["bits_per_local_record"] - 2
    if (codes.shape != (report["sites"],) or not np.issubdtype(codes.dtype, np.integer)
            or np.any(codes < 0) or np.any(codes > count + 1) or len(clock_bits) != 2
            or any(bit not in (0, 1) for bit in clock_bits)):
        raise ValueError("Use every finite local code and both clock bits, including saturation codes.")
    if not np.isfinite(receipt) or receipt < report["all_records_ready_at"]:
        raise ValueError("Wait for every serially routed field and clock bit.")
    decoded_units = np.where(codes == 0, -count, np.where(codes == count + 1, count, 2 * codes - 1 - count))
    signed_sum = sum(int(weight) * int(unit) for weight, unit in zip(report["local_weight_integer_table"], decoded_units))
    return int(signed_sum < 0), *clock_bits


def source_record_schedule(sites, duration, parameters=PARAMETERS):
    _, spacing = field_grid(sites, parameters)
    if not np.isfinite(duration) or duration < 0:
        raise ValueError("Use a finite nonnegative read parameter.")
    sources = [("position", site) for site in range(sites)] + [("clock", 0), ("clock", 1)]
    receipt, schedule = duration, []
    for kind, site in sources:
        hops = min(site, sites - site)
        receipt += spacing * hops
        schedule.append({"kind": kind, "source_site": site, "read_at": duration,
                         "hops": hops, "received_at_site_zero": receipt})
    return schedule


def release_source_record(raw_positions, clock_bits, receipt, duration):
    schedule = source_record_schedule(len(raw_positions), duration)
    if (any(bit not in (0, 1) for bit in raw_positions) or sum(raw_positions) != 1
            or len(clock_bits) != 2 or any(bit not in (0, 1) for bit in clock_bits)):
        raise ValueError("Retain every position bit and both clock results.")
    if not np.isfinite(receipt) or receipt < schedule[-1]["received_at_site_zero"]:
        raise ValueError("Wait for all raw source and clock records.")
    occupied = list(raw_positions).index(1)
    if occupied not in (0, 2):
        raise ValueError("This held-source protocol only prepares locations 0 or 2.")
    return int(occupied == 2), *clock_bits


def run_report():
    source_rows = []
    for link_scale in (0., 1.):
        current = HeldSourceField(link_scale=link_scale)
        records = [current.clock_record(occupied, 2.) for occupied in (0, 2)]
        probabilities = [float(record[..., 0].sum()) for record in records]
        branches = current.sectors(2.)
        source_rows.append({"link_scale": link_scale, "source_occupied_sites": [0, 2],
                            "clock_C_Y_plus_probabilities": probabilities,
                            "absolute_remote_record_difference": abs(probabilities[0] - probabilities[1]),
                            "joint_data_clock_probabilities": [record.tolist() for record in records],
                            "minimum_sector_stiffness_eigenvalue": min(float(np.linalg.eigvalsh(branch["matrix"]).min()) for branch in branches),
                            "maximum_relative_width_bound": max(branch["relative_width_bound"] for branch in branches),
                            "energy_source_at_zero": [dict(duration=time, **current.energy_report(0, time)) for time in (0., 1., 2.)]})
    free_rows = [smooth_record(sites) for sites in (16, 32, 64, 128)]
    scaled_rows = []
    for sites in (16, 32, 64, 128):
        refinement = math.log2(sites / 16)
        scaled_rows.append(smooth_record(sites, bits=16 + math.ceil(3 * refinement),
                           weight_bits=24 + math.ceil(4 * refinement), saturation_target=1e-8 * (16 / sites) ** 2))
    for row in free_rows:
        spacing = PARAMETERS.length / row["sites"]
        probabilities = []
        for location, bare, basis in ((0., 1., "X"), (spacing, 1.2, "Y")):
            phase = bare * BACKGROUND.lapse(location) * row["duration"]
            probability = (1 + (math.cos(phase) if basis == "X" else math.sin(phase))) / 2
            probabilities.append(np.array([probability, 1 - probability]))
        row["ideal_lattice_sign_clock_joint_probabilities"] = np.einsum("i,j,k->ijk",
            [row["ideal_lattice_sign_probability"], 1 - row["ideal_lattice_sign_probability"]], *probabilities).tolist()
    current = HeldSourceField()
    conservative_margin = PARAMETERS.gap ** 2 - 2 * PARAMETERS.coupling * (
        max(abs(current.mass_values)) + max(current.clock_rates) / 2)
    return {"model_version": "C8", "route_round": "08", "research_scope": "sufficiency only",
            "parameters": vars(PARAMETERS),
            "explicit_generator_change": "C7 number-conserving exchange replaced by one-half sum p_j squared + gap squared q_j squared + speed squared over spacing squared times nearest q difference squared",
            "field_equation": "q_ddot = speed^2 discrete_Laplacian q - gap^2 q - 2 coupling S q; closed c-number equation only in frozen source/clock-Z sectors",
            "canonical_scaling": "q_j=sqrt(spacing)*phi_j; p_j=sqrt(spacing)*pi_j",
            "held_source_conservative_stability_margin": float(conservative_margin),
            "free_continuum_records": free_rows, "scaled_precision_continuum_records": scaled_rows,
            "precision_scaling": "detector bits 16+ceil(3 log2(N/16)), weight fractional bits 24+ceil(4 log2(N/16)), saturation target 1e-8*(16/N)^2; finite report N=16..128",
            "held_source_clock_task": source_rows,
            "source_task_schedule": source_record_schedule(8, 2.),
            "source_task_resources": {"chain_qubits": 8, "clock_qubits": 2, "background_modes": 8,
                                      "background_Fock_cutoff": None, "conditional_Gaussian_branches": 8,
                                      "raw_data_clock_bits": 10, "classical_bit_hops": 17,
                                      "read_parameter": 2., "all_records_ready_at": 44.5,
                                      "background_retained_unmeasured": True, "resets": 0},
            "limits": ["new oscillator gradient law, gap, speed and quantum rules are working inputs, not a derivation of Einstein dynamics",
                       "free continuum record proof uses coupling zero and five fixed smooth modes; source task uses held data hopping, not active C7 data propagation",
                       "field modes now occupy every node; source site 2 is no longer uncoupled storage",
                       "same connected-field ground preparation used in both link settings, including its correlations and quench cost",
                       "Gaussian method is exact without Fock cutoff for this frozen-source quadratic slice, not general active quantum matter",
                       "free ground zero-point energy grows with refinement; no ultraviolet renormalization or state preparation cost closure",
                       "finite-bin record interval is a deterministic approximation budget, not a sampling confidence interval or directed-rounding certificate",
                       "no continuous sourced field convergence, dynamic metric constraints, Einstein equation, or lab data"]}


class BackgroundFieldLimitTests(unittest.TestCase):
    def test_periodic_stiffness_has_the_derived_Klein_Gordon_dispersion(self):
        sites = 16
        positions, _ = field_grid(sites)
        matrix = stiffness(sites)
        for mode in range(-7, 8):
            vector = np.exp(2j * np.pi * mode * positions / PARAMETERS.length) / math.sqrt(sites)
            frequency = frequencies(sites, [mode])[0][0]
            np.testing.assert_allclose(matrix @ vector, frequency ** 2 * vector, atol=3e-14)

    def test_symplectic_flow_is_reversible_and_preserves_quadratic_energy(self):
        sites = 8
        matrix = stiffness(sites)
        flow = harmonic_flow(matrix, 1.3)
        identity = np.eye(sites)
        symplectic = np.block([[np.zeros_like(identity), identity], [-identity, np.zeros_like(identity)]])
        energy = np.block([[matrix, np.zeros_like(identity)], [np.zeros_like(identity), identity]])
        np.testing.assert_allclose(flow.T @ symplectic @ flow, symplectic, atol=3e-14)
        np.testing.assert_allclose(flow.T @ energy @ flow, energy, atol=3e-14)
        np.testing.assert_allclose(harmonic_flow(matrix, -1.3) @ flow, np.eye(2 * sites), atol=3e-14)

    def test_fixed_low_modes_obey_explicit_second_order_frequency_bound(self):
        previous_error = None
        for sites in (16, 32, 64, 128):
            lattice, continuum, bound = frequencies(sites, np.arange(-2, 3))
            self.assertTrue(np.all(abs(lattice - continuum) <= bound + 1e-15))
            error = np.max(abs(lattice - continuum))
            if previous_error is not None:
                self.assertLess(error, previous_error / 3.8)
            previous_error = error

    def test_original_C7_exchange_is_not_the_new_field_generator(self):
        sites = 16
        phase = 2 * np.pi * np.arange(sites) / sites
        old_frequency = 1. + 0.6 * np.cos(phase)
        new_frequency = frequencies(sites, np.arange(sites))[0]
        self.assertGreater(np.max(abs(old_frequency - new_frequency)), 0.2)

    def test_gaussian_width_and_ground_phase_agree_with_exact_free_evolution(self):
        matrix = stiffness(8)
        width = positive_root(matrix)
        branch = gaussian_branch(matrix, width, 12.)
        np.testing.assert_allclose(branch["width"], width, atol=3e-14)
        self.assertAlmostEqual(abs(branch["log_factor"] + 6j * np.trace(width)), 0)
        self.assertAlmostEqual(gaussian_overlap(branch, branch, width).real, 1)

    def test_conditional_gaussians_retain_normalization_covariance_and_energy(self):
        current = HeldSourceField()
        for branch in current.sectors(2.):
            self.assertLess(branch["relative_width_bound"], 1)
            self.assertAlmostEqual(abs(gaussian_overlap(branch, branch, current.width) - 1), 0)
            flow = harmonic_flow(branch["matrix"], 2.)
            covariance = flow @ current.covariance @ flow.T
            np.testing.assert_allclose(np.linalg.inv(branch["width"].real) / 2, covariance[:8, :8], atol=3e-14)
        self.assertAlmostEqual(current.energy_report(0, 0.)["total"], current.energy_report(0, 2.)["total"])

    def test_joint_clock_records_are_positive_and_source_changes_remote_probability(self):
        current = HeldSourceField()
        records = [current.clock_record(occupied, 2.) for occupied in (0, 2)]
        for record in records:
            self.assertAlmostEqual(record.sum(), 1)
            self.assertGreaterEqual(record.min(), -1e-14)
        difference = abs(records[0][..., 0].sum() - records[1][..., 0].sum())
        self.assertGreater(difference, 1e-10)

    def test_switching_off_gradient_links_removes_source_dependent_remote_signal(self):
        current = HeldSourceField(link_scale=0.)
        probabilities = [current.clock_record(occupied, 2.)[..., 0].sum() for occupied in (0, 2)]
        self.assertAlmostEqual(probabilities[0], probabilities[1])

    def test_preparation_mixture_is_affine_with_all_clock_coherences(self):
        current = HeldSourceField()
        direct = sum(current.reduced_matter(np.diag([int(label == 0), int(label == 1)]), 1.) / 2 for label in (0, 1))
        alternate = sum(current.reduced_matter(np.outer([1, sign], [1, sign]) / 2, 1.) / 2 for sign in (1, -1))
        np.testing.assert_allclose(direct, alternate, atol=3e-14)
        self.assertGreaterEqual(np.linalg.eigvalsh(direct).min(), -1e-14)

    def test_smooth_quantum_record_has_spatial_and_finite_bin_budgets(self):
        for sites in (16, 32, 64, 128):
            report = smooth_record(sites)
            self.assertLessEqual(abs(report["mean_lattice"] - report["mean_continuum"]), report["mean_error_bound"])
            self.assertLessEqual(abs(report["variance_lattice"] - report["variance_continuum"]), report["variance_error_bound"])
            self.assertLessEqual(abs(report["ideal_lattice_sign_probability"] - report["continuum_sign_probability"]), report["spatial_sign_probability_error_bound"])
            self.assertLess(report["saturation_union_bound"], 1e-6)

    def test_finite_quadrature_bins_keep_saturation_and_wait_for_complete_records(self):
        report = smooth_record(8, bits=8)
        radius = report["quadrature_radius"]
        values = np.linspace(-1.2 * radius, 1.2 * radius, 8)
        codes = quadrature_bins(values, radius, 8)
        self.assertEqual(codes[0], 0)
        self.assertEqual(codes[-1], 255)
        with self.assertRaises(ValueError):
            release_field_record(codes, [0, 1], report["all_records_ready_at"] - 0.1, report)
        self.assertEqual(release_field_record(codes, [0, 1], report["all_records_ready_at"], report)[1:], (0, 1))

    def test_invalid_stiffness_and_record_inputs_are_rejected(self):
        with self.assertRaises(ValueError):
            harmonic_flow(-np.eye(8), 1.)
        with self.assertRaises(ValueError):
            gaussian_branch(np.eye(8), 3 * np.eye(8), 1.)
        with self.assertRaises(ValueError):
            smooth_record(8, bits=1)

    def test_gaussian_relative_phase_matches_independent_two_mode_Fock_evolution(self):
        cutoff = 14
        operators = oscillator_operators(cutoff)
        position, quadratic = operators["position"], operators["position_squared_projected"]
        momentum_squared = 2 * operators["number"] + np.eye(cutoff) - quadratic
        width = 0.8 * np.eye(2)
        coefficients = np.zeros(cutoff)
        for level in range(0, cutoff, 2):
            half = level // 2
            coefficients[level] = math.sqrt(2 * math.sqrt(0.8) / 1.8) * math.sqrt(math.factorial(level)) / (2 ** half * math.factorial(half)) * (0.2 / 1.8) ** half
        initial = np.kron(coefficients, coefficients)
        initial /= np.linalg.norm(initial)
        matrices = (np.array([[1.1, -0.12], [-0.12, 0.7]]), np.array([[0.75, -0.08], [-0.08, 0.9]]))
        constants = (0.31, -0.23)
        states, branches = [], []
        for matrix, constant in zip(matrices, constants):
            hamiltonian = (np.kron(momentum_squared + matrix[0, 0] * quadratic, np.eye(cutoff))
                           + np.kron(np.eye(cutoff), momentum_squared + matrix[1, 1] * quadratic)) / 2
            hamiltonian += matrix[0, 1] * np.kron(position, position) + constant * np.eye(cutoff ** 2)
            values, vectors = np.linalg.eigh(hamiltonian)
            states.append(vectors @ (np.exp(-2.3j * values) * (vectors.T @ initial)))
            branches.append(gaussian_branch(matrix, width, 2.3, constant))
        self.assertLess(abs(np.vdot(states[0], states[1]) - gaussian_overlap(branches[0], branches[1], width)), 2e-8)

    def test_full_joint_preparation_mixture_and_partial_trace_match(self):
        current = HeldSourceField()
        direct = sum(current.joint_density_on_span(np.diag([int(label == 0), int(label == 1)]), 1.) / 2 for label in (0, 1))
        alternate = sum(current.joint_density_on_span(np.outer([1, sign], [1, sign]) / 2, 1.) / 2 for sign in (1, -1))
        np.testing.assert_allclose(direct, alternate, atol=3e-14)
        reduced = np.trace(direct.reshape(8, 8, 8, 8), axis1=1, axis2=3)
        np.testing.assert_allclose(reduced, current.reduced_matter(np.eye(2) / 2, 1.), atol=3e-14)

    def test_smooth_mode_mean_and_variance_match_the_finite_quantum_covariance(self):
        sites = 16
        positions, _ = field_grid(sites)
        modes, preparation = smooth_modes(0.35)
        _, weight = smooth_modes(0.5)
        wave = np.exp(2j * np.pi * np.outer(positions / PARAMETERS.length, modes)) / math.sqrt(sites)
        initial = np.concatenate([(wave @ (1.5 * preparation)).real, np.zeros(sites)])
        flow = harmonic_flow(stiffness(sites), 1.)
        mean = flow @ initial
        covariance = flow @ ground_covariance(stiffness(sites)) @ flow.T
        observable = (wave @ weight).real
        report = smooth_record(sites)
        self.assertAlmostEqual(observable @ mean[:sites], report["mean_lattice"])
        self.assertAlmostEqual(observable @ covariance[:sites, :sites] @ observable, report["variance_lattice"])

    def test_coordinated_detector_precision_restores_second_order_record_convergence(self):
        previous = None
        for exponent in range(4):
            sites = 16 * 2 ** exponent
            target = 1e-8 * (16 / sites) ** 2
            report = smooth_record(sites, bits=16 + 3 * exponent, weight_bits=24 + 4 * exponent, saturation_target=target)
            self.assertLess(report["saturation_union_bound"], target)
            error = report["complete_sign_probability_error_bound"]
            if previous is not None:
                self.assertLess(error, previous / 3.5)
            previous = error

    def test_held_source_records_wait_for_all_position_and_clock_bits(self):
        positions = [0] * 8
        positions[2] = 1
        schedule = source_record_schedule(8, 2.)
        self.assertEqual(sum(event["hops"] for event in schedule), 17)
        self.assertEqual(schedule[-1]["received_at_site_zero"], 44.5)
        with self.assertRaises(ValueError):
            release_source_record(positions, [1, 0], 44.4, 2.)
        self.assertEqual(release_source_record(positions, [1, 0], 44.5, 2.), (1, 1, 0))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    outcome = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(BackgroundFieldLimitTests))
    if not outcome.wasSuccessful():
        raise SystemExit(1)
    if arguments.write_results:
        report = run_report()
        Path(__file__).with_name("background_field_limit_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"tests": outcome.testsRun, "free_record_64": report["free_continuum_records"][2],
                          "source_clock_task": report["held_source_clock_task"]}, indent=2))