"""Construction 13: finite radial-mode enlargement and omission budgets.

The first two overtones begin in vacuum. Bounds compare this three-mode
linear model with C12, not with the entire fluid or nonlinear gravity.
"""

import argparse
import json
import math
import unittest
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from radial_fluid_modes import RadialProblem, RadialMode, dynamic_probe_terms, record_schedule, release_record
from quantum_radial_backreaction import (canonical_inertia, Parameters as SingleParameters, model as single_model,
                                        total_record_budget)
from quantum_chain_background import oscillator_operators
from massive_envelope import postprocessing_table
from clock_calibration import clock_effect
from common_orientation_structure import encode_state
from baseline import real_lift


@lru_cache(maxsize=1)
def radial_modes():
    problem = RadialProblem()
    return tuple(RadialMode(problem, problem.find_mode(index)) for index in range(3))


def kinetic_gram(quadrature_order=256):
    modes = radial_modes()
    nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
    radii = modes[0].radius * (nodes + 1) / 2
    weights *= modes[0].radius / 2
    matrix = np.zeros((3, 3))
    for radius, weight in zip(radii, weights):
        row = modes[0].problem.equilibrium(float(radius))
        displacements = np.array([mode.local(float(radius))["displacement"] for mode in modes])
        matrix += weight * 4 * math.pi * radius ** 2 * (row["density"] + row["pressure"]) * row["radial_factor"] ** 1.5 / row["lapse"] * np.outer(displacements, displacements)
    return matrix


@lru_cache(maxsize=4)
def mode_data(sites=8, scale=1024.):
    if not np.isfinite(scale) or scale < 1:
        raise ValueError("Use a finite physical homology scale >=1.")
    rows = []
    for index, mode in enumerate(radial_modes()):
        inertia = canonical_inertia(mode) * scale ** 3
        frequency = mode.omega / scale
        sigma = 1 / math.sqrt(2 * inertia * frequency)
        terms = dynamic_probe_terms(mode, sites, 0.35, 0.01)
        clock_source = terms["clock_rates"] * terms["clock_lapse_fractional_derivatives"] / scale
        sources = []
        for first, second in ((0, 0), (0, 1), (1, 0), (1, 1)):
            signs = np.array([1 - 2 * first, 1 - 2 * second])
            sources.append(terms["mode_hamiltonian_derivative"] / scale + (signs @ clock_source / 2) * np.eye(sites))
        source_norm = max(float(np.linalg.norm(source, ord=2)) for source in sources)
        rows.append({"index": index, "mode": mode, "reference_frequency": mode.omega,
                     "physical_frequency": frequency, "physical_inertia": inertia, "sigma": sigma,
                     "terms": terms, "sources": tuple(sources), "source_norm": source_norm,
                     "vacuum_coupling_norm": sigma * source_norm,
                     "source_to_frequency_ratio": sigma * source_norm / frequency})
    return tuple(rows)


def omission_budget(retained=(0,), reference_duration=2., sites=8, scale=1024.):
    if tuple(retained) not in ((0,), (0, 1), (0, 1, 2)) or not np.isfinite(reference_duration) or reference_duration < 0:
        raise ValueError("Retain an initial prefix of the three modes and use a nonnegative time.")
    omitted = [row for row in mode_data(sites, scale) if row["index"] not in retained]
    coupling_squared = sum(row["vacuum_coupling_norm"] ** 2 for row in omitted)
    duration = reference_duration * scale
    vector_bound = duration * math.sqrt(coupling_squared)
    squared_time_coupling = vector_bound ** 2
    reduced_bound = math.expm1(squared_time_coupling / 2) + squared_time_coupling / 2 if squared_time_coupling < 1 else 1.
    return {"retained_modes": list(retained), "omitted_modes": [row["index"] for row in omitted],
            "physical_duration": duration, "vacuum_coupling_squared_sum": coupling_squared,
            "full_state_vector_error_upper_bound": vector_bound,
            "same_full_record_TV_upper_bound": min(1., vector_bound),
            "retained_record_TV_upper_bound_after_omitted_trace": min(1., reduced_bound),
            "scope": "finite omitted set initially independent vacuum; only omitted interactions are removed, free vacuum phases retained; not an infinite-mode tail bound"}


@dataclass(frozen=True)
class Parameters:
    cutoffs: tuple = (20, 4, 4)
    interaction_weights: tuple = (1., 1., 1.)
    sites: int = 8
    scale: float = 1024.

    def __post_init__(self):
        if (len(self.cutoffs) != 3 or any(not isinstance(value, int) or value < 2 for value in self.cutoffs)
                or len(self.interaction_weights) != 3
                or any(not np.isfinite(value) or not 0 <= value <= 1 for value in self.interaction_weights)
                or not isinstance(self.sites, int) or self.sites < 8 or self.sites % 4
                or not np.isfinite(self.scale) or self.scale < 1):
            raise ValueError("Use three finite cutoffs >=2, three coupling weights in [0,1], and a valid inherited ring and scale.")


def mode_operator(local, index, cutoffs):
    result = sparse.csr_matrix([[1.]])
    for axis, cutoff in enumerate(cutoffs):
        result = sparse.kron(result, sparse.csr_matrix(local) if axis == index else sparse.eye(cutoff), format="csr")
    return result


class ThreeRadialModes:
    def __init__(self, parameters=Parameters()):
        self.parameters = parameters
        self.rows = mode_data(parameters.sites, parameters.scale)
        self.shape = (parameters.sites, *parameters.cutoffs, 4)
        self.background_dimension = math.prod(parameters.cutoffs)
        self.positions, self.projected_squares, mode_terms = [], [], []
        for index, (row, cutoff) in enumerate(zip(self.rows, parameters.cutoffs)):
            operators = oscillator_operators(cutoff)
            self.positions.append(mode_operator(math.sqrt(2) * row["sigma"] * operators["position"], index, parameters.cutoffs))
            self.projected_squares.append(mode_operator(2 * row["sigma"] ** 2 * operators["position_squared_projected"], index, parameters.cutoffs))
            mode_terms.append(mode_operator(row["physical_frequency"] * (operators["number"] + np.eye(cutoff) / 2), index, parameters.cutoffs))
        fundamental = single_model(SingleParameters(sites=parameters.sites, scale=parameters.scale,
                                   cutoff=max(4, parameters.cutoffs[0])))
        self.initial_fundamental = fundamental.coherent_vector()[:parameters.cutoffs[0]].copy()
        self.initial_fundamental /= np.linalg.norm(self.initial_fundamental)
        self.coherent_initial_error = fundamental.coherent_initial_budget() if parameters.cutoffs[0] >= 4 else None
        self.blocks = []
        for sector, first_block in enumerate(fundamental.blocks):
            matter = sparse.kron(sparse.csr_matrix(first_block["matter_small"]), sparse.eye(self.background_dimension), format="csr")
            background = sparse.kron(sparse.eye(parameters.sites), sum(mode_terms), format="csr")
            interactions = [weight * sparse.kron(sparse.csr_matrix(row["sources"][sector]), position, format="csr")
                            for row, position, weight in zip(self.rows, self.positions, parameters.interaction_weights)]
            hamiltonian = matter + background + sum(interactions)
            generator = parameters.scale * hamiltonian
            self.blocks.append({"matter": matter, "background": background, "interaction": sum(interactions),
                                "mode_terms": [sparse.kron(sparse.eye(parameters.sites), term, format="csr") for term in mode_terms],
                                "hamiltonian": hamiltonian, "generator": generator,
                                "generator_trace": float(generator.diagonal().sum())})

    def initial_vector(self, data=None):
        from quantum_chain_background import initial_data
        data = initial_data(self.parameters.sites) if data is None else np.asarray(data, dtype=complex)
        if data.shape != (self.parameters.sites,) or not np.isclose(np.linalg.norm(data), 1., rtol=0, atol=1e-12):
            raise ValueError("Use normalized one-excitation data; overtones are explicitly fresh vacua.")
        state = np.zeros(self.shape, dtype=complex)
        state[:, :, 0, 0, :] = data[:, None, None] * self.initial_fundamental[None, :, None] / 2
        return state

    def evolve(self, initial, reference_duration):
        if initial.shape != self.shape or not np.isfinite(reference_duration):
            raise ValueError("Use the full three-mode joint vector and finite reference time t/scale.")
        output = np.empty(self.shape, dtype=complex)
        for sector, block in enumerate(self.blocks):
            vector = initial[..., sector].reshape(-1)
            output[..., sector] = expm_multiply(-1j * reference_duration * block["generator"], vector,
                                                 traceA=-1j * reference_duration * block["generator_trace"]).reshape(self.shape[:-1])
        return output

    def diagnostics(self, state):
        energies = {term: float(sum(np.vdot(state[..., sector].reshape(-1), block[term] @ state[..., sector].reshape(-1)).real
                                   for sector, block in enumerate(self.blocks)))
                    for term in ("matter", "background", "interaction", "hamiltonian")}
        mode_energies = [float(sum(np.vdot(state[..., sector].reshape(-1), block["mode_terms"][index] @ state[..., sector].reshape(-1)).real
                                    for sector, block in enumerate(self.blocks))) for index in range(3)]
        reduced_amplitudes = state.transpose(1, 2, 3, 0, 4).reshape(self.background_dimension, -1)
        means = [float(np.vdot(reduced_amplitudes, position @ reduced_amplitudes).real) for position in self.positions]
        total_square = sum(self.projected_squares)
        for first in range(3):
            for second in range(first + 1, 3):
                total_square = total_square + 2 * self.positions[first] @ self.positions[second]
        total_second = float(np.vdot(reduced_amplitudes, total_square @ reduced_amplitudes).real)
        return {"norm": float(np.linalg.norm(state)), "energies": energies, "individual_mode_energies": mode_energies,
                "mode_mean_amplitudes": means, "total_surface_fraction_second_moment": total_second,
                "overtone_nonvacuum_probabilities": [float(np.sum(abs(np.take(state, np.arange(1, self.parameters.cutoffs[index]), axis=index + 1)) ** 2)) for index in (1, 2)]}

    def record_probabilities(self, state):
        probabilities = np.zeros((2, 2, 2, 2, 2, 2))
        thresholds, denominator = postprocessing_table(self.parameters.sites, center_fraction=0.5)
        mode_weights = []
        for cutoff in self.parameters.cutoffs:
            vacuum = np.zeros(cutoff)
            vacuum[0] = 1.
            mode_weights.append((vacuum, 1 - vacuum))
        for data_bit, first_mode, second_mode, third_mode, first_clock, second_clock in np.ndindex(probabilities.shape):
            data_weight = thresholds / denominator if data_bit == 0 else 1 - thresholds / denominator
            clock = np.kron(clock_effect("X", first_clock), clock_effect("Y", second_clock))
            probabilities[data_bit, first_mode, second_mode, third_mode, first_clock, second_clock] = np.einsum(
                "jklma,ab,jklmb,j,k,l,m->", state.conj(), clock, state, data_weight,
                mode_weights[0][first_mode], mode_weights[1][second_mode], mode_weights[2][third_mode], optimize=True).real
        if probabilities.min() < -1e-12 or abs(probabilities.sum() - 1) > 1e-11:
            raise ValueError("The complete six-bit law must remain positive and normalized.")
        return probabilities

    def cutoff_residual(self, sector, vector):
        state = vector.reshape(self.shape[:-1])
        pieces = []
        for index, (row, cutoff, weight) in enumerate(zip(self.rows, self.parameters.cutoffs, self.parameters.interaction_weights)):
            top = np.take(state, cutoff - 1, axis=index + 1)
            exterior = self.parameters.scale * row["sigma"] * math.sqrt(cutoff) * weight * row["sources"][sector] @ top.reshape(self.parameters.sites, -1)
            pieces.append(exterior.reshape(-1))
        return np.concatenate(pieces)

    def cutoff_budget(self, initial, reference_duration=2., samples=129):
        if (initial.shape != self.shape or not np.isclose(np.linalg.norm(initial), 1., rtol=0, atol=1e-12)
                or not np.isfinite(reference_duration) or reference_duration <= 0
                or not isinstance(samples, int) or samples < 3):
            raise ValueError("Use a normalized joint state, positive finite horizon and at least three quadrature nodes.")
        if self.coherent_initial_error is None or not np.allclose(initial, self.initial_vector(), rtol=0, atol=1e-13):
            raise ValueError("This preparation budget only covers the declared coherent fundamental and two vacuum overtones with the fixed data preparation.")
        times = np.linspace(0., reference_duration, samples)
        squared = np.zeros(samples)
        second_derivative_norms = []
        for sector, block in enumerate(self.blocks):
            vector = initial[..., sector].reshape(-1)
            weight = float(np.vdot(vector, vector).real)
            shift = float(np.vdot(vector, block["generator"] @ vector).real / weight) if weight else 0.
            shifted = block["generator"] - shift * sparse.eye(len(vector), format="csr")
            trajectories = expm_multiply(-1j * shifted, vector, start=0., stop=reference_duration, num=samples,
                                          traceA=-1j * (block["generator_trace"] - shift * len(vector)))
            for position, trajectory in enumerate(trajectories):
                residual = self.cutoff_residual(sector, trajectory)
                squared[position] += float(np.vdot(residual, residual).real)
            exterior_norm_bound = math.sqrt(sum((self.parameters.scale * row["sigma"] * math.sqrt(cutoff) * weight * np.linalg.norm(row["sources"][sector], ord=2)) ** 2
                                               for row, cutoff, weight in zip(self.rows, self.parameters.cutoffs, self.parameters.interaction_weights)))
            second_derivative_norms.append(exterior_norm_bound * float(np.linalg.norm(shifted @ (shifted @ vector))))
        sampled = reference_duration / (samples - 1) * (np.sqrt(squared).sum() - (math.sqrt(squared[0]) + math.sqrt(squared[-1])) / 2)
        second_bound = float(np.linalg.norm(second_derivative_norms))
        quadrature = second_bound * reference_duration ** 3 / (12 * (samples - 1) ** 2)
        preparation = self.coherent_initial_error["initial_vector_error_bound"]
        return {"sampled_residual_norm_trapezoid": float(sampled), "vector_interpolation_remainder_bound": quadrature,
                "residual_vector_second_derivative_bound": second_bound,
                "evolution_vector_error_upper_bound": float(sampled + quadrature),
                "coherent_initial_vector_error": preparation,
                "same_final_record_TV_upper_bound": min(1., float(sampled + quadrature + preparation)),
                "quadrature_nodes": samples,
                "scope": "finite three-mode linear completion; exact residual/interpolation inequalities evaluated with floating sparse exponential action, not a rounding certificate"}

    def validity_budget(self):
        from quantum_chain_background import initial_data
        data = initial_data(self.parameters.sites)
        scale = self.parameters.scale
        base = single_model(SingleParameters(sites=self.parameters.sites, scale=scale))
        mean = 0.00025
        restoring = np.array([row["physical_inertia"] * row["physical_frequency"] ** 2 for row in self.rows])
        sources = np.array([row["source_norm"] * weight for row, weight in zip(self.rows, self.parameters.interaction_weights)])
        energy = sum(row["physical_frequency"] / 2 for row in self.rows) + restoring[0] * mean ** 2 / 2
        energy += sum(np.vdot(data, (block["matter_small"] + mean * self.parameters.interaction_weights[0] * self.rows[0]["sources"][sector]) @ data).real
                      for sector, block in enumerate(base.blocks)) / 4
        matter_floor = min(float(np.linalg.eigvalsh(block["matter_small"]).min()) for block in base.blocks)
        moment = 4 * (energy - matter_floor + np.sum(sources ** 2 / restoring)) * np.sum(1 / restoring)
        return {"ideal_three_mode_total_energy": float(energy), "sum_absolute_amplitudes_squared_expectation_upper_bound": float(moment),
                "pointwise_probability_sum_absolute_amplitudes_above_001": min(1., float(moment / 0.01 ** 2)),
                "scope": "finite three-mode energy/Cauchy/Markov bound for ideal coherent fundamental and two independent vacua; not an infinite-mode or nonlinear error bound"}


@lru_cache(maxsize=3)
def model(parameters=Parameters()):
    return ThreeRadialModes(parameters)


def embed_state(state, cutoffs):
    if len(cutoffs) != 3 or any(larger < smaller for larger, smaller in zip(cutoffs, state.shape[1:4])):
        raise ValueError("The comparison embedding cannot discard a mode or Fock coordinate.")
    output = np.zeros((state.shape[0], *cutoffs, 4), dtype=complex)
    output[:, :state.shape[1], :state.shape[2], :state.shape[3], :] = state
    return output


def terminal_schedule(parameters=Parameters(), reference_duration=2.):
    result = record_schedule(parameters.sites, reference_duration * parameters.scale, parameters.scale)
    for index in range(3):
        result.append({"kind": "collective_mode_vacuum_bit", "mode_index": index, "source_site": 0,
                       "hops": 0, "read_at": reference_duration * parameters.scale,
                       "received_at": reference_duration * parameters.scale})
    return sorted(result, key=lambda event: event["received_at"])


def release_joint_record(positions, clocks, mode_bits, random_integer, receipt, parameters=Parameters()):
    if len(positions) != parameters.sites or len(mode_bits) != 3 or any(bit not in (0, 1) for bit in mode_bits):
        raise ValueError("Retain every declared site and all three collective-mode bits.")
    data_bit, first, second = release_record(positions, clocks, random_integer, receipt, 2 * parameters.scale, parameters.scale)
    return data_bit, *mode_bits, first, second


def run_report():
    current = model()
    initial = current.initial_vector()
    final = current.evolve(initial, 2.)
    records = current.record_probabilities(final)
    inherited = single_model()
    inherited_initial = inherited.initial_vector()
    inherited_records = inherited.record_probabilities(inherited.evolve(inherited_initial, 2.))
    marginal = records.sum(axis=(2, 3))
    current_budget = current.cutoff_budget(initial, samples=257)
    inherited_budget = total_record_budget(inherited, inherited_initial)
    reference = model(Parameters(cutoffs=(24, 5, 5)))
    reference_initial = reference.initial_vector()
    reference_final = reference.evolve(reference_initial, 2.)
    reference_records = reference.record_probabilities(reference_final)
    preliminary_reference_budget = reference.cutoff_budget(reference_initial, samples=513)
    reference_budget = reference.cutoff_budget(reference_initial, samples=1025)
    refined_single = single_model(SingleParameters(cutoff=24))
    refined_single_initial = refined_single.initial_vector()
    refined_single_records = refined_single.record_probabilities(refined_single.evolve(refined_single_initial, 2.))
    refined_single_budget = total_record_budget(refined_single, refined_single_initial)
    refined_marginal_difference = float(abs(reference_records.sum(axis=(2, 3)) - refined_single_records).sum() / 2)
    omitted = omission_budget()
    matrix = kinetic_gram()
    summaries = []
    for row in current.rows:
        summaries.append({key: value for key, value in row.items() if key not in ("mode", "terms", "sources")})
        summaries[-1]["source_sector_commutator_with_fundamental_max_norm"] = max(
            float(np.linalg.norm(source @ fundamental - fundamental @ source, ord=2))
            for source, fundamental in zip(row["sources"], current.rows[0]["sources"]))
    cutoff_comparisons = []
    for cutoffs in ((20, 2, 2), (20, 3, 3), (20, 4, 4)):
        candidate = model(Parameters(cutoffs=cutoffs))
        prepared = candidate.initial_vector()
        evolved = candidate.evolve(prepared, 2.)
        cutoff_comparisons.append({"cutoffs": list(cutoffs),
                                  "vector_difference_from_refined": float(np.linalg.norm(embed_state(evolved, reference.parameters.cutoffs) - reference_final)),
                                  "record_TV_difference_from_refined": float(abs(candidate.record_probabilities(evolved) - reference_records).sum() / 2),
                                  "budget": candidate.cutoff_budget(prepared, samples=129)})
    retained_two = model(Parameters(interaction_weights=(1., 1., 0.)))
    retained_two_records = retained_two.record_probabilities(retained_two.evolve(retained_two.initial_vector(), 2.))
    return {"model_version": "C13", "route_round": "13", "research_scope": "sufficiency only",
            "parameters": {"cutoffs": list(current.parameters.cutoffs), "scale": current.parameters.scale, "sites": current.parameters.sites,
                           "initial_mode_states": ["C12 coherent amplitude 0.00025", "vacuum", "vacuum"]},
            "physical_mode_parameters": summaries, "physical_reference_kinetic_gram": matrix.tolist(),
            "normalized_kinetic_gram": (matrix / np.sqrt(np.outer(np.diag(matrix), np.diag(matrix)))).tolist(),
            "omission_budgets": [omitted, omission_budget((0, 1)), omission_budget((0, 1, 2))],
            "current_Fock_budget": current_budget, "C12_Fock_budget": inherited_budget,
            "preliminary_513_node_refined_budget": preliminary_reference_budget,
            "refined_three_mode_budget": reference_budget, "refined_C12_budget": refined_single_budget,
            "cutoff_comparisons": cutoff_comparisons,
            "three_mode_records": records.tolist(), "retained_C12_record_marginal": marginal.tolist(),
            "single_mode_comparison": {"current_retained_record_TV_difference": float(abs(marginal - inherited_records).sum() / 2),
                                       "refined_retained_record_TV_difference": refined_marginal_difference,
                                       "preliminary_513_node_numerical_allowance_sum": preliminary_reference_budget["same_final_record_TV_upper_bound"] + refined_single_budget["same_final_record_TV_upper_bound"],
                                       "refined_numerical_allowance_sum": reference_budget["same_final_record_TV_upper_bound"] + refined_single_budget["same_final_record_TV_upper_bound"],
                                       "ideal_retained_record_difference_lower_with_Fock_allowance": max(0., refined_marginal_difference - reference_budget["same_final_record_TV_upper_bound"] - refined_single_budget["same_final_record_TV_upper_bound"]),
                                       "ideal_retained_record_analytic_upper_bound": omitted["retained_record_TV_upper_bound_after_omitted_trace"],
                                       "current_three_mode_vs_two_mode_record_TV_after_third_trace": float(abs(records.sum(axis=3) - retained_two_records.sum(axis=3)).sum() / 2),
                                       "not_a_laboratory_detection_or_rounding_certificate": True},
            "time_rows": [dict(reference_time=time, physical_time=time * current.parameters.scale, **current.diagnostics(current.evolve(initial, time))) for time in (0., 1., 2.)],
            "three_mode_linear_window_budget": current.validity_budget(),
            "resources": {"data_qubits": 8, "clock_qubits": 2, "retained_collective_modes": 3,
                          "default_cutoffs": list(current.parameters.cutoffs), "occupied_complex_sector_dimension": math.prod(current.shape),
                          "refined_cutoffs": [24, 5, 5], "refined_occupied_complex_sector_dimension": math.prod(reference.shape),
                          "raw_bits": 13, "coarse_joint_outcomes": 64, "independent_fair_bits": 12, "threshold_table_bits": 104,
                          "bit_hops": 20, "supplied_hop_interval": 1024., "all_records_ready_at": 22528., "resets": 0},
            "terminal_schedule": terminal_schedule(),
            "limits": ["only the first three inherited radial modes, not a bound on all higher radial or nonradial modes",
                       "independent vacuum preparation of the two added modes is a new explicit input needed by the reduced-record bound",
                       "higher frequency alone does not remove virtual or vacuum effects; no rotating-wave or adiabatic elimination is assumed",
                       "each finite-mode linear Hamiltonian has a lower-bounded oscillator completion; no infinite-mode source summability or renormalization proved",
                       "three collective vacuum readouts are supplied instruments, not derived local instantaneous detectors",
                       "energy closed only for the finite linear multimode model; nonlinear constraints, higher-amplitude terms and material interactions remain outside scope",
                       "kinetic off-diagonal numerical residuals, sparse exponential roundoff and background-mode errors are not interval certified",
                       "no experiments, samples, necessity proofs, or derivation of adopted gravity and quantum rules from cognition"]}


class RadialMultimodeTests(unittest.TestCase):
    def test_physical_radial_modes_are_kinetically_orthogonal(self):
        matrix = kinetic_gram()
        normal = matrix / np.sqrt(np.outer(np.diag(matrix), np.diag(matrix)))
        np.testing.assert_allclose(normal, np.eye(3), atol=2e-8)
        for index, mode in enumerate(radial_modes()):
            self.assertAlmostEqual(matrix[index, index], canonical_inertia(mode), places=7)

    def test_fundamental_normalization_and_source_equal_C12(self):
        first = mode_data()[0]
        inherited = single_model()
        self.assertAlmostEqual(first["sigma"], inherited.sigma)
        self.assertAlmostEqual(first["physical_frequency"], inherited.frequency)
        for source, block in zip(first["sources"], inherited.blocks):
            np.testing.assert_allclose(source, block["source"], atol=2e-12)

    def test_overtones_have_nonzero_vacuum_coupling_despite_zero_mean(self):
        data = mode_data()
        self.assertLess(data[0]["physical_frequency"], data[1]["physical_frequency"])
        self.assertLess(data[1]["physical_frequency"], data[2]["physical_frequency"])
        for row in data[1:]:
            self.assertGreater(row["vacuum_coupling_norm"], 0.)
        budget = omission_budget()
        self.assertGreater(budget["same_full_record_TV_upper_bound"], 1e-9)
        self.assertLess(budget["same_full_record_TV_upper_bound"], 0.01)

    def test_omission_budget_vanishes_for_full_retention_or_zero_time(self):
        self.assertEqual(omission_budget((0, 1, 2))["same_full_record_TV_upper_bound"], 0.)
        self.assertEqual(omission_budget(reference_duration=0.)["same_full_record_TV_upper_bound"], 0.)
        self.assertLess(omission_budget((0, 1))["same_full_record_TV_upper_bound"], omission_budget()["same_full_record_TV_upper_bound"])

    def test_zero_overtone_interactions_reproduce_C12_up_to_the_free_vacuum_phase(self):
        parameters = Parameters(cutoffs=(8, 2, 2), interaction_weights=(1., 0., 0.))
        current = model(parameters)
        state = current.evolve(current.initial_vector(), 2.)
        reference = single_model(SingleParameters(cutoff=8)).evolve(single_model(SingleParameters(cutoff=8)).initial_vector(), 2.)
        phase = np.exp(-1j * sum(row["reference_frequency"] for row in current.rows[1:]))
        np.testing.assert_allclose(state[:, :, 0, 0, :], reference * phase, atol=4e-13)
        self.assertLess(np.linalg.norm(state[:, :, 1:, :, :]), 1e-14)

    def test_three_mode_joint_energy_and_inverse_evolution_are_consistent(self):
        current = model(Parameters(cutoffs=(8, 2, 2)))
        initial = current.initial_vector()
        final = current.evolve(initial, 2.)
        np.testing.assert_allclose(current.evolve(final, -2.), initial, atol=3e-13)
        first, last = current.diagnostics(initial), current.diagnostics(final)
        self.assertAlmostEqual(first["energies"]["hamiltonian"], last["energies"]["hamiltonian"], places=13)
        self.assertAlmostEqual(last["energies"]["hamiltonian"], sum(last["energies"][term] for term in ("matter", "background", "interaction")), places=13)

    def test_full_joint_records_and_original_C12_marginal_remain_normalized(self):
        current = model(Parameters(cutoffs=(8, 2, 2)))
        records = current.record_probabilities(current.evolve(current.initial_vector(), 2.))
        self.assertEqual(records.shape, (2,) * 6)
        self.assertAlmostEqual(records.sum(), 1.)
        self.assertEqual(records.sum(axis=(2, 3)).shape, (2,) * 4)
        self.assertGreater(records[:, :, 1, :, :, :].sum(), 0.)

    def test_rectangular_exterior_residual_matches_larger_tensor_Hamiltonian(self):
        current = model(Parameters(cutoffs=(4, 2, 2)))
        larger = model(Parameters(cutoffs=(5, 3, 3)))
        generator = np.random.default_rng(13)
        initial = generator.normal(size=current.shape) + 1j * generator.normal(size=current.shape)
        initial /= np.linalg.norm(initial)
        extended = embed_state(initial, larger.parameters.cutoffs)
        for sector, block in enumerate(current.blocks):
            actual = (larger.blocks[sector]["generator"] @ extended[..., sector].reshape(-1)).reshape(larger.shape[:-1])
            actual[:, :4, :2, :2] -= (block["generator"] @ initial[..., sector].reshape(-1)).reshape(current.shape[:-1])
            self.assertAlmostEqual(np.linalg.norm(actual), np.linalg.norm(current.cutoff_residual(sector, initial[..., sector].reshape(-1))), places=12)

    def test_finite_cutoff_difference_is_covered_by_independent_residual_budgets(self):
        first = model(Parameters(cutoffs=(12, 2, 2)))
        second = model(Parameters(cutoffs=(16, 3, 3)))
        first_initial, second_initial = first.initial_vector(), second.initial_vector()
        difference = np.linalg.norm(embed_state(first.evolve(first_initial, 2.), second.parameters.cutoffs) - second.evolve(second_initial, 2.))
        allowance = first.cutoff_budget(first_initial)["same_final_record_TV_upper_bound"] + second.cutoff_budget(second_initial)["same_final_record_TV_upper_bound"]
        self.assertLess(difference, allowance)

    def test_zero_mean_vacuum_reduced_bound_matches_exact_forced_oscillator_dephasing(self):
        for coupling in (0.01, 0.1):
            for duration in (0.2, 1., 2.):
                frequency = 0.7
                exact_reduced_distance = -math.expm1(-4 * (coupling / frequency) ** 2 * (1 - math.cos(frequency * duration))) / 2
                squared = (coupling * duration) ** 2
                bound = math.expm1(squared / 2) + squared / 2
                self.assertLessEqual(exact_reduced_distance, bound)

    def test_omitted_vacuum_effect_on_C12_records_obeys_the_reduced_bound(self):
        current = model(Parameters(cutoffs=(12, 3, 3)))
        reference = model(Parameters(cutoffs=(12, 3, 3), interaction_weights=(1., 0., 0.)))
        full = current.record_probabilities(current.evolve(current.initial_vector(), 2.)).sum(axis=(2, 3))
        reduced = reference.record_probabilities(reference.evolve(reference.initial_vector(), 2.)).sum(axis=(2, 3))
        self.assertLess(abs(full - reduced).sum() / 2, omission_budget()["retained_record_TV_upper_bound_after_omitted_trace"])

    def test_total_surface_second_moment_is_not_the_sum_of_individual_variances(self):
        current = model(Parameters(cutoffs=(8, 3, 3)))
        diagnostics = current.diagnostics(current.evolve(current.initial_vector(), 2.))
        self.assertGreater(diagnostics["total_surface_fraction_second_moment"], 0.)
        budget = current.validity_budget()
        self.assertLess(diagnostics["total_surface_fraction_second_moment"], budget["sum_absolute_amplitudes_squared_expectation_upper_bound"])
        self.assertLess(budget["pointwise_probability_sum_absolute_amplitudes_above_001"], 0.1)

    def test_equivalent_ensembles_give_the_same_full_three_mode_density(self):
        current = model(Parameters(cutoffs=(2, 2, 2)))
        first, second = np.eye(8)[:, 0], np.eye(8)[:, 2]
        def output(data):
            state = current.evolve(current.initial_vector(data), 0.7).reshape(-1)
            return np.outer(state, state.conj())
        np.testing.assert_allclose((output(first) + output(second)) / 2,
                                   sum(output((first + sign * second) / math.sqrt(2)) for sign in (1, -1)) / 2, atol=4e-13)

    def test_common_J_encoding_matches_a_complete_small_multimode_gate(self):
        current = model(Parameters(cutoffs=(2, 2, 2)))
        initial = current.initial_vector().reshape(-1)
        gate = np.zeros((len(initial), len(initial)), dtype=complex)
        for sector, block in enumerate(current.blocks):
            values, vectors = np.linalg.eigh(block["generator"].toarray())
            indices = np.arange(sector, len(initial), 4)
            gate[np.ix_(indices, indices)] = (vectors * np.exp(-0.7j * values)) @ vectors.T
        final = current.evolve(current.initial_vector(), 0.7).reshape(-1)
        np.testing.assert_allclose(gate @ initial, final, atol=4e-13)
        lifted = real_lift(gate)
        np.testing.assert_allclose(lifted @ encode_state(np.outer(initial, initial.conj())) @ lifted.T,
                                   encode_state(np.outer(final, final.conj())), atol=4e-13)

    def test_raw_records_keep_all_mode_flags_and_wait_for_complete_delivery(self):
        schedule = terminal_schedule()
        self.assertEqual(len(schedule), 13)
        self.assertEqual(sum(event["hops"] for event in schedule), 20)
        positions = [0] * 8
        positions[3] = 1
        with self.assertRaises(ValueError):
            release_joint_record(positions, [0, 1], [1, 0, 0], 0, 22527.)
        self.assertEqual(release_joint_record(positions, [0, 1], [1, 0, 0], 0, 22528.)[1:], (1, 0, 0, 0, 1))

    def test_invalid_cutoff_and_retention_contracts_are_rejected(self):
        with self.assertRaises(ValueError):
            Parameters(cutoffs=(4, 1, 3))
        with self.assertRaises(ValueError):
            omission_budget((1,))
        with self.assertRaises(ValueError):
            omission_budget(reference_duration=-1.)
        current = model(Parameters(cutoffs=(4, 2, 2)))
        with self.assertRaises(ValueError):
            current.cutoff_budget(current.initial_vector(np.eye(8)[:, 0]))

    def test_vector_interpolation_remainder_reduces_quadratically_on_refinement(self):
        current = model(Parameters(cutoffs=(8, 2, 2)))
        initial = current.initial_vector()
        coarse = current.cutoff_budget(initial, samples=33)
        refined = current.cutoff_budget(initial, samples=65)
        self.assertAlmostEqual(coarse["vector_interpolation_remainder_bound"], 4 * refined["vector_interpolation_remainder_bound"])
        self.assertGreater(coarse["same_final_record_TV_upper_bound"], refined["same_final_record_TV_upper_bound"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    outcome = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RadialMultimodeTests))
    if not outcome.wasSuccessful():
        raise SystemExit(1)
    if arguments.write_results:
        report = run_report()
        Path(__file__).with_name("radial_multimode_reduction_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"tests": outcome.testsRun, "modes": report["physical_mode_parameters"],
                          "omission": report["omission_budgets"], "comparison": report["single_mode_comparison"],
                          "Fock_budget": report["current_Fock_budget"], "refined_Fock_budget": report["refined_three_mode_budget"],
                          "window": report["three_mode_linear_window_budget"]}, indent=2))