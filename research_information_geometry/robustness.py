"""Study identifiability and finite-count robustness of local phase responses."""

import argparse
import json
import math
import platform
import unittest
from pathlib import Path

import numpy as np

from experiment import CHAIN_EDGES, STAR_EDGES, graph_laplacian, propagator


SITE_COUNT = 4
PAIRS = tuple((source, target) for source in range(SITE_COUNT) for target in range(source + 1, SITE_COUNT))
LOCAL_OFFSETS = np.array([0.3, -0.2, 0.1, -0.2])
PROBE_TIMES = (0.1, 0.25, 0.5, 0.9)
GENERIC_PARAMETERS = np.array([1.0, 0.25, 0.15, 0.8, 0.2, 1.1, 0.3, -0.2, 0.1, 0.8])
SHOT_COUNTS = (1_000, 10_000, 100_000, 1_000_000)
SHOT_TIME_GRID = tuple(np.geomspace(0.001, 0.6, 49))
MONTE_CARLO_REPETITIONS = 2_000


def partially_coherent_state(coherence, site_count=SITE_COUNT):
    if not 0 <= coherence <= 1:
        raise ValueError("The coherence parameter must lie in [0, 1].")
    return (coherence * np.ones((site_count, site_count)) + (1 - coherence) * np.eye(site_count)) / site_count


def phase_kicked_state(state, source, phase):
    phases = np.ones(state.shape[0], dtype=complex)
    phases[source] = np.exp(1j * phase)
    return phases[:, None] * state * phases.conj()[None, :]


def occupation_probabilities(hamiltonian, state, source, phase, time):
    kicked = phase_kicked_state(state, source, phase)
    evolution = propagator(hamiltonian, time)
    probabilities = np.diag(evolution @ kicked @ evolution.conj().T).real
    if probabilities.min() < -1e-12 or not np.isclose(probabilities.sum(), 1):
        raise ValueError("Invalid occupation probabilities.")
    probabilities = np.maximum(probabilities, 0)
    return probabilities / probabilities.sum()


def response(hamiltonian, state, source, time):
    return (
        occupation_probabilities(hamiltonian, state, source, math.pi / 2, time)
        - occupation_probabilities(hamiltonian, state, source, -math.pi / 2, time)
    )


def analytic_finite_time_response(hamiltonian, coherence, source, time):
    evolution = propagator(hamiltonian, time)
    evolved_uniform_amplitudes = evolution.sum(axis=1)
    return -4 * coherence * (evolved_uniform_amplitudes.conj() * evolution[:, source]).imag / SITE_COUNT


def commutator_response_slope(hamiltonian, state, source):
    slopes = []
    for phase in (math.pi / 2, -math.pi / 2):
        kicked = phase_kicked_state(state, source, phase)
        derivative = -1j * (hamiltonian @ kicked - kicked @ hamiltonian)
        slopes.append(np.diag(derivative).real)
    return slopes[0] - slopes[1]


def finite_time_edge_estimate(hamiltonian, state, time, assumed_coherence):
    if time <= 0 or assumed_coherence <= 0:
        raise ValueError("Time and assumed coherence must be positive.")
    responses = [response(hamiltonian, state, source, time) for source in range(SITE_COUNT - 1)]
    return np.array([
        -SITE_COUNT * responses[source][target] / (4 * time * assumed_coherence)
        for source, target in PAIRS
    ])


def parameter_model(parameters):
    parameters = np.asarray(parameters, dtype=float)
    if parameters.shape != (10,):
        raise ValueError("Use six edge weights, three offset coordinates, and coherence.")
    edges = tuple((source, target, weight) for (source, target), weight in zip(PAIRS, parameters[:6]) if weight != 0)
    offsets = np.r_[parameters[6:9], -parameters[6:9].sum()]
    return graph_laplacian(SITE_COUNT, edges) + np.diag(offsets), partially_coherent_state(parameters[9])


def initial_slope_features(parameters):
    hamiltonian, state = parameter_model(parameters)
    slopes = [commutator_response_slope(hamiltonian, state, source) for source in range(SITE_COUNT - 1)]
    return np.array([slopes[source][target] for source, target in PAIRS])


def full_probability_features(parameters, times=PROBE_TIMES):
    hamiltonian, state = parameter_model(parameters)
    return np.concatenate([
        occupation_probabilities(hamiltonian, state, source, phase, time)
        for time in times
        for source in range(SITE_COUNT)
        for phase in (math.pi / 2, -math.pi / 2)
    ])


def finite_difference_jacobian(function, parameters, step=1e-5):
    columns = []
    for index in range(len(parameters)):
        displacement = np.zeros_like(parameters)
        displacement[index] = step
        columns.append((function(parameters + displacement) - function(parameters - displacement)) / (2 * step))
    return np.stack(columns, axis=1)


def jacobian_summary(jacobian):
    singular_values = np.linalg.svd(jacobian, compute_uv=False)
    rank = int(np.sum(singular_values > 1e-7 * singular_values[0])) if singular_values[0] else 0
    return {
        "shape": list(jacobian.shape),
        "singular_values": singular_values.tolist(),
        "rank_relative_threshold_1e_minus7": rank,
        "parameter_nullity": jacobian.shape[1] - rank,
        "coordinate_dependent_condition_number": float(singular_values[0] / singular_values[-1]) if rank == jacobian.shape[1] else None,
    }


def probability_settings(hamiltonian, state, time):
    return np.array([
        [occupation_probabilities(hamiltonian, state, source, phase, time)
         for phase in (math.pi / 2, -math.pi / 2)]
        for source in range(SITE_COUNT - 1)
    ])


def edge_moments(probabilities, time, assumed_coherence, shots):
    scale = -SITE_COUNT / (4 * time * assumed_coherence)
    positive = np.array([probabilities[source, 0, target] for source, target in PAIRS])
    negative = np.array([probabilities[source, 1, target] for source, target in PAIRS])
    mean = scale * (positive - negative)
    variance = scale**2 * (positive * (1 - positive) + negative * (1 - negative)) / shots
    return mean, variance


def sample_edges(probabilities, time, assumed_coherence, shots, repetitions, generator):
    counts = np.stack([
        np.stack([generator.multinomial(shots, phase_probabilities, size=repetitions)
                  for phase_probabilities in source_probabilities], axis=1)
        for source_probabilities in probabilities
    ], axis=1)
    scale = -SITE_COUNT / (4 * time * assumed_coherence * shots)
    return np.stack([
        scale * (counts[:, source, 0, target] - counts[:, source, 1, target])
        for source, target in PAIRS
    ], axis=1)


def wilson_interval(successes, repetitions):
    if repetitions <= 0 or not 0 <= successes <= repetitions:
        raise ValueError("Use a positive repetition count and a valid success count.")
    quantile = 1.959963984540054
    fraction = successes / repetitions
    denominator = 1 + quantile**2 / repetitions
    center = (fraction + quantile**2 / (2 * repetitions)) / denominator
    half_width = quantile * math.sqrt(fraction * (1 - fraction) / repetitions + quantile**2 / (4 * repetitions**2)) / denominator
    return [max(0.0, center - half_width), min(1.0, center + half_width)]


def count_study():
    cases = (
        ("chain_pure", CHAIN_EDGES, 1.0, 1.0, np.zeros(SITE_COUNT)),
        ("star_pure", STAR_EDGES, 1.0, 1.0, np.zeros(SITE_COUNT)),
        ("chain_mixed_known_coherence", CHAIN_EDGES, 0.8, 0.8, np.zeros(SITE_COUNT)),
        ("chain_mixed_wrong_pure_assumption", CHAIN_EDGES, 0.8, 1.0, np.zeros(SITE_COUNT)),
        ("chain_mixed_with_offsets", CHAIN_EDGES, 0.8, 0.8, LOCAL_OFFSETS),
    )
    case_results = {}
    for case_index, (name, edges, actual_coherence, assumed_coherence, offsets) in enumerate(cases):
        hamiltonian = graph_laplacian(SITE_COUNT, edges) + np.diag(offsets)
        state = partially_coherent_state(actual_coherence)
        actual_weights = np.array([-hamiltonian[source, target] for source, target in PAIRS])
        probabilities_by_time = [probability_settings(hamiltonian, state, time) for time in SHOT_TIME_GRID]
        count_rows = []
        for count_index, shots in enumerate(SHOT_COUNTS):
            moments = [edge_moments(probabilities, time, assumed_coherence, shots)
                       for time, probabilities in zip(SHOT_TIME_GRID, probabilities_by_time)]
            errors = [float(np.mean((mean - actual_weights)**2 + variance)) for mean, variance in moments]
            best_index = int(np.argmin(errors))
            time = SHOT_TIME_GRID[best_index]
            mean, variance = moments[best_index]
            generator = np.random.default_rng(np.random.SeedSequence([20260914, case_index, count_index]))
            samples = sample_edges(probabilities_by_time[best_index], time, assumed_coherence,
                                   shots, MONTE_CARLO_REPETITIONS, generator)
            graph_correct = np.all((samples > 0.5) == (actual_weights > 0.5), axis=1)
            success_rate = float(np.mean(graph_correct))
            count_rows.append({
                "shots_per_setting": shots,
                "total_shots_per_repetition": 2 * (SITE_COUNT - 1) * shots,
                "oracle_grid_time": float(time),
                "oracle_choice_on_grid_boundary": best_index in (0, len(SHOT_TIME_GRID) - 1),
                "mean_edge_estimates": mean.tolist(),
                "bias_rms": float(np.sqrt(np.mean((mean - actual_weights)**2))),
                "standard_deviation_rms": float(np.sqrt(np.mean(variance))),
                "predicted_edge_rmse": math.sqrt(errors[best_index]),
                "sampled_edge_rmse": float(np.sqrt(np.mean((samples - actual_weights)**2))),
                "all_six_edges_classified_correct_fraction": success_rate,
                "graph_fraction_monte_carlo_standard_error": math.sqrt(success_rate * (1 - success_rate) / MONTE_CARLO_REPETITIONS),
                "graph_fraction_95_percent_wilson_interval": wilson_interval(int(np.sum(graph_correct)), MONTE_CARLO_REPETITIONS),
                "expected_rmse_at_time_0_001": math.sqrt(errors[0]),
            })
        log_counts = np.log(np.array(SHOT_COUNTS, dtype=float))
        log_errors = np.log([row["predicted_edge_rmse"] for row in count_rows])
        case_results[name] = {
            "actual_coherence": actual_coherence,
            "assumed_coherence": assumed_coherence,
            "local_offsets": offsets.tolist(),
            "actual_weights": actual_weights.tolist(),
            "initial_slope_bias_floor_rms": float(np.sqrt(np.mean(((actual_coherence / assumed_coherence - 1) * actual_weights)**2))),
            "fitted_log_rmse_log_shots_slope": float(np.polyfit(log_counts, log_errors, 1)[0]),
            "rows": count_rows,
        }
    return {
        "evidence": "Synthetic multinomial counts from the specified quantum model, not measured laboratory data.",
        "shots_per_setting": SHOT_COUNTS,
        "time_grid": [float(time) for time in SHOT_TIME_GRID],
        "independent_repetitions": MONTE_CARLO_REPETITIONS,
        "seed_construction": "SeedSequence([20260914, case_index, count_index])",
        "sampling": "Independent source/phase settings; joint N-outcome multinomial counts preserve within-setting target correlations.",
        "time_choice": "Oracle minimization of the exact mean squared error on a fixed grid, using known simulation truth. Not an experimental data-driven design algorithm.",
        "edge_threshold": "0.5 for this zero-or-one edge family only; estimates are not clipped or rounded.",
        "cases": case_results,
    }


def reference_measurement_probabilities(state):
    mixer = np.eye(SITE_COUNT)
    mixer[:2, :2] = np.array([[1, 1], [1, -1]]) / math.sqrt(2)
    probabilities = np.diag(mixer @ state @ mixer.T).real
    return probabilities / probabilities.sum()


def calibration_study():
    shots, repetitions, time, coherence = 10_000, 4_000, 0.08, 0.8
    generator = np.random.default_rng(20260915)
    state = partially_coherent_state(coherence)
    hamiltonian = graph_laplacian(SITE_COUNT, CHAIN_EDGES)
    probabilities = probability_settings(hamiltonian, state, time)
    uncalibrated_samples = sample_edges(probabilities, time, 1.0, shots, repetitions, generator)
    calibration_counts = generator.multinomial(shots, reference_measurement_probabilities(state), size=repetitions)
    coherence_samples = SITE_COUNT * (calibration_counts[:, 0] - calibration_counts[:, 1]) / (2 * shots)
    usable = coherence_samples > 0
    calibrated_samples = uncalibrated_samples[usable] / coherence_samples[usable, None]
    actual_weights = np.array([-hamiltonian[source, target] for source, target in PAIRS])
    return {
        "evidence": "Synthetic counts with a trusted independent pair-mixing measurement. This is additional calibration information.",
        "shots_per_setting": shots,
        "time": time,
        "independent_repetitions": repetitions,
        "actual_coherence": coherence,
        "seed": 20260915,
        "settings_for_edges": 6,
        "extra_settings_for_calibration": 1,
        "total_shots_per_repetition": 7 * shots,
        "reference_probabilities": reference_measurement_probabilities(state).tolist(),
        "coherence_sample_mean": float(coherence_samples.mean()),
        "coherence_sample_standard_deviation": float(coherence_samples.std(ddof=1)),
        "coherence_exact_standard_deviation": math.sqrt((SITE_COUNT / 2 - coherence**2) / shots),
        "nonpositive_calibration_samples": int(np.sum(~usable)),
        "calibration_samples_above_one": int(np.sum(coherence_samples > 1)),
        "wrong_pure_assumption_edge_rmse": float(np.sqrt(np.mean((uncalibrated_samples - actual_weights)**2))),
        "known_exact_coherence_edge_rmse": float(np.sqrt(np.mean((uncalibrated_samples / coherence - actual_weights)**2))),
        "measured_coherence_edge_rmse": float(np.sqrt(np.mean((calibrated_samples - actual_weights)**2))) if np.any(usable) else None,
        "boundary_policy": "No clipping; nonpositive calibration estimates are counted and omitted from ratio statistics, whose usable sample count is reported.",
        "usable_ratio_samples": int(np.sum(usable)),
        "limitations": "The preparation-error family and reference mixer are assumed correct; readout and control errors are not included.",
    }


def local_fisher_analysis(parameters, times=PROBE_TIMES):
    function = lambda candidate: full_probability_features(candidate, times)
    probabilities = function(parameters)
    jacobian = finite_difference_jacobian(function, parameters)
    weighted = jacobian / np.sqrt(probabilities)[:, None]
    information = weighted.T @ weighted
    eigenvalues = np.linalg.eigvalsh(information)
    return {
        "rank": int(np.linalg.matrix_rank(information, tol=1e-8 * eigenvalues[-1])),
        "eigenvalues_per_one_shot_per_setting": eigenvalues.tolist(),
        "settings": len(times) * SITE_COUNT * 2,
        "local_cramer_rao_standard_deviation_at_1000_shots_per_setting": np.sqrt(np.diag(np.linalg.inv(1000 * information))).tolist(),
        "interpretation": "Coordinate-dependent local Cramer-Rao lower bounds under a correct model and regular unbiased estimation. Not achieved estimator errors or global uniqueness.",
    }


def structural_results():
    first_jacobian = finite_difference_jacobian(initial_slope_features, GENERIC_PARAMETERS)
    time_jacobian = finite_difference_jacobian(full_probability_features, GENERIC_PARAMETERS)
    reference_hamiltonian = graph_laplacian(SITE_COUNT, CHAIN_EDGES)
    alternatives = (
        (reference_hamiltonian, partially_coherent_state(1.0)),
        (2 * reference_hamiltonian, partially_coherent_state(0.5)),
    )
    slopes = [commutator_response_slope(hamiltonian, state, 0) for hamiltonian, state in alternatives]
    finite_responses = [response(hamiltonian, state, 0, 0.4) for hamiltonian, state in alternatives]
    wrong_uniform_formula_residual = np.max(abs(
        response(reference_hamiltonian + np.diag(LOCAL_OFFSETS), partially_coherent_state(1), 0, 0.5)
        + 4 * propagator(reference_hamiltonian + np.diag(LOCAL_OFFSETS), 0.5)[:, 0].imag / SITE_COUNT
    ))
    return {
        "parameters_for_local_rank_check": GENERIC_PARAMETERS.tolist(),
        "parameter_order": ["J01", "J02", "J03", "J12", "J13", "J23", "v0", "v1", "v2", "eta"],
        "offset_gauge": "v3=-(v0+v1+v2); a uniform energy shift is not identifiable from populations.",
        "probe_times": PROBE_TIMES,
        "first_slope_jacobian": jacobian_summary(first_jacobian),
        "finite_time_probability_jacobian": jacobian_summary(time_jacobian),
        "finite_time_multinomial_information": local_fisher_analysis(GENERIC_PARAMETERS),
        "equal_slope_example": {
            "models": ["unit chain with eta=1", "double-strength chain with eta=0.5"],
            "source": 0,
            "first_slopes": [slope.tolist() for slope in slopes],
            "max_first_slope_difference": float(np.max(abs(slopes[0] - slopes[1]))),
            "time_0_4_responses": [values.tolist() for values in finite_responses],
            "max_finite_time_difference": float(np.max(abs(finite_responses[0] - finite_responses[1]))),
        },
        "invalid_old_finite_time_identity_with_offsets_max_residual": float(wrong_uniform_formula_residual),
        "claim_scope": "Initial-slope non-identifiability is analytic; finite-time full rank is a local check at the stated generic parameters, not global uniqueness.",
    }


class RobustnessTests(unittest.TestCase):
    def test_correct_finite_time_formula_includes_evolved_uniform_amplitude(self):
        for edges in (CHAIN_EDGES, STAR_EDGES):
            hamiltonian = graph_laplacian(SITE_COUNT, edges) + np.diag(LOCAL_OFFSETS)
            for coherence in (0.3, 1.0):
                for source in range(SITE_COUNT):
                    for time in (0.0, 0.2, 0.5, 1.1):
                        direct = response(hamiltonian, partially_coherent_state(coherence), source, time)
                        analytic = analytic_finite_time_response(hamiltonian, coherence, source, time)
                        np.testing.assert_allclose(direct, analytic, atol=1e-12)

    def test_success_interval_does_not_claim_certainty_at_sample_endpoints(self):
        lower, upper = wilson_interval(2000, 2000)
        self.assertGreater(lower, 0.998)
        self.assertLess(lower, 1)
        self.assertAlmostEqual(upper, 1)
        lower, upper = wilson_interval(0, 2000)
        self.assertAlmostEqual(lower, 0)
        self.assertGreater(upper, 0)
        self.assertLess(upper, 0.002)

    def test_probability_settings_preserve_normalization(self):
        hamiltonian, state = parameter_model(GENERIC_PARAMETERS)
        probabilities = probability_settings(hamiltonian, state, 0.2)
        self.assertEqual(probabilities.shape, (3, 2, 4))
        np.testing.assert_allclose(probabilities.sum(axis=2), 1, atol=1e-12)

    def test_multinomial_simulation_agrees_with_exact_edge_moments(self):
        time, shots, repetitions = 0.08, 10_000, 12_000
        hamiltonian = graph_laplacian(SITE_COUNT, CHAIN_EDGES)
        probabilities = probability_settings(hamiltonian, partially_coherent_state(0.8), time)
        mean, variance = edge_moments(probabilities, time, 0.8, shots)
        samples = sample_edges(probabilities, time, 0.8, shots, repetitions, np.random.default_rng(8102))
        self.assertTrue(np.all(abs(samples.mean(axis=0) - mean) < 6 * np.sqrt(variance / repetitions)))
        np.testing.assert_allclose(samples.var(axis=0, ddof=1), variance, rtol=0.06)

    def test_finite_sample_variance_has_correct_inverse_shot_scaling(self):
        hamiltonian = graph_laplacian(SITE_COUNT, STAR_EDGES)
        probabilities = probability_settings(hamiltonian, partially_coherent_state(1), 0.1)
        _, variance = edge_moments(probabilities, 0.1, 1, 10_000)
        _, more_shots_variance = edge_moments(probabilities, 0.1, 1, 40_000)
        np.testing.assert_allclose(more_shots_variance, variance / 4, atol=1e-14)

    def test_two_time_coherence_confusion_has_nonzero_weight_bias_floor(self):
        hamiltonian = graph_laplacian(SITE_COUNT, CHAIN_EDGES)
        actual = np.array([-hamiltonian[source, target] for source, target in PAIRS])
        wrong = finite_time_edge_estimate(hamiltonian, partially_coherent_state(0.8), 1e-4, 1)
        np.testing.assert_allclose(wrong, 0.8 * actual, atol=1e-7)
        self.assertAlmostEqual(float(np.sqrt(np.mean((wrong - actual)**2))), math.sqrt(0.02), places=6)

    def test_second_order_bias_survives_real_local_offsets(self):
        hamiltonian = graph_laplacian(SITE_COUNT, CHAIN_EDGES) + np.diag(LOCAL_OFFSETS)
        state = partially_coherent_state(0.8)
        actual = np.array([-hamiltonian[source, target] for source, target in PAIRS])
        coarse = np.linalg.norm(finite_time_edge_estimate(hamiltonian, state, 0.02, 0.8) - actual)
        fine = np.linalg.norm(finite_time_edge_estimate(hamiltonian, state, 0.01, 0.8) - actual)
        self.assertAlmostEqual(float(coarse / fine), 4, delta=0.01)

    def test_reference_interference_calibrates_coherence(self):
        for coherence in (0.0, 0.3, 0.8, 1.0):
            probabilities = reference_measurement_probabilities(partially_coherent_state(coherence))
            expected = np.array([(1 + coherence) / 4, (1 - coherence) / 4, 0.25, 0.25])
            np.testing.assert_allclose(probabilities, expected, atol=1e-12)
            self.assertAlmostEqual(SITE_COUNT * (probabilities[0] - probabilities[1]) / 2, coherence)
            measured_variance = (SITE_COUNT / 2)**2 * (probabilities[0] + probabilities[1] - (probabilities[0] - probabilities[1])**2)
            self.assertAlmostEqual(measured_variance, SITE_COUNT / 2 - coherence**2)

    def test_local_information_is_positive_definite_at_selected_parameters(self):
        information = local_fisher_analysis(GENERIC_PARAMETERS)
        self.assertEqual(information["rank"], 10)
        self.assertGreater(information["eigenvalues_per_one_shot_per_setting"][0], 0.1)

    def test_initial_states_are_valid_and_populations_do_not_calibrate_coherence(self):
        for coherence in (0, 0.25, 0.8, 1):
            state = partially_coherent_state(coherence)
            np.testing.assert_allclose(np.diag(state), 0.25)
            self.assertAlmostEqual(float(np.trace(state)), 1)
            self.assertGreaterEqual(float(np.linalg.eigvalsh(state).min()), -1e-12)

    def test_general_slope_is_coherence_times_coupling_independent_of_offsets(self):
        for edges in (CHAIN_EDGES, STAR_EDGES):
            base = graph_laplacian(SITE_COUNT, edges)
            for coherence in (0.3, 0.8, 1):
                state = partially_coherent_state(coherence)
                for hamiltonian in (base, base + np.diag(LOCAL_OFFSETS)):
                    for source in range(SITE_COUNT):
                        slope = commutator_response_slope(hamiltonian, state, source)
                        expected = 4 * coherence * hamiltonian[:, source] / SITE_COUNT
                        expected[source] = 0
                        expected[source] = -expected.sum()
                        np.testing.assert_allclose(slope, expected, atol=1e-12)

    def test_finite_differences_match_commutator_slope_with_offsets(self):
        hamiltonian = graph_laplacian(SITE_COUNT, STAR_EDGES) + np.diag(LOCAL_OFFSETS)
        state = partially_coherent_state(0.7)
        for source in range(SITE_COUNT):
            expected = commutator_response_slope(hamiltonian, state, source)
            np.testing.assert_allclose(response(hamiltonian, state, source, 1e-4) / 1e-4, expected, atol=1e-6)

    def test_more_general_slope_depends_on_prepared_pair_coherence(self):
        amplitudes = np.array([0.4, 0.5j, 0.6, -0.3 + 0.2j])
        amplitudes /= np.linalg.norm(amplitudes)
        state = np.outer(amplitudes, amplitudes.conj())
        hamiltonian = graph_laplacian(SITE_COUNT, CHAIN_EDGES) + np.diag(LOCAL_OFFSETS)
        for source, target in PAIRS:
            expected = 4 * hamiltonian[target, source] * state[source, target].real
            self.assertAlmostEqual(commutator_response_slope(hamiltonian, state, source)[target], expected)

    def test_equal_initial_slopes_do_not_imply_equal_finite_time_data(self):
        hamiltonian = graph_laplacian(SITE_COUNT, CHAIN_EDGES)
        pure = partially_coherent_state(1)
        mixed = partially_coherent_state(0.5)
        np.testing.assert_allclose(commutator_response_slope(hamiltonian, pure, 0), commutator_response_slope(2 * hamiltonian, mixed, 0))
        difference = response(hamiltonian, pure, 0, 0.4) - response(2 * hamiltonian, mixed, 0, 0.4)
        self.assertGreater(float(np.max(abs(difference))), 0.01)

    def test_maximally_mixed_sector_state_has_no_dynamic_population_information(self):
        state = partially_coherent_state(0)
        for edges in (CHAIN_EDGES, STAR_EDGES):
            hamiltonian = graph_laplacian(SITE_COUNT, edges) + np.diag(LOCAL_OFFSETS)
            np.testing.assert_allclose(occupation_probabilities(hamiltonian, state, 0, 0.7, 0.9), 0.25, atol=1e-12)

    def test_uniform_energy_offset_is_exactly_unobservable(self):
        hamiltonian, state = parameter_model(GENERIC_PARAMETERS)
        for time in PROBE_TIMES:
            expected = occupation_probabilities(hamiltonian, state, 1, math.pi / 2, time)
            shifted = occupation_probabilities(hamiltonian + 3.7 * np.eye(SITE_COUNT), state, 1, math.pi / 2, time)
            np.testing.assert_allclose(shifted, expected, atol=1e-12)

    def test_first_slope_rank_matches_analytic_jacobian(self):
        expected = np.zeros((6, 10))
        expected[:, :6] = -4 * GENERIC_PARAMETERS[9] * np.eye(6) / SITE_COUNT
        expected[:, 9] = -4 * GENERIC_PARAMETERS[:6] / SITE_COUNT
        numeric = finite_difference_jacobian(initial_slope_features, GENERIC_PARAMETERS)
        np.testing.assert_allclose(numeric, expected, atol=1e-9)
        self.assertEqual(jacobian_summary(numeric)["rank_relative_threshold_1e_minus7"], 6)
        null_direction = np.r_[GENERIC_PARAMETERS[:6], np.zeros(3), -GENERIC_PARAMETERS[9]]
        np.testing.assert_allclose(numeric @ null_direction, 0, atol=1e-9)

    def test_finite_time_rank_check_is_stable_under_difference_step(self):
        coarse = finite_difference_jacobian(full_probability_features, GENERIC_PARAMETERS, 1e-5)
        fine = finite_difference_jacobian(full_probability_features, GENERIC_PARAMETERS, 2.5e-6)
        np.testing.assert_allclose(coarse, fine, atol=1e-8)
        self.assertEqual(jacobian_summary(coarse)["rank_relative_threshold_1e_minus7"], 10)

    def test_unmodeled_offsets_invalidate_old_finite_time_identity(self):
        hamiltonian = graph_laplacian(SITE_COUNT, CHAIN_EDGES) + np.diag(LOCAL_OFFSETS)
        actual = response(hamiltonian, partially_coherent_state(1), 0, 0.5)
        old_formula = -4 * propagator(hamiltonian, 0.5)[:, 0].imag / SITE_COUNT
        self.assertGreater(float(np.max(abs(actual - old_formula))), 0.01)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RobustnessTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    results = {
        "scope": "Finite quantum model calculations, not laboratory data or evidence for quantum gravity.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "structural_identifiability": structural_results(),
        "finite_count_study": count_study(),
        "independent_calibration_example": calibration_study(),
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    encoded = json.dumps(results, indent=2)
    if arguments.write_results:
        output = Path(__file__).resolve().parent / "robustness_results.json"
        output.write_text(encoded + "\n", encoding="utf-8")
        print(f"Results written to {output}")
    compact = {
        "automated_checks": results["automated_checks"],
        "first_slope_rank": results["structural_identifiability"]["first_slope_jacobian"],
        "finite_time_rank": results["structural_identifiability"]["finite_time_probability_jacobian"],
        "shot_study": {
            name: {
                "rmse_log_slope": case["fitted_log_rmse_log_shots_slope"],
                "rows": [{key: row[key] for key in ("shots_per_setting", "oracle_grid_time", "predicted_edge_rmse", "sampled_edge_rmse", "all_six_edges_classified_correct_fraction", "expected_rmse_at_time_0_001")}
                         for row in case["rows"]],
            } for name, case in results["finite_count_study"]["cases"].items()
        },
        "calibration_example": results["independent_calibration_example"],
        "scope": results["scope"],
    }
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()