"""Compare realizable updates using a fixed, one-future-read distinction task.

The square-root bound follows from classical branch-kernel positivity. Its
attainment is proved only for the explicitly narrower fresh-window class.
No quantum state, quantum instrument, or global optimality claim is used.
"""

import argparse
import json
import math
import platform
import unittest
from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path

import numpy as np

from bounded_responses import (
    harmonic_features,
    locally_equivalent_populations,
    population_features,
    probability_effect,
    translation_matrix,
)
from operational_capacity import PREPARATION_HALF_WIDTH, prepared_score
from outcome_updates import WINDOW_ATTENUATION, partial_branch_matrix, validate_unit_interval
from predictive_closure import ScorePopulation


def fresh_branch_matrix(setting, outcome, contrast, longitudinal=1.0):
    """First moments of a kernel that emits mixtures of original-width windows."""
    validate_unit_interval(contrast)
    validate_unit_interval(longitudinal)
    probability_effect(setting, outcome)  # Validate the shared interface.
    sign = 2 * outcome - 1
    alpha = WINDOW_ATTENUATION
    root = math.sqrt(1 - contrast**2)
    local = np.array([
        [0.5, sign * contrast / 2, 0],
        [alpha * longitudinal * sign * contrast / 2, alpha * longitudinal / 2, 0],
        [0, 0, alpha * root / 2],
    ])
    rotation = translation_matrix(setting, (1.0,))
    return rotation @ local @ rotation.T


def fresh_kernel_modes(scores, setting, outcome, contrast, longitudinal=1.0):
    """Return two window centers and unnormalized weights for every input x.

    A branch of zero probability has arbitrary centers and zero weights.
    All output windows retain the original half-width; delta outputs are not
    preparations in this construction. No matrix update is used here.
    """
    validate_unit_interval(contrast)
    validate_unit_interval(longitudinal)
    probability_effect(setting, outcome)
    scores = np.asarray(scores, dtype=float)
    if not np.all(np.isfinite(scores)):
        raise ValueError("Scores must be finite.")
    relative = scores - setting
    sign = 2 * outcome - 1
    horizontal = np.cos(relative) + sign * contrast
    vertical = math.sqrt(1 - contrast**2) * np.sin(relative)
    likelihood = (1 + sign * contrast * np.cos(relative)) / 2
    centers = setting + np.stack((np.arctan2(vertical, horizontal),
                                  np.arctan2(vertical, -horizontal)), axis=-1)
    weights = likelihood[..., None] * np.array([(1 + longitudinal) / 2,
                                                (1 - longitudinal) / 2])
    return centers, weights


def fixed_tangent_kernel_modes(scores, setting, outcome, contrast):
    """A finite two-window witness with the same moments as longitudinal=0.

    Unlike the larger optimal family, its centers do not depend on the input.
    It therefore preserves finite mixtures as full distributions, too.
    """
    validate_unit_interval(contrast)
    probability_effect(setting, outcome)
    scores = np.asarray(scores, dtype=float)
    if not np.all(np.isfinite(scores)):
        raise ValueError("Scores must be finite.")
    relative = scores - setting
    denominator = 1 + (2 * outcome - 1) * contrast * np.cos(relative)
    transverse = math.sqrt(1 - contrast**2) * np.sin(relative)
    weights = np.stack((denominator + transverse, denominator - transverse), axis=-1) / 4
    centers = np.broadcast_to(np.array([setting + math.pi / 2, setting - math.pi / 2]), weights.shape)
    return centers, weights


@lru_cache(maxsize=None)
def quadrature(order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return nodes, weights / 2


def initial_cloud(population, order):
    nodes, weights = quadrature(order)
    points, masses = [], []
    for window, weight in zip(population.windows, population.weights):
        center = (window.lower + window.upper) / 2
        width = (window.upper - window.lower) / 2
        points.append(center + width * nodes)
        masses.append(float(weight) * weights)
    return np.concatenate(points), np.concatenate(masses)


def direct_kernel_summary(population, test, contrast, longitudinal=1.0, order=20):
    """Independent quadrature of the actual classical window-kernel paths.

    Intermediate clouds grow exponentially, so this is for short verification
    protocols, not the production predictor. Final window moments are also
    integrated numerically, without using the derived branch matrix.
    """
    points, masses = initial_cloud(population, order)
    nodes, weights = quadrature(order)
    for index, (setting, outcome) in enumerate(test):
        centers, mode_weights = fresh_kernel_modes(points, setting, outcome,
                                                   contrast, longitudinal)
        window_masses = (masses[:, None] * mode_weights).reshape(-1)
        window_centers = centers.reshape(-1)
        emitted = window_centers[:, None] + PREPARATION_HALF_WIDTH * nodes
        if index == len(test) - 1:
            moments = np.einsum("ijk,j->ik", harmonic_features(emitted, (1.0,)), weights)
            return window_masses @ moments
        points = emitted.reshape(-1)
        masses = (window_masses[:, None] * weights).reshape(-1)
    return masses @ harmonic_features(points, (1.0,))


def protocol_summary(initial, test, contrast, longitudinal=1.0):
    summary = np.asarray(initial, dtype=float)
    for setting, outcome in test:
        summary = fresh_branch_matrix(setting, outcome, contrast, longitudinal) @ summary
    return summary


def one_read_distance(first, second, contrast):
    """Supremum binary-record TV over all settings, for normalized summaries."""
    validate_unit_interval(contrast)
    first, second = np.asarray(first), np.asarray(second)
    if not (np.isclose(first[0], 1) and np.isclose(second[0], 1)):
        raise ValueError("Conditional summaries must be normalized.")
    return float(contrast * np.linalg.norm(first[1:] - second[1:]) / 2)


def retention_rows():
    alpha = WINDOW_ATTENUATION
    rows = []
    for contrast in (0.0, 0.005, 0.02, 0.25, 0.5, 0.75, 1.0):
        root = math.sqrt(1 - contrast**2)
        rows.append({
            "contrast": contrast,
            "baseline_tangent_gain": 1 - contrast,
            "fresh_window_tangent_gain": alpha * root,
            "general_periodic_kernel_upper_bound": root,
            "baseline_next_read_gap": contrast * alpha * (1 - contrast),
            "fresh_window_next_read_gap": contrast * alpha**2 * root,
            "general_next_read_upper_bound": contrast * alpha * root,
            "fresh_strictly_improves_baseline": 0 < contrast < 1 and alpha * root > 1 - contrast,
            "scope": "Gains are coefficients; at contrast=0 the read-distance ratio is undefined.",
        })
    return rows


def nonuniqueness_rows():
    contrast = 0.5
    initial = np.array([1.0, 0.0, 0.0])
    rows = []
    for longitudinal in (0.0, 0.5, 1.0):
        branch = fresh_branch_matrix(0, 1, contrast, longitudinal)
        after = branch @ initial
        probability_11 = (branch @ after)[0]
        rows.append({
            "longitudinal_parameter": longitudinal,
            "plus_branch_matrix": branch.tolist(),
            "branch_rank": int(np.linalg.matrix_rank(branch)),
            "conditional_summary_after_one": (after / after[0]).tolist(),
            "same_setting_next_one_probability": float(probability_11 / after[0]),
            "joint_probability_11": float(probability_11),
            "tangent_gain": WINDOW_ATTENUATION * math.sqrt(1 - contrast**2),
        })
    return rows


def joint_diagnostics():
    contrast = 0.5
    first = np.array([1.0, 0.0, WINDOW_ATTENUATION])
    second = np.array([1.0, 0.0, -WINDOW_ATTENUATION])
    initial = (np.outer(first, first) + np.outer(second, second)) / 2
    left = fresh_branch_matrix(0.2, 1, contrast, 0.0)
    right = fresh_branch_matrix(-0.4, 0, contrast, 1.0)
    nonselective = sum(fresh_branch_matrix(0.2, outcome, contrast) for outcome in (0, 1))
    direct = sum(float((left @ state)[0] * (right @ state)[0]) for state in (first, second)) / 2
    return {
        "matrix_joint_probability": float((left @ initial @ right.T)[0, 0]),
        "classical_component_probability": direct,
        "remote_marginal_error": float(np.max(abs((nonselective @ initial)[0] - initial[0]))),
        "swap_composition_error": float(np.max(abs((left @ initial @ right.T).T
                                                    - right @ initial.T @ left.T))),
    }


def verification_diagnostics():
    population = ScorePopulation((prepared_score(-0.6), prepared_score(1.1)),
                                 (Fraction(1, 3), Fraction(2, 3)))
    initial = population_features(population, (1.0,))
    protocol = ((0.3, 1), (-0.7, 0), (0.8, 1))
    predicted = protocol_summary(initial, protocol, 0.5, 0.7)
    errors = []
    for order in (12, 20, 28):
        actual = direct_kernel_summary(population, protocol, 0.5, 0.7, order)
        errors.append({"quadrature_order": order, "maximum_summary_error": float(np.max(abs(actual - predicted)))})
    return {"three_record_protocol": protocol, "matrix_summary": predicted.tolist(),
            "independent_kernel_quadrature": errors}


class RetentionTradeoffTests(unittest.TestCase):
    def test_finite_two_window_witness_has_nonnegative_correct_weights(self):
        scores = np.linspace(-math.pi, math.pi, 257)
        for contrast, outcome in product((0, 0.5, 1), (0, 1)):
            centers, weights = fixed_tangent_kernel_modes(scores, 0.3, outcome, contrast)
            self.assertGreaterEqual(float(weights.min()), -1e-15)
            np.testing.assert_allclose(weights.sum(axis=1),
                                       (1 + (2 * outcome - 1) * contrast * np.cos(scores - 0.3)) / 2)
            np.testing.assert_allclose(centers, np.broadcast_to([0.3 + math.pi / 2, 0.3 - math.pi / 2], centers.shape))

    def test_finite_two_window_witness_realizes_the_same_optimal_summary(self):
        scores = np.linspace(-math.pi, math.pi, 101)
        for contrast, outcome in product((0, 0.5, 1), (0, 1)):
            centers, weights = fixed_tangent_kernel_modes(scores, 0.2, outcome, contrast)
            features = harmonic_features(centers, (1.0,))
            features[..., 1:] *= WINDOW_ATTENUATION
            direct = np.einsum("ij,ijk->ik", weights, features)
            predicted = harmonic_features(scores, (1.0,)) @ fresh_branch_matrix(0.2, outcome, contrast, 0).T
            np.testing.assert_allclose(direct, predicted, atol=1e-14)

    def test_fixed_preparation_width_and_effect(self):
        self.assertEqual(PREPARATION_HALF_WIDTH, 0.25)
        for setting, outcome, contrast in product((0, 0.7), (0, 1), (0, 0.5, 1)):
            effect = probability_effect(setting, outcome)
            effect[1:] *= contrast
            np.testing.assert_allclose(fresh_branch_matrix(setting, outcome, contrast)[0], effect, atol=1e-14)

    def test_actual_kernel_has_nonnegative_correct_branch_weights(self):
        scores = np.linspace(-math.pi, math.pi, 257)
        for contrast, longitudinal, outcome in product((0, 0.5, 1), (0, 0.4, 1), (0, 1)):
            centers, weights = fresh_kernel_modes(scores, 0, outcome, contrast, longitudinal)
            self.assertTrue(np.all(np.isfinite(centers)))
            self.assertGreaterEqual(float(weights.min()), 0)
            np.testing.assert_allclose(weights.sum(axis=1), (1 + (2 * outcome - 1) * contrast * np.cos(scores)) / 2)

    def test_output_window_integration_matches_pointwise_first_moments(self):
        population = ScorePopulation((prepared_score(0.7),), (Fraction(1),))
        for contrast, longitudinal, outcome in product((0, 0.5, 1), (0, 0.7, 1), (0, 1)):
            direct = direct_kernel_summary(population, ((0.2, outcome),), contrast, longitudinal)
            matrix = fresh_branch_matrix(0.2, outcome, contrast, longitudinal)
            np.testing.assert_allclose(direct, matrix @ population_features(population, (1.0,)), atol=2e-14)

    def test_zero_probability_points_are_safe(self):
        for score, outcome in ((math.pi, 1), (0, 0)):
            centers, weights = fresh_kernel_modes(score, 0, outcome, 1, 0.5)
            self.assertTrue(np.all(np.isfinite(centers)))
            np.testing.assert_array_equal(weights, [0, 0])

    def test_kernel_respects_fresh_window_moment_bound(self):
        scores = np.linspace(-math.pi, math.pi, 1201)
        features = harmonic_features(scores, (1.0,))
        for contrast, longitudinal, outcome in product((0, 0.5, 0.99, 1), (0, 0.3, 1), (0, 1)):
            outputs = features @ fresh_branch_matrix(0, outcome, contrast, longitudinal).T
            self.assertLessEqual(float(np.max(np.linalg.norm(outputs[:, 1:], axis=1)
                                              - WINDOW_ATTENUATION * outputs[:, 0])), 2e-14)

    def test_allowed_summary_cone_is_preserved(self):
        summaries = harmonic_features(np.linspace(-math.pi, math.pi, 301), (1.0,))
        summaries[:, 1:] *= WINDOW_ATTENUATION
        for contrast, longitudinal, outcome in product((0, 0.5, 1), (0, 1), (0, 1)):
            outputs = summaries @ fresh_branch_matrix(0.3, outcome, contrast, longitudinal).T
            self.assertGreater(float(outputs[:, 0].min()), 0)
            self.assertLessEqual(float(np.max(np.linalg.norm(outputs[:, 1:], axis=1)
                                              - WINDOW_ATTENUATION * outputs[:, 0])), 2e-14)

    def test_one_read_distance_is_attained_by_an_allowed_setting(self):
        first = np.array([1, 0.2, 0.5])
        second = np.array([1, -0.4, 0.1])
        difference = first[1:] - second[1:]
        setting = math.atan2(difference[1], difference[0])
        effect = fresh_branch_matrix(setting, 1, 0.5)[0]
        self.assertAlmostEqual(one_read_distance(first, second, 0.5), float(effect @ (first - second)))

    def test_tangent_pair_has_equal_first_outcome_weights(self):
        alpha = WINDOW_ATTENUATION
        for outcome, contrast, longitudinal in product((0, 1), (0, 0.5, 1), (0, 1)):
            matrix = fresh_branch_matrix(0, outcome, contrast, longitudinal)
            self.assertAlmostEqual(float((matrix @ [1, 0, alpha])[0]), 0.5)
            self.assertAlmostEqual(float((matrix @ [1, 0, -alpha])[0]), 0.5)

    def test_transverse_gap_matches_prediction_for_both_outcomes(self):
        alpha = WINDOW_ATTENUATION
        for outcome, contrast, longitudinal in product((0, 1), (0.25, 0.5, 0.75), (0, 0.6, 1)):
            matrix = fresh_branch_matrix(0, outcome, contrast, longitudinal)
            first, second = matrix @ [1, 0, alpha], matrix @ [1, 0, -alpha]
            gap = one_read_distance(first / first[0], second / second[0], contrast)
            self.assertAlmostEqual(gap, contrast * alpha**2 * math.sqrt(1 - contrast**2))

    def test_actual_four_record_probabilities_include_all_first_outcomes(self):
        alpha, contrast = WINDOW_ATTENUATION, 0.5
        records = []
        for initial in ([1, 0, alpha], [1, 0, -alpha]):
            records.append(np.array([protocol_summary(initial, ((0, first), (math.pi / 2, second)), contrast)[0]
                                     for first, second in product((0, 1), repeat=2)]))
        self.assertAlmostEqual(float(records[0].sum()), 1)
        self.assertAlmostEqual(float(records[1].sum()), 1)
        self.assertAlmostEqual(float(abs(records[0] - records[1]).sum() / 2),
                               contrast * alpha**2 * math.sqrt(1 - contrast**2))

    def test_analytic_bound_is_saturated_in_fresh_window_class(self):
        alpha = WINDOW_ATTENUATION
        for contrast in (0, 0.25, 0.5, 0.99):
            angle = math.acos(-contrast)
            root = math.sqrt(1 - contrast**2)
            for longitudinal in (0, 0.5, 1):
                matrix = fresh_branch_matrix(0, 1, contrast, longitudinal)
                outputs = harmonic_features(np.array([angle, -angle]), (1.0,)) @ matrix.T
                self.assertAlmostEqual(float(np.linalg.norm(matrix[1:, 2])), alpha * root / 2)
                np.testing.assert_allclose(outputs[:, 1], 0, atol=1e-14)
                np.testing.assert_allclose(abs(outputs[:, 2]), alpha * outputs[:, 0], atol=1e-14)

    def test_exceeding_fresh_bound_has_a_specific_invalid_kernel_witness(self):
        contrast = 0.5
        matrix = fresh_branch_matrix(0, 1, contrast, 0)
        matrix[2, 2] *= 1.01
        output = matrix @ harmonic_features(math.acos(-contrast), (1.0,))
        self.assertGreater(float(np.linalg.norm(output[1:])), WINDOW_ATTENUATION * output[0])

    def test_general_bound_is_not_silently_replaced_by_fresh_bound(self):
        contrast = 0.5
        coefficient = math.sqrt(1 - contrast**2) / 2
        self.assertGreater(coefficient, fresh_branch_matrix(0, 1, contrast)[2, 2])
        # Saturating the looser pointwise bound with delta outputs would not
        # establish preservation of the allowed preparation set.
        self.assertLess(WINDOW_ATTENUATION, 1)

    def test_baseline_comparison_and_low_contrast_exception(self):
        rows = {row["contrast"]: row for row in retention_rows()}
        self.assertFalse(rows[0.005]["fresh_strictly_improves_baseline"])
        self.assertTrue(rows[0.5]["fresh_strictly_improves_baseline"])
        threshold = (1 - WINDOW_ATTENUATION**2) / (1 + WINDOW_ATTENUATION**2)
        self.assertAlmostEqual(WINDOW_ATTENUATION * math.sqrt(1 - threshold**2), 1 - threshold)

    def test_optimal_family_is_operationally_nonunique(self):
        first, _, last = nonuniqueness_rows()
        self.assertEqual(first["tangent_gain"], last["tangent_gain"])
        self.assertEqual(first["branch_rank"], 2)
        self.assertEqual(last["branch_rank"], 3)
        self.assertGreater(last["same_setting_next_one_probability"], first["same_setting_next_one_probability"])

    def test_improved_tangent_retention_does_not_mean_all_direction_dominance(self):
        contrast, alpha = 0.5, WINDOW_ATTENUATION
        baseline = sum(partial_branch_matrix(0, outcome, contrast) for outcome in (0, 1))
        fresh = sum(fresh_branch_matrix(0, outcome, contrast) for outcome in (0, 1))
        self.assertGreater(fresh[2, 2], baseline[2, 2])
        self.assertLess(fresh[1, 1], baseline[1, 1])
        np.testing.assert_allclose(np.diag(fresh), [1, alpha, alpha * math.sqrt(1 - contrast**2)])

    def test_short_sequences_match_independent_classical_path_integration(self):
        population = ScorePopulation((prepared_score(-0.4), prepared_score(1.3)),
                                     (Fraction(1, 3), Fraction(2, 3)))
        initial = population_features(population, (1.0,))
        for outcomes, longitudinal in product(product((0, 1), repeat=2), (0, 0.4, 1)):
            test = tuple(zip((0.2, -0.7), outcomes))
            direct = direct_kernel_summary(population, test, 0.5, longitudinal)
            np.testing.assert_allclose(direct, protocol_summary(initial, test, 0.5, longitudinal), atol=3e-14)

    def test_adaptive_protocol_normalizes_all_eight_records(self):
        initial = np.array([1, 0.2, -0.3])
        probabilities = []
        for outcomes in product((0, 1), repeat=3):
            settings = (0.3, -0.6 if outcomes[0] else 0.8, 0.4 if outcomes[1] else -1.1)
            probabilities.append(protocol_summary(initial, tuple(zip(settings, outcomes)), 0.6, 0.7)[0])
        self.assertGreaterEqual(min(probabilities), 0)
        self.assertAlmostEqual(sum(probabilities), 1)

    def test_previous_local_equivalences_survive_new_sequential_updates(self):
        first, second = locally_equivalent_populations()
        test = ((0.3, 1), (-0.4, 0), (1.2, 1))
        np.testing.assert_allclose(protocol_summary(population_features(first, (1.0,)), test, 0.5),
                                   protocol_summary(population_features(second, (1.0,)), test, 0.5), atol=1e-14)

    def test_covariance_and_outcome_relabeling(self):
        rotation = translation_matrix(0.7, (1.0,))
        for outcome in (0, 1):
            np.testing.assert_allclose(fresh_branch_matrix(0.9, outcome, 0.4) @ rotation,
                                       rotation @ fresh_branch_matrix(0.2, outcome, 0.4), atol=1e-14)
            np.testing.assert_allclose(fresh_branch_matrix(0.2 + math.pi, 1 - outcome, 0.4),
                                       fresh_branch_matrix(0.2, outcome, 0.4), atol=1e-14)

    def test_joint_moments_and_remote_marginal(self):
        data = joint_diagnostics()
        self.assertAlmostEqual(data["matrix_joint_probability"], data["classical_component_probability"])
        self.assertLess(data["remote_marginal_error"], 1e-14)
        self.assertLess(data["swap_composition_error"], 1e-14)

    def test_repeatability_and_latent_copy_obstructions_remain(self):
        contrast = 0.5
        self.assertLess((1 + contrast * WINDOW_ATTENUATION) / 2, 1)
        beta = math.sin(2 * PREPARATION_HALF_WIDTH) / (2 * PREPARATION_HALF_WIDTH)
        self.assertGreater(contrast**2 * beta / 4, 0)

    def test_invalid_parameters_are_rejected(self):
        for invalid in (-0.1, 1.1, math.nan, math.inf):
            with self.assertRaises(ValueError):
                fresh_branch_matrix(0, 1, invalid)
            with self.assertRaises(ValueError):
                fresh_kernel_modes(0, 0, 1, 0.5, invalid)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(RetentionTradeoffTests)
    outcome = unittest.TextTestRunner(verbosity=2).run(suite)
    if not outcome.wasSuccessful():
        raise SystemExit(1)
    alpha = WINDOW_ATTENUATION
    results = {
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "automated_checks": {"run": outcome.testsRun, "failures": len(outcome.failures), "errors": len(outcome.errors)},
        "preparation_half_width": PREPARATION_HALF_WIDTH,
        "window_attenuation": alpha,
        "strict_improvement_contrast_threshold": (1 - alpha**2) / (1 + alpha**2),
        "comparison": retention_rows(),
        "equally_optimal_fresh_window_family": nonuniqueness_rows(),
        "kernel_crosscheck": verification_diagnostics(),
        "joint_checks": joint_diagnostics(),
        "scope": {
            "metric": "Supremum one-future-read TV at the original fixed contrast, for a pair with identical first-outcome probabilities.",
            "proved_optimum": "Fresh-window kernels only: every pointwise conditional output is a mixture of the original half-width windows; periodic output moments and exact first-harmonic closure.",
            "broader_class": "For general periodic classical kernels preserving preparations only after input averaging, a looser upper bound is proved; attainment is unresolved.",
            "mixtures": "General probability mixtures of original windows. The lambda=0 witness also has an exactly finite two-window output for every input distribution.",
            "not_established": "Global optimum for all kernels or all future adaptive protocols; a unique update; quantum statistics; gravity.",
        },
    }
    if arguments.write_results:
        output = Path(__file__).with_name("retention_tradeoff_results.json")
        output.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        print(f"Results written to {output}")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
