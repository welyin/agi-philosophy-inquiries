"""Test realizable partial-memory updates of the bounded cognitive response."""

import argparse
import json
import math
import platform
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from bounded_responses import (
    MINIMAL_RESPONSE,
    harmonic_features,
    locally_equivalent_populations,
    population_expectation,
    population_features,
    probability_effect,
    reset_branch_matrix,
    reset_population,
    reset_test_probability,
    retain_score_test_probability,
    translation_matrix,
)
from operational_capacity import PREPARATION_HALF_WIDTH, prepared_score
from predictive_closure import ScorePopulation


WINDOW_ATTENUATION = math.sin(PREPARATION_HALF_WIDTH) / PREPARATION_HALF_WIDTH


def validate_unit_interval(value):
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError("The parameter must be a finite probability.")


def naive_retention_probability(population, retention):
    validate_unit_interval(retention)
    test = ((0.0, 1), (0.0, 1))
    return (retention * retain_score_test_probability(population, test)
            + (1 - retention) * reset_test_probability(population, test))


def naive_retention_results():
    first, second = locally_equivalent_populations()
    second_attenuation = math.sin(2 * PREPARATION_HALF_WIDTH) / (2 * PREPARATION_HALF_WIDTH)
    rows = []
    for retention in (0.0, 0.25, 0.5, 1.0):
        first_probability = naive_retention_probability(first, retention)
        second_probability = naive_retention_probability(second, retention)
        rows.append({
            "old_score_retention_probability": retention,
            "probability_11_first": first_probability,
            "probability_11_second": second_probability,
            "gap": first_probability - second_probability,
            "analytic_gap": retention * second_attenuation / 4,
        })
    return {
        "rows": rows,
        "scope": "After the first sharp-response read, a hidden coin either retains the old score or performs the preceding study's reset. The two source preparations have the same first-harmonic summary. Any nonzero constant retention probability exposes their different second harmonics.",
    }


def apparent_tangent_matrix(setting, outcome):
    effect = probability_effect(setting, outcome)
    tangent = np.array([-math.sin(setting), math.cos(setting)])
    scale = math.sqrt(1 - WINDOW_ATTENUATION**2) / 2
    matrix = np.zeros((3, 3))
    matrix[0] = effect
    matrix[1:, 1:] = scale * np.outer(tangent, tangent)
    return matrix


def cone_only_diagnostics():
    attenuation = WINDOW_ATTENUATION
    matrix = apparent_tangent_matrix(0.0, 1)
    angles = np.linspace(-math.pi, math.pi, 1001)
    summaries = harmonic_features(angles, (1.0,))
    summaries[:, 1:] *= attenuation
    outputs = summaries @ matrix.T
    slack = attenuation * outputs[:, 0] - np.linalg.norm(outputs[:, 1:], axis=1)
    boundary_input = np.array([1.0, -attenuation**2,
                               attenuation * math.sqrt(1 - attenuation**2)])
    boundary_output = matrix @ boundary_input
    latent_score = math.pi - 0.01
    latent_output = matrix @ harmonic_features(latent_score, (1.0,))
    return {
        "apparent_positive_branch_matrix": matrix.tolist(),
        "sampled_minimum_allowed_summary_slack": float(slack.min()),
        "saturating_input_summary": boundary_input.tolist(),
        "saturating_output_slack": float(attenuation * boundary_output[0]
                                         - np.linalg.norm(boundary_output[1:])),
        "latent_score_witness": latent_score,
        "latent_branch_probability": float(latent_output[0]),
        "required_conditional_first_harmonic_norm": float(np.linalg.norm(latent_output[1:]) / latent_output[0]),
        "universal_conditional_first_harmonic_norm_bound": 1.0,
        "allowed_preparation_summary_radius": attenuation,
        "scope": "This affine branch map is positive on the disk of permitted window summaries. It nevertheless violates the universal pointwise moment bound of any classical score kernel with periodic first-harmonic output moments realizing the same sharp effect and exact window-averaged closure. The bound on any latent conditional distribution is one, not the smaller radius of the allowed window preparation summaries. The implication from averaged closure to pointwise form uses the nonzero window Fourier multipliers and is proved in the research note.",
    }


def sharp_boundary_results():
    constraints = np.array([[1.0, -1.0, 0.0], [0.0, 0.0, -1.0]])
    coefficients = np.array([1.0, 1.0, 0.0])
    return {
        "constraints_on_each_output_moment_coefficients": constraints.tolist(),
        "constraint_rank": int(np.linalg.matrix_rank(constraints)),
        "surviving_coefficient_direction": coefficients.tolist(),
        "residual": (constraints @ coefficients).tolist(),
        "scope": "For a periodic classical kernel, exact closure after averaging over every fixed-width h=.25 window forces the unnormalized output first harmonic into span(1,cos,sin). Its norm is bounded by the sharp branch probability, which has a quadratic zero at pi. Value and first derivative therefore vanish there: a=b and c=0. The surviving harmonic is proportional to the branch effect, so only the observable first harmonic resets. Hidden higher moments need not reset.",
    }


def partial_branch_matrix(setting, outcome, strength):
    validate_unit_interval(strength)
    return strength * reset_branch_matrix(setting, outcome) + (1 - strength) * np.eye(3) / 2


def partial_test_summary(population, test, strength):
    summary = population_features(population, (1.0,))
    for setting, outcome in test:
        summary = partial_branch_matrix(setting, outcome, strength) @ summary
    return summary


def explicit_partial_paths(population, test, strength):
    validate_unit_interval(strength)
    branches = [(1.0, population)]
    for setting, outcome in test:
        probability_effect(setting, outcome)
        updated = []
        for weight, source in branches:
            if strength < 1:
                updated.append((weight * (1 - strength) / 2, source))
            if strength > 0:
                def sharp_likelihood(values):
                    probability_one = MINIMAL_RESPONSE.probability(values - setting)
                    return probability_one if outcome else 1 - probability_one

                probability = population_expectation(source, sharp_likelihood)
                updated.append((weight * strength * probability, reset_population(setting, outcome)))
        branches = updated
    return branches


def partial_retention_results():
    first = ScorePopulation((prepared_score(math.pi / 2),), (Fraction(1),))
    second = ScorePopulation((prepared_score(-math.pi / 2),), (Fraction(1),))
    read_one = ((0.0, 1),)
    read_two = ((0.0, 1), (math.pi / 2, 1))
    rows = []
    for strength in (0.0, 0.25, 0.5, 0.75, 1.0):
        first_branch = partial_test_summary(first, read_one, strength)
        second_branch = partial_test_summary(second, read_one, strength)
        first_conditional = partial_test_summary(first, read_two, strength)[0] / first_branch[0]
        second_conditional = partial_test_summary(second, read_two, strength)[0] / second_branch[0]
        rows.append({
            "readout_contrast": strength,
            "unconditional_old_population_retention_weight": 1 - strength,
            "first_record_one_probability_each": float(first_branch[0]),
            "first_conditional_summary": (first_branch / first_branch[0]).tolist(),
            "second_conditional_summary": (second_branch / second_branch[0]).tolist(),
            "next_record_one_probability_first": float(first_conditional),
            "next_record_one_probability_second": float(second_conditional),
            "next_record_probability_gap": float(first_conditional - second_conditional),
            "analytic_gap": strength * (1 - strength) * WINDOW_ATTENUATION,
            "maximum_same_setting_repeat_probability": (1 + strength * WINDOW_ATTENUATION) / 2,
        })
    equivalent_first, equivalent_second = locally_equivalent_populations()
    test = ((0.3, 1), (-0.4, 0), (1.2, 1))
    first_summary = partial_test_summary(equivalent_first, test, 0.5)
    second_summary = partial_test_summary(equivalent_second, test, 0.5)
    return {
        "rows": rows,
        "example_plus_branch_matrix_strength_0_5": partial_branch_matrix(0.0, 1, 0.5).tolist(),
        "max_equivalent_preparation_branch_difference": float(np.max(abs(first_summary - second_summary))),
        "scope": "An explicitly changed readout f=(1+strength*cos(x-setting))/2. A hidden classical mode selects sharp read-and-reset with probability strength, or a fair recorded bit and unchanged old population with probability 1-strength. The outcome is never overwritten. Each resulting branch is a mixture of previously allowed window preparations and the reset reference. The tradeoff is a constructive example, not an optimal or uniquely derived cognitive update.",
    }


def joint_window_components(aligned):
    first = ScorePopulation((prepared_score(math.pi / 2),), (Fraction(1),))
    second = ScorePopulation((prepared_score(-math.pi / 2),), (Fraction(1),))
    if aligned:
        return ((0.5, first, first), (0.5, second, second))
    return ((0.5, first, second), (0.5, second, first))


def joint_summary(components):
    return sum(weight * np.outer(population_features(first, (1.0,)), population_features(second, (1.0,)))
               for weight, first, second in components)


def joint_operation_results():
    strength = 0.5
    matrix = partial_branch_matrix(math.pi / 2, 1, strength)
    nonselective = partial_branch_matrix(0.4, 0, strength) + partial_branch_matrix(0.4, 1, strength)
    rows = []
    for name, aligned in (("aligned", True), ("opposed", False)):
        components = joint_window_components(aligned)
        initial = joint_summary(components)
        updated = matrix @ initial @ matrix.T
        direct = sum(weight * sum(branch_weight for branch_weight, _ in explicit_partial_paths(first, ((math.pi / 2, 1),), strength))
                     * sum(branch_weight for branch_weight, _ in explicit_partial_paths(second, ((math.pi / 2, 1),), strength))
                     for weight, first, second in components)
        remote_after = nonselective @ initial
        rows.append({
            "joint_preparation": name,
            "joint_summary": initial.tolist(),
            "first_marginal_summary": initial[:, 0].tolist(),
            "second_marginal_summary": initial[0].tolist(),
            "joint_record_11_probability": float(updated[0, 0]),
            "direct_classical_probability": float(direct),
            "analytic_probability": (1 + (1 if aligned else -1) * strength**2 * WINDOW_ATTENUATION**2) / 4,
            "max_remote_marginal_change_after_unconditional_local_read": float(np.max(abs(remote_after[0] - initial[0]))),
        })
    first_source, second_source = locally_equivalent_populations()
    copied_rows = []
    for source in (first_source, second_source):
        probability = population_expectation(source,
                                             lambda values: ((1 + strength * np.cos(values)) / 2)**2)
        copied_rows.append(float(probability))
    return {
        "readout_contrast": strength,
        "rows": rows,
        "same_latent_score_copy_probabilities_11": copied_rows,
        "same_latent_score_copy_gap": copied_rows[0] - copied_rows[1],
        "analytic_copy_gap": strength**2 * math.sin(2 * PREPARATION_HALF_WIDTH) / (8 * PREPARATION_HALF_WIDTH),
        "scope": "Nine joint moments W=E[b(X)b(Y)^T] suffice for the declared local shifts, partial instruments and carrier swap (W transposes). Independent mixtures of window products are used for these joint tests. Local outcomes update as M_A W M_B^T and cannot be predicted from marginal summaries alone. A separately added exact pre-read latent-score copy still exposes source second harmonics and is not determined by a three-coordinate source summary. No arbitrary interacting composite theory, quantum tensor product or Bell violation is constructed.",
    }


class OutcomeUpdateTests(unittest.TestCase):
    def test_no_retention_preserves_previous_local_equivalence(self):
        row = naive_retention_results()["rows"][0]
        self.assertAlmostEqual(row["probability_11_first"], row["probability_11_second"])
        self.assertAlmostEqual(row["probability_11_first"], (1 + WINDOW_ATTENUATION) / 4)

    def test_any_tested_nonzero_naive_retention_breaks_equivalence(self):
        for row in naive_retention_results()["rows"]:
            self.assertAlmostEqual(row["gap"], row["analytic_gap"], places=13)
            if row["old_score_retention_probability"] > 0:
                self.assertGreater(row["gap"], 0)

    def test_naive_retention_is_a_convex_mixture_of_valid_probabilities(self):
        for population in locally_equivalent_populations():
            reset = naive_retention_probability(population, 0)
            retained = naive_retention_probability(population, 1)
            actual = naive_retention_probability(population, 0.3)
            self.assertAlmostEqual(actual, 0.7 * reset + 0.3 * retained)
            self.assertGreaterEqual(actual, 0)
            self.assertLessEqual(actual, 1)

    def test_apparent_tangent_map_is_positive_on_allowed_summary_boundary(self):
        diagnostics = cone_only_diagnostics()
        self.assertGreaterEqual(diagnostics["sampled_minimum_allowed_summary_slack"], -1e-13)
        self.assertLess(abs(diagnostics["saturating_output_slack"]), 1e-13)

    def test_apparent_map_uses_the_existing_sharp_probability_effect(self):
        for setting in (0.0, 0.4, -1.1):
            for outcome in (0, 1):
                np.testing.assert_array_equal(apparent_tangent_matrix(setting, outcome)[0],
                                              probability_effect(setting, outcome))

    def test_apparent_map_is_covariant_under_simultaneous_relabeling(self):
        rotation = translation_matrix(0.7, (1.0,))
        original = apparent_tangent_matrix(0.2, 1)
        shifted = apparent_tangent_matrix(0.9, 1)
        np.testing.assert_allclose(shifted @ rotation, rotation @ original, atol=1e-14)

    def test_positive_summary_map_can_fail_classical_kernel_bound(self):
        diagnostics = cone_only_diagnostics()
        self.assertGreater(diagnostics["latent_branch_probability"], 0)
        self.assertGreater(diagnostics["required_conditional_first_harmonic_norm"], 20)
        self.assertEqual(diagnostics["universal_conditional_first_harmonic_norm_bound"], 1)
        self.assertLess(diagnostics["allowed_preparation_summary_radius"], 1)

    def test_sharp_boundary_leaves_only_an_effect_proportional_output(self):
        diagnostics = sharp_boundary_results()
        self.assertEqual(diagnostics["constraint_rank"], 2)
        np.testing.assert_array_equal(diagnostics["residual"], [0, 0])

    def test_partial_update_effect_has_the_declared_reduced_contrast(self):
        for strength in (0.0, 0.25, 0.5, 1.0):
            for outcome in (0, 1):
                effect = probability_effect(0.7, outcome)
                effect[1:] *= strength
                np.testing.assert_allclose(partial_branch_matrix(0.7, outcome, strength)[0], effect, atol=1e-14)

    def test_closed_branch_summary_matches_explicit_classical_kernel_paths(self):
        population = ScorePopulation((prepared_score(-0.4), prepared_score(1.3)), (Fraction(1, 3), Fraction(2, 3)))
        for strength in (0.0, 0.3, 1.0):
            for outcomes in product((0, 1), repeat=3):
                test = tuple(zip((0.1, -0.3, 1.2), outcomes))
                branches = explicit_partial_paths(population, test, strength)
                direct = sum(weight * population_features(source, (1.0,)) for weight, source in branches)
                np.testing.assert_allclose(partial_test_summary(population, test, strength), direct, atol=1e-13)

    def test_partial_branches_remain_in_the_allowed_window_mixture_cone(self):
        for angle in np.linspace(-math.pi, math.pi, 65):
            summary = harmonic_features(angle, (1.0,))
            summary[1:] *= WINDOW_ATTENUATION
            for strength in (0.0, 0.5, 1.0):
                outputs = [partial_branch_matrix(0.3, outcome, strength) @ summary for outcome in (0, 1)]
                self.assertAlmostEqual(sum(output[0] for output in outputs), 1)
                for output in outputs:
                    self.assertGreater(output[0], 0)
                    self.assertLessEqual(float(np.linalg.norm(output[1:])), WINDOW_ATTENUATION * output[0] + 1e-13)

    def test_partial_update_retains_observable_old_relations(self):
        for row in partial_retention_results()["rows"]:
            self.assertAlmostEqual(row["next_record_probability_gap"], row["analytic_gap"], places=13)
            if 0 < row["readout_contrast"] < 1:
                self.assertGreater(row["next_record_probability_gap"], 0)
        self.assertLess(partial_retention_results()["max_equivalent_preparation_branch_difference"], 1e-13)

    def test_partial_update_endpoints_and_relabeling_are_consistent(self):
        for outcome in (0, 1):
            np.testing.assert_array_equal(partial_branch_matrix(0.2, outcome, 0), np.eye(3) / 2)
            np.testing.assert_allclose(partial_branch_matrix(0.2, outcome, 1), reset_branch_matrix(0.2, outcome))
            rotation = translation_matrix(0.7, (1.0,))
            np.testing.assert_allclose(partial_branch_matrix(0.9, outcome, 0.4) @ rotation,
                                       rotation @ partial_branch_matrix(0.2, outcome, 0.4), atol=1e-14)

    def test_window_width_prevents_perfect_same_setting_repeatability(self):
        for row in partial_retention_results()["rows"]:
            strength = row["readout_contrast"]
            source = reset_population(0, 1)
            first = partial_test_summary(source, ((0.0, 1),), strength)[0]
            repeated = partial_test_summary(source, ((0.0, 1), (0.0, 1)), strength)[0] / first
            self.assertAlmostEqual(repeated, row["maximum_same_setting_repeat_probability"], places=13)
            self.assertLess(repeated, 1)

    def test_joint_moments_match_correlated_classical_protocols(self):
        for row in joint_operation_results()["rows"]:
            self.assertAlmostEqual(row["joint_record_11_probability"], row["direct_classical_probability"], places=13)
            self.assertAlmostEqual(row["joint_record_11_probability"], row["analytic_probability"], places=13)
            np.testing.assert_allclose(row["first_marginal_summary"], [1, 0, 0], atol=1e-13)
            np.testing.assert_allclose(row["second_marginal_summary"], [1, 0, 0], atol=1e-13)

    def test_ignored_local_outcome_leaves_remote_marginal_unchanged(self):
        for row in joint_operation_results()["rows"]:
            self.assertLess(row["max_remote_marginal_change_after_unconditional_local_read"], 1e-13)

    def test_joint_summary_swap_and_local_operations_compose(self):
        initial = joint_summary(joint_window_components(True))
        first = partial_branch_matrix(0.2, 1, 0.3)
        second = partial_branch_matrix(-0.4, 0, 0.7)
        np.testing.assert_allclose((first @ initial @ second.T).T, second @ initial.T @ first.T, atol=1e-14)

    def test_softening_does_not_remove_the_latent_copy_obstruction(self):
        diagnostics = joint_operation_results()
        self.assertAlmostEqual(diagnostics["same_latent_score_copy_gap"], diagnostics["analytic_copy_gap"], places=13)
        self.assertGreater(diagnostics["same_latent_score_copy_gap"], 0.05)

    def test_invalid_retention_is_rejected(self):
        population = locally_equivalent_populations()[0]
        for invalid in (-0.1, 1.1, float("nan")):
            with self.assertRaises(ValueError):
                naive_retention_probability(population, invalid)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(OutcomeUpdateTests)
    outcome = unittest.TextTestRunner(verbosity=2).run(suite)
    if not outcome.wasSuccessful():
        raise SystemExit(1)
    results = {
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "automated_checks": {"run": outcome.testsRun, "failures": len(outcome.failures),
                             "errors": len(outcome.errors)},
        "constant_retention_control": naive_retention_results(),
        "summary_positive_vs_kernel_positive": cone_only_diagnostics(),
        "sharp_boundary_constraint": sharp_boundary_results(),
        "closed_partial_retention": partial_retention_results(),
        "joint_operation_checks": joint_operation_results(),
    }
    if arguments.write_results:
        output = Path(__file__).with_name("outcome_update_results.json")
        output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
        print(f"Results written to {output}")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()