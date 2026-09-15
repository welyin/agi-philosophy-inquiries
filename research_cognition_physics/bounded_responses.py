"""Check bounded finite translation spans as a new cognitive readout branch."""

import argparse
import json
import math
import platform
import unittest
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

from operational_capacity import PREPARATION_HALF_WIDTH, prepared_score
from predictive_closure import ScorePopulation


def harmonic_features(values, frequencies):
    values = np.asarray(values, dtype=float)
    columns = [np.ones_like(values)]
    for frequency in frequencies:
        columns.extend((np.cos(frequency * values), np.sin(frequency * values)))
    return np.stack(columns, axis=-1)


def translation_matrix(offset, frequencies):
    if not math.isfinite(offset):
        raise ValueError("The score translation must be finite.")
    matrix = np.eye(1 + 2 * len(frequencies))
    for index, frequency in enumerate(frequencies):
        cosine, sine = math.cos(frequency * offset), math.sin(frequency * offset)
        start = 1 + 2 * index
        matrix[start:start + 2, start:start + 2] = ((cosine, -sine), (sine, cosine))
    return matrix


def window_features(window, frequencies):
    center = (window.lower + window.upper) / 2
    half_width = (window.upper - window.lower) / 2
    features = harmonic_features(center, frequencies)
    for index, frequency in enumerate(frequencies):
        attenuation = float(np.sinc(frequency * half_width / math.pi))
        features[1 + 2 * index:3 + 2 * index] *= attenuation
    return features


def population_features(population, frequencies):
    return sum(float(weight) * window_features(window, frequencies)
               for window, weight in zip(population.windows, population.weights))


@dataclass(frozen=True)
class BoundedResponse:
    constant: float
    frequencies: tuple
    cosine_weights: tuple
    sine_weights: tuple

    def __post_init__(self):
        if not (len(self.frequencies) == len(self.cosine_weights) == len(self.sine_weights)):
            raise ValueError("Each frequency needs two real response coefficients.")
        parameters = (self.constant,) + self.frequencies + self.cosine_weights + self.sine_weights
        if not all(math.isfinite(value) for value in parameters):
            raise ValueError("Response parameters must be finite.")
        if any(frequency <= 0 for frequency in self.frequencies) or len(set(self.frequencies)) != len(self.frequencies):
            raise ValueError("Distinct positive frequencies are required.")
        radius_bound = sum(math.hypot(cosine, sine)
                           for cosine, sine in zip(self.cosine_weights, self.sine_weights))
        if not radius_bound <= min(self.constant, 1 - self.constant):
            raise ValueError("This control requires a sufficient global probability bound.")

    @property
    def coefficients(self):
        return np.array([self.constant] + [coefficient for pair in zip(self.cosine_weights, self.sine_weights)
                                          for coefficient in pair])

    def probability(self, values):
        return harmonic_features(values, self.frequencies) @ self.coefficients

    def population_probability(self, population):
        return float(self.coefficients @ population_features(population, self.frequencies))


MINIMAL_RESPONSE = BoundedResponse(0.5, (1.0,), (0.5,), (0.0,))
MULTIFREQUENCY_RESPONSE = BoundedResponse(0.5, (1.0, math.sqrt(2)), (0.2, 0.1), (0.0, 0.0))


def probability_effect(setting, outcome):
    if not math.isfinite(setting) or outcome not in (0, 1):
        raise ValueError("A finite readout setting and binary outcome are required.")
    sign = 2 * outcome - 1
    return np.array([0.5, 0.5 * sign * math.cos(setting), 0.5 * sign * math.sin(setting)])


def reset_population(setting, outcome):
    probability_effect(setting, outcome)
    center = setting if outcome == 1 else setting + math.pi
    return ScorePopulation((prepared_score(center),), (Fraction(1),))


def reset_branch_matrix(setting, outcome):
    after = population_features(reset_population(setting, outcome), (1.0,))
    return np.outer(after, probability_effect(setting, outcome))


def reset_test_probability(population, test):
    summary = population_features(population, (1.0,))
    for setting, outcome in test:
        summary = reset_branch_matrix(setting, outcome) @ summary
    return float(summary[0])


def population_expectation(population, function):
    nodes, weights = np.polynomial.legendre.leggauss(64)
    result = 0.0
    for window, probability in zip(population.windows, population.weights):
        center = (window.lower + window.upper) / 2
        half_width = (window.upper - window.lower) / 2
        result += float(probability) * float(np.dot(weights, function(center + half_width * nodes))) / 2
    return result


def retain_score_test_probability(population, test):
    def likelihood(values):
        product = np.ones_like(values)
        for setting, outcome in test:
            probability_effect(setting, outcome)
            probability_one = MINIMAL_RESPONSE.probability(values - setting)
            product *= probability_one if outcome else 1 - probability_one
        return product

    return population_expectation(population, likelihood)


def direct_reset_test_probability(population, test):
    probability = 1.0
    for setting, outcome in test:
        branch = retain_score_test_probability(population, ((setting, outcome),))
        probability *= branch
        population = reset_population(setting, outcome)
    return probability


def repeated_one_probability(population, count):
    if not isinstance(count, int) or count < 0:
        raise ValueError("A nonnegative repetition count is required.")
    frequencies = tuple(range(1, count + 1))
    summary = population_features(population, frequencies)
    probability = math.comb(2 * count, count) / 4**count
    for frequency in frequencies:
        probability += 2 * math.comb(2 * count, count - frequency) / 4**count * summary[2 * frequency - 1]
    return float(probability)


def locally_equivalent_populations():
    return (
        ScorePopulation((prepared_score(0), prepared_score(math.pi)), (Fraction(1, 2), Fraction(1, 2))),
        ScorePopulation((prepared_score(math.pi / 2), prepared_score(3 * math.pi / 2)),
                        (Fraction(1, 2), Fraction(1, 2))),
    )


def record_update_results():
    first, second = locally_equivalent_populations()
    first_summary = population_features(first, (1.0,))
    second_summary = population_features(second, (1.0,))
    test = ((0.0, 1), (0.0, 1))
    first_retained = retain_score_test_probability(first, test)
    second_retained = retain_score_test_probability(second, test)
    second_frequency_attenuation = float(np.sinc(2 * PREPARATION_HALF_WIDTH / math.pi))
    return {
        "first_preparation_centers": [0.0, math.pi],
        "second_preparation_centers": [math.pi / 2, 3 * math.pi / 2],
        "weights_each": [0.5, 0.5],
        "first_local_summary": first_summary.tolist(),
        "second_local_summary": second_summary.tolist(),
        "max_local_summary_difference": float(np.max(abs(first_summary - second_summary))),
        "second_frequency_attenuation": second_frequency_attenuation,
        "retain_score_probability_11_first": first_retained,
        "retain_score_probability_11_second": second_retained,
        "retain_score_probability_gap": abs(first_retained - second_retained),
        "analytic_retain_score_probability_11_first": 3 / 8 + second_frequency_attenuation / 8,
        "analytic_retain_score_probability_11_second": 3 / 8 - second_frequency_attenuation / 8,
        "reset_probability_11_first": reset_test_probability(first, test),
        "reset_probability_11_second": reset_test_probability(second, test),
        "same_setting_repeat_after_reset": reset_test_probability(reset_population(0, 1), ((0.0, 1),)),
        "scope": "The new response alone does not specify a state update. Retain-score reads draw conditionally independent Bernoulli outcomes without changing the latent score; they reveal higher harmonics. The explicitly chosen reset branch discards old details and prepares a fixed-width window at the selected extremum. Its rank-one real matrix gives closed local statistics, not a derived quantum measurement rule.",
    }


def repeated_read_rank_results():
    rows = []
    for depth in (1, 2, 4, 8):
        centers = np.linspace(0.0, math.pi, depth + 1)
        populations = [ScorePopulation((prepared_score(float(center)),), (Fraction(1),)) for center in centers]
        tests = np.array([[repeated_one_probability(population, count) for population in populations]
                          for count in range(depth + 1)])
        rows.append({
            "maximum_read_count": depth,
            "number_of_tested_preparations": len(populations),
            "sampled_rank_including_normalization": int(np.linalg.matrix_rank(tests, tol=1e-11)),
            "smallest_sampled_singular_value": float(np.linalg.svd(tests, compute_uv=False)[-1]),
        })
    return {
        "rows": rows,
        "scope": "Retain-score repeated-one likelihoods are powers of one nonconstant response. Their ranks are unbounded analytically, even without varying the readout setting. Finite-window averaging leaves each highest harmonic nonzero at the declared half-width .25. Sample ranks are checks, not the proof for arbitrary depth.",
    }


def copying_results():
    first, second = locally_equivalent_populations()
    test = ((0.0, 1), (0.0, 1))
    rows = []
    for name, population in (("first", first), ("second", second)):
        local_probability = MINIMAL_RESPONSE.population_probability(population)
        copied_probability = retain_score_test_probability(population, test)
        rows.append({
            "source": name,
            "local_probability_one": local_probability,
            "two_independently_prepared_carriers_probability_11": local_probability**2,
            "copied_same_latent_score_probability_11": copied_probability,
        })
    return {
        "rows": rows,
        "scope": "A separately added exact pre-read copy operation gives two carriers sharing the same old score, with independent readout randomness. Each carrier may reset after its own read. Copying therefore distinguishes source preparations identical under the restricted local reset interface; it cannot be represented by the local three-coordinate summary. This is a compatibility obstruction for this copying rule, not a universal quantum no-cloning theorem or Bell test.",
    }


def translation_span_results():
    values = np.linspace(-10, 10, 257)
    offsets = np.linspace(-2, 2, 21)
    rows = []
    for name, response in (("minimal_full_contrast", MINIMAL_RESPONSE),
                           ("two_incommensurate_frequencies", MULTIFREQUENCY_RESPONSE)):
        tests = np.column_stack([np.ones_like(values)] + [response.probability(values + offset) for offset in offsets])
        features = harmonic_features(values, response.frequencies)
        maximum_error = max(float(np.max(abs(harmonic_features(values + offset, response.frequencies)
                                              - features @ translation_matrix(float(offset), response.frequencies).T)))
                            for offset in offsets)
        rows.append({
            "response": name,
            "frequencies": list(response.frequencies),
            "constant": response.constant,
            "cosine_weights": list(response.cosine_weights),
            "sine_weights": list(response.sine_weights),
            "analytic_span_dimension": 1 + 2 * len(response.frequencies),
            "sampled_translation_rank": int(np.linalg.matrix_rank(tests, tol=1e-11)),
            "sampled_singular_values": np.linalg.svd(tests, compute_uv=False).tolist(),
            "maximum_shift_intertwining_error": maximum_error,
        })
    return {
        "rows": rows,
        "rank_tolerance": 1e-11,
        "scope": "Numerical checks of forms classified analytically under continuous bounded responses, every real shift and exact finite-dimensional closure. No finite sample of shifts proves the classification or excludes every other response. Minimality, full contrast, frequencies and the restricted interface remain declared assumptions.",
    }


def calibrated_response_results():
    centers = (0.0, math.pi / 2, math.pi, 3 * math.pi / 2, 2 * math.pi)
    width_factor = float(np.sinc(PREPARATION_HALF_WIDTH / math.pi))
    rows = []
    for center in centers:
        population = ScorePopulation((prepared_score(center),), (Fraction(1),))
        rows.append({
            "center": center,
            "pointwise_response": float(MINIMAL_RESPONSE.probability(center)),
            "half_angle_square": math.cos(center / 2)**2,
            "finite_window_probability": MINIMAL_RESPONSE.population_probability(population),
        })
    return {
        "frequency": 1.0,
        "preparation_half_width": PREPARATION_HALF_WIDTH,
        "first_frequency_attenuation": width_factor,
        "maximum_finite_window_probability": (1 + width_factor) / 2,
        "minimum_finite_window_probability": (1 - width_factor) / 2,
        "rows": rows,
        "nonmonotonic_control": [float(MINIMAL_RESPONSE.probability(center))
                                  for center in (0.0, math.pi, 2 * math.pi)],
        "scope": "The half-angle square is an identity for the calibrated pointwise response, not a Born-rule derivation. Allowed preparations here remain fixed-width uniform windows and their mixtures, so extrema zero and one are not exactly attained by their observed probabilities. The pointwise score dependence is not globally monotonic and cannot silently replace an ordered priority rule.",
    }


class BoundedResponseTests(unittest.TestCase):
    def test_translation_representation_matches_direct_evaluation(self):
        for row in translation_span_results()["rows"]:
            self.assertEqual(row["sampled_translation_rank"], row["analytic_span_dimension"])
            self.assertLess(row["maximum_shift_intertwining_error"], 1e-13)

    def test_translation_has_composition_and_a_legal_inverse_on_features(self):
        frequencies = (1.0, math.sqrt(2))
        first, second = translation_matrix(0.3, frequencies), translation_matrix(-0.7, frequencies)
        np.testing.assert_allclose(first @ second, translation_matrix(-0.4, frequencies), atol=1e-14)
        np.testing.assert_allclose(first @ translation_matrix(-0.3, frequencies), np.eye(5), atol=1e-14)
        np.testing.assert_allclose(first.T @ first, np.eye(5), atol=1e-14)

    def test_multiple_frequency_bounded_response_is_also_allowed(self):
        values = MULTIFREQUENCY_RESPONSE.probability(np.linspace(-100, 100, 2001))
        self.assertGreaterEqual(float(values.min()), 0.2 - 1e-14)
        self.assertLessEqual(float(values.max()), 0.8 + 1e-14)
        self.assertEqual(translation_span_results()["rows"][1]["analytic_span_dimension"], 5)

    def test_full_contrast_needs_an_added_calibration(self):
        response = BoundedResponse(0.4, (1.0,), (0.2,), (0.0,))
        self.assertAlmostEqual(float(response.probability(0)), 0.6)
        self.assertAlmostEqual(float(response.probability(math.pi)), 0.2)
        values = np.linspace(-4, 4, 101)
        np.testing.assert_allclose(MINIMAL_RESPONSE.probability(values), np.cos(values / 2)**2, atol=1e-14)

    def test_original_window_preparations_have_analytic_harmonic_averages(self):
        results = calibrated_response_results()
        attenuation = math.sin(PREPARATION_HALF_WIDTH) / PREPARATION_HALF_WIDTH
        self.assertAlmostEqual(results["first_frequency_attenuation"], attenuation)
        for row in results["rows"]:
            self.assertAlmostEqual(row["finite_window_probability"], (1 + attenuation * math.cos(row["center"])) / 2)
        self.assertGreater(results["minimum_finite_window_probability"], 0)

    def test_harmonic_summaries_preserve_preparation_mixing(self):
        population = ScorePopulation((prepared_score(-0.4), prepared_score(1.2)), (Fraction(1, 3), Fraction(2, 3)))
        features = population_features(population, MINIMAL_RESPONSE.frequencies)
        expected = window_features(population.windows[0], (1.0,)) / 3 + 2 * window_features(population.windows[1], (1.0,)) / 3
        np.testing.assert_allclose(features, expected, atol=1e-14)
        shifted = population.shifted(0.7)
        np.testing.assert_allclose(population_features(shifted, (1.0,)), translation_matrix(0.7, (1.0,)) @ features, atol=1e-14)

    def test_calibrated_minimal_response_is_not_an_ordered_priority_rule(self):
        np.testing.assert_allclose(calibrated_response_results()["nonmonotonic_control"], [1, 0, 1], atol=1e-14)

    def test_linear_and_real_exponential_modes_fail_global_probability_bounds(self):
        self.assertGreater(0.5 + 0.1 * 6, 1)
        self.assertGreater(0.5 + 0.1 * math.exp(0.4 * 10), 1)
        self.assertLess(0.5 + 0.1 * (-6), 0)

    def test_analytic_window_probabilities_match_independent_quadrature(self):
        population = locally_equivalent_populations()[0]
        for setting in (0.0, 0.4, 1.3):
            summary = population_features(population, (1.0,))
            actual = retain_score_test_probability(population, ((setting, 1),))
            self.assertAlmostEqual(actual, float(probability_effect(setting, 1) @ summary), places=13)

    def test_reset_branches_preserve_normalization_and_positive_probabilities(self):
        for center in (-1.1, 0.2, 3.0):
            population = ScorePopulation((prepared_score(center),), (Fraction(1),))
            summary = population_features(population, (1.0,))
            for setting in (0.0, 0.7):
                outputs = [reset_branch_matrix(setting, outcome) @ summary for outcome in (0, 1)]
                self.assertAlmostEqual(sum(output[0] for output in outputs), 1)
                for output in outputs:
                    self.assertGreater(output[0], 0)
                    self.assertLessEqual(float(np.linalg.norm(output[1:])), output[0] + 1e-14)

    def test_serial_reset_summary_matches_direct_classical_protocol(self):
        population = ScorePopulation((prepared_score(-0.4), prepared_score(1.2)), (Fraction(1, 3), Fraction(2, 3)))
        for outcomes in ((0, 0, 1), (1, 0, 0), (1, 1, 1)):
            test = tuple(zip((0.3, 1.1, -0.4), outcomes))
            self.assertAlmostEqual(reset_test_probability(population, test),
                                   direct_reset_test_probability(population, test), places=13)

    def test_same_local_summary_does_not_suffice_for_retained_score_reads(self):
        diagnostics = record_update_results()
        self.assertLess(diagnostics["max_local_summary_difference"], 1e-13)
        self.assertAlmostEqual(diagnostics["retain_score_probability_11_first"],
                               diagnostics["analytic_retain_score_probability_11_first"], places=13)
        self.assertAlmostEqual(diagnostics["retain_score_probability_11_second"],
                               diagnostics["analytic_retain_score_probability_11_second"], places=13)
        self.assertGreater(diagnostics["retain_score_probability_gap"], 0.2)

    def test_resetting_restores_local_equivalence_but_is_an_extra_rule(self):
        diagnostics = record_update_results()
        self.assertAlmostEqual(diagnostics["reset_probability_11_first"], diagnostics["reset_probability_11_second"])
        self.assertAlmostEqual(diagnostics["reset_probability_11_first"],
                               (1 + math.sin(0.25) / 0.25) / 4)

    def test_finite_window_reset_is_not_perfectly_repeatable(self):
        probability = record_update_results()["same_setting_repeat_after_reset"]
        self.assertAlmostEqual(probability, (1 + math.sin(0.25) / 0.25) / 2)
        self.assertLess(probability, 1)

    def test_repeated_query_formula_matches_direct_integration(self):
        for population in locally_equivalent_populations():
            for count in (0, 1, 2, 4, 8):
                direct = retain_score_test_probability(population, ((0.0, 1),) * count)
                self.assertAlmostEqual(direct, repeated_one_probability(population, count), places=13)

    def test_repeated_reads_require_more_than_a_fixed_first_harmonic_summary(self):
        for row in repeated_read_rank_results()["rows"]:
            self.assertEqual(row["sampled_rank_including_normalization"], row["maximum_read_count"] + 1)

    def test_copying_context_reveals_locally_invisible_preparation_detail(self):
        first, second = copying_results()["rows"]
        self.assertAlmostEqual(first["local_probability_one"], second["local_probability_one"])
        self.assertAlmostEqual(first["two_independently_prepared_carriers_probability_11"], 0.25)
        self.assertAlmostEqual(second["two_independently_prepared_carriers_probability_11"], 0.25)
        self.assertGreater(first["copied_same_latent_score_probability_11"]
                           - second["copied_same_latent_score_probability_11"], 0.2)

    def test_invalid_response_parameters_are_rejected(self):
        with self.assertRaises(ValueError):
            BoundedResponse(0.5, (1.0,), (0.6,), (0.0,))
        with self.assertRaises(ValueError):
            BoundedResponse(0.5, (float("nan"),), (0.2,), (0.0,))
        with self.assertRaises(ValueError):
            BoundedResponse(0.5, (1.0, 1.0), (0.1, 0.1), (0.0, 0.0))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(BoundedResponseTests)
    outcome = unittest.TextTestRunner(verbosity=2).run(suite)
    if not outcome.wasSuccessful():
        raise SystemExit(1)
    results = {
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "automated_checks": {"run": outcome.testsRun, "failures": len(outcome.failures),
                             "errors": len(outcome.errors)},
        "finite_translation_spans": translation_span_results(),
        "minimal_calibrated_response": calibrated_response_results(),
        "record_updates": record_update_results(),
        "repeated_read_ranks": repeated_read_rank_results(),
        "pre_read_copying": copying_results(),
    }
    if arguments.write_results:
        output = Path(__file__).with_name("bounded_response_results.json")
        output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
        print(f"Results written to {output}")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()