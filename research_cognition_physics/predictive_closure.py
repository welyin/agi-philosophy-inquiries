"""Test exact and approximate finite summaries of continuous cognitive scores."""

import argparse
import json
import math
import platform
import unittest
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

from operational_capacity import Query, UniformScore, prepared_score, transcript_distribution


@dataclass(frozen=True)
class ScorePopulation:
    windows: tuple
    weights: tuple

    def __post_init__(self):
        if not self.windows or len(self.windows) != len(self.weights):
            raise ValueError("Every score window needs a weight.")
        if not all(isinstance(weight, Fraction) and weight >= 0 for weight in self.weights):
            raise ValueError("Exact nonnegative rational weights are required.")
        if sum(self.weights) != 1:
            raise ValueError("The population weights must sum to one.")

    def moment(self, order):
        if not isinstance(order, int) or order < 0:
            raise ValueError("Moment order must be a nonnegative integer.")
        result = Fraction(0)
        for window, weight in zip(self.windows, self.weights):
            lower, upper = Fraction(window.lower), Fraction(window.upper)
            result += weight * (upper**(order + 1) - lower**(order + 1)) / ((order + 1) * (upper - lower))
        return result

    def moments(self, order):
        return tuple(self.moment(degree) for degree in range(order + 1))

    def threshold_one(self, threshold):
        threshold = Fraction(threshold)
        result = Fraction(0)
        for window, weight in zip(self.windows, self.weights):
            lower, upper = Fraction(window.lower), Fraction(window.upper)
            fraction_above = min(Fraction(1), max(Fraction(0), (upper - threshold) / (upper - lower)))
            result += weight * fraction_above
        return result

    def shifted(self, offset):
        return ScorePopulation(tuple(window.shifted(float(offset)) for window in self.windows), self.weights)


def matched_moment_populations(order):
    if not isinstance(order, int) or order < 0:
        raise ValueError("The matched moment order must be nonnegative.")
    difference_order = order + 1
    windows = [[], []]
    weights = [[], []]
    for index in range(difference_order + 1):
        parity = index % 2
        center = 2 * index - difference_order
        windows[parity].append(prepared_score(center))
        weights[parity].append(Fraction(math.comb(difference_order, index), 2**order))
    return tuple(ScorePopulation(tuple(branch_windows), tuple(branch_weights))
                 for branch_windows, branch_weights in zip(windows, weights))


def witness_threshold(order):
    split = order // 2 + 1
    return 2 * split - 1 - (order + 1)


def translated_moments(moments, offset):
    offset = Fraction(offset)
    return tuple(sum(Fraction(math.comb(degree, lower_degree)) * offset**(degree - lower_degree)
                     * moments[lower_degree] for lower_degree in range(degree + 1))
                 for degree in range(len(moments)))


def mixed_population(first, second, first_weight=Fraction(1, 2)):
    first_weight = Fraction(first_weight)
    if not 0 <= first_weight <= 1:
        raise ValueError("The mixture weight must be a probability.")
    return ScorePopulation(first.windows + second.windows,
                           tuple(first_weight * weight for weight in first.weights)
                           + tuple((1 - first_weight) * weight for weight in second.weights))


def population_transcripts(population, query):
    distribution = {}
    for window, weight in zip(population.windows, population.weights):
        for record, probability in transcript_distribution(window, query).items():
            distribution[record] = distribution.get(record, 0.0) + float(weight) * probability
    return distribution


def histogram_population(population, edges):
    edges = np.asarray(edges, dtype=float)
    if (edges.ndim != 1 or len(edges) < 2 or not np.all(np.isfinite(edges))
            or np.any(np.diff(edges) <= 0)):
        raise ValueError("Histogram edges must be finite and strictly increasing.")
    if edges[0] > min(window.lower for window in population.windows) or edges[-1] < max(
            window.upper for window in population.windows):
        raise ValueError("The histogram domain must cover the entire score population.")
    windows, weights = [], []
    for lower, upper in zip(edges[:-1], edges[1:]):
        weight = population.threshold_one(float(lower)) - population.threshold_one(float(upper))
        if weight > 0:
            windows.append(UniformScore(float(lower), float(upper)))
            weights.append(weight)
    return ScorePopulation(tuple(windows), tuple(weights))


def maximum_threshold_gap(first, second):
    breakpoints = sorted({endpoint for population in (first, second) for window in population.windows
                          for endpoint in (window.lower, window.upper)})
    differences = [abs(first.threshold_one(threshold) - second.threshold_one(threshold))
                   for threshold in breakpoints]
    index = max(range(len(breakpoints)), key=lambda candidate: differences[candidate])
    return differences[index], breakpoints[index]


def transcript_distance(first, second, query):
    first_records = population_transcripts(first, query)
    second_records = population_transcripts(second, query)
    return sum(abs(first_records.get(record, 0.0) - second_records.get(record, 0.0))
               for record in set(first_records) | set(second_records)) / 2


def approximation_results():
    first, _ = matched_moment_populations(2)
    original = first.shifted(Fraction(1, 3))
    density_bound = 1 / min(window.upper - window.lower for window in original.windows)
    rows = []
    for bins in (8, 16, 32, 64, 128):
        edges = np.linspace(-4, 4, bins + 1)
        compressed = histogram_population(original, edges)
        error, threshold = maximum_threshold_gap(original, compressed)
        query = Query(float(threshold), Query(float(threshold) - 0.7), Query(float(threshold) + 0.9))
        rows.append({
            "grid_bins": bins,
            "grid_spacing": float(np.max(np.diff(edges))),
            "nonempty_bins": len(compressed.windows),
            "max_threshold_probability_error": float(error),
            "witness_threshold": float(threshold),
            "density_based_uniform_bound": float(density_bound * np.max(np.diff(edges)) / 4),
            "adaptive_record_total_variation": transcript_distance(original, compressed, query),
        })
    return {
        "domain": [-4.0, 4.0],
        "initial_population": "The first mean/variance witness translated by 1/3, represented with floating endpoints and exact rational weights.",
        "density_upper_bound": float(density_bound),
        "rows": rows,
        "scope": "Finite-bin memory reconstructs a uniform density within each bin. The true and reconstructed CDFs are piecewise linear; all breakpoints are checked, not just a sampled threshold grid. The L*spacing/4 bound assumes a density bounded by L and complete support coverage. Only the initial state is compressed in the readout comparison; later result-dependent reset states are handled exactly.",
    }


def recompression_results():
    original = ScorePopulation((prepared_score(0),), (Fraction(1),))
    edges = np.linspace(-4, 4, 33)
    offset = Fraction(1, 8)
    initial_summary = histogram_population(original, edges)
    moved_summary = histogram_population(initial_summary.shifted(offset), edges)
    returned_summary = histogram_population(moved_summary.shifted(-offset), edges)
    exact_return = original.shifted(offset).shifted(-offset)
    error, threshold = maximum_threshold_gap(original, returned_summary)
    return {
        "grid_spacing": float(edges[1] - edges[0]),
        "forward_offset": float(offset),
        "exact_round_trip_error": float(maximum_threshold_gap(original, exact_return)[0]),
        "recompressed_round_trip_error": float(error),
        "witness_threshold": float(threshold),
        "initial_variance_exact": str(original.moment(2) - original.moment(1)**2),
        "recompressed_variance_exact": str(returned_summary.moment(2) - returned_summary.moment(1)**2),
        "variance_increase_exact": str(returned_summary.moment(2) - original.moment(2)),
        "scope": "The underlying score translations are exactly inverse; projecting to a fixed histogram after each one loses detail. The difference belongs to the memory approximation, not a physical arrow of time or fundamental irreversibility.",
    }


def independent_sum_moments(first, second, order):
    first_moments, second_moments = first.moments(order), second.moments(order)
    return tuple(sum(Fraction(math.comb(degree, first_degree)) * first_moments[first_degree]
                     * second_moments[degree - first_degree] for first_degree in range(degree + 1))
                 for degree in range(order + 1))


def aggregation_results():
    local = ScorePopulation((prepared_score(-1), prepared_score(1)), (Fraction(1, 2), Fraction(1, 2)))
    window_noise_variance = Fraction(1, 48)
    rows = []
    for name, alignment in (("aligned_centers", 1), ("opposed_centers", -1)):
        center_pairs = [(-1, -alignment), (1, alignment)]
        cross_moment = sum(Fraction(first * second, 2) for first, second in center_pairs)
        sum_variance = sum(Fraction((first + second)**2, 2) for first, second in center_pairs)
        sum_variance += 2 * window_noise_variance
        rows.append({
            "preparation": name,
            "center_pairs": [list(pair) for pair in center_pairs],
            "cross_moment_exact": str(cross_moment),
            "aggregate_variance_exact": str(sum_variance),
            "aggregate_variance": float(sum_variance),
        })
    return {
        "same_local_moments_exact": [str(value) for value in local.moments(2)],
        "independent_aggregate_variance_exact": str(independent_sum_moments(local, local, 2)[2]),
        "rows": rows,
        "scope": "An added score-fusion control: X=S+U and Y=alignment*S+V, with a shared fair sign S and independent uniform noises U,V of half-width .25. Both local distributions agree exactly in both cases, but X+Y depends on their correlation. Sum fusion is an additional declared operation, not derived from the previous interface; its full output distribution need not be a finite uniform-window mixture.",
    }


def moment_results():
    rows = []
    for order in (1, 2, 4, 8, 12):
        first, second = matched_moment_populations(order)
        threshold = witness_threshold(order)
        first_probability = first.threshold_one(threshold)
        second_probability = second.threshold_one(threshold)
        difference = abs(first_probability - second_probability)
        rows.append({
            "highest_matched_moment": order,
            "moments_match_exactly": first.moments(order) == second.moments(order),
            "threshold": threshold,
            "first_probability_one": float(first_probability),
            "second_probability_one": float(second_probability),
            "probability_difference_exact": str(difference),
            "probability_difference": float(difference),
            "summary_only_worst_case_error_lower_bound": float(difference / 2),
            "first_unmatched_moment_difference_exact": str(first.moment(order + 1) - second.moment(order + 1)),
        })
    first, second = matched_moment_populations(2)
    query = Query(0, Query(-0.9), Query(0.9))
    first_records = population_transcripts(first, query)
    second_records = population_transcripts(second, query)
    records = sorted(set(first_records) | set(second_records))
    return {
        "mean_variance_witness": {
            "first_centers": [(window.lower + window.upper) / 2 for window in first.windows],
            "first_weights": [str(weight) for weight in first.weights],
            "second_centers": [(window.lower + window.upper) / 2 for window in second.windows],
            "second_weights": [str(weight) for weight in second.weights],
            "common_moments_exact": [str(value) for value in first.moments(2)],
            "common_variance": float(first.moment(2) - first.moment(1)**2),
            "adaptive_transcripts": [{"record": list(record), "first_probability": first_records.get(record, 0.0),
                                      "second_probability": second_records.get(record, 0.0)} for record in records],
            "adaptive_record_total_variation": sum(abs(first_records.get(record, 0.0)
                                                       - second_records.get(record, 0.0)) for record in records) / 2,
        },
        "matched_order_rows": rows,
        "scope": "Exact rational equalities for specified mixtures of the preceding study's uniform score windows. Finite raw moments update exactly under translation and mixing, yet do not determine threshold selection probabilities. This is a classical summary counterexample, not a quantum reconstruction.",
    }


class PredictiveClosureTests(unittest.TestCase):
    def test_matched_populations_are_normalized_positive_mixtures(self):
        for order in range(6):
            for population in matched_moment_populations(order):
                self.assertEqual(population.moment(0), 1)
                self.assertTrue(all(weight > 0 for weight in population.weights))

    def test_mean_and_variance_match_exactly(self):
        first, second = matched_moment_populations(2)
        self.assertEqual(first.moments(2), second.moments(2))
        self.assertEqual(first.moments(2), (Fraction(1), Fraction(0), Fraction(145, 48)))

    def test_same_mean_and_variance_do_not_determine_selection_probability(self):
        first, second = matched_moment_populations(2)
        self.assertEqual(first.threshold_one(0), Fraction(3, 4))
        self.assertEqual(second.threshold_one(0), Fraction(1, 4))

    def test_finite_difference_construction_matches_all_requested_moments(self):
        for order in (0, 1, 2, 4, 8, 12):
            first, second = matched_moment_populations(order)
            self.assertEqual(first.moments(order), second.moments(order))
            expected_difference = 2 * (-1)**(order + 1) * math.factorial(order + 1)
            self.assertEqual(first.moment(order + 1) - second.moment(order + 1), expected_difference)

    def test_threshold_witness_has_the_analytic_binomial_gap(self):
        for order in (0, 1, 2, 3, 4, 8, 12):
            first, second = matched_moment_populations(order)
            threshold = witness_threshold(order)
            difference = abs(first.threshold_one(threshold) - second.threshold_one(threshold))
            self.assertEqual(difference, Fraction(math.comb(order, order // 2), 2**order))

    def test_summary_dynamics_close_exactly_under_translation(self):
        for population in matched_moment_populations(4):
            shifted = population.shifted(Fraction(3, 4))
            self.assertEqual(shifted.moments(4), translated_moments(population.moments(4), Fraction(3, 4)))
            self.assertEqual(shifted.threshold_one(Fraction(3, 4)), population.threshold_one(0))

    def test_moment_and_read_probabilities_respect_random_mixtures(self):
        first, second = matched_moment_populations(2)
        mixture = mixed_population(first, second, Fraction(1, 3))
        self.assertEqual(mixture.moments(2), first.moments(2))
        expected = Fraction(1, 3) * first.threshold_one(0) + Fraction(2, 3) * second.threshold_one(0)
        self.assertEqual(mixture.threshold_one(0), expected)

    def test_existing_complete_read_protocol_detects_the_difference(self):
        diagnostics = moment_results()["mean_variance_witness"]
        self.assertAlmostEqual(diagnostics["adaptive_record_total_variation"], 0.5)
        self.assertAlmostEqual(sum(row["first_probability"] for row in diagnostics["adaptive_transcripts"]), 1)
        self.assertAlmostEqual(sum(row["second_probability"] for row in diagnostics["adaptive_transcripts"]), 1)

    def test_histogram_preserves_bin_masses_and_grid_readouts(self):
        original = matched_moment_populations(2)[0].shifted(Fraction(1, 3))
        edges = np.linspace(-4, 4, 33)
        compressed = histogram_population(original, edges)
        self.assertEqual(compressed.moment(0), 1)
        for threshold in edges:
            self.assertEqual(compressed.threshold_one(float(threshold)), original.threshold_one(float(threshold)))

    def test_histogram_projection_respects_random_mixing(self):
        first, second = matched_moment_populations(2)
        edges = np.linspace(-4, 4, 17)
        first, second = first.shifted(Fraction(1, 3)), second.shifted(Fraction(1, 3))
        mixed_first = histogram_population(mixed_population(first, second, Fraction(1, 3)), edges)
        projected_first = mixed_population(histogram_population(first, edges),
                                           histogram_population(second, edges), Fraction(1, 3))
        self.assertEqual(maximum_threshold_gap(mixed_first, projected_first)[0], 0)

    def test_density_based_interpolation_bound_can_be_attained(self):
        original = ScorePopulation((UniformScore(0.125, 0.625),), (Fraction(1),))
        compressed = histogram_population(original, np.linspace(-1, 1, 9))
        self.assertEqual(maximum_threshold_gap(original, compressed)[0], Fraction(1, 8))
        self.assertEqual(Fraction(2) * Fraction(1, 4) / 4, Fraction(1, 8))

    def test_all_threshold_and_complete_record_errors_respect_the_bound(self):
        for row in approximation_results()["rows"]:
            self.assertLessEqual(row["max_threshold_probability_error"], row["density_based_uniform_bound"] + 1e-13)
            self.assertAlmostEqual(row["adaptive_record_total_variation"], row["max_threshold_probability_error"], places=12)

    def test_recompression_breaks_an_exact_translation_round_trip(self):
        diagnostics = recompression_results()
        self.assertEqual(diagnostics["exact_round_trip_error"], 0)
        self.assertAlmostEqual(diagnostics["recompressed_round_trip_error"], 0.125)
        self.assertEqual(diagnostics["initial_variance_exact"], "1/48")
        self.assertEqual(diagnostics["recompressed_variance_exact"], "5/96")
        self.assertEqual(diagnostics["variance_increase_exact"], "1/32")

    def test_histogram_rejects_uncovered_mass_instead_of_silently_discarding_it(self):
        original = ScorePopulation((prepared_score(5),), (Fraction(1),))
        with self.assertRaises(ValueError):
            histogram_population(original, np.linspace(-4, 4, 9))

    def test_independent_score_fusion_has_closed_finite_moments(self):
        local = ScorePopulation((prepared_score(-1), prepared_score(1)), (Fraction(1, 2), Fraction(1, 2)))
        self.assertEqual(independent_sum_moments(local, local, 2), (Fraction(1), Fraction(0), Fraction(49, 24)))

    def test_identical_local_statistics_do_not_fix_correlated_score_fusion(self):
        diagnostics = aggregation_results()
        self.assertEqual(diagnostics["same_local_moments_exact"], ["1", "0", "49/48"])
        self.assertEqual([row["aggregate_variance_exact"] for row in diagnostics["rows"]], ["97/24", "1/24"])
        self.assertEqual([row["cross_moment_exact"] for row in diagnostics["rows"]], ["1", "-1"])

    def test_invalid_summary_inputs_are_rejected(self):
        with self.assertRaises(ValueError):
            matched_moment_populations(-1)
        with self.assertRaises(ValueError):
            ScorePopulation((prepared_score(0),), (Fraction(1, 2),))
        with self.assertRaises(ValueError):
            ScorePopulation((prepared_score(0),), (-1,))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(PredictiveClosureTests)
    outcome = unittest.TextTestRunner(verbosity=2).run(suite)
    if not outcome.wasSuccessful():
        raise SystemExit(1)
    results = {
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "automated_checks": {"run": outcome.testsRun, "failures": len(outcome.failures),
                             "errors": len(outcome.errors)},
        "finite_moment_summaries": moment_results(),
        "finite_histogram_approximation": approximation_results(),
        "recompression_round_trip": recompression_results(),
        "inter_agent_aggregation": aggregation_results(),
    }
    if arguments.write_results:
        output = Path(__file__).with_name("predictive_closure_results.json")
        output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
        print(f"Results written to {output}")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()