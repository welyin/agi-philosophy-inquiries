"""Test finite readout capacity of a classical continuous cognitive score."""

import argparse
import json
import math
import platform
import unittest
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np


PREPARATION_HALF_WIDTH = 0.25


@dataclass(frozen=True)
class UniformScore:
    lower: float
    upper: float

    def __post_init__(self):
        if not (math.isfinite(self.lower) and math.isfinite(self.upper)
                and self.lower < self.upper):
            raise ValueError("A score window needs finite ordered endpoints.")

    def shifted(self, offset):
        if not math.isfinite(offset):
            raise ValueError("The silent score offset must be finite.")
        return UniformScore(self.lower + offset, self.upper + offset)


@dataclass(frozen=True)
class Query:
    threshold: float
    below: "Query | None" = None
    above: "Query | None" = None

    def __post_init__(self):
        if not math.isfinite(self.threshold):
            raise ValueError("The threshold must be finite.")


def prepared_score(center, half_width=PREPARATION_HALF_WIDTH):
    if not math.isfinite(half_width) or half_width <= 0:
        raise ValueError("The preparation half-width must be positive.")
    return UniformScore(center - half_width, center + half_width)


def threshold_probability(score, threshold):
    if not math.isfinite(threshold):
        raise ValueError("The threshold must be finite.")
    return min(1.0, max(0.0, (score.upper - threshold) / (score.upper - score.lower)))


def committed_score(threshold, outcome):
    if outcome not in (0, 1):
        raise ValueError("The committed outcome must be binary.")
    return prepared_score(threshold + 2 * outcome - 1)


def read_branches(score, threshold, erase=True, bit_error=0.0):
    if not math.isfinite(bit_error) or not 0 <= bit_error < 0.5:
        raise ValueError("The bit error must lie in [0, 0.5).")
    if not erase and bit_error != 0:
        raise ValueError("The nondestructive control uses noiseless reads only.")
    true_one = threshold_probability(score, threshold)
    reported_one = bit_error + (1 - 2 * bit_error) * true_one
    branches = {}
    for outcome, probability in ((0, 1 - reported_one), (1, reported_one)):
        if probability == 0:
            continue
        if erase:
            after = committed_score(threshold, outcome)
        elif outcome == 0:
            after = UniformScore(score.lower, min(score.upper, threshold))
        else:
            after = UniformScore(max(score.lower, threshold), score.upper)
        branches[outcome] = probability, after
    return branches


def transcript_distribution(score, query, erase=True, bit_error=0.0):
    if query is None:
        return {(): 1.0}
    distribution = {}
    for outcome, (probability, after) in read_branches(score, query.threshold, erase, bit_error).items():
        next_query = query.above if outcome else query.below
        for suffix, suffix_probability in transcript_distribution(after, next_query, erase, bit_error).items():
            distribution[(outcome,) + suffix] = probability * suffix_probability
    return distribution


def balanced_queries(centers):
    centers = tuple(sorted(centers))
    if not centers or not all(math.isfinite(center) for center in centers):
        raise ValueError("A nonempty finite preparation code is required.")
    if len(set(centers)) != len(centers):
        raise ValueError("The code centers must be distinct.")
    if len(centers) == 1:
        return None
    middle = len(centers) // 2
    threshold = (centers[middle - 1] + centers[middle]) / 2
    return Query(threshold, balanced_queries(centers[:middle]), balanced_queries(centers[middle:]))


def matrix_from_distributions(distributions):
    transcripts = sorted(set().union(*(distribution.keys() for distribution in distributions)))
    channel = np.array([[distribution.get(transcript, 0.0) for distribution in distributions]
                       for transcript in transcripts])
    return transcripts, channel


def channel_matrix(preparations, query, erase=True, bit_error=0.0):
    distributions = [transcript_distribution(score, query, erase, bit_error) for score in preparations]
    return matrix_from_distributions(distributions)


def uniform_guess_success(channel):
    return float(np.max(channel, axis=1).sum() / channel.shape[1])


def disjoint_support_code_size(channel):
    count = channel.shape[1]
    if count > 12:
        raise ValueError("This exact code search is restricted to at most twelve preparations.")
    support = channel > 0
    for size in range(count, 0, -1):
        for code in combinations(range(count), size):
            if all(not np.any(support[:, first] & support[:, second])
                   for first, second in combinations(code, 2)):
                return size
    raise AssertionError("A nonempty channel must have a one-preparation code.")


def reset_factorization_error(preparations, query, bit_error=0.0):
    if query is None:
        raise ValueError("At least one query is needed for the first-result factorization.")
    conditional_suffixes = {}
    for outcome in (0, 1):
        next_query = query.above if outcome else query.below
        reset_state = committed_score(query.threshold, outcome)
        conditional_suffixes[outcome] = transcript_distribution(reset_state, next_query, True, bit_error)
    maximum = 0.0
    for score in preparations:
        first_one = bit_error + (1 - 2 * bit_error) * threshold_probability(score, query.threshold)
        factorized = {}
        for outcome, probability in ((0, 1 - first_one), (1, first_one)):
            for suffix, suffix_probability in conditional_suffixes[outcome].items():
                factorized[(outcome,) + suffix] = probability * suffix_probability
        actual = transcript_distribution(score, query, True, bit_error)
        for transcript in set(actual) | set(factorized):
            maximum = max(maximum, abs(actual.get(transcript, 0.0) - factorized.get(transcript, 0.0)))
    return maximum


def uniform_total_variation(first, second):
    overlap = max(0.0, min(first.upper, second.upper) - max(first.lower, second.lower))
    return 1 - overlap / max(first.upper - first.lower, second.upper - second.lower)


def continuity_results():
    reference = prepared_score(0.0)
    rows = []
    for offset in (0.0, 0.025, 0.05, 0.1, 0.25, 0.5):
        shifted = reference.shifted(offset)
        thresholds = (reference.lower, reference.upper, shifted.lower, shifted.upper)
        max_change = max(abs(threshold_probability(reference, threshold)
                             - threshold_probability(shifted, threshold)) for threshold in thresholds)
        restored = shifted.shifted(-offset)
        rows.append({
            "offset": offset,
            "total_variation": uniform_total_variation(reference, shifted),
            "max_threshold_probability_change": max_change,
            "max_inverse_endpoint_error": max(abs(restored.lower - reference.lower),
                                               abs(restored.upper - reference.upper)),
        })
    return {
        "preparation_width": reference.upper - reference.lower,
        "rows": rows,
        "scope": "Silent translations of finite-width classical score distributions are continuous in total variation. Data processing bounds any fixed finite readout protocol by this distance. The coordinate is a candidate priority margin, not physical position, phase, time, or a derived continuous ontic variable.",
    }


def route_score(score, destinations, slot_count=3):
    if not isinstance(slot_count, int) or slot_count < 1:
        raise ValueError("At least one workspace slot is required.")
    slots = [None] * slot_count
    slots[0] = score
    owner = 0
    for destination in destinations:
        if not isinstance(destination, int) or not 0 <= destination < slot_count:
            raise ValueError("The destination must name an existing workspace slot.")
        slots[owner], slots[destination] = slots[destination], slots[owner]
        owner = destination
    return tuple(slots)


def copied_score_transcripts(shared_latent_score, query, carriers):
    if not isinstance(carriers, int) or carriers < 0:
        raise ValueError("The number of available copied carriers must be nonnegative.")
    if query is None or carriers == 0:
        return {(): 1.0}
    distribution = {}
    for outcome, (probability, conditional_latent_score) in read_branches(
            shared_latent_score, query.threshold, erase=False).items():
        next_query = query.above if outcome else query.below
        for suffix, suffix_probability in copied_score_transcripts(
                conditional_latent_score, next_query, carriers - 1).items():
            distribution[(outcome,) + suffix] = probability * suffix_probability
    return distribution


def communication_results():
    centers = np.arange(-7.0, 8.0, 2.0)
    preparations = [prepared_score(float(center)) for center in centers]
    query = balanced_queries(centers)
    route = (1, 2, 0, 2)
    moved = [next(score for score in route_score(preparation, route) if score is not None)
             for preparation in preparations]
    _, original_channel = channel_matrix(preparations, query)
    _, moved_channel = channel_matrix(moved, query)
    rows = []
    for carriers in (1, 2, 3):
        distributions = [copied_score_transcripts(score, query, carriers) for score in preparations]
        transcripts, channel = matrix_from_distributions(distributions)
        rows.append({
            "carriers_with_same_initial_score": carriers,
            "maximum_consumed_carriers": max(map(len, transcripts)),
            "zero_error_code_size_for_this_protocol": disjoint_support_code_size(channel),
            "optimal_uniform_guess_success_for_this_protocol": uniform_guess_success(channel),
        })
    same_question = Query(0, Query(0), Query(0))
    correlated = copied_score_transcripts(UniformScore(-1, 1), same_question, 2)
    return {
        "unique_carrier_transfer_route": list(route),
        "occupied_slots_after_transfer": sum(score is not None for score in route_score(preparations[0], route)),
        "max_channel_change_after_unique_carrier_transfer": float(np.max(abs(original_channel - moved_channel))),
        "copy_before_read_rows": rows,
        "same_score_copy_correlation_control": [{"transcript": list(record), "probability": probability}
                                                 for record, probability in sorted(correlated.items())],
        "scope": "Transfers move one uncommitted score without retaining a backup. The alternative copy operation duplicates the actual latent score before any read; copied registers are correlated, not independent re-preparations. Each queried copy is consumed and unused copies retain the old score. Allowing unlimited such copying removes the two-preparation bound on one original input. Neither prohibition nor physical availability of exact analog copying is derived from SoCA.",
    }


def finite_code_results():
    centers = np.arange(-7.0, 8.0, 2.0)
    preparations = [prepared_score(float(center)) for center in centers]
    query = balanced_queries(centers)
    rows = []
    for name, erase, bit_error in (
        ("commit_and_clear", True, 0.0),
        ("read_without_clearing", False, 0.0),
        ("commit_and_clear_bit_error_0_05", True, 0.05),
    ):
        transcripts, channel = channel_matrix(preparations, query, erase, bit_error)
        rows.append({
            "interface": name,
            "transcripts": [list(transcript) for transcript in transcripts],
            "conditional_channel": channel.tolist(),
            "maximum_query_count": max(map(len, transcripts)),
            "zero_error_code_size_for_this_protocol": disjoint_support_code_size(channel),
            "optimal_uniform_guess_success_for_this_protocol": uniform_guess_success(channel),
            "max_column_normalization_error": float(np.max(abs(channel.sum(axis=0) - 1))),
            "max_first_result_factorization_error": reset_factorization_error(preparations, query, bit_error) if erase else None,
        })
    return {
        "code_centers": centers.tolist(),
        "preparation_half_width": PREPARATION_HALF_WIDTH,
        "rows": rows,
        "scope": "One initially supplied continuous score, a fixed adaptive threshold tree, and fresh randomness independent of the preparation. Clear resets to a known result-dependent window and leaves no accessible backup. Finite table code sizes concern this tree; all-protocol bounds are proved separately in the research note.",
    }


def threshold_rank_results():
    rows = []
    for count in (2, 4, 8, 16):
        centers = np.arange(count, dtype=float) * 2
        preparations = [prepared_score(float(center)) for center in centers]
        thresholds = (centers[:-1] + centers[1:]) / 2
        effects = np.vstack((np.ones(count),
                             [[threshold_probability(score, float(threshold)) for score in preparations]
                              for threshold in thresholds]))
        rows.append({"number_of_preparations": count,
                     "test_rank_including_normalization": int(np.linalg.matrix_rank(effects))})
    return {
        "rows": rows,
        "scope": "Different thresholds are tested on separately prepared inputs, not all on the same consumed score. These exact step rows establish growing predictive rank, not growing single-score destructive-read capacity.",
    }


class OperationalCapacityTests(unittest.TestCase):
    def test_score_translation_is_continuous_and_has_an_inverse(self):
        score = prepared_score(0.3)
        restored = score.shifted(0.7).shifted(-0.7)
        self.assertAlmostEqual(restored.lower, score.lower)
        self.assertAlmostEqual(restored.upper, score.upper)
        direct = score.shifted(0.2 + 0.4)
        composed = score.shifted(0.2).shifted(0.4)
        self.assertAlmostEqual(direct.lower, composed.lower)
        self.assertAlmostEqual(direct.upper, composed.upper)
        self.assertAlmostEqual(threshold_probability(prepared_score(0.0).shifted(0.05), 0), 0.6)

    def test_uniform_window_threshold_probabilities(self):
        score = prepared_score(0.0)
        self.assertEqual(threshold_probability(score, -1), 1)
        self.assertEqual(threshold_probability(score, 1), 0)
        self.assertAlmostEqual(threshold_probability(score, 0), 0.5)

    def test_clearing_erases_every_input_detail_except_the_reported_result(self):
        first = read_branches(prepared_score(-3.0), 0)
        second = read_branches(prepared_score(-1.0), 0)
        self.assertEqual(first, second)
        self.assertEqual(first[0][1], prepared_score(-1.0))

    def test_nondestructive_read_conditions_instead_of_repreparing(self):
        score = UniformScore(-2, 2)
        branches = read_branches(score, 0, erase=False)
        self.assertEqual(branches[0], (0.5, UniformScore(-2, 0)))
        self.assertEqual(branches[1], (0.5, UniformScore(0, 2)))

    def test_complete_adaptive_record_factors_through_first_result(self):
        preparations = [prepared_score(center) for center in (-3.0, -0.1, 0.0, 0.4, 2.0)]
        query = Query(0, Query(-0.9, Query(-0.8), None), Query(0.9, None, Query(1.1)))
        for bit_error in (0.0, 0.05, 0.2):
            self.assertLess(reset_factorization_error(preparations, query, bit_error), 1e-14)
            _, channel = channel_matrix(preparations, query, True, bit_error)
            np.testing.assert_allclose(channel.sum(axis=0), 1, atol=1e-14)

    def test_binary_readouts_do_not_bound_capacity_without_clearing(self):
        rows = finite_code_results()["rows"]
        self.assertEqual(rows[0]["zero_error_code_size_for_this_protocol"], 2)
        self.assertEqual(rows[1]["zero_error_code_size_for_this_protocol"], 8)
        self.assertEqual(rows[1]["maximum_query_count"], 3)
        self.assertAlmostEqual(rows[0]["optimal_uniform_guess_success_for_this_protocol"], 0.25)
        self.assertAlmostEqual(rows[1]["optimal_uniform_guess_success_for_this_protocol"], 1)

    def test_nonzero_read_noise_changes_zero_error_but_not_to_zero_useful_information(self):
        row = finite_code_results()["rows"][2]
        self.assertEqual(row["zero_error_code_size_for_this_protocol"], 1)
        self.assertAlmostEqual(row["optimal_uniform_guess_success_for_this_protocol"], 0.2375)
        self.assertGreater(row["optimal_uniform_guess_success_for_this_protocol"], 1 / 8)

    def test_finite_destructive_capacity_does_not_bound_predictive_rank(self):
        for row in threshold_rank_results()["rows"]:
            self.assertEqual(row["test_rank_including_normalization"], row["number_of_preparations"])

    def test_noiseless_commit_is_repeatable_for_the_same_threshold(self):
        for threshold in (-4.0, 0.0, 6.0):
            for outcome, (_, after) in read_branches(prepared_score(threshold), threshold).items():
                repeated = read_branches(after, threshold)
                self.assertEqual(set(repeated), {outcome})
                self.assertAlmostEqual(repeated[outcome][0], 1)

    def test_read_and_commit_are_covariant_under_score_relabeling(self):
        score = prepared_score(0.1)
        threshold = 0.0
        offset = 0.7
        original = read_branches(score, threshold)
        shifted = read_branches(score.shifted(offset), threshold + offset)
        for outcome, (probability, after) in original.items():
            new_probability, new_after = shifted[outcome]
            self.assertAlmostEqual(probability, new_probability)
            self.assertAlmostEqual(after.lower + offset, new_after.lower)
            self.assertAlmostEqual(after.upper + offset, new_after.upper)

    def test_translation_continuity_bounds_all_threshold_effects(self):
        for row in continuity_results()["rows"]:
            bound = min(row["offset"] / (2 * PREPARATION_HALF_WIDTH), 1)
            self.assertAlmostEqual(row["total_variation"], bound)
            self.assertAlmostEqual(row["max_threshold_probability_change"], bound)
            self.assertLess(row["max_inverse_endpoint_error"], 1e-14)

    def test_no_read_has_only_one_indistinguishable_transcript(self):
        preparations = [prepared_score(center) for center in (-1.0, 0.0, 1.0)]
        _, channel = channel_matrix(preparations, None)
        self.assertEqual(disjoint_support_code_size(channel), 1)
        self.assertAlmostEqual(uniform_guess_success(channel), 1 / 3)

    def test_unique_carrier_transfer_does_not_duplicate_information(self):
        diagnostics = communication_results()
        self.assertEqual(diagnostics["occupied_slots_after_transfer"], 1)
        self.assertEqual(diagnostics["max_channel_change_after_unique_carrier_transfer"], 0)
        score = prepared_score(3)
        self.assertEqual(route_score(score, (1, 2, 0)), (score, None, None))

    def test_copying_before_reading_expands_complete_protocol_capacity(self):
        rows = communication_results()["copy_before_read_rows"]
        self.assertEqual([row["zero_error_code_size_for_this_protocol"] for row in rows], [2, 4, 8])
        np.testing.assert_allclose([row["optimal_uniform_guess_success_for_this_protocol"] for row in rows],
                                   [0.25, 0.5, 1.0])

    def test_copies_share_a_latent_value_instead_of_independent_noise(self):
        records = communication_results()["same_score_copy_correlation_control"]
        self.assertEqual(records, [{"transcript": [0, 0], "probability": 0.5},
                                   {"transcript": [1, 1], "probability": 0.5}])

    def test_irregular_adaptive_tree_cannot_beat_the_binary_first_result_bound(self):
        preparations = [prepared_score(center) for center in (-0.4, -0.1, 0.2, 0.6)]
        query = Query(0, Query(-0.8, None, Query(0.3)), Query(0.9, Query(1.1), None))
        for bit_error in (0.0, 0.05, 0.2):
            _, channel = channel_matrix(preparations, query, True, bit_error)
            self.assertLessEqual(uniform_guess_success(channel), 2 / len(preparations) + 1e-14)
            self.assertLess(reset_factorization_error(preparations, query, bit_error), 1e-14)

    def test_invalid_parameters_are_rejected(self):
        with self.assertRaises(ValueError):
            UniformScore(1, 1)
        with self.assertRaises(ValueError):
            prepared_score(0, 0)
        with self.assertRaises(ValueError):
            balanced_queries((0, 0))
        with self.assertRaises(ValueError):
            read_branches(prepared_score(0), 0, bit_error=-0.1)
        with self.assertRaises(ValueError):
            read_branches(prepared_score(0), 0, erase=False, bit_error=0.1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(OperationalCapacityTests)
    outcome = unittest.TextTestRunner(verbosity=2).run(suite)
    if not outcome.wasSuccessful():
        raise SystemExit(1)
    results = {
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "automated_checks": {"run": outcome.testsRun, "failures": len(outcome.failures),
                             "errors": len(outcome.errors)},
        "finite_preparation_code": finite_code_results(),
        "predictive_rank": threshold_rank_results(),
        "continuous_reorganization": continuity_results(),
        "communication": communication_results(),
    }
    if arguments.write_results:
        output = Path(__file__).with_name("operational_capacity_results.json")
        output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
        print(f"Results written to {output}")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()