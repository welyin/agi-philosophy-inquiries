"""Rounds 200/202: finite memory for predicting a later-selected own action.

The encoder knows a uniform classical (x,m); the later decoder sees only a
finite message and the selected query. There is no post-encoding state access.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import unittest

import history_process_audit as loop


INPUTS = loop.STATES


def code_statistics(code, first_query_weight=F(1, 2)):
    if len(code) != 4 or not 0 <= first_query_weight <= 1:
        raise ValueError("Four input codewords and a probability are required.")
    weights = (first_query_weight, 1-first_query_weight)
    cells = defaultdict(list)
    for state, message in zip(INPUTS, code):
        cells[message].append(state)
    correct, brier, forecasts = F(0), F(0), {}
    for message, states in cells.items():
        for query, weight in enumerate(weights):
            ones = sum(state[query] for state in states)
            probability = F(ones, len(states))
            forecasts[(message, query)] = probability
            correct += weight*F(max(ones, len(states)-ones), 4)
            brier += weight*F(len(states), 4)*probability*(1-probability)
    return correct, brier, forecasts


@lru_cache(None)
def optimize_codes(messages, first_query_weight=F(1, 2)):
    if messages < 1:
        raise ValueError("At least one memory label is required.")
    best_success, best_brier = F(-1), F(2)
    success_code = brier_code = None
    for code in product(range(messages), repeat=4):
        success, brier, _ = code_statistics(code, first_query_weight)
        if success > best_success:
            best_success, success_code = success, code
        if brier < best_brier:
            best_brier, brier_code = brier, code
    return best_success, best_brier, success_code, brier_code


def whole_pair_success(code):
    # Each occupied message supports at most one correct full-pair guess.
    return F(len(set(code)), 4)


def reread_with_restoration(state, query):
    word = ("read",) if query == 0 else ("swap", "read", "swap")
    final, records = loop.execute(state, word)
    answer = next(record[1] for record in records if record[0] == "read")
    return answer, final, len(word)


class InternalPredictionMemoryTests(unittest.TestCase):
    def test_optimal_classical_memory_curve_for_success_and_probability_loss(self):
        expected = {1: (F(1, 2), F(1, 4)), 2: (F(3, 4), F(1, 8)),
                    3: (F(7, 8), F(1, 16)), 4: (F(1), F(0))}
        for messages, pair in expected.items():
            self.assertEqual(optimize_codes(messages)[:2], pair)

    def test_one_bit_upper_bound_matches_independent_encoder_decoder_enumeration(self):
        for weight in (F(0), F(1, 5), F(1, 2), F(4, 5), F(1)):
            best = F(0)
            for code in product((0, 1), repeat=4):
                for decoder in product((0, 1), repeat=4):
                    score = sum((F(1, 4)*query_weight
                                 for state, message in zip(INPUTS, code)
                                 for query, query_weight in enumerate((weight, 1-weight))
                                 if decoder[2*message+query] == state[query]), F(0))
                    best = max(best, score)
            self.assertEqual(best, optimize_codes(2, weight)[0])
            self.assertEqual(best, F(1, 2)+max(weight, 1-weight)/2)

    def test_weighted_probability_loss_and_three_label_accuracy(self):
        for numerator in range(11):
            weight = F(numerator, 10)
            minor = min(weight, 1-weight)
            self.assertEqual(optimize_codes(2, weight)[1], minor/4)
            three = optimize_codes(3, weight)
            self.assertEqual(three[0], 1-minor/4)
            self.assertEqual(three[1], minor/8)

    def test_exact_all_query_prediction_requires_four_distinct_messages(self):
        for messages in range(1, 5):
            for code in product(range(messages), repeat=4):
                exact = code_statistics(code)[0] == 1
                self.assertEqual(exact, len(set(code)) == 4)
        self.assertEqual(max(whole_pair_success(code) for code in product((0, 1), repeat=4)), F(1, 2))

    def test_all_bayesian_forecasts_are_calibrated_under_declared_source(self):
        for code in product((0, 1), repeat=4):
            _, _, forecasts = code_statistics(code)
            bins = defaultdict(lambda: [F(0), F(0)])
            for state, message in zip(INPUTS, code):
                for query in (0, 1):
                    forecast = forecasts[(message, query)]
                    bins[forecast][0] += F(1, 8)
                    bins[forecast][1] += F(state[query], 8)
            for forecast, (mass, one_mass) in bins.items():
                self.assertEqual(one_mass/mass, forecast)

    def test_early_query_information_or_post_query_read_changes_the_task(self):
        for query in (0, 1):
            early_code = {state: state[query] for state in INPUTS}
            self.assertEqual(len(set(early_code.values())), 2)
            for state in INPUTS:
                self.assertEqual(early_code[state], state[query])
                answer, restored, cost = reread_with_restoration(state, query)
                self.assertEqual(answer, state[query])
                self.assertEqual(restored, state)
                self.assertEqual(cost, 1 if query == 0 else 3)

    def test_accessible_audit_log_contains_the_missing_input_information(self):
        for state in INPUTS:
            final, log = loop.execute(state, ("read", "swap", "read", "swap"))
            self.assertEqual(final, state)
            logged_state = (log[0][1], log[2][1])
            self.assertEqual(logged_state, state)
        self.assertEqual(len({(s[0], s[1]) for s in INPUTS}), 4)

    def test_equal_guess_accuracy_can_hide_different_probability_forecast_quality(self):
        stored_first = (0, 0, 1, 1)
        singleton = (0, 1, 1, 1)
        first, second = code_statistics(stored_first), code_statistics(singleton)
        self.assertEqual(first[0], F(3, 4))
        self.assertEqual(second[0], F(3, 4))
        self.assertEqual(first[1], F(1, 8))
        self.assertEqual(second[1], F(1, 6))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checked = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(InternalPredictionMemoryTests))
    if not checked.wasSuccessful():
        raise SystemExit(1)
    result = {
        "rounds": [200, 202],
        "scope": "Uniform classical two-bit input, query selected independently after compression, one answer; no later source/log access",
        "internal_memory_is_distinct_from_original_plant_registers": True,
        "classical_curve": [
            {"memory_labels": m, "optimal_mean_success": str(optimize_codes(m)[0]),
             "optimal_mean_brier": str(optimize_codes(m)[1]),
             "success_encoder": optimize_codes(m)[2], "brier_encoder": optimize_codes(m)[3]}
            for m in range(1, 5)],
        "weighted_one_bit_success": "1/2 + max(w,1-w)/2",
        "weighted_one_bit_brier": "min(w,1-w)/4",
        "exact_prediction_required_memory_labels": 4,
        "one_bit_full_pair_retrieval_max_success": "1/2",
        "early_choice_or_fresh_read_can_make_prediction_exact": True,
        "audit_log_is_free_input_dependent_memory": False,
        "calibration_implies_per_trial_certainty": False,
        "quantum_structure_derived": False,
        "automated_tests": checked.testsRun,
    }
    if args.write_results:
        Path(__file__).with_name("internal_prediction_memory_results.json").write_text(
            json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"tests": checked.testsRun, "rounds": result["rounds"]}))


if __name__ == "__main__":
    main()
