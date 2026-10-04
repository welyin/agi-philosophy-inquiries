"""Round 199: sharp all-policy distance when memory is kept and noise is rounded.

Budget counts tick uses; other ideal controls are permitted. The analytical
upper bound follows from pre-sampling iid noise and data processing.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
from math import comb
from pathlib import Path
import unittest

import history_process_audit as process


def innovation_distance(ticks, first, second):
    if ticks < 0 or not (0 <= first <= 1 and 0 <= second <= 1):
        raise ValueError("Invalid tick count or probabilities.")
    return sum((comb(ticks, k)*abs(first**k*(1-first)**(ticks-k)
                -second**k*(1-second)**(ticks-k)) for k in range(ticks+1)), F(0))/2


def extraction_word(ticks, known_zero=False):
    if ticks < 0:
        raise ValueError("Tick budget cannot be negative.")
    return (() if known_zero else ("read", "swap", "read"))+("tick", "swap", "read")*ticks


def decode_innovations(records, known_zero=False):
    if known_zero:
        x, memory, offset = 0, 0, 0
    else:
        if records[0][0] != "read" or records[2][0] != "read":
            raise ValueError("Missing initial tomography records.")
        # After read, swap, read: current x=old m, current m=old x.
        x, memory, offset = records[2][1], records[0][1], 3
    result = []
    for start in range(offset, len(records), 3):
        expected = (("done", "tick"), ("done", "swap"))
        if records[start:start+2] != expected or records[start+2][0] != "read":
            raise ValueError("Records do not implement the extraction protocol.")
        following_x = records[start+2][1]
        result.append(following_x ^ x ^ memory)
        x = following_x
    return tuple(result)


def observed_distance(initial, word, first, second):
    policy = lambda history: word[len(history)]
    records = [process.transcript_distribution(initial, policy, len(word), eta)
               for eta in (first, second)]
    return process.total_variation(*records)


def apply_vector(weights, action, noise):
    results = {}
    for state, mass in zip(process.STATES, weights):
        if mass:
            for (successor, record), probability in process.branches(state, action, noise).items():
                target = results.setdefault(record, [F(0)]*4)
                target[process.STATES.index(successor)] += mass*probability
    return {record: tuple(values) for record, values in results.items()}


def optimal_adaptive_distance(initial, commands, ticks, first, second):
    """Exact small finite decision-tree audit, with all outcome-conditioned choices.

    Unnormalized weights keep branch probability differences. Randomizing action
    choices cannot improve a maximum of this linear expected discrimination score.
    """
    zero = (F(0),)*4
    @lru_cache(None)
    def solve(p, q, remaining, budget):
        if remaining == 0:
            return abs(sum(p)-sum(q))/2
        best = F(0)
        for action in process.ACTIONS:
            if action == "tick" and budget == 0:
                continue
            p_out, q_out = apply_vector(p, action, first), apply_vector(q, action, second)
            value = sum((solve(p_out.get(r, zero), q_out.get(r, zero), remaining-1,
                               budget-int(action == "tick"))
                         for r in p_out.keys() | q_out.keys()), F(0))
            best = max(best, value)
        return best
    weights = tuple(initial.get(s, F(0)) for s in process.STATES)
    result = solve(weights, weights, commands, ticks)
    return result, solve.cache_info().currsize


class AdaptiveNoiseDistanceTests(unittest.TestCase):
    def test_binomial_grouping_matches_full_noise_word_total_variation(self):
        for ticks in range(7):
            first, second = F(1, 5), F(49, 100)
            distributions = []
            for eta in (first, second):
                distributions.append({word: eta**sum(word)*(1-eta)**(ticks-sum(word))
                                      for word in product((0, 1), repeat=ticks)})
            self.assertEqual(process.total_variation(*distributions),
                             innovation_distance(ticks, first, second))

    def test_extraction_recovers_all_hidden_noise_words_without_noise_access(self):
        for state in process.STATES:
            for bits in product((0, 1), repeat=3):
                word = extraction_word(3)
                supplied_bits = (0, 0, 0)+tuple(bit for e in bits for bit in (e, 0, 0))
                _, records = process.execute(state, word, supplied_bits)
                self.assertEqual(decode_innovations(records), bits)

    def test_common_initial_distribution_does_not_prevent_sharp_bound_attainment(self):
        preparations = ({s: F(1, 4) for s in process.STATES},
                        {s: F(i+1, 10) for i, s in enumerate(process.STATES)})
        for initial in preparations:
            for ticks in range(5):
                for first, second in ((F(0), F(1, 2)), (F(1, 5), F(1, 2)),
                                      (F(49, 100), F(1, 2))):
                    actual = observed_distance(initial, extraction_word(ticks), first, second)
                    self.assertEqual(actual, innovation_distance(ticks, first, second))

    def test_bounded_adaptive_decision_tree_agrees_with_analytic_optimum(self):
        initial = {(0, 0): F(1)}
        for ticks in (1, 2):
            actual, states_checked = optimal_adaptive_distance(initial, 3*ticks, ticks,
                                                               F(1, 5), F(1, 2))
            self.assertEqual(actual, innovation_distance(ticks, F(1, 5), F(1, 2)))
            self.assertGreater(states_checked, 1)

    def test_record_dependent_tick_allocation_respects_noise_budget(self):
        initial = {s: F(1, 4) for s in process.STATES}
        def policy(history):
            used = sum(record == ("done", "tick") for record in history)
            if len(history) % 2 == 0:
                return "read"
            return "tick" if used < 2 and history[-1][-1] else "swap"
        records = [process.transcript_distribution(initial, policy, 9, eta)
                   for eta in (F(1, 5), F(1, 2))]
        self.assertLessEqual(process.total_variation(*records),
                             innovation_distance(2, F(1, 5), F(1, 2)))
        for distribution in records:
            self.assertTrue(all(sum(r == ("done", "tick") for r in h) <= 2
                                for h in distribution))

    def test_noiseless_echo_cannot_distinguish_actual_noise_parameters(self):
        initial = {s: F(1, 4) for s in process.STATES}
        echo = ("read", "swap", "read", "swap", "read")
        self.assertEqual(observed_distance(initial, echo, F(0), F(1, 2)), 0)
        for ticks in range(12):
            self.assertEqual(innovation_distance(ticks, F(1, 5), F(1, 5)), 0)

    def test_sharp_noise_distance_is_below_uniform_adaptive_coupling_bound(self):
        for ticks in range(21):
            for first, second in ((F(0), F(1, 2)), (F(1, 5), F(1, 2)),
                                  (F(49, 100), F(1, 2))):
                sharp = innovation_distance(ticks, first, second)
                delta = abs(first-second)
                self.assertLessEqual(sharp, 1-(1-delta)**ticks)
                self.assertLessEqual(1-(1-delta)**ticks, min(1, ticks*delta))
        self.assertEqual(innovation_distance(4, F(0), F(1, 2)), F(15, 16))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checked = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(AdaptiveNoiseDistanceTests))
    if not checked.wasSuccessful():
        raise SystemExit(1)
    result = {
        "round": 199,
        "scope": "Two full four-state processes with same preparation and controls; only iid tick-noise probabilities differ",
        "resource_budget": "At most k ticks, with finite noiseless read/swap/flip controls; no separate total-command optimum claimed",
        "all_policy_sharp_tv": "0.5 sum_j binom(k,j)|eta^j(1-eta)^(k-j)-zeta^j(1-zeta)^(k-j)|",
        "uniform_coupling_upper_bound": "1-(1-|eta-zeta|)^k <= min(1,k|eta-zeta|)",
        "attaining_protocol_uniform_or_arbitrary_common_preparation": "read,swap,read; then (tick,swap,read)^k",
        "total_commands_in_attaining_protocol": "3+3*k, or 3*k for known (0,0)",
        "environment_noise_directly_readable": False,
        "memory_erasure_equivalent_to_noise_rounding": False,
        "geometry_changed": False,
        "samples": [{"eta": str(eta), "zeta": "1/2", "tick_budget": ticks,
                     "tv_exact": str(innovation_distance(ticks, eta, F(1, 2))),
                     "tv_float": float(innovation_distance(ticks, eta, F(1, 2)))}
                    for eta in (F(0), F(1, 5), F(49, 100)) for ticks in (1, 2, 10, 100)],
        "automated_tests": checked.testsRun,
    }
    if args.write_results:
        Path(__file__).with_name("adaptive_noise_distance_results.json").write_text(
            json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"tests": checked.testsRun, "round": 199}))


if __name__ == "__main__":
    main()
