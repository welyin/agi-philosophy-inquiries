"""Rounds 197-198: active access invalidates a passive memoryless projection.

All probabilities are exact Fractions. A reduced state is only x; giving its
transition kernel an old transcript would supply additional memory.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import unittest

import history_process_audit as full


UNIFORM = {s: F(1, 4) for s in full.STATES}


def project(initial):
    out = defaultdict(F)
    for (x, _), probability in initial.items():
        out[x] += probability
    return dict(out)


def reduced_branches(x, action):
    if x not in (0, 1) or action not in full.ACTIONS:
        raise ValueError("Invalid reduced state or action.")
    if action == "read":
        return {(x, ("read", x)): F(1)}
    record = ("done", action)
    if action == "flip_memory":
        return {(x, record): F(1)}
    # Stationary averaging forgets m before each swap/tick; it resamples a fair x.
    return {(successor, record): F(1, 2) for successor in (0, 1)}


def reduced_transcripts(initial, policy, horizon):
    active = {(x, ()): p for x, p in initial.items() if p}
    for _ in range(horizon):
        following = defaultdict(F)
        for (x, history), probability in active.items():
            for (successor, record), weight in reduced_branches(x, policy(history)).items():
                following[(successor, history+(record,))] += probability*weight
        active = following
    result = defaultdict(F)
    for (_, history), probability in active.items():
        result[history] += probability
    return dict(result)


def word_policy(word):
    return lambda history: word[len(history)]


def word_error(initial, word, noise=F(1, 2)):
    policy = word_policy(word)
    original = full.transcript_distribution(initial, policy, len(word), noise)
    reduced = reduced_transcripts(project(initial), policy, len(word))
    return full.total_variation(original, reduced)


def projected_kernel(state, action, noise):
    out = defaultdict(F)
    for (successor, record), p in full.branches(state, action, noise).items():
        out[(successor[0], record)] += p
    return dict(out)


def local_defects(noise):
    return {a: max(full.total_variation(projected_kernel(s, a, noise),
                                        reduced_branches(s[0], a)) for s in full.STATES)
            for a in full.ACTIONS}


def coupling_bound(initial_error, defects):
    survival = 1-initial_error
    for defect in defects:
        if not 0 <= defect <= 1:
            raise ValueError("Defects must be probabilities.")
        survival *= 1-defect
    return 1-survival


class ActiveMemoryProjectionTests(unittest.TestCase):
    def test_passive_complete_records_match_at_half_noise_but_echo_does_not(self):
        for readings in range(1, 7):
            word = ("read",)+("tick", "read")*(readings-1)
            self.assertEqual(word_error(UNIFORM, word), 0)
        echo = ("read", "swap", "read", "swap", "read")
        for noise in (F(0), F(1, 5), F(49, 100), F(1, 2)):
            self.assertEqual(word_error(UNIFORM, echo, noise), F(1, 2))

    def test_stationary_one_step_averaging_is_exact_for_every_action(self):
        for action in full.ACTIONS:
            for x in (0, 1):
                averaged = defaultdict(F)
                for memory in (0, 1):
                    for key, probability in projected_kernel((x, memory), action, F(1, 5)).items():
                        averaged[key] += probability/2
                self.assertEqual(dict(averaged), reduced_branches(x, action))

    def test_echo_error_grows_exactly_and_known_preparation_saturates_coupling_bound(self):
        for swaps in range(1, 8):
            word = ("read",)+("swap", "read")*swaps
            self.assertEqual(word_error(UNIFORM, word), 1-F(1, 2**(swaps-1)))
            error = word_error({(0, 0): F(1)}, word)
            self.assertEqual(error, 1-F(1, 2**swaps))
            self.assertEqual(error, coupling_bound(F(0), [F(1, 2)]*swaps))

    def test_same_current_state_has_future_distance_one_for_every_noise(self):
        word = ("swap", "read")
        for x in (0, 1):
            for noise in (F(0), F(49, 100), F(1, 2)):
                distributions = [full.transcript_distribution({(x, m): F(1)},
                                 word_policy(word), len(word), noise) for m in (0, 1)]
                self.assertEqual(full.total_variation(*distributions), 1)

    def test_read_swap_read_identifies_every_preparation_and_mixture(self):
        word = ("read", "swap", "read")
        images = {}
        for state in full.STATES:
            records = full.transcript_distribution({state: F(1)}, word_policy(word), 3)
            self.assertEqual(len(records), 1)
            images[state] = next(iter(records))
        self.assertEqual(len(set(images.values())), 4)
        mixture = {state: F(i+1, 10) for i, state in enumerate(full.STATES)}
        records = full.transcript_distribution(mixture, word_policy(word), 3)
        self.assertEqual(records, {images[state]: p for state, p in mixture.items()})

    def test_uniform_pointwise_kernel_defects_do_not_shrink_with_tick_noise(self):
        expected = {"read": F(0), "swap": F(1, 2), "flip_memory": F(0), "tick": F(1, 2)}
        for noise in (F(0), F(1, 5), F(49, 100), F(1, 2)):
            self.assertEqual(local_defects(noise), expected)
        # Minimax over any binary successor law q is >= 1/2, attained by q=1/2.
        for numerator in range(11):
            q = F(numerator, 10)
            self.assertGreaterEqual(max(q, 1-q), F(1, 2))

    def test_short_words_and_feedback_respect_uniform_coupling_contract(self):
        noise = F(1, 5)
        defects = local_defects(noise)
        initial = {state: F(i+1, 10) for i, state in enumerate(full.STATES)}
        for size in range(1, 5):
            for word in product(full.ACTIONS, repeat=size):
                self.assertLessEqual(word_error(initial, word, noise),
                                     coupling_bound(F(0), [defects[a] for a in word]))
        def policy(history):
            if len(history) % 2 == 0:
                return "read"
            return "swap" if history[-1][-1] else "flip_memory"
        original = full.transcript_distribution(initial, policy, 9, noise)
        reduced = reduced_transcripts(project(initial), policy, 9)
        # At most four potentially defective actions, including either branch.
        self.assertLessEqual(full.total_variation(original, reduced), 1-F(1, 16))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checked = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(ActiveMemoryProjectionTests))
    if not checked.wasSuccessful():
        raise SystemExit(1)
    result = {
        "rounds": [197, 198],
        "scope": "Same two-bit classical interface, ideal swap/read, noise only on tick",
        "passive_half_noise_uniform_source_error": "0 at all finite passive horizons",
        "five_command_echo_error": "1/2 for all tick noise parameters; no tick used",
        "uniform_echo_with_k_swaps_error": "1-2^(-(k-1)), k>=1",
        "known_state_echo_with_k_swaps_error": "1-2^(-k), k>=0",
        "exact_active_deterministic_labels": 4,
        "exact_active_ensemble_affine_dimension": 3,
        "hilbert_dimension_lower_bound_for_exact_four_label_task": 4,
        "one_visible_bit_all_preparation_worst_error_lower_bound": "1/2",
        "per_action_pointwise_defects": {a: str(d) for a, d in local_defects(F(1, 2)).items()},
        "adaptive_contract": "TV <= 1-(1-delta_initial)*product_j(1-delta_j), for uniform conditional projected outcome-successor defects",
        "stationary_average_error_is_sufficient": False,
        "erased_information_recreated": False,
        "quantum_necessity_derived": False,
        "automated_tests": checked.testsRun,
    }
    if args.write_results:
        Path(__file__).with_name("active_memory_projection_results.json").write_text(
            json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"tests": checked.testsRun, "rounds": result["rounds"]}))


if __name__ == "__main__":
    main()
