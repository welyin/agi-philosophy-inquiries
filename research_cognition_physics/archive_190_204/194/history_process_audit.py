"""Rounds 194/196: reversible classical memory and stationary temporal records.

Selected history/rollback primitives, not a complete implementation of SoCA.
Exact probabilities use Fraction. The environment noise record is not public.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import json
from math import comb
from pathlib import Path
import unittest


STATES = tuple(product((0, 1), repeat=2))  # visible register x, retained register m
ACTIONS = ("read", "swap", "flip_memory", "tick")


def step(state, action, noise_bit=0):
    x, memory = state
    if state not in STATES or action not in ACTIONS or noise_bit not in (0, 1):
        raise ValueError("Invalid state, action or binary noise.")
    if action == "read":
        return state, ("read", x)
    if action == "swap":
        successor = (memory, x)
    elif action == "flip_memory":
        successor = (x, 1-memory)
    else:
        successor = (memory, x ^ memory ^ noise_bit)
    return successor, ("done", action)


def branches(state, action, noise=F(0)):
    if not 0 <= noise <= 1:
        raise ValueError("Noise must be a probability.")
    output = defaultdict(F)
    for bit, weight in ((0, 1-noise), (1, noise)):
        if weight:
            successor, record = step(state, action, bit)
            output[(successor, record)] += weight
    return dict(output)


def execute(state, actions, noise_bits=None):
    actions = tuple(actions)
    noise_bits = (0,)*len(actions) if noise_bits is None else tuple(noise_bits)
    if len(noise_bits) != len(actions):
        raise ValueError("One supplied noise bit per command is required.")
    history = []
    for action, bit in zip(actions, noise_bits):
        state, record = step(state, action, bit)
        history.append(record)
    return state, tuple(history)


def inverse_word(actions):
    # This reverses internal noiseless states; an append-only transcript remains.
    return tuple(part for action in reversed(actions)
                 for part in (("tick", "tick") if action == "tick" else (action,)))


def undo_noisy_tick(successor, noise_bit):
    x_new, memory_new = successor
    return (memory_new ^ x_new ^ noise_bit, x_new)


def minimal_partition(actions):
    labels = {state: 0 for state in STATES}
    while True:
        identities, updated = {}, {}
        for state in STATES:
            signature = (labels[state], tuple((record, labels[successor])
                         for successor, record in (step(state, a) for a in actions)))
            updated[state] = identities.setdefault(signature, len(identities))
        if updated == labels:
            return updated
        labels = updated


def transcript_distribution(initial, policy, horizon, noise=F(0)):
    active = {(state, ()): weight for state, weight in initial.items() if weight}
    for _ in range(horizon):
        following = defaultdict(F)
        for (state, history), probability in active.items():
            for (successor, record), weight in branches(state, policy(history), noise).items():
                following[(successor, history+(record,))] += probability*weight
        active = following
    result = defaultdict(F)
    for (_, history), probability in active.items():
        result[history] += probability
    return dict(result)


def passive_records(length, noise=F(0)):
    """Uniform initial state; read x, then tick between consecutive readings."""
    if length < 1:
        raise ValueError("At least one reading is required.")
    active = {(state, (state[0],)): F(1, 4) for state in STATES}
    for _ in range(length-1):
        following = defaultdict(F)
        for (state, history), probability in active.items():
            for (successor, _), weight in branches(state, "tick", noise).items():
                following[(successor, history+(successor[0],))] += probability*weight
        active = following
    result = defaultdict(F)
    for (_, history), probability in active.items():
        result[history] += probability
    return dict(result)


def total_variation(first, second):
    return sum((abs(first.get(key, 0)-second.get(key, 0))
                for key in first.keys() | second.keys()), F(0))/2


def memoryless_error(length, noise):
    """Exact TV from the unique stationary one-step fitted fair-bit chain."""
    if length < 1 or not 0 <= noise <= 1:
        raise ValueError("Invalid length or probability.")
    innovations = max(0, length-2)
    return sum((comb(innovations, k)*abs(noise**k*(1-noise)**(innovations-k)
                -F(1, 2**innovations)) for k in range(innovations+1)), F(0))/2


class HistoryProcessTests(unittest.TestCase):
    def test_future_queries_require_four_classes_and_read_only_requires_two(self):
        self.assertEqual(len(set(minimal_partition(("read", "flip_memory")).values())), 2)
        self.assertEqual(len(set(minimal_partition(ACTIONS).values())), 4)
        readouts = {execute(s, ("read", "swap", "read"))[1] for s in STATES}
        self.assertEqual(len(readouts), 4)

    def test_same_current_reading_different_history_and_order_are_classical(self):
        base = (0, 0)
        changed = execute(base, ("flip_memory",))[0]
        self.assertEqual(step(base, "read")[1], step(changed, "read")[1])
        self.assertNotEqual(execute(base, ("swap", "read"))[1],
                            execute(changed, ("swap", "read"))[1])
        self.assertEqual(execute(base, ("swap", "flip_memory"))[0], (0, 1))
        self.assertEqual(execute(base, ("flip_memory", "swap"))[0], (1, 0))

    def test_all_short_noiseless_words_undo_and_replay_internal_states(self):
        for size in range(5):
            for word in product(ACTIONS, repeat=size):
                for state in STATES:
                    final, history = execute(state, word)
                    restored, undo_history = execute(final, inverse_word(word))
                    self.assertEqual(restored, state)
                    self.assertEqual(execute(restored, word), (final, history))
                    # Undo does not erase the accumulated external audit trail.
                    self.assertGreaterEqual(len(history+undo_history), len(history))

    def test_unknown_noise_requires_side_information_for_universal_undo(self):
        for state in STATES:
            for bit in (0, 1):
                successor, _ = step(state, "tick", bit)
                self.assertEqual(undo_noisy_tick(successor, bit), state)
        for successor in STATES:
            # Two different input/noise pairs have the same output and public ack.
            first, second = (undo_noisy_tick(successor, b) for b in (0, 1))
            self.assertNotEqual(first, second)
            self.assertEqual(step(first, "tick", 0), step(second, "tick", 1))
        word, seed = ("tick", "read", "tick", "read"), (1, 0, 0, 0)
        initial_run = execute((0, 0), word, seed)
        self.assertEqual(initial_run, execute((0, 0), word, seed))
        self.assertNotEqual(initial_run, execute((0, 0), word, (0, 0, 0, 0)))

    def test_uniform_full_state_is_stationary_but_triples_retain_memory(self):
        for noise in (F(0), F(1, 5), F(49, 100), F(1, 2)):
            stationary = defaultdict(F)
            for state in STATES:
                for (successor, _), weight in branches(state, "tick", noise).items():
                    stationary[successor] += weight/4
            self.assertEqual(dict(stationary), {s: F(1, 4) for s in STATES})
            self.assertEqual(passive_records(2, noise), {r: F(1, 4) for r in STATES})
            for record, probability in passive_records(3, noise).items():
                parity = record[0] ^ record[1] ^ record[2]
                self.assertEqual(probability, (noise if parity else 1-noise)/4)

    def test_record_formula_and_exact_tv_agree_with_independent_path_enumeration(self):
        for noise in (F(0), F(1, 5), F(49, 100), F(1, 2)):
            for length in range(2, 9):
                actual = passive_records(length, noise)
                fair = {r: F(1, 2**length) for r in product((0, 1), repeat=length)}
                for record in fair:
                    k = sum(record[j] ^ record[j+1] ^ record[j+2]
                            for j in range(length-2))
                    predicted = F(1, 4)*noise**k*(1-noise)**(length-2-k)
                    self.assertEqual(actual.get(record, 0), predicted)
                self.assertEqual(total_variation(actual, fair), memoryless_error(length, noise))

    def test_memoryless_error_endpoints_and_horizon_bounds(self):
        for length in range(2, 35):
            self.assertEqual(memoryless_error(length, F(0)), 1-F(1, 2**(length-2)))
            self.assertEqual(memoryless_error(length, F(1, 2)), 0)
            for noise in (F(1, 5), F(49, 100)):
                error = memoryless_error(length, noise)
                self.assertLessEqual(error, min(1, (length-2)*abs(F(1, 2)-noise)))
                self.assertLessEqual(error, memoryless_error(length+1, noise))

    def test_exact_second_order_prediction_and_failure_of_first_order_when_biased(self):
        noise = F(1, 5)
        triple = passive_records(3, noise)
        for first, current in STATES:
            predicted_next = first ^ current
            self.assertEqual(triple[(first, current, predicted_next)]/F(1, 4), 1-noise)
        for current in (0, 1):
            next_one = sum(p for r, p in triple.items() if r[1:] == (current, 1))
            self.assertEqual(next_one/F(1, 2), F(1, 2))
        self.assertEqual(memoryless_error(3, noise), F(3, 10))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checked = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(HistoryProcessTests))
    if not checked.wasSuccessful():
        raise SystemExit(1)
    result = {
        "rounds": [194, 196],
        "scope": "Chosen two-bit classical memory interface; not the full SoCA architecture",
        "minimal_deterministic_labels": {"read_and_flip_memory": 2, "full_interface": 4},
        "noise_is_public": False,
        "internal_undo_erases_external_history": False,
        "stationary_process": "Uniform (x,m); tick=(m,x xor m xor e); e iid Bernoulli(eta)",
        "memoryless_comparison": "Unique stationary first-order chain matching all adjacent pairs; passive read/tick task only",
        "exact_tv": "0.5 sum_{k=0}^{n-2} binom(n-2,k) |eta^k(1-eta)^(n-2-k) - 2^(-(n-2))|, n>=2",
        "fixed_finite_horizon_eta_to_half_error": "0",
        "fixed_eta_not_half_infinite_horizon_error": "1 (analytic proof in note 196)",
        "horizon_samples": [
            {"eta": str(eta), "readings": n, "tv_exact": str(memoryless_error(n, eta)),
             "tv_float": float(memoryless_error(n, eta))}
            for eta in (F(0), F(1, 5), F(49, 100), F(1, 2)) for n in (3, 8, 32, 102)
        ],
        "automated_tests": checked.testsRun,
        "quantum_structure_derived": False,
    }
    if args.write_results:
        Path(__file__).with_name("history_process_audit_results.json").write_text(
            json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"tests": checked.testsRun, "rounds": result["rounds"]}))


if __name__ == "__main__":
    main()
