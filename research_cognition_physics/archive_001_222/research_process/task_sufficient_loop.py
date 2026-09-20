"""Round 191: derive an exact task quotient of a small classical cognitive loop.

This implements selected SoCA primitives, not the full SoCA architecture.
The external transcript is retained; the finite state governs future tasks.
"""

import argparse
from collections import defaultdict
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import unittest


@dataclass(frozen=True)
class TaskState:
    target: int
    sensor_flip: int
    memory: int
    goal: int
    calibrated: int


STATES = tuple(TaskState(*bits) for bits in product((0, 1), repeat=5))
UNANCHORED = tuple(state for state in STATES if not state.calibrated)
RELATIVE_ACTIONS = ("audit", "sense", "flip", "seek", "goal_0", "goal_1", "seek_absolute")
ANCHORED_ACTIONS = RELATIVE_ACTIONS + ("calibrate",)


def step(state, command, noise_bit=0):
    """The controller uses the true sensor flip only with a retained certificate.

    Calibration uses a separately supplied known-zero probe, leaves the target
    untouched, and retains the measured sensor flip. The sensor does not drift,
    so the retained calibration value equals sensor_flip when calibrated=1.
    """
    if command not in ANCHORED_ACTIONS or noise_bit not in (0, 1):
        raise ValueError("Unknown command or nonbinary readout noise.")
    if command == "audit":
        certificate = state.sensor_flip if state.calibrated else None
        return state, ("audit", state.memory, state.goal, state.calibrated, certificate)
    if command.startswith("goal_"):
        value = int(command[-1])
        return replace(state, goal=value), ("goal", value)
    if command == "calibrate":
        return replace(state, calibrated=1), ("calibrate", "known_zero_probe", state.sensor_flip)
    if command == "seek_absolute" and not state.calibrated:
        return state, ("blocked", "no_signed_calibration")
    action = {"sense": 0, "flip": 1, "seek": state.memory ^ state.goal,
              "seek_absolute": state.memory ^ state.sensor_flip ^ state.goal}[command]
    prediction = state.memory ^ action
    target = state.target ^ action
    observation = target ^ state.sensor_flip ^ noise_bit
    record = (command, action, prediction, observation, prediction ^ observation)
    return replace(state, target=target, memory=observation), record


def branches(state, command, noise=F(0)):
    if not 0 <= noise <= 1:
        raise ValueError("Noise must be a probability.")
    result = defaultdict(F)
    for bit, probability in ((0, 1 - noise), (1, noise)):
        if probability:
            successor, record = step(state, command, bit)
            result[(successor, record)] += probability
    return dict(result)


def predictive_key(state):
    if state.calibrated:
        return (1, state.target, state.sensor_flip, state.memory, state.goal)
    return (0, state.target ^ state.sensor_flip, state.memory, state.goal)


def minimal_partition(states, commands):
    """Deterministic output/successor refinement; no geometric assumptions."""
    labels = {state: 0 for state in states}
    refinements = []
    while True:
        ids = {}
        updated = {}
        for state in states:
            signature = (labels[state], tuple((record, labels[successor])
                         for successor, record in (step(state, command) for command in commands)))
            updated[state] = ids.setdefault(signature, len(ids))
        refinements.append(len(ids))
        if updated == labels:
            return updated, refinements
        labels = updated


def projected_branches(state, command, noise=F(0)):
    result = defaultdict(F)
    for (successor, record), probability in branches(state, command, noise).items():
        result[(predictive_key(successor), record)] += probability
    return dict(result)


def transcript_distribution(initial, policy, horizon, noise=F(0), quotient=False):
    representatives = {predictive_key(state): state for state in STATES}
    active = defaultdict(F)
    for state, probability in initial.items():
        key = predictive_key(state) if quotient else state
        active[(key, ())] += probability
    for _ in range(horizon):
        following = defaultdict(F)
        for (state_or_key, history), probability in active.items():
            state = representatives[state_or_key] if quotient else state_or_key
            command = policy(history)
            if quotient and command not in RELATIVE_ACTIONS:
                raise ValueError("This quotient has not been proved for signed calibration tasks.")
            for (successor, record), weight in branches(state, command, noise).items():
                key = predictive_key(successor) if quotient else successor
                following[(key, history + (record,))] += probability * weight
        active = following
    records = defaultdict(F)
    for (_, history), probability in active.items():
        records[history] += probability
    return dict(records)


def readout_word(state, commands):
    records = []
    for command in commands:
        state, record = step(state, command)
        records.append(record)
    return tuple(records)


class TaskSufficientLoopTests(unittest.TestCase):
    def test_partition_refinement_derives_eight_twenty_four_and_thirty_two_classes(self):
        for states, commands, expected in ((UNANCHORED, RELATIVE_ACTIONS, 8),
                                           (STATES, RELATIVE_ACTIONS, 24),
                                           (STATES, ANCHORED_ACTIONS, 32)):
            labels, _ = minimal_partition(states, commands)
            self.assertEqual(len(set(labels.values())), expected)
            if commands == RELATIVE_ACTIONS:
                for first in states:
                    for second in states:
                        self.assertEqual(labels[first] == labels[second],
                                         predictive_key(first) == predictive_key(second))

    def test_exact_outcome_and_successor_lumping_survives_fixed_sensor_noise(self):
        for first in UNANCHORED:
            second = replace(first, target=1-first.target, sensor_flip=1-first.sensor_flip)
            for command in RELATIVE_ACTIONS:
                for noise in (F(0), F(1, 5), F(1, 2)):
                    self.assertEqual(projected_branches(first, command, noise),
                                     projected_branches(second, command, noise))

    def test_prediction_error_is_recorded_before_feedback_overwrites_memory(self):
        state = TaskState(0, 0, 1, 1, 0)
        successor, record = step(state, "sense")
        self.assertEqual(record, ("sense", 0, 1, 0, 1))
        self.assertEqual(successor.memory, 0)
        self.assertEqual(state.memory, 1)
        self.assertEqual(step(successor, "seek")[1][1], 1)

    def test_parity_alone_discards_memory_and_goal_needed_for_future_actions(self):
        base = TaskState(0, 0, 0, 0, 0)
        changed_memory = replace(base, memory=1)
        changed_goal = replace(base, goal=1)
        self.assertEqual(step(base, "seek")[1][1], 0)
        self.assertEqual(step(changed_memory, "seek")[1][1], 1)
        self.assertEqual(step(changed_goal, "seek")[1][1], 1)

    def test_calibration_splits_old_classes_and_absolute_action_requires_provenance(self):
        first = TaskState(0, 0, 0, 0, 0)
        second = TaskState(1, 1, 0, 0, 0)
        self.assertEqual(predictive_key(first), predictive_key(second))
        self.assertEqual(step(first, "seek_absolute")[1], ("blocked", "no_signed_calibration"))
        a, record_a = step(first, "calibrate")
        b, record_b = step(second, "calibrate")
        self.assertNotEqual(record_a, record_b)
        self.assertNotEqual(predictive_key(a), predictive_key(b))
        self.assertEqual((a.target, b.target), (first.target, second.target))
        self.assertEqual(step(a, "seek_absolute")[0].target, 0)
        self.assertEqual(step(b, "seek_absolute")[0].target, 0)

    def test_short_readout_words_distinguish_all_quotient_labels_and_their_mixtures(self):
        relative_records = {predictive_key(s): readout_word(s, ("audit", "sense")) for s in STATES}
        self.assertEqual(len(set(relative_records.values())), 24)
        unanchored_records = {predictive_key(s): readout_word(s, ("audit", "sense")) for s in UNANCHORED}
        self.assertEqual(len(set(unanchored_records.values())), 8)
        anchored_records = [readout_word(s, ("audit", "calibrate", "sense")) for s in STATES]
        self.assertEqual(len(set(anchored_records)), 32)

    def test_full_and_quotient_models_match_complete_adaptive_history_probabilities(self):
        initial = {state: F(index + 1, 528) for index, state in enumerate(STATES)}
        policies = (
            lambda history: ("audit", "sense", "seek", "goal_1", "seek_absolute")[len(history)],
            lambda history: "sense" if not history else
                ("flip" if history[-1][-1] else "seek"),
        )
        for policy in policies:
            for noise in (F(0), F(1, 5)):
                full = transcript_distribution(initial, policy, 5, noise)
                projected = transcript_distribution(initial, policy, 5, noise, quotient=True)
                self.assertEqual(full, projected)
                self.assertEqual(sum(full.values()), 1)

    def test_individual_quotient_is_not_certified_for_new_cross_subject_reference_task(self):
        a = TaskState(0, 0, 0, 0, 0)
        a_flipped = TaskState(1, 1, 0, 0, 0)
        b = TaskState(0, 0, 0, 0, 0)
        self.assertEqual((predictive_key(a), predictive_key(b)),
                         (predictive_key(a_flipped), predictive_key(b)))
        # New declared resource: equal probe bits z reach both classical sensors.
        # Their output parity measures c_A xor c_B, even if z is unknown.
        for z in (0, 1):
            first = (z ^ a.sensor_flip) ^ (z ^ b.sensor_flip)
            second = (z ^ a_flipped.sensor_flip) ^ (z ^ b.sensor_flip)
            self.assertEqual((first, second), (0, 1))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(TaskSufficientLoopTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = {
        "round": 191,
        "scope": "Selected SoCA primitives in a chosen finite classical controlled loop, not a complete SoCA implementation",
        "state_bits": ["target", "sensor_flip", "last_raw_prediction_memory", "goal", "calibration_certificate_present"],
        "uncalibrated_states": 16,
        "relative_task_minimal_labels_uncalibrated": 8,
        "all_states": 32,
        "relative_task_minimal_labels_all_states": 24,
        "known_probe_task_minimal_labels_all_states": 32,
        "minimality_counts_assume_noiseless_target_readout": True,
        "noiseless_uncalibrated_operational_belief_space": "full classical simplex with 8 vertices and normalized dimension 7",
        "finite_tasks_and_fixed_noise_have_exact_outcome_conditioned_quotient": True,
        "quotient_remains_sufficient_after_arbitrary_new_joint_access": False,
        "kahler_or_quantum_structure_derived": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("task_sufficient_loop_results.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
