"""A classical baseline for selected SoCA primitives, without quantum inputs."""

import argparse
import json
import math
import platform
import unittest
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class LoopState:
    peer_value: int
    hypothesis: int

    def __post_init__(self):
        if self.peer_value not in (0, 1) or self.hypothesis not in (0, 1):
            raise ValueError("The peer value and retained hypothesis must be binary.")


@dataclass(frozen=True)
class Candidate:
    action: int
    predicted_observation: int
    predicted_target_loss: int


@dataclass(frozen=True)
class StepRecord:
    state_before: LoopState
    candidates: tuple
    requested_action: int
    executed_action: int
    vetoed: bool
    prediction: int
    observation: int
    prediction_error: int
    state_after: LoopState


STATES = tuple(LoopState(peer_value, hypothesis)
               for peer_value in (0, 1) for hypothesis in (0, 1))


def candidates_for(state):
    return tuple(Candidate(action, state.hypothesis ^ action,
                           int((state.hypothesis ^ action) != 1))
                 for action in (0, 1))


def loop_step(state, imposed_action=None, allow_toggle=True, observation_flip=False):
    if imposed_action is not None and imposed_action not in (0, 1):
        raise ValueError("An imposed action must be zero or one.")
    candidates = candidates_for(state)
    selected = min(candidates, key=lambda candidate: candidate.predicted_target_loss)
    requested = selected.action if imposed_action is None else imposed_action
    executed = 0 if requested == 1 and not allow_toggle else requested
    prediction = state.hypothesis ^ executed
    peer_after = state.peer_value ^ executed
    observation = peer_after ^ int(observation_flip)
    updated = LoopState(peer_after, observation)
    return StepRecord(state, candidates, requested, executed, requested != executed,
                      prediction, observation, int(prediction != observation), updated)


def state_index(state):
    return 2 * state.peer_value + state.hypothesis


def transition_kernel(imposed_action=None, observation_noise=0.0, allow_toggle=True):
    if not math.isfinite(observation_noise) or not 0 <= observation_noise <= 1:
        raise ValueError("Observation noise must be a probability.")
    kernel = np.zeros((len(STATES), len(STATES)))
    for column, state in enumerate(STATES):
        for observation_flip, probability in ((False, 1 - observation_noise),
                                               (True, observation_noise)):
            record = loop_step(state, imposed_action, allow_toggle, observation_flip)
            kernel[state_index(record.state_after), column] += probability
    return kernel


def entropy_nats(probabilities):
    probabilities = np.asarray(probabilities, dtype=float)
    if (probabilities.ndim != 1 or not np.all(np.isfinite(probabilities))
            or np.any(probabilities < 0) or not np.isclose(probabilities.sum(), 1)):
        raise ValueError("A normalized nonnegative distribution is required.")
    positive = probabilities[probabilities > 0]
    return float(-np.dot(positive, np.log(positive)))


def reference_trace():
    state = LoopState(0, 1)
    records = []
    for _ in range(3):
        record = loop_step(state)
        records.append(asdict(record))
        state = record.state_after
    return records


def baseline_results():
    autonomous = transition_kernel()
    noisy = transition_kernel(observation_noise=0.2)
    first = np.array([0.1, 0.2, 0.3, 0.4])
    second = np.array([0.4, 0.1, 0.2, 0.3])
    mixture_weight = 0.35
    mixture = mixture_weight * first + (1 - mixture_weight) * second
    mixture_error = noisy @ mixture - (mixture_weight * (noisy @ first)
                                       + (1 - mixture_weight) * (noisy @ second))
    before_update = np.array([0.5, 0.5, 0.0, 0.0])
    after_update = transition_kernel(imposed_action=0) @ before_update
    permutation = np.array([[0.0, 1.0], [1.0, 0.0]])
    mixed_update = 0.9 * np.eye(2) + 0.1 * permutation
    inverse = np.linalg.inv(mixed_update)
    return {
        "state_order": [asdict(state) for state in STATES],
        "autonomous_trace": reference_trace(),
        "veto_control": asdict(loop_step(LoopState(0, 0), allow_toggle=False)),
        "autonomous_kernel": autonomous.tolist(),
        "autonomous_kernel_rank": int(np.linalg.matrix_rank(autonomous)),
        "noisy_kernel": noisy.tolist(),
        "max_column_sum_error": float(np.max(abs(noisy.sum(axis=0) - 1))),
        "max_mixture_affinity_error": float(np.max(abs(mixture_error))),
        "memory_overwrite": {
            "before_distribution": before_update.tolist(),
            "after_distribution": after_update.tolist(),
            "before_retained_state_entropy_nats": entropy_nats(before_update),
            "after_retained_state_entropy_nats": entropy_nats(after_update),
            "scope": "Entropy of the retained peer-memory state only; an external audit record is excluded. No thermodynamic bath, heat, or entropy production is modeled.",
        },
        "finite_simplex_reversibility_control": {
            "forward_stochastic_matrix": mixed_update.tolist(),
            "algebraic_inverse": inverse.tolist(),
            "minimum_inverse_entry": float(inverse.min()),
            "scope": "This invertible forward stochastic matrix has no stochastic inverse. The finite-simplex permutation theorem is stated analytically in research_direction.md, not proved by this numerical example.",
        },
        "scope": "A chosen classical realization of prediction, two-candidate arbitration, a deterministic veto, action-feedback and memory overwrite. Not a complete SoCA implementation, a quantum reconstruction, or a proof that every cognitive model is classical.",
    }


class CognitiveLoopTests(unittest.TestCase):
    def test_candidates_are_scored_using_retained_hypothesis(self):
        candidates = candidates_for(LoopState(0, 0))
        self.assertEqual([candidate.predicted_observation for candidate in candidates], [0, 1])
        self.assertEqual([candidate.predicted_target_loss for candidate in candidates], [1, 0])
        self.assertEqual(candidates, candidates_for(LoopState(1, 0)))

    def test_prediction_is_compared_with_observation_before_memory_update(self):
        record = loop_step(LoopState(0, 1))
        self.assertEqual((record.prediction, record.observation, record.prediction_error), (1, 0, 1))
        self.assertEqual(record.state_after, LoopState(0, 0))
        self.assertEqual(record.state_before.hypothesis, 1)

    def test_feedback_changes_the_next_action(self):
        records = reference_trace()
        self.assertEqual([record["executed_action"] for record in records], [0, 1, 0])
        self.assertEqual([record["prediction_error"] for record in records], [1, 0, 0])
        self.assertEqual(records[-1]["state_after"], {"peer_value": 1, "hypothesis": 1})

    def test_veto_overrides_arbitration_and_prediction_uses_actual_action(self):
        record = loop_step(LoopState(0, 0), allow_toggle=False)
        self.assertEqual((record.requested_action, record.executed_action), (1, 0))
        self.assertTrue(record.vetoed)
        self.assertEqual((record.prediction, record.observation), (0, 0))

    def test_imposed_action_separates_protocol_from_autonomous_policy(self):
        state = LoopState(0, 0)
        self.assertEqual(loop_step(state).observation, 1)
        self.assertEqual(loop_step(state, imposed_action=0).observation, 0)

    def test_observation_noise_does_not_rewrite_peer_value(self):
        record = loop_step(LoopState(0, 0), observation_flip=True)
        self.assertEqual(record.state_after.peer_value, 1)
        self.assertEqual((record.observation, record.state_after.hypothesis), (0, 0))

    def test_all_protocol_kernels_are_stochastic(self):
        for protocol in (None, 0, 1):
            for noise in (0.0, 0.2, 0.5, 1.0):
                for allowed in (False, True):
                    kernel = transition_kernel(protocol, noise, allowed)
                    self.assertGreaterEqual(float(kernel.min()), 0)
                    np.testing.assert_allclose(kernel.sum(axis=0), 1, atol=1e-15)

    def test_noiseless_kernel_matches_executed_steps(self):
        kernel = transition_kernel()
        for column, state in enumerate(STATES):
            expected = np.zeros(len(STATES))
            expected[state_index(loop_step(state).state_after)] = 1
            np.testing.assert_array_equal(kernel[:, column], expected)

    def test_fixed_protocol_respects_mixture_affinity(self):
        self.assertLess(baseline_results()["max_mixture_affinity_error"], 1e-14)

    def test_overwriting_memory_can_decrease_retained_state_entropy(self):
        diagnostics = baseline_results()["memory_overwrite"]
        self.assertAlmostEqual(diagnostics["before_retained_state_entropy_nats"], math.log(2))
        self.assertAlmostEqual(diagnostics["after_retained_state_entropy_nats"], 0)
        np.testing.assert_array_equal(diagnostics["after_distribution"], [1, 0, 0, 0])

    def test_forward_stochastic_invertibility_is_not_physical_reversibility(self):
        control = baseline_results()["finite_simplex_reversibility_control"]
        forward = np.asarray(control["forward_stochastic_matrix"])
        inverse = np.asarray(control["algebraic_inverse"])
        np.testing.assert_allclose(inverse @ forward, np.eye(2), atol=1e-14)
        self.assertAlmostEqual(control["minimum_inverse_entry"], -0.125)

    def test_permutations_have_stochastic_inverses(self):
        permutation = np.array([[0.0, 1.0], [1.0, 0.0]])
        np.testing.assert_array_equal(permutation @ permutation, np.eye(2))
        self.assertGreaterEqual(float(permutation.min()), 0)

    def test_invalid_inputs_are_rejected(self):
        with self.assertRaises(ValueError):
            LoopState(2, 0)
        with self.assertRaises(ValueError):
            loop_step(LoopState(0, 0), imposed_action=3)
        for invalid_noise in (-0.1, 1.1, float("nan")):
            with self.assertRaises(ValueError):
                transition_kernel(observation_noise=invalid_noise)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(CognitiveLoopTests)
    outcome = unittest.TextTestRunner(verbosity=2).run(suite)
    if not outcome.wasSuccessful():
        raise SystemExit(1)
    results = {
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "automated_checks": {"run": outcome.testsRun, "failures": len(outcome.failures),
                             "errors": len(outcome.errors)},
        "baseline": baseline_results(),
    }
    if arguments.write_results:
        output = Path(__file__).with_name("cognitive_loop_results.json")
        output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
        print(f"Results written to {output}")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()