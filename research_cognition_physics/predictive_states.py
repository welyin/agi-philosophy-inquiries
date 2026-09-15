"""Find predictive states of the cognitive loop from allowed future protocols."""

import argparse
import json
import math
import platform
import unittest
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from cognitive_loop import LoopState, STATES, loop_step, state_index, transition_kernel


PROTOCOLS = {"keep": 0, "toggle": 1, "autonomous": None}
COORDINATES = np.array([
    [1.0, 1.0, 1.0, 1.0],
    [0.0, 0.0, 1.0, 1.0],
    [1.0, 0.0, 0.0, 1.0],
])


@dataclass(frozen=True)
class PredictiveBasis:
    effects: np.ndarray
    orthonormal_columns: np.ndarray
    witness_tests: tuple
    closure_residual: float

    @property
    def rank(self):
        return len(self.effects)


def outcome_kernels(protocols=tuple(PROTOCOLS), observation_noise=0.0, expose_action=False):
    if not math.isfinite(observation_noise) or not 0 <= observation_noise <= 1:
        raise ValueError("Observation noise must be a probability.")
    kernels = {}
    for protocol in protocols:
        if protocol not in PROTOCOLS:
            raise ValueError("Unknown cognitive protocol.")
        outcomes = {}
        for column, state in enumerate(STATES):
            for flipped, probability in ((False, 1 - observation_noise),
                                         (True, observation_noise)):
                if probability == 0:
                    continue
                record = loop_step(state, imposed_action=PROTOCOLS[protocol],
                                   observation_flip=flipped)
                label = ((record.executed_action, record.observation) if expose_action
                         else (record.observation,))
                if label not in outcomes:
                    outcomes[label] = np.zeros((len(STATES), len(STATES)))
                outcomes[label][state_index(record.state_after), column] += probability
        kernels[protocol] = outcomes
    return kernels


def test_effect(kernels, test):
    dimension = next(iter(next(iter(kernels.values())).values())).shape[0]
    effect = np.ones(dimension)
    for protocol, label in reversed(test):
        effect = kernels[protocol][label].T @ effect
    return effect


def test_probability(kernels, probabilities, test):
    return float(test_effect(kernels, test) @ np.asarray(probabilities))


def predictive_basis(kernels, tolerance=1e-12):
    dimension = next(iter(next(iter(kernels.values())).values())).shape[0]
    effects = [np.ones(dimension)]
    orthonormal = [effects[0] / np.linalg.norm(effects[0])]
    tests = [()]
    frontier = 0
    while frontier < len(effects):
        for protocol, outcomes in kernels.items():
            for label, kernel in outcomes.items():
                candidate = kernel.T @ effects[frontier]
                residual = candidate.copy()
                for _ in range(2):
                    for basis_vector in orthonormal:
                        residual -= np.dot(basis_vector, residual) * basis_vector
                residual_norm = float(np.linalg.norm(residual))
                if residual_norm > tolerance * max(1.0, float(np.linalg.norm(candidate))):
                    effects.append(candidate)
                    orthonormal.append(residual / residual_norm)
                    tests.append(((protocol, label),) + tests[frontier])
        frontier += 1
    columns = np.column_stack(orthonormal)
    projector = columns @ columns.T
    closure_error = 0.0
    for outcomes in kernels.values():
        for kernel in outcomes.values():
            transformed = kernel.T @ columns
            closure_error = max(closure_error,
                                float(np.max(abs(transformed - projector @ transformed))))
    return PredictiveBasis(np.asarray(effects), columns, tuple(tests), closure_error)


def basis_summary(basis):
    return {
        "linear_predictive_rank_including_normalization": basis.rank,
        "normalized_affine_dimension": basis.rank - 1,
        "max_selected_test_length": max(map(len, basis.witness_tests)),
        "invariant_span_residual": basis.closure_residual,
        "selected_effects": basis.effects.tolist(),
        "selected_tests": [[{"protocol": protocol, "outcome": list(label)}
                            for protocol, label in test] for test in basis.witness_tests],
    }


def distinguishing_test(kernels, first, second, tolerance=1e-12):
    basis = predictive_basis(kernels)
    difference = np.asarray(first) - np.asarray(second)
    for effect, test in zip(basis.effects, basis.witness_tests):
        if abs(float(effect @ difference)) > tolerance:
            return {
                "test": [{"protocol": protocol, "outcome": list(label)}
                         for protocol, label in test],
                "first_probability": float(effect @ first),
                "second_probability": float(effect @ second),
            }
    return None


def exchange_hypotheses():
    dimension = len(STATES)
    exchange = np.zeros((dimension**2, dimension**2))
    for first_index, first_state in enumerate(STATES):
        for second_index, second_state in enumerate(STATES):
            first_after = LoopState(first_state.peer_value, second_state.hypothesis)
            second_after = LoopState(second_state.peer_value, first_state.hypothesis)
            row = dimension * state_index(first_after) + state_index(second_after)
            column = dimension * first_index + second_index
            exchange[row, column] = 1.0
    return exchange


def sharing_context_kernels(observation_noise=0.0):
    dimension = len(STATES)
    receiver = np.zeros(dimension)
    receiver[state_index(LoopState(0, 0))] = 1.0
    exchange = exchange_hypotheses()
    receiver_operations = outcome_kernels(("autonomous",), observation_noise)["autonomous"]
    outcomes = {}
    for label, receiver_kernel in receiver_operations.items():
        joint_operation = np.kron(np.eye(dimension), receiver_kernel) @ exchange
        source_kernel = np.zeros((dimension, dimension))
        for column in range(dimension):
            initial_joint = np.kron(np.eye(dimension)[:, column], receiver)
            after_joint = joint_operation @ initial_joint
            source_kernel[:, column] = after_joint.reshape(dimension, dimension).sum(axis=1)
        outcomes[label] = source_kernel
    return outcomes


def reduced_local_results():
    kernels = outcome_kernels()
    right_inverse = np.linalg.pinv(COORDINATES)
    closure_error = 0.0
    for outcomes in kernels.values():
        for kernel in outcomes.values():
            reduced = COORDINATES @ kernel @ right_inverse
            closure_error = max(closure_error,
                                float(np.max(abs(reduced @ COORDINATES - COORDINATES @ kernel))))
    keep_zero = COORDINATES @ kernels["keep"][(0,)] @ right_inverse
    return {
        "normalized_state_space_vertices": COORDINATES[1:].T.tolist(),
        "keep_zero_unnormalized_coordinate_update": keep_zero.tolist(),
        "minimum_reduced_entry": float(keep_zero.min()),
        "max_local_intertwining_error": closure_error,
        "interpretation": "The local predictive state space is a square, a projection of a classical simplex under restricted access. Negative entries in a compressed-coordinate update are not negative probabilities or evidence of quantum amplitudes.",
    }


def composition_results():
    dimension = len(STATES)
    first = np.array([0.5, 0.0, 0.5, 0.0])
    second = np.array([0.0, 0.5, 0.0, 0.5])
    receiver = np.zeros(dimension)
    receiver[state_index(LoopState(0, 0))] = 1.0
    local = outcome_kernels()
    exchange = exchange_hypotheses()
    receiver_one = test_effect(local, (("autonomous", (1,)),))
    joint_effect = exchange.T @ np.kron(np.ones(dimension), receiver_one)
    first_joint = np.kron(first, receiver)
    second_joint = np.kron(second, receiver)
    joint_coordinates = np.kron(COORDINATES, COORDINATES)
    reduced_exchange = joint_coordinates @ exchange @ np.linalg.pinv(joint_coordinates)
    exchange_error = reduced_exchange @ joint_coordinates - joint_coordinates @ exchange
    extended = dict(local)
    extended["share_then_query"] = sharing_context_kernels()
    complete_coordinates = np.vstack((COORDINATES, [0.0, 1.0, 0.0, 1.0]))
    return {
        "receiver_preparation": receiver.tolist(),
        "operation": "Exchange the two retained hypotheses, leaving both peer values unchanged; then execute only the receiver's autonomous loop and read its feedback observation.",
        "max_initial_product_coordinate_difference": float(np.max(abs(joint_coordinates @ (first_joint - second_joint)))),
        "receiver_one_probability_first": float(joint_effect @ first_joint),
        "receiver_one_probability_second": float(joint_effect @ second_joint),
        "context_effect_on_source": (joint_effect.reshape(dimension, dimension) @ receiver).tolist(),
        "max_exchange_intertwining_error": float(np.max(abs(exchange_error))),
        "composition_aware_basis": basis_summary(predictive_basis(extended)),
        "context_distinguishing_test": distinguishing_test(extended, first, second),
        "execution_audit_basis": basis_summary(predictive_basis(outcome_kernels(expose_action=True))),
        "complete_coordinate_map": complete_coordinates.tolist(),
        "complete_coordinate_inverse": np.linalg.inv(complete_coordinates).tolist(),
        "scope": "An explicitly added classical memory-exchange operation motivated by SoCA's inter-agent sharing. Cartesian joint probabilities and independent prepared inputs are assumed for the witness. No Bell nonlocality, quantum tensor product, or general quantum composition theorem is derived.",
    }


def local_results():
    restricted = outcome_kernels(("keep", "toggle"))
    complete = outcome_kernels()
    noisy = outcome_kernels(observation_noise=0.2)
    blind = outcome_kernels(observation_noise=0.5)
    aligned = np.array([0.5, 0.0, 0.0, 0.5])
    opposed = np.array([0.0, 0.5, 0.5, 0.0])
    hidden_zero = np.array([0.5, 0.0, 0.5, 0.0])
    hidden_one = np.array([0.0, 0.5, 0.0, 0.5])
    basis = predictive_basis(complete)
    return {
        "interfaces": {
            "forced_actions_observation_only": basis_summary(predictive_basis(restricted)),
            "with_autonomy_observation_only": basis_summary(basis),
            "with_autonomy_noise_0_2": basis_summary(predictive_basis(noisy)),
            "with_autonomy_noise_0_5": basis_summary(predictive_basis(blind)),
        },
        "correlation_witness": {
            "aligned_preparation": aligned.tolist(),
            "opposed_preparation": opposed.tolist(),
            "peer_one_probability_each": float(COORDINATES[1] @ aligned),
            "hypothesis_one_probability_each": 0.5,
            "forced_interface_witness": distinguishing_test(restricted, aligned, opposed),
            "autonomous_interface_witness": distinguishing_test(complete, aligned, opposed),
            "noisy_autonomous_interface_witness": distinguishing_test(noisy, aligned, opposed),
        },
        "local_equivalence": {
            "hypothesis_zero_preparation": hidden_zero.tolist(),
            "hypothesis_one_preparation": hidden_one.tolist(),
            "first_coordinates": (COORDINATES @ hidden_zero).tolist(),
            "second_coordinates": (COORDINATES @ hidden_one).tolist(),
            "max_basis_probability_difference": float(np.max(abs(basis.effects @ (hidden_zero - hidden_one)))),
            "distinguishing_test": distinguishing_test(complete, hidden_zero, hidden_one),
        },
        "coordinate_map": COORDINATES.tolist(),
        "coordinate_names": ["normalization", "probability_peer_one", "probability_peer_equals_hypothesis"],
        "scope": "Prepared distributions over the existing four classical states; only feedback observations are exposed and executed actions are hidden. Future-test span closure covers arbitrary finite adaptive local protocols under the declared interface. This is not a quantum state reconstruction or a claim about all reachable post-feedback histories.",
    }


class PredictiveStateTests(unittest.TestCase):
    def test_outcome_kernels_sum_to_existing_unconditional_update(self):
        for noise in (0.0, 0.2, 0.5):
            kernels = outcome_kernels(observation_noise=noise)
            for protocol, outcomes in kernels.items():
                total = sum(outcomes.values())
                np.testing.assert_allclose(total, transition_kernel(PROTOCOLS[protocol], noise))
                np.testing.assert_allclose(total.sum(axis=0), 1)
                self.assertTrue(all(np.min(kernel) >= 0 for kernel in outcomes.values()))

    def test_ordered_test_matches_direct_noiseless_feedback(self):
        kernels = outcome_kernels()
        preparation = np.array([0.0, 1.0, 0.0, 0.0])
        test = (("autonomous", (0,)), ("autonomous", (1,)), ("keep", (1,)))
        self.assertAlmostEqual(test_probability(kernels, preparation, test), 1)
        self.assertAlmostEqual(test_probability(kernels, preparation,
                                               (("autonomous", (1,)),)), 0)

    def test_forced_actions_need_only_the_peer_marginal(self):
        basis = predictive_basis(outcome_kernels(("keep", "toggle")))
        self.assertEqual(basis.rank, 2)
        self.assertLess(basis.closure_residual, 1e-13)
        self.assertIsNone(local_results()["correlation_witness"]["forced_interface_witness"])

    def test_autonomy_exposes_correlation_not_fixed_by_two_marginals(self):
        witness = local_results()["correlation_witness"]["autonomous_interface_witness"]
        self.assertEqual(witness["test"][0]["protocol"], "autonomous")
        self.assertAlmostEqual(abs(witness["first_probability"] - witness["second_probability"]), 1)
        self.assertEqual(predictive_basis(outcome_kernels()).rank, 3)

    def test_future_effect_span_is_closed_without_a_horizon_cutoff(self):
        kernels = outcome_kernels()
        basis = predictive_basis(kernels)
        self.assertEqual(max(map(len, basis.witness_tests)), 1)
        self.assertLess(basis.closure_residual, 1e-13)
        coordinate_projector = np.linalg.pinv(COORDINATES) @ COORDINATES
        np.testing.assert_allclose(basis.effects @ coordinate_projector, basis.effects, atol=1e-13)

    def test_distinct_preparations_can_be_locally_predictively_equivalent(self):
        diagnostics = local_results()["local_equivalence"]
        self.assertNotEqual(diagnostics["hypothesis_zero_preparation"],
                            diagnostics["hypothesis_one_preparation"])
        self.assertEqual(diagnostics["first_coordinates"], diagnostics["second_coordinates"])
        self.assertLess(diagnostics["max_basis_probability_difference"], 1e-13)
        self.assertIsNone(diagnostics["distinguishing_test"])

    def test_equivalence_survives_local_conditioned_updates(self):
        kernels = outcome_kernels(observation_noise=0.2)
        basis = predictive_basis(kernels)
        first = np.array([0.5, 0.0, 0.5, 0.0])
        second = np.array([0.0, 0.5, 0.0, 0.5])
        for outcomes in kernels.values():
            for kernel in outcomes.values():
                first_branch = kernel @ first
                second_branch = kernel @ second
                probability = float(first_branch.sum())
                self.assertAlmostEqual(probability, float(second_branch.sum()))
                if probability > 0:
                    np.testing.assert_allclose(basis.effects @ (first_branch / probability),
                                               basis.effects @ (second_branch / probability), atol=1e-13)

    def test_observation_noise_controls_identifiability(self):
        self.assertEqual(predictive_basis(outcome_kernels(observation_noise=0.2)).rank, 3)
        self.assertEqual(predictive_basis(outcome_kernels(observation_noise=0.5)).rank, 1)
        witness = local_results()["correlation_witness"]["noisy_autonomous_interface_witness"]
        self.assertAlmostEqual(abs(witness["first_probability"] - witness["second_probability"]), 0.6)


    def test_local_compressed_dynamics_close_despite_signed_coefficients(self):
        diagnostics = reduced_local_results()
        expected = np.array([[1.0, -1.0, 0.0], [0.0, 0.0, 0.0], [1.0, -1.0, 0.0]])
        np.testing.assert_allclose(diagnostics["keep_zero_unnormalized_coordinate_update"], expected, atol=1e-13)
        self.assertLess(diagnostics["max_local_intertwining_error"], 1e-13)
        self.assertAlmostEqual(diagnostics["minimum_reduced_entry"], -1)

    def test_square_center_has_distinct_classical_preimages(self):
        first = np.array([0.5, 0.0, 0.5, 0.0])
        second = np.array([0.0, 0.5, 0.0, 0.5])
        np.testing.assert_array_equal(COORDINATES @ first, [1.0, 0.5, 0.5])
        np.testing.assert_array_equal(COORDINATES @ second, [1.0, 0.5, 0.5])
        vertices = {tuple(vertex) for vertex in COORDINATES[1:].T}
        self.assertEqual(vertices, {(0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0)})

    def test_exchange_is_reversible_and_does_not_change_peer_records(self):
        dimension = len(STATES)
        exchange = exchange_hypotheses()
        np.testing.assert_array_equal(exchange @ exchange, np.eye(dimension**2))
        np.testing.assert_array_equal(exchange.sum(axis=0), 1)
        first_peer = np.array([state.peer_value for state in STATES for _ in STATES])
        second_peer = np.array([state.peer_value for _ in STATES for state in STATES])
        np.testing.assert_array_equal(first_peer @ exchange, first_peer)
        np.testing.assert_array_equal(second_peer @ exchange, second_peer)

    def test_local_equivalence_fails_after_allowed_memory_sharing(self):
        diagnostics = composition_results()
        self.assertLess(diagnostics["max_initial_product_coordinate_difference"], 1e-13)
        self.assertAlmostEqual(diagnostics["receiver_one_probability_first"], 1)
        self.assertAlmostEqual(diagnostics["receiver_one_probability_second"], 0)
        self.assertGreater(diagnostics["max_exchange_intertwining_error"], 0.1)

    def test_sharing_context_is_a_valid_classical_operation(self):
        for noise in (0.0, 0.2):
            kernels = sharing_context_kernels(noise)
            np.testing.assert_allclose(sum(kernels.values()).sum(axis=0), 1)
            self.assertTrue(all(np.min(kernel) >= 0 for kernel in kernels.values()))
        effect = np.ones(len(STATES)) @ sharing_context_kernels()[(1,)]
        np.testing.assert_array_equal(effect, [1.0, 0.0, 1.0, 0.0])

    def test_sharing_closes_only_after_restoring_missing_hypothesis_statistic(self):
        diagnostics = composition_results()
        self.assertEqual(diagnostics["composition_aware_basis"]["linear_predictive_rank_including_normalization"], 4)
        self.assertLess(diagnostics["composition_aware_basis"]["invariant_span_residual"], 1e-13)
        witness = diagnostics["context_distinguishing_test"]
        self.assertEqual(witness["test"][0]["protocol"], "share_then_query")
        complete = np.asarray(diagnostics["complete_coordinate_map"])
        inverse = np.asarray(diagnostics["complete_coordinate_inverse"])
        np.testing.assert_allclose(inverse @ complete, np.eye(len(STATES)), atol=1e-13)

    def test_exposing_execution_audit_already_restores_full_rank(self):
        basis = predictive_basis(outcome_kernels(expose_action=True))
        self.assertEqual(basis.rank, 4)
        self.assertLess(basis.closure_residual, 1e-13)

    def test_conditional_compressed_updates_predict_the_same_future_tests(self):
        kernels = outcome_kernels()
        preparation = np.array([0.1, 0.2, 0.3, 0.4])
        right_inverse = np.linalg.pinv(COORDINATES)
        for outcomes in kernels.values():
            for kernel in outcomes.values():
                reduced = COORDINATES @ kernel @ right_inverse
                predicted_branch = reduced @ (COORDINATES @ preparation)
                full_branch = kernel @ preparation
                np.testing.assert_allclose(predicted_branch, COORDINATES @ full_branch, atol=1e-13)
                probability = float(predicted_branch[0])
                self.assertGreater(probability, 0)
                np.testing.assert_allclose(predicted_branch / probability,
                                           COORDINATES @ (full_branch / full_branch.sum()), atol=1e-13)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(PredictiveStateTests)
    outcome = unittest.TextTestRunner(verbosity=2).run(suite)
    if not outcome.wasSuccessful():
        raise SystemExit(1)
    results = {
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "automated_checks": {"run": outcome.testsRun, "failures": len(outcome.failures),
                             "errors": len(outcome.errors)},
        "local_predictive_states": local_results(),
        "reduced_local_dynamics": reduced_local_results(),
        "composition_consistency": composition_results(),
    }
    if arguments.write_results:
        output = Path(__file__).with_name("predictive_state_results.json")
        output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
        print(f"Results written to {output}")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()