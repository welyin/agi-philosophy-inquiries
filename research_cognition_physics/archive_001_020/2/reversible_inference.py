"""Separate silent cognitive updates, statistical injectivity, and reversibility."""

import argparse
import json
import math
import platform
import unittest
from itertools import permutations
from pathlib import Path

import numpy as np

from cognitive_loop import LoopState, STATES, state_index
from predictive_states import outcome_kernels, test_probability


def hypothesis_update(outputs):
    if len(outputs) != 2 or any(value not in (0, 1) for value in outputs):
        raise ValueError("A binary output is required for each retained hypothesis.")
    update = np.zeros((len(STATES), len(STATES)))
    for column, state in enumerate(STATES):
        after = LoopState(state.peer_value, outputs[state.hypothesis])
        update[state_index(after), column] = 1.0
    return update


def noisy_reorganization(flip_probability=0.2):
    if not math.isfinite(flip_probability) or not 0 <= flip_probability <= 1:
        raise ValueError("The flip probability must lie between zero and one.")
    return ((1 - flip_probability) * np.eye(len(STATES))
            + flip_probability * hypothesis_update((1, 0)))


def silent_protocol(update):
    update = np.asarray(update, dtype=float)
    if (update.shape != (len(STATES), len(STATES))
            or not np.all(np.isfinite(update)) or np.any(update < 0)
            or not np.allclose(update.sum(axis=0), 1)):
        raise ValueError("The silent update must be a stochastic map on the loop states.")
    protocols = outcome_kernels()
    protocols["silent"] = {(): update}
    return protocols


def total_variation(first, second):
    return float(np.sum(abs(np.asarray(first) - np.asarray(second))) / 2)


def distinguishable_preparations():
    first = np.zeros(len(STATES))
    second = np.zeros(len(STATES))
    first[state_index(LoopState(0, 0))] = 1.0
    second[state_index(LoopState(0, 1))] = 1.0
    return first, second


def permutation_matrix(destinations):
    destinations = tuple(destinations)
    dimension = len(destinations)
    if sorted(destinations) != list(range(dimension)):
        raise ValueError("Destinations must specify a permutation.")
    matrix = np.zeros((dimension, dimension))
    matrix[list(destinations), np.arange(dimension)] = 1.0
    return matrix


def buffer_exchange():
    destinations = []
    for state in STATES:
        for buffer_value in (0, 1):
            after = LoopState(state.peer_value, buffer_value)
            destinations.append(2 * state_index(after) + state.hypothesis)
    return permutation_matrix(destinations)


def seed_controlled_flip():
    destinations = []
    for state in STATES:
        for seed in (0, 1):
            after = LoopState(state.peer_value, state.hypothesis ^ seed)
            destinations.append(2 * state_index(after) + seed)
    return permutation_matrix(destinations)


def retained_storage_results():
    dimension = len(STATES)
    marginal = np.kron(np.eye(dimension), np.ones((1, 2)))
    first, second = distinguishable_preparations()
    rows = []
    for name, operation, auxiliary, expected in (
        ("clear_via_blank_buffer", buffer_exchange(), np.array([1.0, 0.0]),
         hypothesis_update((0, 0))),
        ("flip_with_retained_seed", seed_controlled_flip(), np.array([0.8, 0.2]),
         noisy_reorganization(0.2)),
    ):
        embedding = np.kron(np.eye(dimension), auxiliary[:, None])
        first_joint = operation @ (embedding @ first)
        second_joint = operation @ (embedding @ second)
        effective = marginal @ operation @ embedding
        recovered = marginal @ operation.T @ operation @ embedding
        rows.append({
            "case": name,
            "auxiliary_preparation": auxiliary.tolist(),
            "max_reduced_map_error": float(np.max(abs(effective - expected))),
            "retained_loop_total_variation": total_variation(marginal @ first_joint, marginal @ second_joint),
            "joint_total_variation": total_variation(first_joint, second_joint),
            "max_recovery_error": float(np.max(abs(recovered - np.eye(dimension)))),
        })
    return {
        "rows": rows,
        "scope": "The auxiliary bit is reversible scratch storage or a retained random seed, not an immutable SoCA audit log. The clear construction requires an initially blank buffer. Access to the full joint register permits recovery; tracing it out gives the original irreversible reduced maps. No heat or thermodynamic energy is modeled.",
    }


def continuous_random_update(parameter):
    if not math.isfinite(parameter) or parameter < 0:
        raise ValueError("The interpolation parameter must be finite and nonnegative.")
    return noisy_reorganization((1 - math.exp(-2 * parameter)) / 2)


def continuous_update_results():
    first, second = distinguishable_preparations()
    rows = []
    for parameter in (0.0, 0.25, 0.5, 1.0):
        update = continuous_random_update(parameter)
        rows.append({
            "parameter": parameter,
            "matrix_rank": int(np.linalg.matrix_rank(update)),
            "total_variation": total_variation(update @ first, update @ second),
            "minimum_algebraic_inverse_entry": float(np.linalg.inv(update).min()),
        })
    return {
        "rows": rows,
        "scope": "A continuous classical Markov semigroup with a dimensionless update parameter, not physical time. It is injective at finite parameter but contracts distinctions and has no positive inverse away from zero.",
    }


def candidate_group():
    return tuple(permutation_matrix(destinations) for destinations in permutations(range(4)))


def invariant_metric_data():
    encoding = np.column_stack((np.eye(3), -np.ones(3)))
    embedding = np.linalg.pinv(encoding)
    representations = np.asarray([encoding @ operation @ embedding for operation in candidate_group()])
    metric = np.mean([representation.T @ representation for representation in representations], axis=0)
    whitening = np.linalg.cholesky(metric).T
    whitened = np.asarray([whitening @ representation @ np.linalg.inv(whitening)
                           for representation in representations])
    return encoding, representations, metric, whitened


def invariant_metric_results():
    encoding, representations, metric, whitened = invariant_metric_data()
    metric_error = max(float(np.max(abs(representation.T @ metric @ representation - metric)))
                       for representation in representations)
    orthogonality_error = max(float(np.max(abs(representation.T @ representation - np.eye(3))))
                             for representation in whitened)
    displacements = [float(np.linalg.norm(representation - np.eye(3))) for representation in whitened]
    flip_first = permutation_matrix((2, 3, 0, 1))
    copy_parity = permutation_matrix((0, 1, 3, 2))
    initial = np.array([1.0, 0.0, 0.0, 0.0])
    first_order = copy_parity @ flip_first @ initial
    second_order = flip_first @ copy_parity @ initial
    return {
        "internal_configurations": ["00", "01", "10", "11"],
        "number_of_reversible_operations": len(representations),
        "coordinate_encoding": encoding.tolist(),
        "invariant_metric": metric.tolist(),
        "metric_eigenvalues": np.linalg.eigvalsh(metric).tolist(),
        "max_metric_invariance_error": metric_error,
        "max_whitened_orthogonality_error": orthogonality_error,
        "minimum_nonidentity_whitened_frobenius_distance": min(value for value in displacements if value > 1e-10),
        "noncommuting_control": {
            "flip_then_controlled_parity_distribution": first_order.tolist(),
            "controlled_parity_then_flip_distribution": second_order.tolist(),
            "output_total_variation": total_variation(first_order, second_order),
        },
        "scope": "All 24 permutations of four configurations of two internal working bits are explicitly allowed in this control. These are not the four peer-hypothesis states of the earlier loop; the peer is an unchanged spectator. Averaging the finite group derives an invariant positive quadratic form, but the state region is still a classical tetrahedron and the group is discrete, not quantum dynamics or a spacetime metric.",
    }


def silent_update_results():
    first, second = distinguishable_preparations()
    operations = {
        "hypothesis_flip": hypothesis_update((1, 0)),
        "hypothesis_clear": hypothesis_update((0, 0)),
        "unrecorded_random_flip": noisy_reorganization(0.2),
    }
    rows = []
    for name, operation in operations.items():
        protocols = silent_protocol(operation)
        rank = int(np.linalg.matrix_rank(operation))
        inverse = np.linalg.inv(operation) if rank == len(STATES) else None
        future_test = (("silent", ()), ("autonomous", (1,)))
        rows.append({
            "operation": name,
            "update": operation.tolist(),
            "public_silent_outcome_probability_first": test_probability(protocols, first, (("silent", ()),)),
            "public_silent_outcome_probability_second": test_probability(protocols, second, (("silent", ()),)),
            "matrix_rank": rank,
            "first_output_distribution": (operation @ first).tolist(),
            "second_output_distribution": (operation @ second).tolist(),
            "total_variation_after_update": total_variation(operation @ first, operation @ second),
            "future_observation_one_probability_first": test_probability(protocols, first, future_test),
            "future_observation_one_probability_second": test_probability(protocols, second, future_test),
            "minimum_algebraic_inverse_entry": None if inverse is None else float(inverse.min()),
            "inverse_is_stochastic": bool(inverse is not None and np.min(inverse) >= -1e-13
                                          and np.allclose(inverse.sum(axis=0), 1)),
        })
    return {
        "initial_first": first.tolist(),
        "initial_second": second.tolist(),
        "initial_total_variation": total_variation(first, second),
        "rows": rows,
        "scope": "Silent means a single public result label, not absence of internal state change or hidden environmental records. The random flip seed and old audit state are inaccessible in this reduced model. No quantum evolution is assumed.",
    }


class ReversibleInferenceTests(unittest.TestCase):
    def test_silent_result_is_identical_for_distinct_inputs(self):
        for row in silent_update_results()["rows"]:
            self.assertAlmostEqual(row["public_silent_outcome_probability_first"], 1)
            self.assertAlmostEqual(row["public_silent_outcome_probability_second"], 1)

    def test_silent_clear_can_destroy_future_distinguishability(self):
        first, second = distinguishable_preparations()
        clear = hypothesis_update((0, 0))
        np.testing.assert_array_equal(clear @ first, clear @ second)
        self.assertEqual(np.linalg.matrix_rank(clear), 2)
        self.assertAlmostEqual(total_variation(clear @ first, clear @ second), 0)

    def test_hypothesis_flip_is_reversible_and_preserves_peer(self):
        flip = hypothesis_update((1, 0))
        np.testing.assert_array_equal(flip @ flip, np.eye(len(STATES)))
        peer_effect = np.array([state.peer_value for state in STATES])
        np.testing.assert_array_equal(peer_effect @ flip, peer_effect)
        first, second = distinguishable_preparations()
        self.assertAlmostEqual(total_variation(flip @ first, flip @ second), 1)

    def test_injective_random_update_has_no_stochastic_inverse(self):
        update = noisy_reorganization(0.2)
        self.assertEqual(np.linalg.matrix_rank(update), len(STATES))
        inverse = np.linalg.inv(update)
        np.testing.assert_allclose(inverse @ update, np.eye(len(STATES)), atol=1e-14)
        self.assertAlmostEqual(float(inverse.min()), -1 / 3)

    def test_random_update_reduces_readable_distinction(self):
        first, second = distinguishable_preparations()
        update = noisy_reorganization(0.2)
        self.assertAlmostEqual(total_variation(update @ first, update @ second), 0.6)
        protocols = silent_protocol(update)
        future_test = (("silent", ()), ("autonomous", (1,)))
        self.assertAlmostEqual(test_probability(protocols, first, future_test), 0.8)
        self.assertAlmostEqual(test_probability(protocols, second, future_test), 0.2)

    def test_reversible_update_can_change_predictions_without_losing_distinctions(self):
        first, second = distinguishable_preparations()
        protocols = silent_protocol(hypothesis_update((1, 0)))
        future_test = (("silent", ()), ("autonomous", (1,)))
        self.assertAlmostEqual(test_probability(protocols, first, future_test), 0)
        self.assertAlmostEqual(test_probability(protocols, second, future_test), 1)

    def test_unrecorded_choice_is_an_affine_mixture_of_operations(self):
        preparation = np.array([0.1, 0.2, 0.3, 0.4])
        expected = 0.8 * preparation + 0.2 * (hypothesis_update((1, 0)) @ preparation)
        np.testing.assert_allclose(noisy_reorganization(0.2) @ preparation, expected)

    def test_invalid_operations_are_rejected(self):
        for invalid in (-0.1, 1.1, float("nan")):
            with self.assertRaises(ValueError):
                noisy_reorganization(invalid)
        with self.assertRaises(ValueError):
            silent_protocol(np.zeros((len(STATES), len(STATES))))
        with self.assertRaises(ValueError):
            hypothesis_update((0, 2))


    def test_retained_auxiliary_bits_reproduce_reduced_maps_and_recover_inputs(self):
        diagnostics = retained_storage_results()
        for row in diagnostics["rows"]:
            self.assertLess(row["max_reduced_map_error"], 1e-14)
            self.assertLess(row["max_recovery_error"], 1e-14)
            self.assertAlmostEqual(row["joint_total_variation"], 1)
        self.assertAlmostEqual(diagnostics["rows"][0]["retained_loop_total_variation"], 0)
        self.assertAlmostEqual(diagnostics["rows"][1]["retained_loop_total_variation"], 0.6)

    def test_buffer_and_seed_updates_are_global_reversible_permutations(self):
        for operation in (buffer_exchange(), seed_controlled_flip()):
            np.testing.assert_array_equal(operation.T @ operation, np.eye(2 * len(STATES)))
            peer_effect = np.repeat([state.peer_value for state in STATES], 2)
            np.testing.assert_array_equal(peer_effect @ operation, peer_effect)

    def test_reusing_nonblank_buffer_does_not_perform_another_clear(self):
        initial = np.zeros(2 * len(STATES))
        initial[2 * state_index(LoopState(0, 1))] = 1
        operation = buffer_exchange()
        marginal = np.kron(np.eye(len(STATES)), np.ones((1, 2)))
        once = marginal @ operation @ initial
        twice = marginal @ operation @ operation @ initial
        self.assertAlmostEqual(once[state_index(LoopState(0, 0))], 1)
        self.assertAlmostEqual(twice[state_index(LoopState(0, 1))], 1)

    def test_continuous_update_has_semigroup_composition(self):
        combined = continuous_random_update(0.3) @ continuous_random_update(0.7)
        np.testing.assert_allclose(combined, continuous_random_update(1.0), atol=1e-14)
        np.testing.assert_array_equal(continuous_random_update(0.0), np.eye(len(STATES)))

    def test_continuous_injectivity_does_not_preserve_distinguishability(self):
        for row in continuous_update_results()["rows"]:
            self.assertEqual(row["matrix_rank"], len(STATES))
            self.assertAlmostEqual(row["total_variation"], math.exp(-2 * row["parameter"]))
            if row["parameter"] > 0:
                self.assertLess(row["minimum_algebraic_inverse_entry"], 0)

    def test_candidate_permutations_preserve_classical_probabilities(self):
        group = candidate_group()
        self.assertEqual(len(group), 24)
        for operation in group:
            np.testing.assert_array_equal(operation.sum(axis=0), np.ones(4))
            np.testing.assert_array_equal(operation.T @ operation, np.eye(4))

    def test_group_averaging_derives_positive_invariant_metric(self):
        _, representations, metric, _ = invariant_metric_data()
        expected = 2 * np.eye(3) - 0.5 * np.ones((3, 3))
        np.testing.assert_allclose(metric, expected, atol=1e-13)
        self.assertGreater(float(np.linalg.eigvalsh(metric).min()), 0)
        for representation in representations:
            np.testing.assert_allclose(representation.T @ metric @ representation, metric, atol=1e-13)

    def test_whitening_produces_orthogonal_maps_without_quantum_assumptions(self):
        diagnostics = invariant_metric_results()
        self.assertLess(diagnostics["max_metric_invariance_error"], 1e-13)
        self.assertLess(diagnostics["max_whitened_orthogonality_error"], 1e-13)

    def test_finite_reversible_group_has_a_nonzero_identity_gap(self):
        diagnostics = invariant_metric_results()
        self.assertAlmostEqual(diagnostics["minimum_nonidentity_whitened_frobenius_distance"], 2)

    def test_reversible_classical_internal_operations_can_fail_to_commute(self):
        control = invariant_metric_results()["noncommuting_control"]
        np.testing.assert_array_equal(control["flip_then_controlled_parity_distribution"], [0, 0, 0, 1])
        np.testing.assert_array_equal(control["controlled_parity_then_flip_distribution"], [0, 0, 1, 0])
        self.assertAlmostEqual(control["output_total_variation"], 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(ReversibleInferenceTests)
    outcome = unittest.TextTestRunner(verbosity=2).run(suite)
    if not outcome.wasSuccessful():
        raise SystemExit(1)
    results = {
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "automated_checks": {"run": outcome.testsRun, "failures": len(outcome.failures),
                             "errors": len(outcome.errors)},
        "silent_updates": silent_update_results(),
        "retained_storage": retained_storage_results(),
        "continuous_updates": continuous_update_results(),
        "invariant_relation_metric": invariant_metric_results(),
    }
    if arguments.write_results:
        output = Path(__file__).with_name("reversible_inference_results.json")
        output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
        print(f"Results written to {output}")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()