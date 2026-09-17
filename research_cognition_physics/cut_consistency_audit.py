"""Round 175: regrouping a fixed experiment differs from losing its relations."""

import argparse
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import unittest

import numpy as np
from cognitive_loop import LoopState, loop_step


def global_records(depth, noise=F(1, 4), overrides=None):
    """Enumerate the initial world and the whole retained random seed tape."""
    overrides = overrides or {}
    result = defaultdict(F)
    for x0 in (0, 1):
        for tape in product((0, 1), repeat=depth):
            mass = F(1, 3) if x0 else F(2, 3)
            state = LoopState(x0, 0)
            record = ()
            for step, flip in enumerate(tape):
                mass *= noise if flip else 1 - noise
                event = loop_step(state, overrides.get(step), step % 3 != 1, flip)
                record += (event.executed_action, event.observation)
                state = event.state_after
            if mass:
                result[record] += mass
    return dict(result)


def local_factor_records(depth, noise=F(1, 4), overrides=None, boundaries=()):
    """Contract controller, actuator and sensor factors, retaining boundary state.

    No controller reads x. Grouping blocks changes the contraction schedule,
    not the physically available controller inputs or action rule.
    """
    overrides = overrides or {}
    nodes = {((), 0, 0): F(2, 3), ((), 1, 0): F(1, 3)}
    cuts = (0,) + tuple(sorted(set(boundaries) - {0, depth})) + (depth,)
    for start, stop in zip(cuts, cuts[1:]):
        for step in range(start, stop):
            following = defaultdict(F)
            for (history, x, memory), weight in nodes.items():
                requested = overrides.get(step, 1 - memory)
                action = 0 if step % 3 == 1 else requested
                world = x ^ action
                for observation in (0, 1):
                    likelihood = 1 - noise if observation == world else noise
                    if likelihood:
                        following[(history + (action, observation), world, observation)] += weight * likelihood
            nodes = dict(following)
    result = defaultdict(F)
    for (history, _, _), weight in nodes.items():
        result[history] += weight
    return dict(result)


def trace_keep(matrix, dims, keep):
    keep = tuple(keep)
    discard = tuple(i for i in range(len(dims)) if i not in keep)
    order = keep + discard
    tensor = matrix.reshape(tuple(dims) * 2).transpose(order + tuple(i + len(dims) for i in order))
    k = int(np.prod([dims[i] for i in keep]))
    e = int(np.prod([dims[i] for i in discard]))
    return np.trace(tensor.reshape(k, e, k, e), axis1=1, axis2=3)


class CutConsistencyAuditTests(unittest.TestCase):
    def test_all_regroupings_of_six_step_feedback_agree_exactly(self):
        for overrides in ({}, {0: 0, 3: 1}, {0: 1, 2: 0, 5: 0}):
            explicit = global_records(6, overrides=overrides)
            for mask in product((0, 1), repeat=5):
                boundaries = tuple(i + 1 for i, enabled in enumerate(mask) if enabled)
                self.assertEqual(local_factor_records(6, overrides=overrides, boundaries=boundaries), explicit)
            self.assertEqual(sum(explicit.values()), 1)

    def test_source_probabilities_zero_half_one_are_not_hidden_assumptions(self):
        for noise in (F(0), F(1, 2), F(1)):
            self.assertEqual(global_records(4, noise), local_factor_records(4, noise, boundaries=(1, 3)))

    def test_conditioning_on_action_is_not_intervening_on_action(self):
        # Initially memory=x, so the policy a=1-memory correlates a with x.
        conditional_weights = defaultdict(F)
        intervention_weights = defaultdict(F)
        for x, mass in ((0, F(2, 3)), (1, F(1, 3))):
            observed = loop_step(LoopState(x, x), observation_flip=False)
            if observed.executed_action == 1:
                conditional_weights[observed.observation] += mass
            imposed = loop_step(LoopState(x, x), imposed_action=1, observation_flip=False)
            intervention_weights[imposed.observation] += mass
        conditional = conditional_weights[1] / sum(conditional_weights.values())
        self.assertEqual(conditional, 1)
        self.assertEqual(intervention_weights[1], F(2, 3))

    def test_real_and_complex_nested_partial_traces_agree(self):
        rng = np.random.default_rng(175)
        for complex_model in (False, True):
            raw = rng.normal(size=(12, 5))
            if complex_model:
                raw = raw + 1j * rng.normal(size=raw.shape)
            rho = raw @ raw.conj().T
            rho /= np.trace(rho)
            direct = trace_keep(rho, (2, 3, 2), (0,))
            ab = trace_keep(rho, (2, 3, 2), (0, 1))
            ac = trace_keep(rho, (2, 3, 2), (0, 2))
            np.testing.assert_allclose(direct, trace_keep(ab, (2, 3), (0,)), atol=2e-16)
            np.testing.assert_allclose(direct, trace_keep(ac, (2, 2), (0,)), atol=2e-16)

    def test_same_intervention_in_schrodinger_and_heisenberg_descriptions(self):
        rng = np.random.default_rng(2175)
        for complex_model in (False, True):
            raw = rng.normal(size=(4, 4))
            if complex_model:
                raw = raw + 1j * rng.normal(size=raw.shape)
            rho = raw @ raw.conj().T
            rho /= np.trace(rho)
            cnot = np.eye(4)[[0, 1, 3, 2]]
            branch = np.kron(np.diag([.8, .6]), np.eye(2))
            effect = np.kron(np.eye(2), np.diag([1., 0.]))
            after = cnot @ branch @ rho @ branch.T @ cnot.T
            reduced = trace_keep(after, (2, 2), (1,))
            direct = np.trace(effect @ after)
            pulled = np.trace(branch.T @ cnot.T @ effect @ cnot @ branch @ rho)
            self.assertAlmostEqual(direct, pulled)
            self.assertAlmostEqual(direct, reduced[0, 0])

    def test_dephasing_a_record_changes_experiment_instead_of_only_moving_cut(self):
        cnot = np.eye(4)[[0, 1, 3, 2]]
        initial = np.outer([1., 0., 1., 0.], [1., 0., 1., 0.]) / 2
        encoded = cnot @ initial @ cnot.T
        dephased = sum(np.kron(np.eye(2), np.diag([int(i == j) for i in (0, 1)]))
                       @ encoded @ np.kron(np.eye(2), np.diag([int(i == j) for i in (0, 1)]))
                       for j in (0, 1))
        plus = np.ones((2, 2)) / 2
        effect = np.kron(plus, np.eye(2))
        self.assertEqual(np.trace(effect @ cnot.T @ encoded @ cnot), 1)
        self.assertEqual(np.trace(effect @ cnot.T @ dephased @ cnot), .5)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(CutConsistencyAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 175,
        "classical_feedback_horizon": 6,
        "groupings_per_protocol": 32,
        "intervention_protocols_checked": 3,
        "full_record_distribution_max_gap_exact": "0",
        "intervention_vs_conditioning_probabilities": ["2/3", "1"],
        "coherent_record_undo_success": "1",
        "dephased_record_undo_success": "1/2",
        "same_experiment_regrouping_selects_purification_or_complex_structure": False,
        "controller_receives_hidden_world_state": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("cut_consistency_audit_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
