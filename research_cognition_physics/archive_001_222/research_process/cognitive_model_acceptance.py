"""Round 58: the same finite predictive controller in real and complex candidates.

This is a calibrated finite-hypothesis controller, not general intelligence or
a derivation of either probability theory. Each trial uses a fresh preparation.
"""

import argparse
import json
import math
import unittest
from dataclasses import dataclass
from itertools import product
from pathlib import Path

import numpy as np

from bipartite_composition import interaction
from complex_bipartite_closure import local_effect
from quantum_interface_audit import ALPHA, IDENTITY, branch_kraus, two_rebit_states
from role_symmetry_and_swap import pauli_word


def entropy(probabilities):
    return -sum(float(p)*math.log(float(p)) for p in probabilities if p > 0)


def predictive_update(prior, likelihood, outcome):
    mass = np.asarray(prior)*likelihood[:, outcome]
    probability = float(mass.sum())
    if probability <= 0:
        raise ValueError("Cannot condition on an impossible record.")
    return probability, mass/probability


def information_gain(prior, likelihood):
    prediction = np.asarray(prior) @ likelihood
    return entropy(prediction)-sum(p*entropy(row) for p, row in zip(prior, likelihood))


def action_tables(candidate, visibility=.2):
    if candidate not in ("real", "complex"):
        raise ValueError("Use real or complex.")
    states = two_rebit_states()
    effects = {}
    effects["local_z"] = [np.kron(IDENTITY, local_effect("Z", s)) for s in (-1, 1)]
    gate = interaction(math.pi/2)
    effects["joint_parity"] = [gate.conj().T @ e @ gate for e in effects["local_z"]]
    if candidate == "complex":
        effects["separated_y"] = [sum(np.kron(local_effect("Y", s, visibility), local_effect("Y", t, visibility))
                                      for s, t in product((-1, 1), repeat=2) if s*t == r)
                                   for r in (-1, 1)]
    return {name: np.array([[np.trace(rho @ e).real for e in pair] for rho in states])
            for name, pair in effects.items()}


def choose_action(prior, tables, joint_cost=.1):
    scores = {name: information_gain(prior, table)-(joint_cost if name == "joint_parity" else 0)
              for name, table in tables.items()}
    return max(scores, key=scores.get)


@dataclass(frozen=True)
class Record:
    action: str
    prediction: tuple
    outcome: int
    posterior: tuple


def record_tree(tables, prior=(.5, .5), depth=3):
    """Enumerate all records; the controller uses its actual observed history."""
    leaves = []
    def visit(belief, mass, history, per_hypothesis):
        if len(history) == depth:
            leaves.append((mass, history, per_hypothesis))
            return
        action = choose_action(belief, tables)
        likelihood = tables[action]
        prediction = tuple(np.asarray(belief) @ likelihood)
        for outcome in (0, 1):
            probability, posterior = predictive_update(belief, likelihood, outcome)
            record = Record(action, prediction, outcome, tuple(posterior))
            visit(posterior, mass*probability, history+(record,), per_hypothesis*likelihood[:, outcome])
    visit(np.asarray(prior), 1., (), np.ones(len(prior)))
    return leaves


class CognitiveModelAcceptanceTests(unittest.TestCase):
    def test_both_candidates_have_normalized_positive_prediction_kernels(self):
        for model in ("real", "complex"):
            for table in action_tables(model).values():
                self.assertGreaterEqual(table.min(), 0)
                np.testing.assert_allclose(table.sum(axis=1), 1, atol=3e-16)

    def test_old_interaction_is_informative_in_both_models(self):
        for model in ("real", "complex"):
            tables = action_tables(model)
            np.testing.assert_allclose(tables["joint_parity"], [[(1-ALPHA)/2, (1+ALPHA)/2], [(1+ALPHA)/2, (1-ALPHA)/2]], atol=3e-16)
            self.assertGreater(information_gain([.5, .5], tables["joint_parity"]), .65)
            self.assertAlmostEqual(information_gain([.5, .5], tables["local_z"]), 0)

    def test_new_separated_information_has_the_predicted_visibility(self):
        table = action_tables("complex")["separated_y"]
        np.testing.assert_allclose(table, [[.48, .52], [.52, .48]], atol=3e-16)

    def test_complete_adaptive_record_distributions_and_bayes_are_consistent(self):
        prior = np.array([.4, .6])
        for model in ("real", "complex"):
            leaves = record_tree(action_tables(model), prior, 5)
            self.assertAlmostEqual(sum(mass for mass, _, _ in leaves), 1, places=13)
            np.testing.assert_allclose(sum(row for _, _, row in leaves), np.ones(2), atol=2e-15)
            for mass, history, likelihood in leaves:
                self.assertAlmostEqual(mass, prior @ likelihood, places=14)
                np.testing.assert_allclose(history[-1].posterior, prior*likelihood/mass, atol=1e-15)

    def test_feedback_changes_the_selected_action_in_both_models(self):
        for model in ("real", "complex"):
            for _, history, _ in record_tree(action_tables(model)):
                self.assertEqual(history[0].action, "joint_parity")
                self.assertNotEqual(history[1].action, "joint_parity")

    def test_audit_records_are_immutable_and_prediction_precedes_observation(self):
        from dataclasses import FrozenInstanceError
        history = record_tree(action_tables("real"))[0][1]
        with self.assertRaises(FrozenInstanceError):
            history[0].outcome = 1
        np.testing.assert_allclose(history[0].prediction, (.5, .5), atol=2e-16, rtol=0)
        self.assertNotEqual(history[0].posterior, (.5, .5))

    def test_impossible_evidence_is_rejected(self):
        with self.assertRaises(ValueError):
            predictive_update([1., 0.], np.eye(2), 1)

    def test_mixture_linearity_holds_for_one_fixed_record_policy(self):
        for model in ("real", "complex"):
            leaves = record_tree(action_tables(model), depth=4)
            for mixing in (.1, .7):
                probabilities = [mixing*l[0]+(1-mixing)*l[1] for _, _, l in leaves]
                self.assertAlmostEqual(sum(probabilities), 1., places=14)
                for (_, history, likelihood), probability in zip(leaves, probabilities):
                    direct = np.ones(2)
                    for record in history:
                        direct *= action_tables(model)[record.action][:, record.outcome]
                    self.assertAlmostEqual(probability, np.array([mixing, 1-mixing]) @ direct, places=14)

    def test_physical_branch_updates_preserve_each_candidate_positive_cone(self):
        for complex_input in (False, True):
            rng = np.random.default_rng(580+complex_input)
            matrix = rng.normal(size=(4, 4))
            if complex_input:
                matrix = matrix+1j*rng.normal(size=(4, 4))
            rho = matrix @ matrix.conj().T
            rho /= np.trace(rho)
            branches = []
            for outcome in (0, 1):
                output = sum(np.kron(k, IDENTITY) @ rho @ np.kron(k, IDENTITY).conj().T
                             for k in branch_kraus(.31, outcome, .73))
                self.assertGreaterEqual(np.linalg.eigvalsh(output).min(), -2e-16)
                if not complex_input:
                    np.testing.assert_allclose(output.imag, 0, atol=0)
                branches.append(output)
            self.assertAlmostEqual(np.trace(sum(branches)).real, 1., places=14)

    def test_role_covariance_and_nonsignaling_hold_in_both_state_cones(self):
        swap = np.eye(4)[[0, 2, 1, 3]]
        for model in ("real", "complex"):
            rho = (np.eye(4)+.3*pauli_word("XZ")+.2*pauli_word("YY"))/4
            if model == "complex":
                rho += .1*pauli_word("YI")/4
            for a in ("X", "Z") if model == "real" else ("X", "Y", "Z"):
                marginals = []
                for b in ("X", "Z") if model == "real" else ("X", "Y", "Z"):
                    ea = local_effect(a, 1)
                    values = []
                    for t in (-1, 1):
                        eb = local_effect(b, t)
                        p = np.trace(rho @ np.kron(ea, eb))
                        exchanged = np.trace((swap @ rho @ swap) @ np.kron(eb, ea))
                        self.assertAlmostEqual(p.real, exchanged.real, places=14)
                        values.append(p.real)
                    marginals.append(sum(values))
                np.testing.assert_allclose(marginals, marginals[0], atol=2e-16)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CognitiveModelAcceptanceTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"round": 58, "scope": "Finite calibrated hypothesis controller; fresh independent trials; candidate matrix rules assumed",
              "both_models_support": ["append_only_records", "prediction_action_feedback", "Bayesian_learning", "fixed_policy_mixtures", "positive_conditional_updates", "role_covariance", "nonsignaling", "joint_calibration"],
              "local_tomography_entailed_by_these_checks": False,
              "controller_examples": {model: {"first_path": [r.action for r in record_tree(action_tables(model))[0][1]],
                                                    "first_posterior": record_tree(action_tables(model))[0][1][0].posterior} for model in ("real", "complex")},
              "quantum_theory_derived_from_cognition": False,
              "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("cognitive_model_acceptance_results.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
