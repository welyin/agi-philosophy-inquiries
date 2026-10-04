"""Round 203: close the different-controller gap explicitly left in round 119.

Reuses the round-118 controller and round-201 code. This is an implementation
audit and an application of old joint-measurement/coupling arguments, not a
derivation of quantum structure or a new optimization programme.
"""

import argparse
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import unittest

import numpy as np
import classical_self_audit as old
import real_prediction_code as code


def record_law(world_prior, controller_prior, noise=F(1, 4), allow_toggle=True):
    """Actual (executed action, sensor outcome) law; reports are not outcomes.

    Both controllers are evaluated on the same actual world prior. The public
    record here deliberately excludes the numerical belief/forecast fields:
    exact equality of such reports is a separate, stronger comparison task.
    """
    world_prior, noise = old.probability(world_prior), old.probability(noise)
    pending = old.commit(controller_prior, noise, allow_toggle)
    propagated = 1-world_prior if pending.executed_action else world_prior
    one = noise+(1-2*noise)*propagated
    return {(a, y): (one if y else 1-one) if a == pending.executed_action else F(0)
            for a, y in product((0, 1), repeat=2)}


def tv(left, right):
    return sum(abs(left.get(k, 0)-right.get(k, 0))
               for k in left.keys() | right.keys())/2


def observed_law(records):
    return {y: sum(mass for (a, b), mass in records.items() if b == y)
            for y in (0, 1)}


def immediate_regret(world_prior, controller_prior):
    """Only the one-step goal x'=1, before sensor noise; not long-run value."""
    p = old.probability(world_prior)
    a = old.commit(controller_prior).executed_action
    achieved = 1-p if a else p
    return max(p, 1-p)-achieved


def soft_toggle_probability(belief, slope=F(2)):
    belief, slope = old.probability(belief), F(slope)
    if slope < 0:
        raise ValueError("Nonnegative slope required.")
    return min(F(1), max(F(0), F(1, 2)+slope*(F(1, 2)-belief)))


def coupling_bound(initial_error, policy_errors, kernel_errors):
    """Requires uniform conditional bounds while histories/states still couple.

    A one-step state approximation does not supply these bounds automatically.
    The common-action kernels must include both outcome and successor state.
    """
    if len(policy_errors) != len(kernel_errors):
        raise ValueError("One policy and one kernel bound per step required.")
    survival = 1-old.probability(initial_error)
    for policy_error, kernel_error in zip(policy_errors, kernel_errors):
        survival *= (1-old.probability(policy_error))*(1-old.probability(kernel_error))
    return 1-survival


def marked_record_law(policy_errors, kernel_errors):
    """Finite witness attaining the product bound against the all-zero record."""
    if len(policy_errors) != len(kernel_errors):
        raise ValueError("Mismatched horizons.")
    law = {(): F(1)}
    for ep, de in zip(policy_errors, kernel_errors):
        ep, de = old.probability(ep), old.probability(de)
        law = {h+(a, y): mass*(ep if a else 1-ep)*(de if y else 1-de)
               for h, mass in law.items() for a, y in product((0, 1), repeat=2)}
    return law


def square_joint_scores():
    """Reuse the old square code with F_ab=rho_ab/2; no new source access."""
    states = code.code_states()
    effects = {s: rho/2 for s, rho in states.items()}
    per_bit = whole_pair = 0.
    for source, rho in states.items():
        for report, effect in effects.items():
            probability = float(np.trace(rho @ effect))/4
            per_bit += probability*sum(a == b for a, b in zip(source, report))/2
            whole_pair += probability*(source == report)
    return per_bit, whole_pair


class FeedbackPolicyAuditTests(unittest.TestCase):
    def test_arbitrarily_small_belief_error_can_change_all_action_records(self):
        for denominator in (10, 10000, 10**9):
            epsilon = F(1, denominator)
            p, q = F(1, 2)-epsilon, F(1, 2)
            self.assertEqual(abs(p-q), epsilon)
            self.assertEqual(tv(record_law(p, p), record_law(p, q)), 1)

    def test_sensor_statistics_and_immediate_value_have_different_errors(self):
        epsilon = F(1, 10000)
        p, q = F(1, 2)-epsilon, F(1, 2)
        for noise in (F(0), F(1, 4), F(1, 2), F(3, 4)):
            left, right = record_law(p, p, noise), record_law(p, q, noise)
            self.assertEqual(sum(left.values()), 1)
            self.assertEqual(sum(right.values()), 1)
            self.assertEqual(tv(observed_law(left), observed_law(right)),
                             2*abs(1-2*noise)*epsilon)
        self.assertEqual(immediate_regret(p, q), 2*epsilon)

    def test_veto_equalizes_executed_actions_but_not_requested_actions(self):
        p, q = F(2, 5), F(1, 2)
        self.assertNotEqual(old.commit(p, allow_toggle=False).requested_action,
                            old.commit(q, allow_toggle=False).requested_action)
        self.assertEqual(tv(record_law(p, p, allow_toggle=False),
                            record_law(p, q, allow_toggle=False)), 0)

    def test_margin_and_regret_for_all_rational_grid_pairs(self):
        for i, j in product(range(33), repeat=2):
            p, q = F(i, 32), F(j, 32)
            epsilon = abs(p-q)
            if abs(p-F(1, 2)) > epsilon:
                self.assertEqual(old.commit(p).executed_action,
                                 old.commit(q).executed_action)
            self.assertLessEqual(immediate_regret(p, q), 2*epsilon)

    def test_soft_policy_is_lipschitz_including_clipping(self):
        for slope in (F(0), F(1, 3), F(2), F(7)):
            for i, j in product(range(17), repeat=2):
                p, q = F(i, 16), F(j, 16)
                error = abs(soft_toggle_probability(p, slope)-
                            soft_toggle_probability(q, slope))
                self.assertLessEqual(error, min(F(1), slope*abs(p-q)))

    def test_full_finite_records_attain_policy_and_kernel_product_bound(self):
        ep, de = [F(1, 5), F(0), F(2, 7)], [F(1, 3), F(1, 4), F(0)]
        law = marked_record_law(ep, de)
        self.assertEqual(sum(law.values()), 1)
        baseline = {(0,)*6: F(1)}
        self.assertEqual(tv(baseline, law), coupling_bound(F(0), ep, de))
        self.assertEqual(tv(baseline, law), F(5, 7))

    def test_zero_policy_error_reduces_to_old_same_policy_contract(self):
        self.assertEqual(coupling_bound(F(1, 10), [F(0), F(0)],
                                        [F(1, 3), F(1, 4)]), F(11, 20))
        self.assertEqual(coupling_bound(F(1, 10), [], []), F(1, 10))

    def test_old_square_code_joint_measurement_attains_three_quarters(self):
        states = code.code_states()
        np.testing.assert_allclose(sum(states.values())/2, np.eye(2), atol=1e-14)
        for rho in states.values():
            self.assertGreaterEqual(float(np.linalg.eigvalsh(rho/2).min()), -1e-14)
        per_bit, pair = square_joint_scores()
        self.assertAlmostEqual(per_bit, .75, places=14)
        self.assertAlmostEqual(pair, .5, places=14)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checked = unittest.TextTestRunner(verbosity=1).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(FeedbackPolicyAuditTests))
    if not checked.wasSuccessful():
        raise SystemExit(1)
    epsilon = F(1, 10000)
    p, q = F(1, 2)-epsilon, F(1, 2)
    results = {
        "round": 203,
        "purpose": "Historical-result reuse and closure of the explicit round-119 different-policy gap",
        "antecedent_rounds": [20, 98, 99, 117, 118, 119, 198, 201, 202],
        "quantum_or_cognitive_necessity_derived": False,
        "task_record_fields": ["executed_action", "sensor_outcome"],
        "exact_numeric_belief_reports_in_record_comparison": False,
        "same_actual_world_and_sensor_for_both_controllers": True,
        "counterexample": {
            "world_prior": str(p), "approximate_belief": str(q),
            "belief_error": str(epsilon),
            "action_observation_tv": str(tv(record_law(p, p), record_law(p, q))),
            "observation_only_tv": str(tv(observed_law(record_law(p, p)),
                                          observed_law(record_law(p, q)))),
            "immediate_goal_regret": str(immediate_regret(p, q)),
        },
        "conditional_feedback_bound": "1-(1-d0)*product_t((1-epsilon_policy_t)*(1-delta_kernel_t))",
        "uniform_future_belief_error_derived_automatically": False,
        "long_run_regret_bound_proved": False,
        "fixed_square_code_joint_mean_bit_success": square_joint_scores()[0],
        "fixed_square_code_joint_whole_pair_success": square_joint_scores()[1],
        "square_joint_bound_scope": "Fixed round-201 encoding; all joint POVMs, including sequential readouts with only input-independent auxiliaries",
        "new_tests": checked.testsRun,
    }
    if args.write_results:
        Path(__file__).with_name("feedback_policy_audit_results.json").write_text(
            json.dumps(results, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(results, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
