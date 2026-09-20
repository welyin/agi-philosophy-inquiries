"""Round 219: finite report experiments, task preservation and source identity.

All scientific arithmetic is rational. The general Blackwell equivalence and
approximation bounds are proved in the note; finite checks do not establish
a universal classification. These channels describe public experimental
reports, not a hidden-variable theory for all physical systems.
"""
import argparse
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import unittest


def channel(rows):
    matrix = tuple(tuple(F(x) for x in row) for row in rows)
    if not matrix or not matrix[0] or any(len(row) != len(matrix[0]) for row in matrix):
        raise ValueError('A nonempty rectangular matrix is required.')
    if any(x < 0 for row in matrix for x in row) or any(sum(row) != 1 for row in matrix):
        raise ValueError('Each channel row must be a probability distribution.')
    return matrix


def multiply(a, b):
    if len(a[0]) != len(b):
        raise ValueError('Incompatible matrix dimensions.')
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(len(b)))
                       for j in range(len(b[0]))) for i in range(len(a)))


def deterministic_map(labels, output_size):
    return channel([[int(label == j) for j in range(output_size)] for label in labels])


def row_tv(a, b):
    if len(a) != len(b) or len(a[0]) != len(b[0]):
        raise ValueError('Compare the same source and report alphabets.')
    return max(sum(abs(x - y) for x, y in zip(ra, rb)) / 2 for ra, rb in zip(a, b))


def binary_channel(error):
    p = F(error)
    if not 0 <= p <= F(1, 2):
        raise ValueError('The calibrated binary error lies in [0,1/2].')
    return channel([[1-p, p], [p, 1-p]])


def binary_simulation(source_error, target_error):
    """Return an optimal source-independent decoder and its exact deficiency."""
    p, q = F(source_error), F(target_error)
    binary_channel(p)
    binary_channel(q)
    decoder_error = (q-p)/(1-2*p) if q >= p and p < F(1, 2) else F(0)
    return binary_channel(decoder_error), max(F(0), p-q)


def best_payoff(experiment, utilities, prior=None):
    n = len(experiment)
    prior = tuple(F(x) for x in prior) if prior is not None else (F(1, n),) * n
    if len(prior) != n or sum(prior) != 1 or any(x < 0 for x in prior):
        raise ValueError('A probability prior on the declared source labels is required.')
    if len(utilities) != n or not utilities[0] or any(len(row) != len(utilities[0]) for row in utilities):
        raise ValueError('Utilities must use the source alphabet and a common action alphabet.')
    return sum(max(sum(prior[x] * experiment[x][y] * utilities[x][a] for x in range(n))
                   for a in range(len(utilities[0]))) for y in range(len(experiment[0])))


def independent_joint(a, b):
    """Conditional independence given the source is an explicit input here."""
    if len(a) != len(b):
        raise ValueError('Joint experiments must have the same source labels.')
    return channel([[pa * pb for pa in ra for pb in rb] for ra, rb in zip(a, b)])


def posterior(experiment, report, prior=None):
    prior = (F(1, len(experiment)),) * len(experiment) if prior is None else prior
    weights = tuple(F(prior[x]) * experiment[x][report] for x in range(len(experiment)))
    total = sum(weights)
    if total == 0:
        raise ValueError('A zero-probability report has no defined posterior here.')
    return tuple(x / total for x in weights)


def complementary_experiments():
    sources = tuple(product((0, 1), repeat=2))
    a = deterministic_map([x[0] for x in sources], 2)
    b = deterministic_map([x[1] for x in sources], 2)
    joint = deterministic_map(range(4), 4)
    u_a = tuple(tuple(F(action == x[0]) for action in range(2)) for x in sources)
    u_b = tuple(tuple(F(action == x[1]) for action in range(2)) for x in sources)
    return a, b, joint, u_a, u_b


def timed_report_audit(error=F(1, 4)):
    """One noisy bit arriving before versus after an irrevocable decision.

    Both complete report records contain exactly the same signal. A physical
    decoder at the decision time has access only to the already arrived prefix.
    """
    signal = binary_channel(error)
    early = channel([[row[0], row[1], 0, 0] for row in signal])
    late = channel([[0, 0, row[0], row[1]] for row in signal])
    swap_time = deterministic_map((2, 3, 0, 1), 4)
    reward = ((F(1), F(0)), (F(0), F(1)))
    before_late_arrival = channel([[1], [1]])
    return early, late, swap_time, best_payoff(signal, reward), best_payoff(before_late_arrival, reward)


class CognitiveTaskOrderTests(unittest.TestCase):
    def test_binary_translation_and_exact_error_certificates(self):
        for p, q in ((F(0), F(1, 3)), (F(1, 4), F(1, 8)),
                     (F(1, 4), F(1, 3)), (F(1, 2), F(0)), (F(1, 2), F(1, 2))):
            decoder, error = binary_simulation(p, q)
            actual = row_tv(multiply(binary_channel(p), decoder), binary_channel(q))
            self.assertEqual(actual, error)
            # Data processing bounds the distance between the two simulated rows.
            lower = max(F(0), ((1-2*q) - (1-2*p)) / 2)
            self.assertEqual(actual, lower)

    def test_garbling_preserves_all_enumerated_decision_inequalities(self):
        source = binary_channel(F(1, 8))
        target = multiply(source, binary_channel(F(1, 6)))
        for entries in product((F(0), F(1)), repeat=4):
            utility = (entries[:2], entries[2:])
            for p in (F(0), F(1, 4), F(1, 2), F(1)):
                self.assertGreaterEqual(best_payoff(source, utility, (p, 1-p)),
                                        best_payoff(target, utility, (p, 1-p)))

    def test_approximation_controls_payoff_loss_and_has_a_tight_case(self):
        source, target = binary_channel(F(1, 4)), binary_channel(F(1, 8))
        for entries in product((F(0), F(1)), repeat=4):
            utility = (entries[:2], entries[2:])
            for p in (F(0), F(1, 4), F(1, 2), F(1)):
                self.assertLessEqual(best_payoff(target, utility, (p, 1-p)) -
                                     best_payoff(source, utility, (p, 1-p)), F(1, 8))
        identity_reward = ((F(1), F(0)), (F(0), F(1)))
        self.assertEqual(best_payoff(target, identity_reward)-best_payoff(source, identity_reward), F(1, 8))

    def test_complementary_agents_are_incomparable_and_joint_retains_both(self):
        a, b, joint, ua, ub = complementary_experiments()
        self.assertEqual((best_payoff(a, ua), best_payoff(b, ua)), (F(1), F(1, 2)))
        self.assertEqual((best_payoff(a, ub), best_payoff(b, ub)), (F(1, 2), F(1)))
        self.assertEqual(multiply(joint, a), a)
        self.assertEqual(multiply(joint, b), b)
        self.assertEqual(best_payoff(joint, ua), 1)
        self.assertEqual(best_payoff(joint, ub), 1)

    def test_incomparability_has_exact_simulation_error_witness(self):
        a, b, _, _, _ = complementary_experiments()
        fair = binary_channel(F(1, 2))
        self.assertEqual(row_tv(multiply(a, fair), b), F(1, 2))
        self.assertEqual(a[0], a[1])
        self.assertEqual(sum(abs(x-y) for x, y in zip(b[0], b[1])) / 2, 1)
        # Equal source rows remain equal after any decoder; two target rows at
        # distance one force at least one row error >= 1/2 by the triangle inequality.

    def test_duplicate_source_is_equivalent_but_not_independent_evidence(self):
        base = binary_channel(F(1, 4))
        duplicate = multiply(base, deterministic_map((0, 3), 4))
        recover = deterministic_map((0, 0, 1, 1), 2)
        self.assertEqual(multiply(duplicate, recover), base)
        independent = independent_joint(base, base)
        self.assertEqual(posterior(duplicate, 3)[1], F(3, 4))
        self.assertEqual(posterior(independent, 3)[1], F(9, 10))
        self.assertNotEqual(duplicate, independent)

    def test_consensus_by_erasure_can_destroy_task_ability(self):
        exact = binary_channel(0)
        erased = multiply(exact, deterministic_map((0, 0), 2))
        reward = ((F(1), F(0)), (F(0), F(1)))
        self.assertEqual(best_payoff(exact, reward), 1)
        self.assertEqual(best_payoff(erased, reward), F(1, 2))
        self.assertEqual(erased[0], erased[1])

    def test_bayes_optimization_matches_independent_policy_enumeration(self):
        experiment = channel([[F(1, 2), F(1, 3), F(1, 6)],
                              [F(1, 5), F(3, 5), F(1, 5)]])
        utility = ((F(1), F(1, 3)), (F(0), F(1)))
        prior = (F(2, 5), F(3, 5))
        values = []
        for choices in product(range(2), repeat=3):
            policy = deterministic_map(choices, 2)
            law = multiply(experiment, policy)
            values.append(sum(prior[x]*law[x][a]*utility[x][a] for x in range(2) for a in range(2)))
        self.assertEqual(best_payoff(experiment, utility, prior), max(values))

    def test_offline_equivalence_does_not_preserve_online_decisions(self):
        for error in (F(0), F(1, 4), F(1, 2)):
            early, late, swap, early_score, late_score = timed_report_audit(error)
            self.assertEqual(multiply(early, swap), late)
            self.assertEqual(multiply(late, swap), early)
            self.assertEqual(early_score, 1-error)
            self.assertEqual(late_score, F(1, 2))
            fair = channel([[F(1, 2), F(1, 2)]] * 2)
            self.assertEqual(row_tv(fair, binary_channel(error)), F(1, 2)-error)
            # A constant prefix cannot produce source-dependent early rows.
            # Their distance 1-2*error gives the matching lower bound / 2.


def report():
    a, b, joint, ua, ub = complementary_experiments()
    base = binary_channel(F(1, 4))
    duplicate = multiply(base, deterministic_map((0, 3), 4))
    return {
        'round': 219, 'date': '2026-09-20',
        'starting_commitments': ['testable prediction and feedback',
                                 'preserve declared old task abilities during expansion',
                                 'account for source identity and internal resources'],
        'extra_model_inputs': ['finite fixed experiments and classical public reports',
                              'shared experimental source labels, inaccessible to a decoder',
                              'randomized report processing independent of the source label',
                              'expected payoff and all finite terminal decision tasks for the equivalence theorem'],
        'known_theorem': 'Finite Blackwell equivalence: task dominance iff a source-independent stochastic decoder exists.',
        'theorem_claimed_original': False,
        'revision': 'r2: scope corrected; numerical definitions and checks preserved',
        'inherited_reconstruction': 'Round189: full framework plus U, C and L implies complex matrix state cones. Ordinary real QM fails L and nontrivial classical theory fails C; report-level examples do not reverse these exclusions.',
        'new_application': 'Report-level task-preserving translations and a partial order, as supporting tools rather than replacements for the full reconstruction premises.',
        'results': {
            'complementary_A_task_scores': [str(best_payoff(a, u)) for u in (ua, ub)],
            'complementary_B_task_scores': [str(best_payoff(b, u)) for u in (ua, ub)],
            'joint_task_scores': [str(best_payoff(joint, u)) for u in (ua, ub)],
            'A_to_B_and_B_to_A_deficiency': '1/2',
            'binary_deficiency_formula': 'max(0, source_error-target_error), errors in [0,1/2]',
            'binary_example_source_1_4_target_1_8': str(binary_simulation(F(1,4),F(1,8))[1]),
            'duplicated_report_correct_posterior': str(posterior(duplicate, 3)[1]),
            'independence_assumption_posterior': str(posterior(independent_joint(base,base),3)[1]),
            'early_and_late_full_report_experiments_equivalent': True,
            'early_and_late_online_scores_at_error_1_4': [str(x) for x in timed_report_audit()[3:]],
            'late_to_early_causal_report_error_at_error_1_4': '1/4',
            'perfect_consensus_after_erasure_guess_success': '1/2'},
        'not_derived': ['continuous reversible pure-state transitivity', 'entanglement',
                        'local tomography', 'Jordan/HSD or complex-qubit existence',
                        'a complete physical theory of cognition', 'gravity'],
        'dynamic_scope': 'The new equivalence concerns fixed report experiments and terminal decisions. Sequential intervention requires causal process simulation; reuse rounds 01, 198 and 203 without treating terminal postprocessing as a dynamic simulator.',
        'next': 'Backchain from the round189 reconstruction theorem; distinguish statistical identification from single-copy stochastic measurement simulation before motivating L. Causal report translation remains a supporting tool.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(CognitiveTaskOrderTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = report()
    data['exact_checks'] = {'run': checks.testsRun, 'failures': 0, 'errors': 0}
    if args.write_results:
        Path(__file__).with_name('cognitive_task_order_results.json').write_text(
            json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
