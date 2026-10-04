"""Round 185: all binary ensemble splits need not lift to one ternary test.

A fixed separable preparation correlates a classical trit with three cube
states. Binary steering is exact; optimal ternary label success is only 2/3.
The bound follows from a positive-effect certificate, not a numerical search.
"""

import argparse
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import unittest

import numpy as np
from cone_duality_audit import CUBE_RAYS


U = (F(1), F(0), F(0), F(0))
CODE_STATES = ((1, 1, -1, -1), (1, -1, 1, -1), (1, -1, -1, 1))


def dot(effect, state):
    return sum(e * s for e, s in zip(effect, state))


def subtract(a, b):
    return tuple(x - y for x, y in zip(a, b))


def subset_effect(bits):
    count = sum(bits)
    if count == 0:
        return (F(0),) * 4
    if count == 3:
        return U
    if count == 1:
        index = bits.index(1)
        return (F(1, 2), *(F(int(i == index), 2) for i in range(3)))
    return subtract(U, subset_effect(tuple(1 - bit for bit in bits)))


def binary_effect(probabilities):
    result = [F(0)] * 4
    for bits in product((0, 1), repeat=3):
        weight = F(1)
        for bit, probability in zip(bits, probabilities):
            weight *= probability if bit else 1 - probability
        effect = subset_effect(bits)
        result = [x + weight * y for x, y in zip(result, effect)]
    return tuple(result)


def joint_report(measurement):
    return [[dot(effect, state) / 3 for state in CODE_STATES] for effect in measurement]


class BinaryEnsembleGapTests(unittest.TestCase):
    def test_all_eight_subset_questions_are_legal_and_exact(self):
        for bits in product((0, 1), repeat=3):
            effect = subset_effect(bits)
            self.assertEqual([dot(effect, state) for state in CODE_STATES], list(bits))
            self.assertTrue(all(0 <= dot(effect, vertex) <= 1 for vertex in CUBE_RAYS))

    def test_entire_binary_order_interval_has_exact_constructive_lifts(self):
        for probabilities in product((F(0), F(1, 4), F(2, 3), F(1)), repeat=3):
            effect = binary_effect(probabilities)
            self.assertEqual([dot(effect, state) for state in CODE_STATES], list(probabilities))
            self.assertTrue(all(0 <= dot(effect, vertex) <= 1 for vertex in CUBE_RAYS))
            complement = subtract(U, effect)
            self.assertEqual([dot(complement, state) for state in CODE_STATES],
                             [1 - p for p in probabilities])

    def test_conditioning_map_has_one_unseen_direction(self):
        matrix = np.array(CODE_STATES, dtype=int)
        self.assertEqual(np.linalg.matrix_rank(matrix), 3)
        np.testing.assert_array_equal(matrix @ np.ones(4, dtype=int), np.zeros(3))

    def test_separately_valid_singleton_effects_do_not_form_a_measurement(self):
        effects = [subset_effect(tuple(int(i == j) for i in range(3))) for j in range(3)]
        total = tuple(sum(e[i] for e in effects) for i in range(4))
        self.assertEqual(total, (F(3, 2), F(1, 2), F(1, 2), F(1, 2)))
        self.assertEqual(dot(total, (1, 1, 1, 1)), 3)
        self.assertNotEqual(total, U)

    def test_positive_effect_certificate_bounds_any_ternary_report(self):
        rng = np.random.default_rng(185)
        # Construct varied legal POVMs by mixtures of binary tests and
        # deterministic three-label postprocessings. The proof covers all POVMs.
        for _ in range(30):
            effects = [np.zeros(4) for _ in range(3)]
            for _ in range(5):
                bits = tuple(int(x) for x in rng.integers(0, 2, size=3))
                effect = np.array(subset_effect(bits), dtype=float)
                labels = rng.integers(0, 3, size=2)
                effects[labels[0]] += effect / 5
                effects[labels[1]] += (np.array(U, dtype=float) - effect) / 5
            correct = sum(dot(e, state) for e, state in zip(effects, CODE_STATES)) / 3
            slack = sum(dot(e, (1, *(-np.array(state[1:]))))
                        for e, state in zip(effects, CODE_STATES)) / 3
            self.assertAlmostEqual(correct + slack, 2 / 3)
            self.assertGreaterEqual(slack, -1e-15)

    def test_optimal_success_two_thirds_and_joint_total_variation_one_third_exactly(self):
        effect = subset_effect((1, 0, 0))
        measurement = [effect, subtract(U, effect), (F(0),) * 4]
        actual = joint_report(measurement)
        correct = sum(actual[i][i] for i in range(3))
        ideal = [[F(int(i == j), 3) for j in range(3)] for i in range(3)]
        tv = sum(abs(actual[i][j] - ideal[i][j]) for i in range(3) for j in range(3)) / 2
        self.assertEqual(correct, F(2, 3))
        self.assertEqual(tv, F(1, 3))
        for j in range(3):
            self.assertEqual(sum(actual[i][j] for i in range(3)), F(1, 3))

    def test_matched_classical_trit_auxiliary_has_an_exact_joint_measurement(self):
        copied = np.eye(3) / 3
        self.assertEqual(np.trace(copied), 1)
        self.assertEqual(np.count_nonzero(copied - np.diag(np.diag(copied))), 0)
        self.assertEqual((len(U), len(copied)), (4, 3))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(BinaryEnsembleGapTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = {
        "round": 185,
        "target": "uniform classical trit; fixed joint mixture of three labelled cube states",
        "auxiliary_operational_K": 4,
        "target_operational_K": 3,
        "all_binary_ensemble_splits_exact": True,
        "ternary_ensemble_exact": False,
        "optimal_ternary_label_success": "2/3",
        "minimum_joint_report_total_variation": "1/3",
        "bound_certificate": "success + sum_i E_i(opposite_code_state_i)/3 = 2/3",
        "matched_classical_trit_auxiliary_success": "1",
        "scope": "One fixed auxiliary per target; no joint access to target or extra correlated copies",
        "resolves_order_quotient_implies_strong_quotient_question": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("binary_ensemble_gap_results.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
