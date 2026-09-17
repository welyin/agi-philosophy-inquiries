"""Round 172: finite classical history retention is weaker than purification.

Exact rational certificates accompany the general construction in note 172.
No quantum state space is assumed in this module.
"""

import argparse
from fractions import Fraction as F
from itertools import product, permutations
import json
from pathlib import Path
import unittest


def function_seed(kernel):
    """Independent random table implements a column-stochastic finite kernel."""
    rows = tuple(tuple(F(x) for x in row) for row in kernel)
    m, n = len(rows), len(rows[0])
    if any(len(row) != n for row in rows):
        raise ValueError("Kernel must be rectangular.")
    if any(x < 0 for row in rows for x in row):
        raise ValueError("Negative probability.")
    if any(sum(rows[y][x] for y in range(m)) != 1 for x in range(n)):
        raise ValueError("Columns must sum to one.")
    seed = {}
    for table in product(range(m), repeat=n):
        probability = F(1)
        for x, y in enumerate(table):
            probability *= rows[y][x]
        seed[table] = probability
    return seed


def reversible_update(x, y, table, output_dimension, inverse=False):
    sign = -1 if inverse else 1
    return x, (y + sign * table[x]) % output_dimension, table


def recovered_kernel(kernel):
    seed = function_seed(kernel)
    m, n = len(kernel), len(kernel[0])
    return [[sum(p for table, p in seed.items()
                 if reversible_update(x, 0, table, m)[1] == y)
             for x in range(n)] for y in range(m)]


class HistoryPurityAuditTests(unittest.TestCase):
    def test_all_columns_reproduced_exactly_including_zeros(self):
        kernels = [
            [[F(2, 3), F(1, 5)], [F(1, 3), F(4, 5)]],
            [[F(1, 2), 0], [F(1, 3), 1], [F(1, 6), 0]],
            [[1, 1], [0, 0]],
        ]
        for kernel in kernels:
            self.assertEqual(recovered_kernel(kernel), kernel)
            self.assertEqual(sum(function_seed(kernel).values()), 1)

    def test_extension_is_permutation_on_full_space(self):
        for m, n in ((2, 3), (3, 2)):
            states = list(product(range(n), range(m), product(range(m), repeat=n)))
            outputs = [reversible_update(*state, m) for state in states]
            self.assertEqual(set(outputs), set(states))
            for state, output in zip(states, outputs):
                self.assertEqual(reversible_update(*output, m, inverse=True), state)

    def test_preserves_external_input_correlations(self):
        kernel = [[F(2, 3), F(1, 5)], [F(1, 3), F(4, 5)]]
        joint = [[F(1, 7), F(2, 7)], [F(3, 7), F(1, 7)]]  # x, external
        observed = [[F(0) for _ in range(2)] for _ in range(2)]
        for table, weight in function_seed(kernel).items():
            for x, external in product(range(2), repeat=2):
                y = reversible_update(x, 0, table, 2)[1]
                observed[y][external] += weight * joint[x][external]
        expected = [[sum(kernel[y][x] * joint[x][e] for x in range(2))
                     for e in range(2)] for y in range(2)]
        self.assertEqual(observed, expected)

    def test_three_step_history_reverses_without_forgetting_seed(self):
        # Includes feedback: each next input is the preceding written result.
        for initial in range(2):
            for tables in product(list(product(range(2), repeat=2)), repeat=3):
                history = [initial]
                for table in tables:
                    history.append(reversible_update(history[-1], 0, table, 2)[1])
                for t in reversed(range(3)):
                    _, blank, _ = reversible_update(history[t], history[t + 1], tables[t], 2, True)
                    self.assertEqual(blank, 0)
                self.assertEqual(history[0], initial)

    def test_every_pure_classical_global_state_has_pure_marginal(self):
        for a, e in product(range(3), range(4)):
            joint = [[F(int((x, y) == (a, e))) for y in range(4)] for x in range(3)]
            self.assertEqual([sum(row) for row in joint], [int(x == a) for x in range(3)])

    def test_pure_seed_permutations_cannot_generate_bitflip_randomness(self):
        target = [F(16, 25), F(9, 25)]
        distances = []
        for permutation in permutations(range(4)):
            for seed in range(2):
                out_x = permutation[seed] // 2  # pure input x=0
                marginal = [F(int(x == out_x)) for x in range(2)]
                distances.append(sum(abs(a - b) for a, b in zip(marginal, target)) / 2)
        self.assertEqual(min(distances), F(9, 25))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(HistoryPurityAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    kernel = [[F(2, 3), F(1, 5)], [F(1, 3), F(4, 5)]]
    report = {
        "round": 172,
        "scope": "ordinary finite classical probability; unrestricted finite auxiliary memory",
        "arbitrary_finite_kernel_has_mixed_seed_permutation_extension": True,
        "general_seed_formula": "w(f)=product_x K(f(x)|x)",
        "example_seed": {str(k): str(v) for k, v in function_seed(kernel).items()},
        "pure_classical_global_state_has_pure_marginal": True,
        "pure_seed_bitflip_counterexample_probability": "9/25",
        "minimum_TV_from_pure_classical_output": "9/25",
        "history_retention_implies_purification": False,
        "memory_or_seed_cost_is_claimed_optimal": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("history_purity_audit_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
