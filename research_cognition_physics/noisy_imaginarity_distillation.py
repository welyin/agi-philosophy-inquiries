"""Purify a nonzero Y resource using the old real gate and noisy Z readout.

Independent preparation of the new conditional states is an explicit added
closure assumption. Success probabilities and resource costs are retained.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from bipartite_composition import interaction
from certified_intervals import Interval as I, SCALE, ceil_div, sin_interval
from imaginary_readout_closure import y_effect
from quantum_interface_audit import ALPHA, IDENTITY, PAULI_Y, branch_kraus


def resource_state(bias):
    if not -1 <= bias <= 1:
        raise ValueError("Use a bias in [-1,1].")
    return (IDENTITY + bias*PAULI_Y)/2


def purify_pair(first_bias, second_bias, readout_contrast=1., outcome=1):
    """Keep system A, read/discard B. Return subnormalized A state."""
    gate = interaction(math.pi/2, "YX")
    density = gate @ np.kron(resource_state(first_bias), resource_state(second_bias)) @ gate.conj().T
    output = np.zeros((4, 4), dtype=complex)
    for k in branch_kraus(0., outcome, readout_contrast):
        joint = np.kron(IDENTITY, k)
        output += joint @ density @ joint.conj().T
    return np.einsum("abcb->ac", output.reshape(2, 2, 2, 2))


def updated_bias(first_bias, second_bias, beta, sign=1):
    probability = (1+sign*beta*first_bias*second_bias)/2
    if probability == 0:
        raise ValueError("The selected branch has zero probability.")
    return (first_bias+sign*beta*second_bias)/(2*probability), probability


def error_step(delta, beta):
    """Exact arithmetic is possible when delta and beta are Fractions."""
    return delta*(1-beta+beta*delta)/(1+beta*(1-delta)**2)


def interval_from_endpoints(lower, upper):
    return I((lower*SCALE).numerator // (lower*SCALE).denominator,
             ceil_div((upper*SCALE).numerator, (upper*SCALE).denominator))


def distillation_certificate(initial_bias=Fraction(1, 5), levels=10):
    initial_bias = Fraction(initial_bias)
    if not 0 < initial_bias < 1 or levels < 0:
        raise ValueError("Use nonzero imperfect initial bias and a nonnegative level count.")
    beta = 4*sin_interval(I.rational(1, 4))
    b_lo, b_hi = Fraction(beta.lo, SCALE), Fraction(beta.hi, SCALE)
    delta = I.exact(1-initial_bias)
    cost = I.exact(1)
    rows = []
    for level in range(levels+1):
        d_lo, d_hi = Fraction(delta.lo, SCALE), Fraction(delta.hi, SCALE)
        error_upper = d_hi/2
        rows.append({"level": level, "one_minus_bias_interval": delta.floats(),
                     "trace_distance_upper_exact": str(error_upper),
                     "trace_distance_interval": (delta/2).floats(),
                     "expected_raw_state_cost_interval": cost.floats(),
                     "expected_new_y_readout_uses_upper_exact": str(Fraction(2*cost.hi, SCALE)),
                     "finite_full_tree_raw_states": 2**level,
                     "raw_state_budget_for_failure_at_most_one_percent": ceil_div(100*cost.hi, SCALE),
                     "certifies_error_below_1e_minus_6": error_upper < Fraction(1, 10**6),
                     "certifies_error_below_1e_minus_12": error_upper < Fraction(1, 10**12)})
        p_lo = (1+b_lo*(1-d_hi)**2)/2
        p_hi = (1+b_hi*(1-d_lo)**2)/2
        probability = interval_from_endpoints(p_lo, p_hi)
        cost = 2*cost/probability
        # Error increases with delta, decreases with beta on the unit square.
        delta = interval_from_endpoints(error_step(d_lo, b_hi), error_step(d_hi, b_lo))
    contraction = (1-b_lo*initial_bias)/(1+b_lo*initial_bias**2)
    return {"initial_bias_exact": str(initial_bias), "old_readout_contrast": 1,
            "beta_interval": beta.floats(), "global_geometric_contraction_upper_exact": str(contraction),
            "rows": rows, "cost_is_expected_with_local_retries": True,
            "raw_state_cost_does_not_count_failed_initial_y_heralds": True,
            "new_y_effect_preparation_success_probability": "1/2"}


class NoisyImaginarityDistillationTests(unittest.TestCase):
    def test_actual_real_gate_and_original_instrument_match_both_branches(self):
        for first, second, eta, outcome in product((.1, .4, .9), (-.3, .2, .7), (.2, .85, 1.), (0, 1)):
            bias, probability = updated_bias(first, second, ALPHA*eta, 2*outcome-1)
            np.testing.assert_allclose(purify_pair(first, second, eta, outcome), probability*resource_state(bias), atol=4e-16)

    def test_failure_records_and_success_records_keep_total_probability(self):
        for first, second in product((.2, .8), repeat=2):
            outputs = [purify_pair(first, second, 1., outcome) for outcome in (0, 1)]
            self.assertAlmostEqual(sum(np.trace(o).real for o in outputs), 1., places=14)
            np.testing.assert_allclose(sum(outputs), resource_state(first), atol=4e-16)

    def test_any_positive_bias_increases_under_any_positive_old_readout(self):
        for numerator, denominator in product(range(1, 10), (10, 11)):
            bias = Fraction(numerator, denominator)
            for beta in (Fraction(1, 10), Fraction(1, 2), Fraction(9, 10)):
                updated, probability = updated_bias(bias, bias, beta)
                self.assertGreater(updated, bias)
                self.assertLess(updated, 1)
                self.assertGreater(probability, Fraction(1, 2))
                self.assertEqual(1-updated, error_step(1-bias, beta))

    def test_zero_imaginary_resource_or_zero_old_readout_gives_no_amplification(self):
        self.assertEqual(updated_bias(0, 0, Fraction(1, 2))[0], 0)
        self.assertEqual(updated_bias(Fraction(1, 5), Fraction(1, 5), 0)[0], Fraction(1, 5))
        self.assertEqual(updated_bias(1, 1, Fraction(1, 2))[0], 1)

    def test_geometric_bound_controls_the_entire_sequence(self):
        initial, beta = Fraction(1, 5), Fraction(1, 2)
        contraction = (1-beta*initial)/(1+beta*initial**2)
        delta = 1-initial
        for level in range(9):
            self.assertLessEqual(delta, (1-initial)*contraction**level)
            delta = error_step(delta, beta)

    def test_monotone_interval_recurrence_contains_exact_rational_iterations(self):
        exact, beta = Fraction(4, 5), Fraction(1, 2)
        interval = I.exact(exact)
        for _ in range(9):
            self.assertLessEqual(Fraction(interval.lo, SCALE), exact)
            self.assertGreaterEqual(Fraction(interval.hi, SCALE), exact)
            exact = error_step(exact, beta)
            interval = interval_from_endpoints(error_step(Fraction(interval.lo, SCALE), beta),
                                               error_step(Fraction(interval.hi, SCALE), beta))

    def test_noisy_readout_certificate_reaches_arbitrary_precision_without_rounding_to_one(self):
        report = distillation_certificate()
        self.assertLess(Fraction(report["global_geometric_contraction_upper_exact"]), 1)
        good = [row for row in report["rows"] if row["certifies_error_below_1e_minus_12"]]
        self.assertTrue(good)
        first = good[0]
        self.assertGreater(Fraction(first["trace_distance_upper_exact"]), 0)
        self.assertGreater(first["one_minus_bias_interval"][0], 0)

    def test_expected_cost_recursion_and_bounded_tree_success_are_distinct(self):
        bias, beta, cost, tree_success = Fraction(1, 5), Fraction(1, 2), Fraction(1), Fraction(1)
        for level in range(5):
            self.assertGreater(tree_success, 0)
            self.assertGreaterEqual(cost, 2**level)
            self.assertLessEqual(cost, 4**level)
            bias, probability = updated_bias(bias, bias, beta)
            cost = 2*cost/probability
            tree_success = tree_success**2 * probability

    def test_four_successive_levels_cross_check_the_state_recurrence(self):
        bias = .2
        for _ in range(4):
            expected, probability = updated_bias(bias, bias, ALPHA)
            output = purify_pair(bias, bias)
            self.assertAlmostEqual(np.trace(output).real, probability, places=14)
            self.assertAlmostEqual(np.trace(output @ PAULI_Y).real/probability, expected, places=14)
            bias = expected

    def test_full_rank_resource_cannot_yield_a_nonreal_rank_one_output_via_real_filters(self):
        density = np.kron(resource_state(.2), resource_state(.2))
        self.assertGreater(np.linalg.eigvalsh(density).min(), 0)
        rng = np.random.default_rng(55)
        for _ in range(10):
            matrix = rng.normal(size=(2, 4))
            output = matrix @ density @ matrix.T
            self.assertGreater(np.linalg.eigvalsh(output).min(), 0)
            real_direction = rng.normal(size=2)
            rank_one = np.outer(real_direction, rng.normal(size=4))
            output = rank_one @ density @ rank_one.T
            rounding_bound = 8*np.finfo(float).eps*np.linalg.norm(rank_one, ord=2)**2*np.linalg.norm(density, ord=2)
            np.testing.assert_allclose(output.imag, 0, atol=rounding_bound)
            np.testing.assert_allclose(output/np.trace(output), np.outer(real_direction, real_direction)/np.dot(real_direction, real_direction), atol=3e-16)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(NoisyImaginarityDistillationTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"scope": "Independent copies of the added conditional Y state, usual tensor composition, old real YX gate and original noisy fresh Z readout",
              "equal_input_success_update": "nu' = nu*(1+beta)/(1+beta*nu^2)",
              "success_probability": "(1+beta*nu^2)/2",
              "error_recurrence": "delta' = delta*(1-beta+beta*delta)/(1+beta*(1-delta)^2)",
              "every_nonzero_bias_purifies_in_closure": True,
              "ideal_parity_readout_required": False,
              "finite_full_rank_resources_yield_exact_nonreal_pure_state_by_real_processing": False,
              "certificate": distillation_certificate(),
              "optimal_resource_cost_proved": False,
              "quantum_theory_derived_from_cognition": False,
              "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("noisy_imaginarity_distillation_results.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
