"""Round 61: account for shared preparation and finite record costs.

Real separability is defined using real local density matrices. The sampling
comparison concerns the specified one-read-per-copy experiments only.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I, SCALE, sin_interval
from role_symmetry_and_swap import pauli_word
from shared_real_reference_tomography import shared_reference


ALPHA_LOWER = Fraction(989615, 1000000)


def trace_distance(first, second):
    return float(np.abs(np.linalg.eigvalsh(first-second)).sum()/2)


def yy_word(systems, first, second):
    return "".join("Y" if site in (first, second) else "I" for site in range(systems))


def pair_sum_witness(systems):
    return sum((pauli_word(yy_word(systems, i, j)) for i, j in combinations(range(systems), 2)), start=np.zeros((2**systems, 2**systems), dtype=complex))


def biseparable_bound(systems):
    if systems < 2:
        raise ValueError("At least two parties are required.")
    return (systems-1)*(systems-2)//2


def majority_error_exact(count, visibility):
    """Optimal equal-prior error for an odd number of binary parity records."""
    visibility = Fraction(visibility)
    if count < 1 or count % 2 != 1 or not 0 <= visibility <= 1:
        raise ValueError("Use positive odd count and visibility in [0,1].")
    error = (1-visibility)/2
    if error == 0:
        return Fraction(0)
    a, denominator = error.numerator, error.denominator
    b = denominator-a
    start = count//2+1
    term = math.comb(count, start)*a**start*b**(count-start)
    total = term
    for k in range(start, count):
        numerator = term*(count-k)*a
        divisor = (k+1)*b
        term, remainder = divmod(numerator, divisor)
        assert remainder == 0
        total += term
    return Fraction(total, denominator**count)


def rounded_upper(value, digits=18):
    scale = 10**digits
    numerator = value.numerator*scale
    return Fraction((numerator+value.denominator-1)//value.denominator, scale)


def minimum_odd_count(visibility, target=Fraction(1, 100)):
    visibility, target = Fraction(visibility), Fraction(target)
    if not 0 < visibility <= 1 or not 0 < target < Fraction(1, 2):
        raise ValueError("Use positive visibility and target error between zero and one half.")
    low, high = -1, 0  # count = 2*index+1
    while majority_error_exact(2*high+1, visibility) > target:
        low, high = high, 2*high+1
    while high-low > 1:
        middle = (low+high)//2
        if majority_error_exact(2*middle+1, visibility) <= target:
            high = middle
        else:
            low = middle
    count = 2*high+1
    error = majority_error_exact(count, visibility)
    previous = majority_error_exact(count-2, visibility) if count > 1 else None
    previous_lower = Fraction(previous.numerator*10**18//previous.denominator, 10**18) if previous is not None else None
    return {"visibility_lower_exact": str(visibility), "count": count,
            "error_upper_exact": str(rounded_upper(error)),
            "previous_odd_error_lower_exact": str(previous_lower) if previous_lower is not None else None,
            "previous_odd_count_fails_for_this_visibility": count == 1 or previous > target,
            "target_exact": str(target)}


def fidelity_necessary_count(visibility, target=Fraction(1, 100)):
    """Exact integer inversion of (1-v^2)^N <= 4 e (1-e)."""
    visibility, target = Fraction(visibility), Fraction(target)
    if not 0 < visibility <= 1 or not 0 < target < Fraction(1, 2):
        raise ValueError("Use nonzero visibility and a nontrivial error target.")
    base, threshold = 1-visibility**2, 4*target*(1-target)
    low, high = 0, 1
    while base**high > threshold:
        low, high = high, 2*high
    while high-low > 1:
        middle = (low+high)//2
        if base**middle <= threshold:
            high = middle
        else:
            low = middle
    return high


class SourceIndependenceAuditTests(unittest.TestCase):
    def test_rational_old_visibility_lower_bound_is_certified(self):
        interval = 4*sin_interval(I.rational(1, 4))
        self.assertLess(ALPHA_LOWER, Fraction(interval.lo, SCALE))
        self.assertLess(Fraction(interval.hi, SCALE)-ALPHA_LOWER, Fraction(1, 1000000))

    def test_every_real_product_across_a_cut_has_zero_cross_yy_witness(self):
        rng = np.random.default_rng(61)
        for n in range(2, 6):
            for cut in range(1, n):
                a = rng.normal(size=(2**cut, 2**cut))
                b = rng.normal(size=(2**(n-cut), 2**(n-cut)))
                a, b = a @ a.T, b @ b.T
                rho = np.kron(a/np.trace(a), b/np.trace(b))
                witness = pauli_word(yy_word(n, 0, cut))
                self.assertAlmostEqual(np.trace(rho @ witness).real, 0, places=15)

    def test_distance_to_real_separable_states_across_each_fixed_cut_is_at_least_half_and_attained(self):
        for n in range(2, 6):
            for cut in range(1, n):
                target = shared_reference(n)
                candidate = np.kron(shared_reference(cut), shared_reference(n-cut))
                witness = pauli_word(yy_word(n, 0, cut))
                self.assertEqual(np.trace(target @ witness).real, 1)
                self.assertEqual(np.trace(candidate @ witness).real, 0)
                self.assertAlmostEqual(trace_distance(target, candidate), .5, places=14)

    def test_biseparable_pair_sum_bound_is_sharp_for_every_tested_size(self):
        for n in range(2, 6):
            witness = pair_sum_witness(n)
            self.assertEqual(np.trace(shared_reference(n) @ witness).real, math.comb(n, 2))
            for cut in range(1, n):
                candidate = np.kron(shared_reference(cut), shared_reference(n-cut))
                value = np.trace(candidate @ witness).real
                self.assertEqual(value, math.comb(cut, 2)+math.comb(n-cut, 2))
                self.assertLessEqual(value, biseparable_bound(n))
            saturating = np.kron(shared_reference(n-1), np.eye(2)/2)
            self.assertEqual(np.trace(saturating @ witness).real, biseparable_bound(n))

    def test_white_noise_witness_threshold_is_exact(self):
        for n in range(2, 6):
            witness = pair_sum_witness(n)
            threshold = Fraction(n-2, n)
            for c in (threshold, (1+threshold)/2):
                rho = float(c)*shared_reference(n)+(1-float(c))*np.eye(2**n)/2**n
                value = np.trace(rho @ witness).real
                self.assertAlmostEqual(value, float(c)*math.comb(n, 2), places=14)
                if c > threshold:
                    self.assertGreater(value, biseparable_bound(n))

    def test_majority_formula_agrees_with_exhaustive_record_likelihood_decisions(self):
        from itertools import product
        for count in (1, 3, 5, 7):
            for visibility in (Fraction(0), Fraction(1, 25), Fraction(2, 5), Fraction(1)):
                p = (1+visibility)/2
                error = Fraction(0)
                for path in product((0, 1), repeat=count):
                    k = sum(path)
                    plus = p**k*(1-p)**(count-k)
                    minus = (1-p)**k*p**(count-k)
                    error += min(plus, minus)/2
                self.assertEqual(error, majority_error_exact(count, visibility))

    def test_exact_majority_error_decreases_with_each_pair_of_added_samples(self):
        for v in (Fraction(1, 25), Fraction(1, 5), Fraction(4, 5)):
            p = (1-v)/2
            for n in (1, 3, 5, 9):
                decrease = majority_error_exact(n, v)-majority_error_exact(n+2, v)
                m = (n-1)//2
                self.assertEqual(decrease, v*math.comb(n, m)*(p*(1-p))**(m+1))
                self.assertGreater(decrease, 0)

    def test_fidelity_necessary_count_is_certified_by_neighboring_integers(self):
        v, e = Fraction(1, 25), Fraction(1, 100)
        count = fidelity_necessary_count(v, e)
        self.assertGreater((1-v*v)**(count-1), 4*e*(1-e))
        self.assertLessEqual((1-v*v)**count, 4*e*(1-e))

    def test_conservative_old_protocol_sample_budgets(self):
        self.assertEqual(minimum_odd_count(ALPHA_LOWER)["count"], 1)
        self.assertEqual(minimum_odd_count(ALPHA_LOWER**2)["count"], 3)
        self.assertGreater(majority_error_exact(1, ALPHA_LOWER**2), Fraction(1, 100))
        weak = minimum_odd_count(Fraction(1, 25))
        self.assertEqual(weak["count"], 3381)
        self.assertLess(Fraction(weak["error_upper_exact"]), Fraction(1, 100))
        self.assertGreater(Fraction(weak["previous_odd_error_lower_exact"]), Fraction(1, 100))

    def test_zero_visibility_has_no_learning_and_invalid_counts_are_rejected(self):
        self.assertEqual(majority_error_exact(101, 0), Fraction(1, 2))
        with self.assertRaises(ValueError):
            minimum_odd_count(0)
        with self.assertRaises(ValueError):
            majority_error_exact(4, Fraction(1, 5))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SourceIndependenceAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"round": 61, "fixed_cut_real_separable_distance_to_reference": "1/2, sharp",
              "biseparable_pair_sum_bound": "(n-1)(n-2)/2, sharp", "reference_pair_sum": "n(n-1)/2",
              "white_noise_sufficient_genuine_real_entanglement_threshold": "c>(n-2)/n",
              "shared_classical_randomness_can_prepare_reference_from_real_local_states": False,
              "sample_budget_scope": "Fresh independent target copies, specified single parity record per copy, equal priors; ancillary resources counted separately",
              "one_percent_error_examples": {
                  "old_joint_gate": minimum_odd_count(ALPHA_LOWER),
                  "weak_local_y_visibility_one_fifth": minimum_odd_count(Fraction(1, 25)),
                  "old_local_operations_with_one_shared_real_reference_pair_per_copy": minimum_odd_count(ALPHA_LOWER**2)},
              "weak_local_y_necessary_sample_count_from_fidelity": fidelity_necessary_count(Fraction(1, 25)),
              "weak_local_y_fixed_protocol_scaling": "Theta(nu^-4) for fixed nontrivial error target as nu approaches zero",
              "global_optimal_cost_over_all_resource_distillation_protocols_claimed": False,
              "source_independence_derived_from_cognition": False,
              "quantum_theory_derived_from_cognition": False,
              "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("source_independence_audit_results.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
