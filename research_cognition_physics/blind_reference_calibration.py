"""Round 188: unknown real reference and target have an unidentifiable sign.

Finite-round real LOCC cannot break the simultaneous partial-transpose
ambiguity, even using arbitrarily many copies of the stated source family.
"""

import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path
import unittest

import numpy as np

from reference_permission_audit import (
    H, complete_product_moments, local_product_in_target_order, partial_transpose,
    relation_state, reorder,
)


def many_pairs(parameters):
    state = np.ones((1, 1))
    for value in parameters:
        state = np.kron(state, relation_state(value))
    return state


def positive_reference_interval(moment_interval, reference_interval):
    lo, hi = moment_interval
    lower_c, upper_c = reference_interval
    if not (lo <= hi and 0 < lower_c <= upper_c <= 1):
        raise ValueError("Ordered moment interval and strictly positive reference bounds required.")
    corners = [m / c for m in (lo, hi) for c in (lower_c, upper_c)]
    lower, upper = max(F(-1), min(corners)), min(F(1), max(corners))
    if lower > upper:
        raise ValueError("Intervals are inconsistent with a physical target parameter.")
    return lower, upper


def iid_sample_bound(c, epsilon, delta):
    if not (0 < abs(c) <= 1 and epsilon > 0 and 0 < delta < 1):
        raise ValueError("Nonzero known c, positive precision, and probability delta required.")
    return math.ceil(2 * math.log(2 / delta) / (c * c * epsilon * epsilon))


class BlindReferenceCalibrationTests(unittest.TestCase):
    def test_single_target_reference_pair_all_product_data_depend_only_on_product(self):
        first = np.kron(relation_state(1), relation_state(.5))
        second = np.kron(relation_state(.5), relation_state(1))
        np.testing.assert_array_equal(complete_product_moments(first), complete_product_moments(second))
        self.assertEqual(np.abs(np.linalg.eigvalsh(relation_state(1) - relation_state(.5))).sum() / 2, .25)

    def test_simultaneous_sign_change_is_partial_transpose_of_one_whole_laboratory(self):
        for q, c in ((1., 1.), (.5, -.25), (-.75, .25)):
            joint = np.kron(relation_state(q), relation_state(c))
            flipped = np.kron(relation_state(-q), relation_state(-c))
            np.testing.assert_array_equal(partial_transpose(joint, (2,) * 4, (1, 3)), flipped)
            np.testing.assert_array_equal(complete_product_moments(joint), complete_product_moments(flipped))
            exact_spectrum = sorted([(1 + a * q) * (1 + b * c) / 16
                                     for a in (-1, 1) for b in (-1, 1)] * 4)
            self.assertGreaterEqual(min(exact_spectrum), 0)
            np.testing.assert_allclose(np.linalg.eigvalsh(flipped), exact_spectrum, atol=2e-16)

    def test_copy_collective_local_effects_do_not_break_global_sign_ambiguity(self):
        rng = np.random.default_rng(188)
        for values in ((.5,), (.5, -.25), (1., .5, .5), (.25, .25, -.5, -.5)):
            count = len(values)
            state = many_pairs(values)
            opposite = many_pairs([-value for value in values])
            np.testing.assert_array_equal(partial_transpose(state, (2,) * (2 * count),
                                                          tuple(range(1, 2 * count, 2))), opposite)
            order = tuple(range(0, 2 * count, 2)) + tuple(range(1, 2 * count, 2))
            difference = reorder(state - opposite, (2,) * (2 * count), order)
            for _ in range(3):
                x, y = rng.normal(size=(2, 2**count))
                x /= np.linalg.norm(x)
                y /= np.linalg.norm(y)
                effect = np.kron(np.outer(x, x), np.outer(y, y))
                self.assertAlmostEqual(np.sum(difference * effect), 0, places=14)

    def test_two_reference_copies_can_reveal_magnitude_but_not_sign(self):
        observable = local_product_in_target_order(H, H)
        for c in (.25, .5, 1.):
            positive = many_pairs((c, c))
            negative = many_pairs((-c, -c))
            self.assertEqual(np.trace(positive @ observable), c * c)
            self.assertEqual(np.trace(negative @ observable), c * c)
        self.assertNotEqual(.25**2, .5**2)

    def test_known_signed_reference_breaks_ambiguity_without_altering_target(self):
        for c in (.25, .5, 1.):
            a = np.kron(relation_state(1), relation_state(c))
            b = np.kron(relation_state(-1), relation_state(c))
            observable = local_product_in_target_order(H, H)
            self.assertEqual(np.trace((a - b) @ observable), 2 * c)

    def test_reference_uncertainty_propagates_to_a_sharp_target_interval(self):
        bounds = positive_reference_interval((F(18, 100), F(22, 100)), (F(2, 5), F(3, 5)))
        self.assertEqual(bounds, (F(3, 10), F(11, 20)))
        self.assertEqual(bounds[0] * F(3, 5), F(18, 100))
        self.assertEqual(bounds[1] * F(2, 5), F(22, 100))
        self.assertEqual(positive_reference_interval((F(-1, 5), F(1, 5)), (F(1, 10), F(1, 2))),
                         (F(-1), F(1)))
        with self.assertRaises(ValueError):
            positive_reference_interval((F(0), F(1)), (F(0), F(1)))

    def test_iid_precision_bound_retains_reference_strength_and_sampling_assumptions(self):
        for c in (1., .5, .1):
            samples = iid_sample_bound(c, .1, .05)
            self.assertLessEqual(2 * math.exp(-samples * c * c * .1**2 / 2), .05)
        self.assertGreaterEqual(iid_sample_bound(.1, .1, .05), 99 * iid_sample_bound(1, .1, .05))
        with self.assertRaises(ValueError):
            iid_sample_bound(0, .1, .05)

    def test_unknown_opposite_references_hide_orthogonal_targets_exactly(self):
        first = np.kron(relation_state(1), relation_state(1))
        second = np.kron(relation_state(-1), relation_state(-1))
        np.testing.assert_array_equal(relation_state(1) @ relation_state(-1), np.zeros((4, 4)))
        np.testing.assert_array_equal(complete_product_moments(first), complete_product_moments(second))
        # For the displayed binary local measurements every decision rule has
        # equal-prior success 1/2; the all-protocol proof uses partial transpose.
        probabilities = [(1 + a * b) / 4 for a in (-1, 1) for b in (-1, 1)]
        for decision_mask in range(16):
            success = sum(p * (((decision_mask >> i) & 1) + (1 - ((decision_mask >> i) & 1))) / 2
                          for i, p in enumerate(probabilities))
            self.assertEqual(success, .5)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(BlindReferenceCalibrationTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = {
        "round": 188,
        "sign_convention": "rho_q=(I+q J tensor J)/4, J real skew; opposite to a YY coefficient",
        "single_copy_pair_identifiable_parameter": "q*c",
        "all_finite_copy_real_LOCC_sign_ambiguity": "(q,c) and (-q,-c)",
        "two_reference_copies_can_reveal_c_squared": True,
        "equal_prior_opposite_sign_best_success_without_signed_resource": "1/2",
        "reference_uncertainty_example_q_interval": ["3/10", "11/20"],
        "iid_samples_sufficient_for_known_c": "ceil(2*ln(2/delta)/(c^2*epsilon^2))",
        "finite_sample_bound_applies_automatically_to_reused_reference": False,
        "partial_transpose_is_claimed_as_a_physical_operation": False,
        "cognitive_necessity_of_unassisted_local_tomography_derived": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("blind_reference_calibration_results.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
