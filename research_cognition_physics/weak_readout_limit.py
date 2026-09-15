"""Public weak-record information can vanish while an admissible update disturbs."""

import argparse
import json
import math
import platform
import unittest
from pathlib import Path

import numpy as np

from outcome_updates import WINDOW_ATTENUATION, validate_unit_interval
from randomized_retention import (
    TANGENCY_CONTRAST,
    mixing_plan,
    nonselective_matrix,
    randomized_branch_matrix,
)
from retention_tradeoff import fresh_branch_matrix


def validate_count(count):
    if isinstance(count, bool) or not isinstance(count, int) or count < 1:
        raise ValueError("A positive integer read count is required.")


def bernoulli_kl(first, second):
    if not 0 < first < 1 or not 0 < second < 1:
        raise ValueError("This check uses strictly interior probabilities.")
    return first * math.log(first / second) + (1 - first) * math.log((1 - first) / (1 - second))


def transcript_bounds(count, contrast):
    validate_count(count)
    validate_unit_interval(contrast)
    radius = WINDOW_ATTENUATION * contrast
    divergence = count * radius * (math.log1p(radius) - math.log1p(-radius))
    return {"kl_upper_bound_nats": divergence,
            "tv_upper_bound": min(1.0, math.sqrt(divergence / 2))}


def weak_loss_rates():
    return np.array([(1 - WINDOW_ATTENUATION) / TANGENCY_CONTRAST, TANGENCY_CONTRAST])


def finite_unconditional_gains(count, total_strength=1.0):
    validate_count(count)
    if not math.isfinite(total_strength) or total_strength < 0:
        raise ValueError("The cumulative dimensionless strength must be finite and nonnegative.")
    contrast = total_strength / count
    if contrast > TANGENCY_CONTRAST:
        raise ValueError("This limit formula is for the weak linear branch only.")
    return np.exp(count * np.log1p(-weak_loss_rates() * contrast))


def limit_unconditional_gains(total_strength=1.0):
    if not math.isfinite(total_strength) or total_strength < 0:
        raise ValueError("The cumulative dimensionless strength must be finite and nonnegative.")
    return np.exp(-weak_loss_rates() * total_strength)


def public_transcripts(initial, count, contrast, adaptive=False):
    """Exact finite binary tree, bounded in size; no trajectory Monte Carlo."""
    validate_count(count)
    if count > 18:
        raise ValueError("Explicit transcript enumeration is limited to 18 reads.")
    states = np.asarray(initial, dtype=float)[None, :]
    previous = np.zeros(1, dtype=int)
    for depth in range(count):
        outputs = []
        for outcome in (0, 1):
            if adaptive and depth:
                updated = np.empty_like(states)
                for last in (0, 1):
                    mask = previous == last
                    setting = 0 if last == 0 else math.pi / 2
                    updated[mask] = states[mask] @ randomized_branch_matrix(setting, outcome, contrast).T
            else:
                updated = states @ randomized_branch_matrix(0, outcome, contrast).T
            outputs.append(updated)
        previous = np.concatenate((np.zeros(len(states), dtype=int), np.ones(len(states), dtype=int)))
        states = np.concatenate(outputs)
    return states[:, 0]


def first_active_distribution(initial, count, contrast):
    """Observable only in the explicitly richer interface that reveals modes."""
    validate_count(count)
    active_contrast, probability = mixing_plan(contrast)
    none = (1 - probability)**count
    first_probabilities = np.array([(fresh_branch_matrix(0, outcome, active_contrast) @ initial)[0]
                                   for outcome in (0, 1)])
    return np.concatenate(([none], (1 - none) * first_probabilities))


def enumerate_first_active(initial, count, contrast):
    """Independent four-event tree to check the mode-revelation control."""
    strength, probability = mixing_plan(contrast)
    states = [(None, np.asarray(initial, dtype=float))]
    for _ in range(count):
        next_states = []
        for recorded, state in states:
            for outcome in (0, 1):
                next_states.append((recorded, (1 - probability) * state / 2))
                next_states.append((outcome if recorded is None else recorded,
                                    probability * fresh_branch_matrix(0, outcome, strength) @ state))
        states = next_states
    distribution = np.zeros(3)
    for recorded, state in states:
        distribution[0 if recorded is None else recorded + 1] += state[0]
    return distribution


def limit_rows():
    rows = []
    alpha = WINDOW_ATTENUATION
    for count in (10, 100, 1000, 10000, 1000000):
        contrast = 1 / count
        parallel, tangent = finite_unconditional_gains(count)
        first = first_active_distribution(np.array([1, alpha, 0]), count, contrast)
        second = first_active_distribution(np.array([1, -alpha, 0]), count, contrast)
        rows.append({
            "read_count": count, "contrast_each": contrast,
            "parallel_gain": float(parallel), "tangent_gain": float(tangent),
            "fixed_probe_contrast": 0.5,
            "fixed_probe_tangent_probability_change": alpha * (1 - tangent) / 4,
            "revealed_first_active_record_tv": float(abs(first - second).sum() / 2),
            **transcript_bounds(count, contrast),
        })
    return rows


def finite_tree_rows():
    rows = []
    alpha = WINDOW_ATTENUATION
    for count in (8, 12, 16):
        contrast = 1 / count
        first = public_transcripts([1, alpha, 0], count, contrast)
        second = public_transcripts([1, -alpha, 0], count, contrast)
        rows.append({
            "read_count": count,
            "public_record_tv": float(abs(first - second).sum() / 2),
            "public_record_kl_nats": float(np.sum(first * np.log(first / second))),
            **transcript_bounds(count, contrast),
        })
    return rows


class WeakReadoutLimitTests(unittest.TestCase):
    def test_extreme_bernoulli_pair_attains_single_step_kl_bound(self):
        contrast = 0.1
        radius = contrast * WINDOW_ATTENUATION
        self.assertAlmostEqual(bernoulli_kl((1 + radius) / 2, (1 - radius) / 2),
                               transcript_bounds(1, contrast)["kl_upper_bound_nats"])

    def test_all_sampled_conditional_probabilities_meet_kl_bound(self):
        radius = 0.3 * WINDOW_ATTENUATION
        bound = transcript_bounds(1, 0.3)["kl_upper_bound_nats"]
        for first in np.linspace((1 - radius) / 2, (1 + radius) / 2, 31):
            for second in np.linspace((1 - radius) / 2, (1 + radius) / 2, 31):
                self.assertLessEqual(bernoulli_kl(float(first), float(second)), bound + 1e-14)

    def test_linear_weak_nonselective_generator_matches_actual_branch_sum(self):
        for contrast in (0, 0.001, 0.05, 0.1):
            expected = np.diag(np.concatenate(([1], 1 - contrast * weak_loss_rates())))
            np.testing.assert_allclose(nonselective_matrix(contrast), expected, atol=1e-14)

    def test_finite_products_match_stable_formula(self):
        for count in (10, 100, 1000):
            actual = np.linalg.matrix_power(nonselective_matrix(1 / count), count)
            np.testing.assert_allclose(np.diag(actual)[1:], finite_unconditional_gains(count), atol=3e-13)

    def test_product_converges_to_nonidentity_limit(self):
        limit = limit_unconditional_gains()
        self.assertTrue(np.all(limit < 1))
        self.assertTrue(np.all(limit > 0))
        errors = [float(np.max(abs(finite_unconditional_gains(count) - limit))) for count in (10, 100, 1000)]
        self.assertTrue(errors[0] > errors[1] > errors[2])

    def test_limit_has_semigroup_composition(self):
        np.testing.assert_allclose(limit_unconditional_gains(0.4) * limit_unconditional_gains(0.7),
                                   limit_unconditional_gains(1.1), atol=1e-14)

    def test_zero_total_strength_preserves_state_and_produces_no_record_information(self):
        np.testing.assert_array_equal(finite_unconditional_gains(10, 0), [1, 1])
        self.assertEqual(transcript_bounds(10, 0)["tv_upper_bound"], 0)

    def test_record_bound_vanishes_as_inverse_square_root_count(self):
        for count in (1000, 10000, 1000000):
            scaled = transcript_bounds(count, 1 / count)["tv_upper_bound"] * math.sqrt(count)
            self.assertAlmostEqual(scaled, WINDOW_ATTENUATION, places=6)

    def test_explicit_public_trees_normalize_and_satisfy_information_bound(self):
        for adaptive in (False, True):
            for count in (4, 8, 12):
                contrast = 0.05
                first = public_transcripts([1, 0.4, 0.5], count, contrast, adaptive)
                second = public_transcripts([1, -0.3, -0.6], count, contrast, adaptive)
                self.assertAlmostEqual(float(first.sum()), 1)
                self.assertAlmostEqual(float(second.sum()), 1)
                divergence = float(np.sum(first * np.log(first / second)))
                tv = float(abs(first - second).sum() / 2)
                bound = transcript_bounds(count, contrast)
                self.assertLessEqual(divergence, bound["kl_upper_bound_nats"] + 1e-14)
                self.assertLessEqual(tv, math.sqrt(max(divergence, 0) / 2) + 1e-14)

    def test_probe_detects_finite_state_change_when_weak_record_bound_is_small(self):
        row = limit_rows()[-1]
        self.assertLess(row["tv_upper_bound"], 0.001)
        self.assertGreater(row["fixed_probe_tangent_probability_change"], 0.03)

    def test_revealed_first_active_record_matches_independent_event_tree(self):
        for initial in (np.array([1, WINDOW_ATTENUATION, 0]), np.array([1, -WINDOW_ATTENUATION, 0])):
            np.testing.assert_allclose(first_active_distribution(initial, 4, 0.02),
                                       enumerate_first_active(initial, 4, 0.02), atol=1e-14)

    def test_revealing_hidden_mode_changes_information_conclusion(self):
        row = limit_rows()[-1]
        self.assertGreater(row["revealed_first_active_record_tv"], 0.14)
        self.assertLess(row["tv_upper_bound"], 0.001)

    def test_first_active_control_has_expected_positive_limit(self):
        expected = (1 - math.exp(-1 / TANGENCY_CONTRAST)) * TANGENCY_CONTRAST * WINDOW_ATTENUATION
        self.assertAlmostEqual(limit_rows()[-1]["revealed_first_active_record_tv"], expected, places=7)

    def test_invalid_count_and_out_of_weak_regime_are_rejected(self):
        for count in (0, -1, 0.5, True):
            with self.assertRaises(ValueError):
                transcript_bounds(count, 0.1)
        with self.assertRaises(ValueError):
            finite_unconditional_gains(2)
        with self.assertRaises(ValueError):
            public_transcripts([1, 0, 0], 19, 0.01)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(WeakReadoutLimitTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    results = {
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "automated_checks": {"run": tests.testsRun, "failures": len(tests.failures), "errors": len(tests.errors)},
        "total_dimensionless_strength": 1.0,
        "longitudinal_parameter": 1.0,
        "loss_rates_parallel_tangent": weak_loss_rates().tolist(),
        "limiting_gains_parallel_tangent": limit_unconditional_gains().tolist(),
        "limit_rows": limit_rows(),
        "exact_finite_transcript_trees": finite_tree_rows(),
        "scope": "The uniform KL/TV bound concerns only the n public weak binary results, with hidden modes unavailable, and also allows adaptive settings. State disturbance uses fixed settings and a separate fixed-strength final probe. Revealing modes is an explicitly different interface. Gamma is not physical time.",
    }
    if args.write_results:
        Path(__file__).with_name("weak_readout_limit_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
