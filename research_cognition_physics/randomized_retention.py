"""Optimal convex randomization of idle preservation and fresh-window reads."""

import argparse
import json
import math
import platform
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from bounded_responses import harmonic_features, probability_effect, translation_matrix
from outcome_updates import WINDOW_ATTENUATION, partial_branch_matrix, validate_unit_interval
from retention_tradeoff import fresh_branch_matrix


TANGENCY_CONTRAST = math.sqrt(1 - WINDOW_ATTENUATION**2)


def envelope(contrast):
    validate_unit_interval(contrast)
    if contrast <= TANGENCY_CONTRAST:
        return 1 - TANGENCY_CONTRAST * contrast
    return WINDOW_ATTENUATION * math.sqrt(1 - contrast**2)


def mixing_plan(contrast):
    validate_unit_interval(contrast)
    active_contrast = max(contrast, TANGENCY_CONTRAST)
    return active_contrast, contrast / active_contrast


def randomized_branch_matrix(setting, outcome, contrast, longitudinal=1.0):
    active_contrast, active_weight = mixing_plan(contrast)
    return ((1 - active_weight) * np.eye(3) / 2
            + active_weight * fresh_branch_matrix(setting, outcome, active_contrast, longitudinal))


def nonselective_matrix(contrast, longitudinal=1.0):
    return sum(randomized_branch_matrix(0, outcome, contrast, longitudinal) for outcome in (0, 1))


def protocol_summary(initial, settings, outcomes, contrast, longitudinal=1.0):
    if len(settings) != len(outcomes):
        raise ValueError("Each result must have a setting.")
    state = np.asarray(initial, dtype=float)
    for setting, outcome in zip(settings, outcomes):
        state = randomized_branch_matrix(setting, outcome, contrast, longitudinal) @ state
    return state


def explicit_mode_summary(initial, settings, outcomes, contrast, longitudinal=1.0):
    """Sum all hidden idle/active mode paths independently of the mixed matrix."""
    active_contrast, probability = mixing_plan(contrast)
    total = np.zeros(3)
    for modes in product((0, 1), repeat=len(outcomes)):
        state = np.asarray(initial, dtype=float)
        for mode, setting, outcome in zip(modes, settings, outcomes):
            if mode:
                state = probability * fresh_branch_matrix(setting, outcome, active_contrast, longitudinal) @ state
            else:
                state = (1 - probability) * state / 2
        total += state
    return total


def comparison_rows():
    rows = []
    alpha = WINDOW_ATTENUATION
    for contrast in (0, 0.005, 0.02, 0.05, 0.1, TANGENCY_CONTRAST, 0.25, 0.5, 1):
        active_contrast, probability = mixing_plan(contrast)
        rows.append({
            "contrast": contrast,
            "active_contrast": active_contrast,
            "active_mode_probability": probability,
            "old_tangent_gain": 1 - contrast,
            "fresh_tangent_gain": alpha * math.sqrt(1 - contrast**2),
            "randomized_tangent_gain": envelope(contrast),
            "general_kernel_upper_bound": math.sqrt(1 - contrast**2),
            "next_read_gap": contrast * alpha * envelope(contrast),
            "unconditional_parallel_gain": float(nonselective_matrix(contrast)[1, 1]),
        })
    return rows


class RandomizedRetentionTests(unittest.TestCase):
    def test_fixed_readout_response_after_randomization(self):
        for contrast, outcome in product((0, 0.02, 0.5, 1), (0, 1)):
            effect = probability_effect(0.7, outcome)
            effect[1:] *= contrast
            np.testing.assert_allclose(randomized_branch_matrix(0.7, outcome, contrast)[0], effect, atol=1e-14)

    def test_zero_readout_is_exactly_idle_with_a_fair_record(self):
        for outcome in (0, 1):
            np.testing.assert_array_equal(randomized_branch_matrix(0, outcome, 0), np.eye(3) / 2)

    def test_active_mode_is_a_valid_mixture(self):
        for contrast in np.linspace(0, 1, 301):
            strength, probability = mixing_plan(float(contrast))
            self.assertTrue(0 < strength <= 1)
            self.assertTrue(0 <= probability <= 1)
            self.assertAlmostEqual(strength * probability, contrast)

    def test_envelope_matches_branch_gain(self):
        for contrast in np.linspace(0, 1, 201):
            for outcome in (0, 1):
                self.assertAlmostEqual(2 * randomized_branch_matrix(0, outcome, float(contrast))[2, 2], envelope(float(contrast)))

    def test_tangency_matches_value_and_slope(self):
        kappa, alpha = TANGENCY_CONTRAST, WINDOW_ATTENUATION
        self.assertAlmostEqual(envelope(kappa), alpha**2)
        self.assertAlmostEqual(-alpha * kappa / math.sqrt(1 - kappa**2), -kappa)

    def test_envelope_majorizes_every_fresh_primitive(self):
        for contrast in np.linspace(0, 1, 2001):
            self.assertGreaterEqual(envelope(float(contrast)) + 2e-15,
                                    WINDOW_ATTENUATION * math.sqrt(1 - contrast**2))

    def test_concavity_bounds_random_mixtures(self):
        rng = np.random.default_rng(7319)
        for _ in range(100):
            strengths = np.concatenate(([0.0], rng.uniform(0, 1, 7)))
            weights = rng.dirichlet(np.ones(8))
            gains = WINDOW_ATTENUATION * np.sqrt(1 - strengths**2)
            gains[0] = 1  # Idle, not the contrast-zero fresh primitive.
            self.assertLessEqual(float(weights @ gains), envelope(float(weights @ strengths)) + 1e-14)

    def test_exact_optimal_two_mode_mixture_attains_bound(self):
        for contrast in (0.005, 0.05, 0.1, 0.3, 0.5):
            strength, probability = mixing_plan(contrast)
            attained = 1 - probability + probability * WINDOW_ATTENUATION * math.sqrt(1 - strength**2)
            self.assertAlmostEqual(attained, envelope(contrast))

    def test_mixed_summary_matches_enumerated_hidden_modes(self):
        initial = np.array([1, 0.2, -0.4])
        settings = (0.1, -0.3, 0.9)
        for contrast in (0, 0.02, 0.1, 0.5):
            for outcomes in product((0, 1), repeat=3):
                np.testing.assert_allclose(protocol_summary(initial, settings, outcomes, contrast),
                                           explicit_mode_summary(initial, settings, outcomes, contrast), atol=1e-14)

    def test_all_branch_probabilities_normalize(self):
        initial = [1, 0.2, -0.4]
        for contrast in (0, 0.02, 0.5, 1):
            total = sum(protocol_summary(initial, (0.1, 0.2, -0.3), outcomes, contrast)[0]
                        for outcomes in product((0, 1), repeat=3))
            self.assertAlmostEqual(float(total), 1)

    def test_full_summary_cone_preserved(self):
        states = harmonic_features(np.linspace(-math.pi, math.pi, 301), (1.0,))
        states[:, 1:] *= WINDOW_ATTENUATION
        for contrast, outcome in product((0, 0.01, 0.1, 0.5, 1), (0, 1)):
            output = states @ randomized_branch_matrix(0.3, outcome, contrast).T
            self.assertLessEqual(float(np.max(np.linalg.norm(output[:, 1:], axis=1)
                                              - WINDOW_ATTENUATION * output[:, 0])), 1e-14)

    def test_original_sharp_reset_is_in_the_primitive_class(self):
        for outcome in (0, 1):
            np.testing.assert_allclose(fresh_branch_matrix(0, outcome, 1),
                                       partial_branch_matrix(0, outcome, 1), atol=1e-14)

    def test_low_contrast_improves_both_old_and_direct_fresh_rules(self):
        contrast = 0.005
        self.assertGreater(envelope(contrast), 1 - contrast)
        self.assertGreater(envelope(contrast), WINDOW_ATTENUATION * math.sqrt(1 - contrast**2))
        self.assertLess(envelope(contrast), math.sqrt(1 - contrast**2))

    def test_family_remains_nonunique_at_fixed_optimum(self):
        contrast = 0.05
        first = randomized_branch_matrix(0, 1, contrast, 0)
        last = randomized_branch_matrix(0, 1, contrast, 1)
        self.assertEqual(first[2, 2], last[2, 2])
        self.assertNotEqual(first[1, 1], last[1, 1])
        for matrix in (first, last):
            self.assertGreaterEqual(np.linalg.matrix_rank(matrix), 2)

    def test_zero_limit_has_a_finite_linear_loss(self):
        for contrast in (1e-3, 1e-4, 1e-5):
            self.assertAlmostEqual((1 - envelope(contrast)) / contrast, TANGENCY_CONTRAST, places=10)

    def test_rotational_covariance(self):
        rotation = translation_matrix(0.7, (1.0,))
        np.testing.assert_allclose(randomized_branch_matrix(0.9, 1, 0.03) @ rotation,
                                   rotation @ randomized_branch_matrix(0.2, 1, 0.03), atol=1e-14)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RandomizedRetentionTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    results = {
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "automated_checks": {"run": tests.testsRun, "failures": len(tests.failures), "errors": len(tests.errors)},
        "tangency_contrast": TANGENCY_CONTRAST,
        "comparison": comparison_rows(),
        "scope": "Exact optimum only among input-independent convex mixtures of the idle fair-record instrument and the declared fresh-window family at arbitrary nonnegative strengths. The hidden mode is unavailable in future reads. This does not close the general-kernel gap.",
    }
    if args.write_results:
        Path(__file__).with_name("randomized_retention_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
