"""Four preparations witness the added preparation-noncontextuality assumption.

The original classical model is allowed to retain inaccessible preparation
information. Violating this inequality therefore does not rule it out.
"""

import argparse
import json
import math
import platform
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from bounded_responses import translation_matrix
from finite_hidden_readout import encode_summary, read_kernel
from outcome_updates import WINDOW_ATTENUATION as ALPHA
from retention_tradeoff import fresh_branch_matrix


SIGNS = np.array(list(product((-1, 1), repeat=2)))


def square_preparations(radius=ALPHA):
    return np.column_stack((np.ones(4), radius * SIGNS / math.sqrt(2)))


def success_table(contrast=1.0, radius=ALPHA):
    preparations = square_preparations(radius)
    return np.array([[fresh_branch_matrix(axis * math.pi / 2, int(sign > 0), contrast)[0] @ preparation
                      for axis, sign in enumerate(signs)]
                     for preparation, signs in zip(preparations, SIGNS)])


def parity_defect(encodings):
    encodings = np.asarray(encodings)
    parity = np.prod(SIGNS, axis=1)
    even, odd = encodings[parity == 1].mean(axis=0), encodings[parity == -1].mean(axis=0)
    return float(np.abs(even - odd).sum() / 2)


def optimal_hidden_success(encodings):
    """Optimize both bounded response functions for a fixed hidden encoding."""
    signed = SIGNS.T @ np.asarray(encodings) / 4
    return float(0.5 + np.abs(signed).sum() / 4)


def saturating_encoding(defect):
    if not 0 <= defect <= 1:
        raise ValueError("Use a defect between zero and one.")
    probabilities = np.zeros((4, 6))
    for index, signs in enumerate(SIGNS):
        probabilities[index, int(signs[0] > 0)] = 1 - defect
        probabilities[index, 2 + index] = defect
    return probabilities


def adaptive_probabilities(initial):
    values = []
    for marks in product((0, 1), repeat=4):
        state = np.asarray(initial, dtype=float)
        for step, mark in enumerate(marks):
            setting = 0.4 + 0.7 * sum(marks[:step])
            state = fresh_branch_matrix(setting, mark, 0.8) @ translation_matrix(0.31, (1.0,)) @ state
        values.append(state[0])
    return np.asarray(values)


def diagnostics():
    encodings = np.array([encode_summary(state, 22) for state in square_preparations()])
    score = float(success_table().mean())
    return {"preparation_radius": ALPHA, "preparations": square_preparations().tolist(),
            "correct_guess_probabilities": success_table().tolist(), "success_probability": score,
            "preparation_noncontextual_bound": 0.75, "violation": score - 0.75,
            "contrast_threshold_for_this_witness": 1 / (math.sqrt(2) * ALPHA),
            "minimum_hidden_parity_tv_from_score": max(0.0, 4 * score - 3),
            "round20_grid_hidden_parity_tv": parity_defect(encodings),
            "original_window_hidden_parity_tv": 1.0,
            "contrasts": [{"contrast": eta, "success": float(success_table(eta).mean())}
                          for eta in (0, 0.5, 0.7, 0.8, 1)]}


class PreparationContextualityTests(unittest.TestCase):
    def test_all_four_preparations_are_allowed_boundary_states(self):
        np.testing.assert_allclose(np.linalg.norm(square_preparations()[:, 1:], axis=1), ALPHA)

    def test_opposite_parity_mixtures_have_identical_predictive_summaries(self):
        states = square_preparations()
        parity = np.prod(SIGNS, axis=1)
        np.testing.assert_allclose(states[parity == 1].mean(axis=0), [1, 0, 0], atol=1e-15)
        np.testing.assert_allclose(states[parity == -1].mean(axis=0), [1, 0, 0], atol=1e-15)

    def test_parity_mixtures_remain_indistinguishable_under_adaptive_protocol(self):
        distributions = np.array([adaptive_probabilities(state) for state in square_preparations()])
        parity = np.prod(SIGNS, axis=1)
        np.testing.assert_allclose(distributions[parity == 1].mean(axis=0),
                                   distributions[parity == -1].mean(axis=0), atol=1e-14)

    def test_success_formula_uses_actual_existing_read_instrument(self):
        for eta in (0, 0.1, 0.5, 1):
            np.testing.assert_allclose(success_table(eta), (1 + eta * ALPHA / math.sqrt(2)) / 2)

    def test_threshold_separates_violation_and_nonviolation(self):
        threshold = 1 / (math.sqrt(2) * ALPHA)
        self.assertAlmostEqual(success_table(threshold).mean(), 0.75)
        self.assertGreater(success_table(threshold + 0.01).mean(), 0.75)
        self.assertLess(success_table(threshold - 0.01).mean(), 0.75)

    def test_robust_bound_for_all_deterministic_four_label_encodings(self):
        for labels in product(range(4), repeat=4):
            encodings = np.eye(4)[list(labels)]
            self.assertLessEqual(optimal_hidden_success(encodings), 0.75 + parity_defect(encodings) / 4 + 1e-14)

    def test_robust_bound_for_stochastic_encodings(self):
        rng = np.random.default_rng(2101)
        for _ in range(100):
            encodings = rng.dirichlet(np.ones(7), size=4)
            self.assertLessEqual(optimal_hidden_success(encodings), 0.75 + parity_defect(encodings) / 4 + 1e-14)

    def test_robust_bound_is_attained_at_every_interpolated_defect(self):
        for defect in np.linspace(0, 1, 11):
            encodings = saturating_encoding(defect)
            np.testing.assert_allclose(encodings.sum(axis=1), 1)
            self.assertAlmostEqual(parity_defect(encodings), defect)
            self.assertAlmostEqual(optimal_hidden_success(encodings), 0.75 + defect / 4)

    def test_existing_grid_model_reproduces_violation_with_nonzero_hidden_defect(self):
        states = square_preparations()
        encodings = np.array([encode_summary(state, 22) for state in states])
        probabilities = []
        for encoding, signs in zip(encodings, SIGNS):
            probabilities.append([np.sum(read_kernel(22, axis * math.pi / 2, int(sign > 0), 1) @ encoding)
                                  for axis, sign in enumerate(signs)])
        np.testing.assert_allclose(probabilities, success_table(), atol=1e-14)
        self.assertGreaterEqual(parity_defect(encodings) + 1e-14, 4 * success_table().mean() - 3)

    def test_original_four_windows_have_disjoint_parity_support(self):
        # The smallest separation is pi/2, while each support has length 1/2.
        self.assertLess(2 * 0.25, math.pi / 2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(PreparationContextualityTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": checks.testsRun, "failures": len(checks.failures), "errors": len(checks.errors)},
               **diagnostics(),
               "scope": "Preparation noncontextuality and convex mixing imply the 3/4 bound, without finite hidden-state or deterministic-response assumptions. The original restricted-access classical model violates preparation noncontextuality, not classical probability. Hidden parity TV is not observable parity leakage."}
    if args.write_results:
        Path(__file__).with_name("preparation_contextuality_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
