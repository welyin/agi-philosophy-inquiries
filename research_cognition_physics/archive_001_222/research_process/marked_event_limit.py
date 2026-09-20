"""Mode-resolved classical jump process and exact finite mark-record trees."""

import argparse
import json
import math
import platform
import unittest
from pathlib import Path

import numpy as np

from outcome_updates import WINDOW_ATTENUATION as ALPHA
from randomized_retention import TANGENCY_CONTRAST as KAPPA, randomized_branch_matrix
from retention_tradeoff import fresh_branch_matrix
from weak_readout_limit import limit_unconditional_gains


def active_matrix(outcome):
    return fresh_branch_matrix(0, outcome, KAPPA)


def conditional_jump(initial, outcome):
    """Normalized prediction after a known active event with its binary mark."""
    if outcome not in (0, 1):
        raise ValueError("A binary mark is required.")
    state = np.asarray(initial, dtype=float)
    if state.shape != (3,) or not np.all(np.isfinite(state)):
        raise ValueError("A finite three-component summary is required.")
    if not np.isclose(state[0], 1) or np.linalg.norm(state[1:]) > ALPHA + 1e-13:
        raise ValueError("A normalized allowed summary is required.")
    sign = 2 * outcome - 1
    denominator = 1 + sign * KAPPA * state[1]
    return np.array([1, ALPHA * (state[1] + sign * KAPPA) / denominator,
                     ALPHA**2 * state[2] / denominator])


def mark_statistics(max_events=18):
    """Enumerate every ordered active-mark word for two opposite preparations."""
    if isinstance(max_events, bool) or not isinstance(max_events, int) or not 0 <= max_events <= 20:
        raise ValueError("Explicit enumeration supports 0 through 20 events.")
    states = np.array([[[1, ALPHA, 0], [1, -ALPHA, 0]]], dtype=float)
    rows = []
    for count in range(max_events + 1):
        probabilities = states[:, :, 0]
        rows.append({"active_events": count,
                     "record_tv": float(np.abs(probabilities[:, 0] - probabilities[:, 1]).sum() / 2),
                     "normalization_error": float(np.max(abs(probabilities.sum(axis=0) - 1)))})
        if count < max_events:
            states = np.concatenate([states @ active_matrix(outcome).T for outcome in (0, 1)])
    return rows


def poisson_weights(total_strength, max_events):
    if not math.isfinite(total_strength) or total_strength < 0:
        raise ValueError("Cumulative dimensionless strength must be finite and nonnegative.")
    mean = total_strength / KAPPA
    if max_events + 2 <= mean:
        raise ValueError("Use a cutoff above the mean for the geometric tail bound.")
    weights = [math.exp(-mean)]
    for count in range(1, max_events + 1):
        weights.append(weights[-1] * mean / count)
    first_omitted = weights[-1] * mean / (max_events + 1)
    tail_upper = min(1.0, first_omitted / (1 - mean / (max_events + 2)))
    return np.array(weights), tail_upper


def continuous_record_rows(statistics):
    count = len(statistics) - 1
    distances = np.array([row["record_tv"] for row in statistics])
    rows = []
    for total in (0, 0.1, 0.3, 1.0):
        weights, tail = poisson_weights(total, count)
        partial = float(weights @ distances)
        rows.append({"total_dimensionless_strength": total,
                     "expected_active_events": total / KAPPA,
                     "first_active_record_tv": -math.expm1(-total / KAPPA) * KAPPA * ALPHA,
                     "all_active_marks_tv_lower": partial,
                     "all_active_marks_tv_upper": min(1.0, partial + tail),
                     "analytic_truncation_error_upper": tail})
    return rows


def ordered_path_density(initial, marks, total_strength):
    """Density on ordered event times; valid only for times inside [0, total]."""
    state = np.asarray(initial, dtype=float)
    for mark in marks:
        state = active_matrix(mark) @ state
    return math.exp(-total_strength / KAPPA) * KAPPA**(-len(marks)) * state[0]


class MarkedEventLimitTests(unittest.TestCase):
    def test_conditional_formula_matches_branch_and_preserves_preparations(self):
        for angle in np.linspace(-math.pi, math.pi, 97):
            initial = np.array([1, ALPHA * math.cos(angle), ALPHA * math.sin(angle)])
            for mark in (0, 1):
                raw = active_matrix(mark) @ initial
                after = conditional_jump(initial, mark)
                np.testing.assert_allclose(after, raw / raw[0], atol=1e-14)
                self.assertLessEqual(np.linalg.norm(after[1:]), ALPHA + 1e-14)

    def test_hiding_modes_recovers_previous_instrument(self):
        for eta in (0, 0.01, 0.1):
            for mark in (0, 1):
                revealed_sum = (1 - eta / KAPPA) * np.eye(3) / 2 + eta / KAPPA * active_matrix(mark)
                np.testing.assert_allclose(revealed_sum, randomized_branch_matrix(0, mark, eta), atol=1e-14)

    def test_event_hazard_sum_is_state_independent(self):
        for u in np.linspace(-ALPHA, ALPHA, 31):
            rates = [(1 + sign * KAPPA * u) / (2 * KAPPA) for sign in (-1, 1)]
            self.assertGreater(min(rates), 0)
            self.assertAlmostEqual(sum(rates), 1 / KAPPA)

    def test_poisson_mixture_matches_previous_mean_state(self):
        weights, tail = poisson_weights(1, 30)
        active_sum = active_matrix(0) + active_matrix(1)
        mean_map = sum(weight * np.linalg.matrix_power(active_sum, k) for k, weight in enumerate(weights))
        target = np.diag(np.concatenate(([1], limit_unconditional_gains())))
        self.assertLessEqual(float(np.max(abs(mean_map - target))), tail + 2e-15)

    def test_exact_mark_records_normalize_and_more_records_cannot_hurt(self):
        rows = mark_statistics(12)
        self.assertLess(max(row["normalization_error"] for row in rows), 1e-13)
        self.assertGreaterEqual(min(np.diff([row["record_tv"] for row in rows])), -1e-14)
        self.assertAlmostEqual(rows[1]["record_tv"], KAPPA * ALPHA)
        self.assertGreater(rows[-1]["record_tv"], rows[1]["record_tv"])

    def test_ordered_time_simplex_integration_gives_poisson_count(self):
        from itertools import product
        total = 0.3
        weights, _ = poisson_weights(total, 8)
        for count in range(6):
            density = sum(ordered_path_density([1, 0.2, 0.3], marks, total)
                          for marks in product((0, 1), repeat=count))
            self.assertAlmostEqual(density * total**count / math.factorial(count), weights[count])

    def test_analytic_tail_bound_controls_omitted_mass(self):
        short, tail = poisson_weights(1, 18)
        long, _ = poisson_weights(1, 50)
        self.assertLessEqual(float(long[19:].sum()), tail)
        self.assertAlmostEqual(float(short.sum() + long[19:].sum()), 1)

    def test_full_record_improves_first_event_and_zero_strength_has_no_information(self):
        rows = continuous_record_rows(mark_statistics(12))
        self.assertEqual(rows[0]["all_active_marks_tv_upper"], 0)
        self.assertGreater(rows[-1]["all_active_marks_tv_lower"], rows[-1]["first_active_record_tv"])

    def test_binomial_no_event_probability_converges_to_exponential(self):
        mean = 1 / KAPPA
        expected = math.exp(-mean)
        errors = [abs(math.exp(n * math.log1p(-mean / n)) - expected) for n in (100, 1000, 10000)]
        self.assertTrue(errors[0] > errors[1] > errors[2])

    def test_illegal_states_and_excessive_enumeration_are_rejected(self):
        with self.assertRaises(ValueError):
            conditional_jump([1, 1, 0], 0)
        with self.assertRaises(ValueError):
            mark_statistics(30)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(MarkedEventLimitTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    statistics = mark_statistics()
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": tests.testsRun, "failures": len(tests.failures), "errors": len(tests.errors)},
               "event_rate_per_dimensionless_strength": 1 / KAPPA,
               "active_mark_statistics": statistics,
               "continuous_record_distances": continuous_record_rows(statistics),
               "scope": "Fixed setting, longitudinal parameter 1, public idle/active flags and binary marks. Poisson truncation has an analytic bound; floating-point sums are not certified interval arithmetic. No latent source copy or physical time is assumed."}
    if args.write_results:
        Path(__file__).with_name("marked_event_limit_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
