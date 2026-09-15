"""Changing preparation resolution permits informative weak classical records.

Every count n defines a separate model, with its own fixed preparation width
within that experiment. This does not solve the original fixed-width problem.
"""

import argparse
import json
import math
import platform
import unittest
from pathlib import Path

import numpy as np


def parameters(count, noise=1.0, signal=1.0):
    if isinstance(count, bool) or not isinstance(count, int) or count < 1:
        raise ValueError("A positive integer count is required.")
    if not math.isfinite(noise) or noise < 0 or not math.isfinite(signal) or signal <= 0:
        raise ValueError("Nonnegative noise and positive signal are required.")
    h, eta = math.sqrt(noise / count), signal / math.sqrt(count)
    if h >= math.pi / 2 or eta > 1:
        raise ValueError("Choose enough steps for h < pi/2 and contrast <= 1.")
    alpha = math.sin(h) / h if h else 1.0
    return h, alpha, eta


def branch_matrix(outcome, count, noise=1.0, signal=1.0):
    _, alpha, eta = parameters(count, noise, signal)
    sign = 2 * outcome - 1
    return np.array([[1, sign * eta, 0], [alpha * sign * eta, alpha, 0],
                     [0, 0, alpha * math.sqrt(max(0, 1 - eta**2))]]) / 2


def finite_mean_gains(count, noise=1.0, signal=1.0):
    _, alpha, eta = parameters(count, noise, signal)
    parallel = math.exp(count * math.log(alpha))
    tangent = parallel * math.exp(count / 2 * math.log1p(-eta**2)) if eta < 1 else 0.0
    return np.array([parallel, tangent])


def limit_mean_gains(noise=1.0, signal=1.0):
    return np.exp([-noise / 6, -noise / 6 - signal**2 / 2])


def geometric_average(alpha, count):
    if alpha == 1:
        return 1.0
    return -math.expm1(count * math.log(alpha)) / (count * (1 - alpha))


def record_mean_and_tv_floor(count, noise=1.0, signal=1.0):
    """Y=sum of binary signs / sqrt(n), for initial mu_0 vs mu_pi."""
    _, alpha, _ = parameters(count, noise, signal)
    mean_plus = signal * alpha * geometric_average(alpha, count)
    second_moment_upper = 1 + signal**2 * (count - 1) / count
    return mean_plus, mean_plus**2 / second_moment_upper


def limit_record_tv_floor(noise=1.0, signal=1.0):
    rate = noise / 6
    average = -math.expm1(-rate) / rate if rate else 1.0
    return (signal * average)**2 / (1 + signal**2)


def enumerate_records(count, noise=1.0, signal=1.0):
    if count > 16:
        raise ValueError("Full binary trees are limited to 16 steps.")
    _, alpha, eta = parameters(count, noise, signal)
    states = np.array([[[1, alpha, 0], [1, -alpha, 0]]], dtype=float)
    signs = np.zeros(1)
    for _ in range(count):
        states = np.concatenate([states @ branch_matrix(mark, count, noise, signal).T for mark in (0, 1)])
        signs = np.concatenate((signs - 1, signs + 1))
    probabilities = states[:, :, 0]
    observations = signs / math.sqrt(count)
    return probabilities, observations


def local_conditional_moments(initial, count, noise=1.0, signal=1.0):
    """Conditional drift and covariance per unit step parameter, directly summed."""
    initial = np.asarray(initial, dtype=float)
    drift, second = np.zeros(2), np.zeros((2, 2))
    for mark in (0, 1):
        raw = branch_matrix(mark, count, noise, signal) @ initial
        delta = raw[1:] / raw[0] - initial[1:]
        drift += raw[0] * delta
        second += raw[0] * np.outer(delta, delta)
    return count * drift, count * second - count * np.outer(drift, drift)


class ShrinkingPreparationsTests(unittest.TestCase):
    def test_actual_window_kernel_has_the_claimed_first_moments(self):
        # Integrate the fresh output window directly at each latent input.
        nodes, weights = np.polynomial.legendre.leggauss(24)
        for count in (4, 10, 100):
            h, alpha, eta = parameters(count)
            for x in np.linspace(-math.pi, math.pi, 33):
                for mark in (0, 1):
                    sign = 2 * mark - 1
                    center = math.atan2(math.sqrt(1 - eta**2) * math.sin(x), math.cos(x) + sign * eta)
                    mass = (1 + sign * eta * math.cos(x)) / 2
                    emitted = center + h * nodes
                    actual = mass * np.array([1, weights @ np.cos(emitted) / 2, weights @ np.sin(emitted) / 2])
                    np.testing.assert_allclose(actual, branch_matrix(mark, count) @ [1, math.cos(x), math.sin(x)], atol=1e-14)

    def test_new_preparation_cone_is_preserved_in_each_separate_model(self):
        for count in (4, 100):
            _, alpha, _ = parameters(count)
            for angle in np.linspace(0, 2 * math.pi, 65):
                initial = [1, alpha * math.cos(angle), alpha * math.sin(angle)]
                for mark in (0, 1):
                    after = branch_matrix(mark, count) @ initial
                    self.assertLessEqual(np.linalg.norm(after[1:]), alpha * after[0] + 1e-14)

    def test_exact_matrix_products_match_finite_gain_formula(self):
        for count in (4, 100, 1000):
            matrix = branch_matrix(0, count) + branch_matrix(1, count)
            actual = np.diag(np.linalg.matrix_power(matrix, count))[1:]
            np.testing.assert_allclose(actual, finite_mean_gains(count), atol=1e-13)

    def test_nontrivial_mean_limit_and_convergence(self):
        errors = [np.max(abs(finite_mean_gains(n) - limit_mean_gains())) for n in (10, 100, 1000)]
        self.assertTrue(errors[0] > errors[1] > errors[2])
        self.assertGreater(float(limit_mean_gains().min()), 0.5)

    def test_record_mean_and_second_moment_match_independent_tree(self):
        for count in (4, 8, 12):
            probabilities, observation = enumerate_records(count)
            np.testing.assert_allclose(probabilities.sum(axis=0), [1, 1], atol=1e-14)
            mean, floor = record_mean_and_tv_floor(count)
            np.testing.assert_allclose(observation @ probabilities, [mean, -mean], atol=1e-14)
            _, alpha, eta = parameters(count)
            second = 1 + 2 * eta**2 / count * sum((count - d) * alpha**d for d in range(1, count))
            np.testing.assert_allclose(observation**2 @ probabilities, [second, second], atol=1e-14)
            distance = np.abs(probabilities[:, 0] - probabilities[:, 1]).sum() / 2
            self.assertGreaterEqual(distance, floor - 1e-14)

    def test_record_lower_bound_stays_positive(self):
        self.assertGreater(limit_record_tv_floor(), 0.42)
        self.assertAlmostEqual(record_mean_and_tv_floor(100000)[1], limit_record_tv_floor(), places=5)

    def test_conditional_local_moments_match_candidate_diffusion_coefficients(self):
        initial = np.array([1, 0.3, 0.4])
        u, v = initial[1:]
        drift_target = np.array([-u / 6, -(1 / 6 + 0.5) * v])
        diffusion = np.array([1 - u**2, -u * v])
        errors = []
        for count in (100, 1000, 10000):
            drift, covariance = local_conditional_moments(initial, count)
            errors.append(max(np.max(abs(drift - drift_target)), np.max(abs(covariance - np.outer(diffusion, diffusion)))))
        self.assertTrue(errors[0] > errors[1] > errors[2])
        self.assertLess(errors[-1], 1e-4)

    def test_fixed_width_fresh_family_instead_loses_mean_memory(self):
        from outcome_updates import WINDOW_ATTENUATION
        self.assertLess(WINDOW_ATTENUATION**10000, 1e-40)
        self.assertGreater(finite_mean_gains(10000)[0], 0.84)

    def test_zero_width_is_an_explicit_different_preparation_set(self):
        h, alpha, _ = parameters(10, noise=0)
        self.assertEqual(h, 0)
        self.assertEqual(alpha, 1)
        self.assertEqual(finite_mean_gains(10, noise=0)[0], 1)

    def test_illegal_scalings_are_rejected(self):
        for args in ((0,), (1, 100), (1, 1, 2), (10, -1)):
            with self.assertRaises(ValueError):
                parameters(*args)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ShrinkingPreparationsTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    rows = []
    for count in (10, 100, 1000, 10000, 100000):
        h, alpha, eta = parameters(count)
        mean, floor = record_mean_and_tv_floor(count)
        rows.append({"step_count": count, "half_width": h, "alpha": alpha, "contrast": eta,
                     "mean_gains": finite_mean_gains(count).tolist(),
                     "scaled_record_mean_mu_zero": mean, "full_binary_record_tv_lower": floor})
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": tests.testsRun, "failures": len(tests.failures), "errors": len(tests.errors)},
               "noise": 1, "signal": 1, "rows": rows,
               "limiting_mean_gains": limit_mean_gains().tolist(),
               "liminf_full_binary_record_tv_lower": limit_record_tv_floor(),
               "scope": "Changed model sequence h_n=sqrt(nu/n), eta_n=sigma/sqrt(n), fixed within each experiment. Exact mean and record moment results with an analytic TV lower bound. One-step diffusion moment consistency is checked, but full path convergence is not proved. This does not remove the original fixed-width obstruction."}
    if args.write_results:
        Path(__file__).with_name("shrinking_preparations_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
