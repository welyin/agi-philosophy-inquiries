"""Exact innovations and joint state-record diffusion characteristics.

The path convergence argument is in note 15. Finite grids and trajectories
audit its algebra; they are not a numerical proof of weak convergence.
"""

import argparse
import json
import math
import platform
import unittest
from pathlib import Path

import numpy as np

from shrinking_preparations import branch_matrix, parameters


def exact_coefficients(state, count, noise=1.0, signal=1.0):
    u, v = np.asarray(state, dtype=float)
    _, alpha, eta = parameters(count, noise, signal)
    root = math.sqrt(1 - eta**2)
    if eta >= 1:
        raise ValueError("The uniform diffusion audit requires contrast < 1.")
    drift = count * np.array([(alpha - 1) * u, (alpha * root - 1) * v])
    coefficient = signal * alpha / (1 - eta**2 * u**2) * np.array([1 - u**2, -root * u * v])
    return drift, coefficient


def limiting_coefficients(state, noise=1.0, signal=1.0):
    u, v = np.asarray(state, dtype=float)
    rate = noise / 6
    drift = np.array([-rate * u, -(rate + signal**2 / 2) * v, signal * u, 0])
    diffusion = np.array([signal * (1 - u**2), -signal * u * v, 1, 1])
    return drift, diffusion


def joint_local_moments(state, count, noise=1.0, signal=1.0):
    state = np.asarray(state, dtype=float)
    _, _, eta = parameters(count, noise, signal)
    deltas, probabilities = [], []
    for outcome in (0, 1):
        sign = 2 * outcome - 1
        after = branch_matrix(outcome, count, noise, signal) @ np.r_[1, state]
        probability = after[0]
        delta = np.r_[after[1:] / probability - state,
                      sign / math.sqrt(count), (sign - eta * state[0]) / math.sqrt(count)]
        deltas.append(delta)
        probabilities.append(probability)
    deltas, probabilities = np.array(deltas), np.array(probabilities)
    mean = probabilities @ deltas
    centered = deltas - mean
    covariance = np.einsum("i,ij,ik->jk", probabilities, centered, centered)
    return count * mean, count * covariance, deltas, probabilities


def radius_defect_coefficients(state, noise=1.0, signal=1.0):
    u, v = state
    defect = 1 - u**2 - v**2
    rate = noise / 6
    return (2 * rate - (2 * rate + signal**2 * (1 - u**2)) * defect,
            -2 * signal * u * defect)


def sample_discrete_path(count=2048, noise=1.0, signal=1.0, seed=1501):
    """Diagnostic path of the exact discrete kernel, without an Euler scheme."""
    _, alpha, eta = parameters(count, noise, signal)
    state = np.array([0.4, 0.3]) * alpha
    rng = np.random.default_rng(seed)
    record = innovation = compensator = quadratic_variation = 0.0
    maximum_defect, max_radius = 0.0, float(np.linalg.norm(state))
    for uniform in rng.random(count):
        drift, coefficient = exact_coefficients(state, count, noise, signal)
        sign = 1 if uniform < (1 + eta * state[0]) / 2 else -1
        d_innovation = (sign - eta * state[0]) / math.sqrt(count)
        raw = branch_matrix((sign + 1) // 2, count, noise, signal) @ np.r_[1, state]
        actual = raw[1:] / raw[0]
        predicted = state + drift / count + coefficient * d_innovation
        maximum_defect = max(maximum_defect, float(np.max(abs(actual - predicted))))
        record += sign / math.sqrt(count)
        innovation += d_innovation
        compensator += signal * state[0] / count
        quadratic_variation += (1 - eta**2 * state[0]**2) / count
        state = actual
        max_radius = max(max_radius, float(np.linalg.norm(state)))
    return {"count": count, "noise": noise, "seed": seed,
            "final_state": state.tolist(), "record": record, "innovation": innovation,
            "record_minus_compensator_minus_innovation": record - compensator - innovation,
            "max_exact_recursion_error": maximum_defect,
            "predictable_innovation_quadratic_variation": quadratic_variation,
            "maximum_radius": max_radius, "allowed_radius": alpha}


def convergence_rows():
    rows = []
    for noise in (0, 1, 6):
        for count in (64, 256, 1024, 4096):
            _, alpha, _ = parameters(count, noise)
            drift_error = covariance_error = 0.0
            for radius in (0, 0.5 * alpha, alpha):
                for angle in np.linspace(0, 2 * math.pi, 25, endpoint=False):
                    state = radius * np.array([math.cos(angle), math.sin(angle)])
                    drift, covariance, _, _ = joint_local_moments(state, count, noise)
                    target, diffusion = limiting_coefficients(state, noise)
                    drift_error = max(drift_error, float(np.max(abs(drift - target))))
                    covariance_error = max(covariance_error, float(np.max(abs(covariance - np.outer(diffusion, diffusion)))))
            rows.append({"noise": noise, "count": count, "sampled_max_drift_error": drift_error,
                         "sampled_max_covariance_error": covariance_error})
    return rows


class JointDiffusionLimitTests(unittest.TestCase):
    def test_exact_innovation_representation_for_both_results(self):
        for noise in (0, 1, 6):
            for state in ([0, 0], [0.3, -0.5], [-0.8, 0.1]):
                count = 128
                _, _, eta = parameters(count, noise)
                drift, coefficient = exact_coefficients(state, count, noise)
                for outcome in (0, 1):
                    sign = 2 * outcome - 1
                    raw = branch_matrix(outcome, count, noise) @ np.r_[1, state]
                    predicted = np.asarray(state) + drift / count + coefficient * (sign - eta * state[0]) / math.sqrt(count)
                    np.testing.assert_allclose(predicted, raw[1:] / raw[0], atol=1e-14)

    def test_innovation_is_centered_and_has_correct_variance(self):
        for state in ([0, 0], [0.4, 0.3], [-0.7, 0.2]):
            drift, covariance, _, _ = joint_local_moments(state, 100)
            self.assertAlmostEqual(drift[3], 0)
            self.assertAlmostEqual(covariance[3, 3], 1 - state[0]**2 / 100)

    def test_record_and_state_share_one_noise_not_independent_noises(self):
        drift, covariance, _, _ = joint_local_moments([0.3, 0.4], 1000)
        self.assertEqual(np.linalg.matrix_rank(covariance, tol=1e-12), 1)
        self.assertGreater(covariance[0, 2], 0.9)
        self.assertLess(covariance[1, 2], -0.11)
        np.testing.assert_allclose(covariance[:, 2], covariance[:, 3], atol=1e-14)

    def test_drift_and_joint_covariance_converge_on_sampled_disk(self):
        rows = convergence_rows()
        for noise in (0, 1, 6):
            selected = [row for row in rows if row["noise"] == noise]
            errors = [max(row["sampled_max_drift_error"], row["sampled_max_covariance_error"]) for row in selected]
            self.assertTrue(all(a > b for a, b in zip(errors, errors[1:])))
            self.assertLess(errors[-1], 0.002)

    def test_joint_generator_cross_term_matches_direct_cubic_increment(self):
        # f(u,v,y,w)=u^2*y tests mixed derivatives involving the record.
        state, y = np.array([0.3, 0.4]), 0.7
        target_drift, diffusion = limiting_coefficients(state)
        u = state[0]
        target = 2 * u * y * target_drift[0] + u**2 * target_drift[2] + y * diffusion[0]**2 + 2 * u * diffusion[0]
        errors = []
        for count in (100, 1000, 10000):
            _, _, deltas, probabilities = joint_local_moments(state, count)
            actual = count * sum(p * ((u + d[0])**2 * (y + d[2]) - u**2 * y)
                                 for p, d in zip(probabilities, deltas))
            errors.append(abs(actual - target))
        self.assertTrue(errors[0] > errors[1] > errors[2])
        self.assertLess(errors[-1], 1e-4)

    def test_radius_defect_identity_matches_ito_formula(self):
        for noise in (0, 1, 6):
            for state in ([0, 0], [0.3, 0.4], [0.6, 0.8]):
                drift, diffusion = limiting_coefficients(state, noise)
                direct = (-2 * np.dot(state, drift[:2]) - np.dot(diffusion[:2], diffusion[:2]),
                          -2 * np.dot(state, diffusion[:2]))
                np.testing.assert_allclose(direct, radius_defect_coefficients(state, noise), atol=1e-14)

    def test_boundary_has_inward_drift_and_zero_normal_noise(self):
        for noise in (0, 1, 6):
            for angle in np.linspace(0, 2 * math.pi, 33):
                drift, diffusion = radius_defect_coefficients([math.cos(angle), math.sin(angle)], noise)
                self.assertAlmostEqual(drift, noise / 3)
                self.assertAlmostEqual(diffusion, 0)

    def test_zero_noise_discrete_updates_preserve_unit_circle(self):
        for angle in np.linspace(0, 2 * math.pi, 33):
            state = np.array([1, math.cos(angle), math.sin(angle)])
            for outcome in (0, 1):
                raw = branch_matrix(outcome, 100, noise=0) @ state
                self.assertAlmostEqual(np.linalg.norm(raw[1:] / raw[0]), 1)

    def test_complete_exact_path_and_predictable_variation(self):
        for noise in (0, 1, 6):
            row = sample_discrete_path(1024, noise)
            self.assertLess(abs(row["record_minus_compensator_minus_innovation"]), 1e-12)
            self.assertLess(row["max_exact_recursion_error"], 1e-14)
            self.assertLessEqual(row["maximum_radius"], row["allowed_radius"] + 1e-13)
            self.assertLessEqual(abs(row["predictable_innovation_quadratic_variation"] - 1), 1 / 1024 + 1e-13)

    def test_jump_sizes_shrink_uniformly_on_the_sampled_disk(self):
        for count in (64, 256, 1024):
            _, alpha, _ = parameters(count)
            for angle in np.linspace(0, 2 * math.pi, 33):
                _, _, deltas, _ = joint_local_moments(alpha * np.array([math.cos(angle), math.sin(angle)]), count)
                self.assertLess(float(np.max(abs(deltas))), 4 / math.sqrt(count))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(JointDiffusionLimitTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": tests.testsRun, "failures": len(tests.failures), "errors": len(tests.errors)},
               "joint_coordinate_order": ["u", "v", "record", "innovation"],
               "characteristic_checks": convergence_rows(),
               "exact_discrete_path_diagnostics": [sample_discrete_path(noise=noise) for noise in (0, 1, 6)],
               "scope": "Fixed setting, finite horizon, h_n=sqrt(nu/n), eta_n=sigma/sqrt(n), convergent legal initial summaries. Note 15 proves weak path convergence using tightness, limiting martingales and uniqueness. Sampled checks do not certify uniform bounds. No convergence in total variation or physical time is claimed."}
    if args.write_results:
        Path(__file__).with_name("joint_diffusion_limit_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
