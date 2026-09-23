"""Round 384: minimal reversible redirection removes a dimension regularity input.

The qualitative circle exclusion uses the maximal-equicontinuous-factor theorem,
not these finite computations. This file independently checks a non-Lipschitz
covariant interface and its explicit Holder calibration modulus. Specified sphere
labels, actual redirection, qubit carrier and calibration remain inputs.
"""
import unittest

import numpy as np

from growing_stream_audit import main
from antipodal_dimension_bound_audit import (
    I, PAULI, coverage_radius_bound, grid, rounded_grid_point,
)


def unit(x):
    x = np.asarray(x, dtype=float)
    return x / np.linalg.norm(x)


def sigma(x):
    return np.einsum("i,ijk->jk", np.asarray(x, dtype=float), PAULI)


def odd_power(x, alpha):
    if alpha <= 0:
        raise ValueError("Positive power required.")
    x = np.asarray(x, dtype=float)
    return np.sign(x) * np.abs(x) ** alpha


def sphere_map(x, alpha=0.5):
    return unit(odd_power(x, alpha))


def effect(x, alpha=0.5):
    return (I + sigma(sphere_map(x, alpha))) / 2


def opnorm(x):
    return float(np.linalg.norm(x, 2))


def probability(e, rho):
    return float(np.trace(e @ rho).real)


def cross_matrix(axis):
    x, y, z = np.asarray(axis, dtype=float)
    return np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])


def rotation(axis, angle):
    n = unit(axis)
    a = cross_matrix(n)
    return np.eye(3) + np.sin(angle) * a + (1 - np.cos(angle)) * (a @ a)


def spin_rotation(axis, angle):
    return np.cos(angle / 2) * I - 1j * np.sin(angle / 2) * sigma(unit(axis))


def redirect(x, r, alpha=0.5):
    return sphere_map(r @ sphere_map(x, alpha), 1 / alpha)


def rotation_between(a, b):
    a, b = unit(a), unit(b)
    axis = np.cross(a, b)
    sine = np.linalg.norm(axis)
    cosine = float(np.clip(np.dot(a, b), -1, 1))
    if sine < 1e-14:
        if cosine > 0:
            return np.eye(3)
        axis = np.cross(a, np.eye(3)[np.argmin(abs(a))])
        return rotation(axis, np.pi)
    return rotation(axis, np.arctan2(sine, cosine))


def holder_constant(alpha):
    if not 0 < alpha < 1:
        raise ValueError("This audit uses 0 < alpha < 1.")
    return float(2 ** (-alpha) * 3 ** ((1 - alpha) / 2))


def modulus(distance, alpha=0.5):
    return min(1., holder_constant(alpha) * distance ** alpha)


def axis_point(epsilon):
    return np.array([epsilon, 0., np.sqrt(1 - epsilon ** 2)])


def axis_case(epsilon, alpha=0.5):
    x, pole = axis_point(epsilon), np.array([0., 0., 1.])
    distance = float(np.linalg.norm(x - pole))
    gap = opnorm(effect(x, alpha) - effect(pole, alpha))
    return {"epsilon": float(epsilon), "chord_distance": distance,
            "effect_gap": gap, "ratio": gap / distance,
            "scaled_asymptotic": gap / distance * epsilon ** (1 - alpha)}


def certificate_case(m=16, n=32, alpha=0.5, eta=0.002):
    points = grid(m, n)
    exact, observed = [], []
    for i, x in enumerate(points):
        e, f = effect(x, alpha), effect(-x, alpha)
        # Same actual state is used for both members of this antipodal pair.
        rho = e
        p, q = probability(e, rho), probability(f, rho)
        # Declared deterministic readout errors, each at most eta; both values
        # remain probabilities. Their worst pair error is attained at i=0.
        noise = eta * (0.5 + 0.5 * np.cos(i))
        p_hat, q_hat = p - noise, q + noise
        exact.append(abs(p - q))
        observed.append(abs(p_hat - q_hat))
    epsilon = coverage_radius_bound(m, n)
    return {"m": m, "n": n, "sample_count": len(points), "alpha": alpha,
            "readout_error_each": eta, "cover_radius_upper": epsilon,
            "holder_constant_effect_opnorm": holder_constant(alpha),
            "minimum_exact_gap": float(min(exact)),
            "minimum_observed_gap": float(min(observed)),
            "modulus_at_cover_radius": modulus(epsilon, alpha),
            "certified_gap_lower": float(min(observed) - 2 * eta
                                          - 2 * modulus(epsilon, alpha))}


def reconstructed_rotation(r, alpha=0.5):
    # Difference basis, not a presumed prior continuous unitary lift.
    base = np.eye(3)
    b = np.column_stack([sphere_map(x, alpha) - sphere_map(-x, alpha)
                         for x in base])
    after = np.column_stack([
        sphere_map(redirect(x, r, alpha), alpha)
        - sphere_map(redirect(-x, r, alpha), alpha) for x in base])
    return after @ np.linalg.inv(b)


def word_case():
    axes = ([1., 0., 0.], [0., 0., 1.])
    angles = (2 * np.pi * np.sqrt(2), 2 * np.pi * np.sqrt(3))
    generators = [rotation(a, t) for a, t in zip(axes, angles)]
    spinors = [spin_rotation(a, t) for a, t in zip(axes, angles)]
    word = [(0, 1), (1, -1), (0, 1), (1, 1), (1, 1), (0, -1)]
    x = unit([0.2, -0.7, 0.4])
    current, r, u = x.copy(), np.eye(3), I.copy()
    for index, power in word:
        step = generators[index] if power == 1 else generators[index].T
        spin = spinors[index] if power == 1 else spinors[index].conj().T
        current = redirect(current, step)
        r, u = step @ r, spin @ u
    return {"action_word_error": float(np.linalg.norm(current - redirect(x, r))),
            "spinor_covariance_error": opnorm(effect(current) - u @ effect(x) @ u.conj().T),
            "rotation_angles_over_2pi": [float(np.sqrt(2)), float(np.sqrt(3))],
            "word_length": len(word)}


def report():
    pole = np.array([0., 0., 1.])
    x = unit([0.21, 0.57, -0.8])
    y = unit([-0.41, 0.32, 0.29])
    r = rotation([0.2, -0.7, 0.5], 1.1)
    gx, gy = redirect(x, r), redirect(y, r)
    no_uniform = []
    for n in (4, 8, 16, 32, 64):
        alpha, epsilon = 1 / n, float(np.exp(-n))
        no_uniform.append({"alpha": alpha, **axis_case(epsilon, alpha)})
    return {
        "round": 384,
        "scope": "Conditional dimension lower bound for minimal reversible full-boundary actions and complete covariant qubit effects; no assumed compact/Lie source group, Lipschitz interface, or derived physical space. Finite certification still needs a known modulus and cover.",
        "scientific_baseline_through": 382,
        "analytic_results": {
            "lower_bound": "minimal continuous G-action on X homeomorphic to S^(d-1), whole-family unitary covariance and C_Delta=Herm_0(2) imply d>=3",
            "with_antipodal_separation": "adding all E_x != E_-x in the declared antipodal sphere identification gives d=3 by round 382",
            "external_theorem": "Hauser-Jager 1711.05672 Theorem 2.4 (MEF existence) and Theorem 2.12 (connected fibers); arbitrary topological group, no abelian/minimal premise for the external theorem",
            "circle_argument": "minimal circle action with a nontrivial equicontinuous factor is equicontinuous, has compact closure conjugate into O(2), and cannot have a rotational S^2 factor",
            "non_lipschitz": "effect_ratio(epsilon) ~ (1/2)*epsilon^(alpha-1), 0<alpha<1",
            "holder_modulus": "omega(r)=min(1,2^(-alpha)*3^((1-alpha)/2)*r^alpha)",
            "no_uniform_modulus_family": "alpha=1/n, epsilon=exp(-n): input separation ->0 but effect gap ->sqrt((1-1/sqrt(1+exp(-2)))/2)>0",
            "minimal_not_transitive": "two irrational-angle rotations about nonparallel axes generate a countable group with SO(3) closure; each orbit is dense and countable, not all S^2",
        },
        "axis_ratio_sequence": [axis_case(float(e)) for e in (1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12)],
        "certificate_cases": [certificate_case(8, 16), certificate_case(16, 32), certificate_case(32, 64)],
        "metric_comparison": {
            "raw_chord_before": float(np.linalg.norm(x - y)),
            "raw_chord_after": float(np.linalg.norm(gx - gy)),
            "pullback_distance_before": float(np.linalg.norm(sphere_map(x) - sphere_map(y))),
            "pullback_distance_after": float(np.linalg.norm(sphere_map(gx) - sphere_map(gy))),
        },
        "countable_control_word": word_case(),
        "no_uniform_modulus_family": no_uniform,
        "limiting_effect_gap": float(np.sqrt((1 - 1 / np.sqrt(1 + np.exp(-2))) / 2)),
        "new_inputs_for_example": ["specified S^2 labels and antipodes", "odd-power interface alpha", "actual conjugated reversible controls", "qubit effects and per-pair preparations", "declared chord-metric cover and known Holder modulus", "readout error bounds"],
        "not_established": ["sphere boundary and physical opposite directions from FUCP", "implementation or cost of all redirections", "known modulus from abstract minimality", "density or the MEF theorem from finite sampling", "a physical spatial metric, time, Lorentz geometry or gravity"],
        "primary_sources": ["https://arxiv.org/pdf/1711.05672", "https://arxiv.org/pdf/1904.12203"],
    }


class Checks(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(384)

    def points(self, count):
        data = self.rng.normal(size=(count, 3))
        return data / np.linalg.norm(data, axis=1)[:, None]

    def test_01_homeomorphism_inverse_oddness_and_effect_spectrum(self):
        for alpha in (0.25, 0.5, 0.8):
            for x in np.vstack((np.eye(3), self.points(25))):
                np.testing.assert_allclose(sphere_map(sphere_map(x, alpha), 1 / alpha), x, atol=1e-14)
                np.testing.assert_allclose(sphere_map(-x, alpha), -sphere_map(x, alpha))
                np.testing.assert_allclose(np.linalg.eigvalsh(effect(x, alpha)), [0, 1], atol=4e-16)
                np.testing.assert_allclose(effect(-x, alpha), I - effect(x, alpha), atol=3e-16)

    def test_02_global_holder_constant_uses_normalization_without_extra_factor(self):
        for alpha in (0.2, 0.5, 0.9):
            for a, b in ((0.3, -0.3), (0., 1.), (-0.9, 0.7)):
                scalar = abs(odd_power(a, alpha) - odd_power(b, alpha))
                self.assertLessEqual(scalar, 2 ** (1 - alpha) * abs(a - b) ** alpha + 1e-14)
            for x, y in zip(self.points(80), self.points(80)):
                a, b = odd_power(x, alpha), odd_power(y, alpha)
                r, s = np.linalg.norm(a), np.linalg.norm(b)
                self.assertGreaterEqual(min(r, s), 1 - 1e-14)
                self.assertAlmostEqual(np.linalg.norm(a - b) ** 2,
                                       (r - s) ** 2 + r * s * np.linalg.norm(a / r - b / s) ** 2)
                self.assertLessEqual(opnorm(effect(x, alpha) - effect(y, alpha)),
                                     modulus(np.linalg.norm(x - y), alpha) + 1e-14)

    def test_03_axis_sequence_matches_analytic_non_lipschitz_asymptotic(self):
        for alpha in (0.25, 0.5, 0.75):
            case = axis_case(1e-12, alpha)
            self.assertAlmostEqual(case["scaled_asymptotic"], 0.5, delta=3e-7)
        sequence = [axis_case(e) for e in (1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12)]
        self.assertTrue(all(a["ratio"] < b["ratio"] for a, b in zip(sequence, sequence[1:])))
        self.assertGreater(sequence[-1]["ratio"], 499999)

    def test_04_conjugated_action_composition_inverse_and_exact_transitivity(self):
        r, s = rotation([1, -2, 3], 0.7), rotation([-2, 1, 1], -1.2)
        for x, y in zip(self.points(20), self.points(20)):
            np.testing.assert_allclose(redirect(redirect(x, s), r), redirect(x, r @ s), atol=2e-14)
            np.testing.assert_allclose(redirect(redirect(x, r), r.T), x, atol=2e-14)
            between = rotation_between(sphere_map(x), sphere_map(y))
            np.testing.assert_allclose(redirect(x, between), y, atol=2e-14)

    def test_05_rodrigues_and_spinor_implement_the_same_complete_born_covariance(self):
        axis, angle = [1, 2, -3], 0.91
        r, u = rotation(axis, angle), spin_rotation(axis, angle)
        for x, v in zip(self.points(30), self.points(30)):
            e = effect(x)
            rho = (I + 0.63 * sigma(v)) / 2
            after = effect(redirect(x, r))
            np.testing.assert_allclose(after, u @ e @ u.conj().T, atol=1e-14)
            self.assertAlmostEqual(probability(after, u @ rho @ u.conj().T), probability(e, rho))

    def test_06_topological_covariance_does_not_mean_raw_metric_isometry(self):
        data = report()["metric_comparison"]
        self.assertGreater(abs(data["raw_chord_before"] - data["raw_chord_after"]), 0.05)
        self.assertAlmostEqual(data["pullback_distance_before"], data["pullback_distance_after"])

    def test_07_discrete_controls_match_words_and_have_the_declared_lie_closure(self):
        data = word_case()
        self.assertLess(data["action_word_error"], 1e-14)
        self.assertLess(data["spinor_covariance_error"], 1e-14)
        x, z = cross_matrix([1, 0, 0]), cross_matrix([0, 0, 1])
        comm = x @ z - z @ x
        np.testing.assert_array_equal(comm, cross_matrix([0, -1, 0]))
        self.assertEqual(np.linalg.matrix_rank(np.column_stack((x.ravel(), z.ravel(), comm.ravel()))), 3)
        # Irrational-angle density and all-orbit density are analytic, not tested by samples.

    def test_08_difference_basis_recovers_the_unique_action_and_its_limit(self):
        r = rotation([1, 2, 3], 0.6)
        np.testing.assert_allclose(reconstructed_rotation(r), r, atol=1e-14)
        errors = []
        for epsilon in (1e-2, 1e-4, 1e-6):
            nearby = rotation([2, -1, 3], epsilon) @ r
            reconstructed = reconstructed_rotation(nearby)
            np.testing.assert_allclose(reconstructed, nearby, atol=1e-14)
            errors.append(np.linalg.norm(reconstructed - r))
        self.assertGreater(errors[0], 90 * errors[1])
        self.assertGreater(errors[1], 90 * errors[2])

    def test_09_known_holder_cover_certifies_even_without_a_lipschitz_constant(self):
        coarse, fine = certificate_case(8, 16), certificate_case(16, 32)
        self.assertEqual(fine["sample_count"], 482)
        self.assertAlmostEqual(fine["minimum_exact_gap"], 1)
        self.assertAlmostEqual(fine["minimum_observed_gap"], 0.996)
        self.assertLess(coarse["certified_gap_lower"], 0)
        self.assertGreater(fine["certified_gap_lower"], 0.16)
        self.assertLess(fine["certified_gap_lower"], 1)

    def test_10_covered_probability_transfer_uses_the_same_probe_at_nearby_points(self):
        epsilon = coverage_radius_bound(16, 32)
        for x in np.vstack((self.points(120), [axis_point(1e-9)])):
            q = rounded_grid_point(x, 16, 32)
            self.assertLessEqual(np.linalg.norm(x - q), epsilon + 1e-14)
            rho = effect(q)
            at_q = probability(effect(q), rho) - probability(effect(-q), rho)
            at_x = probability(effect(x), rho) - probability(effect(-x), rho)
            self.assertLessEqual(abs(at_q - at_x), 2 * modulus(np.linalg.norm(x - q)) + 1e-14)
            self.assertGreaterEqual(opnorm(effect(x) - effect(-x)), abs(at_q) - 2 * modulus(epsilon) - 1e-14)

    def test_11_qualitative_contract_does_not_supply_one_uniform_modulus_for_all_models(self):
        limit = np.sqrt((1 - 1 / np.sqrt(1 + np.exp(-2))) / 2)
        values = [axis_case(float(np.exp(-n)), 1 / n) for n in (4, 8, 16, 32, 64)]
        self.assertLess(values[-1]["chord_distance"], 2e-28)
        self.assertGreater(values[-1]["effect_gap"], 0.17)
        self.assertAlmostEqual(values[-1]["effect_gap"], limit, places=14)
        self.assertLess(abs(values[-1]["effect_gap"] - limit), abs(values[0]["effect_gap"] - limit))


if __name__ == "__main__":
    main(__name__, "equivariant_boundary_regularity_audit", report)
