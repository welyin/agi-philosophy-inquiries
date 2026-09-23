"""Round 382: antipodal dimension bound and certified finite coverage.

Borsuk--Ulam supplies the analytic upper bound; finite calculations do not prove
that theorem or establish an unknown physical boundary. The sphere, its antipodal
map, a declared metric and a Lipschitz bound are separate operational inputs.
"""
import unittest

import numpy as np

from growing_stream_audit import main


I = np.eye(2, dtype=complex)
PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]],
                  [[1, 0], [0, -1]]], dtype=complex)
TETRA = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1],
                  [-1, -1, 1]], dtype=float) / np.sqrt(3)


def sigma(v):
    return np.einsum("i,ijk->jk", np.asarray(v, dtype=float), PAULI)


def effect(w, c=0.5):
    return c * I + sigma(w) / 2


def opnorm(a):
    return float(np.linalg.norm(a, 2))


def projection_effect(x, visibility=0.8):
    """S^2 or S^3 coordinate projection; x[:3] is not always a unit vector."""
    return effect(visibility * np.asarray(x)[:3])


def varying_trace_effect(x, alpha=0.3):
    x = np.asarray(x)
    return (0.5 + alpha * x[0]) * I + alpha * sigma(x[1:])


def sphere_samples(rng, count, dimension):
    a = rng.normal(size=(count, dimension))
    return a / np.linalg.norm(a, axis=1)[:, None]


def point(theta, phi):
    return np.array([np.sin(theta) * np.cos(phi),
                     np.sin(theta) * np.sin(phi), np.cos(theta)])


def grid(m=8, n=16):
    """Latitude/longitude net; poles are included only once."""
    if m < 2 or n < 2:
        raise ValueError("Need at least two latitude and longitude divisions.")
    return np.vstack(([0., 0., 1.],
                      [point(j * np.pi / m, k * 2 * np.pi / n)
                       for j in range(1, m) for k in range(n)],
                      [0., 0., -1.]))


def rounded_grid_point(x, m=8, n=16):
    theta = np.arccos(np.clip(x[2], -1, 1))
    phi = np.arctan2(x[1], x[0]) % (2 * np.pi)
    j = int(np.floor(theta * m / np.pi + 0.5))
    k = int(np.floor(phi * n / (2 * np.pi) + 0.5)) % n
    return point(j * np.pi / m, k * 2 * np.pi / n)


def coverage_radius_bound(m=8, n=16):
    # A meridian segment followed by a parallel segment bounds geodesic distance.
    # The same upper bound holds for chord distance, which is no larger.
    return np.pi / (2 * m) + np.pi / n


def certified_gap(sample_lower, lipschitz, epsilon, pair_error=0.,
                  antipode_lipschitz=1.):
    return float(sample_lower - pair_error
                 - lipschitz * (1 + antipode_lipschitz) * epsilon)


def best_probe(delta):
    eigenvalues, eigenvectors = np.linalg.eigh(delta)
    v = eigenvectors[:, int(np.argmax(np.abs(eigenvalues)))]
    return np.outer(v, v.conj())


def probability(e, state):
    return float(np.trace(e @ state).real)


def certificate_case(m=8, n=16):
    visibility, eta = 0.8, 0.002
    points = grid(m, n)
    exact, observed = [], []
    for index, x in enumerate(points):
        plus, minus = projection_effect(x), projection_effect(-x)
        state = best_probe(plus - minus)
        p, q = probability(plus, state), probability(minus, state)
        # Deterministic bounded readout errors, with each pair allowed to differ.
        e_plus = eta * np.cos(0.41 * index)
        e_minus = eta * np.sin(0.73 * index)
        exact.append(abs(p - q))
        observed.append(abs(p + e_plus - q - e_minus))
    observed = np.asarray(observed)
    epsilon = coverage_radius_bound(m, n)
    return {"m": m, "n": n, "points": points,
            "visibility": visibility, "readout_error_each": eta,
            "exact_gaps": np.asarray(exact), "observed_gaps": observed,
            "epsilon": epsilon, "lipschitz": visibility / 2,
            "lower": certified_gap(float(observed.min()), visibility / 2,
                                    epsilon, 2 * eta)}


def report():
    cal = certificate_case()
    coarse = certificate_case(2, 4)
    points = np.column_stack((TETRA, np.zeros(4)))
    pole = np.array([0., 0., 0., 1.])
    hidden_sample_gaps = [opnorm(projection_effect(x) - projection_effect(-x))
                          for x in points]
    rng = np.random.default_rng(3822026)
    trial = sphere_samples(rng, 2000, 3)
    distances = [np.linalg.norm(x - rounded_grid_point(x)) for x in trial]
    alpha = 0.3
    return {
        "round": 382,
        "scope": {
            "analytic_result": "Continuous fixed-trace qubit effects on a complete antipodal S^(d-1), with every antipodal pair distinct, imply d<=3 by Borsuk-Ulam. Full-family transitive covariance and contrast completeness, plus a locally Lipschitz map on a metric boundary of Hausdorff dimension d-1, additionally imply d>=3.",
            "inputs": ["Complete physical position boundary and its standard antipodal identification", "Selected complete qubit interface", "Fixed trace for the upper bound", "Actual whole-family covariant redirections and contrast completeness for the lower bound", "Metric regularity, certified coverage and Lipschitz constants for the finite certificate"],
            "not_proved": ["A physical boundary, its dimension or complete coverage from finite samples", "Any of these extra inputs from FUCP alone", "A new proof or numerical verification of Borsuk-Ulam", "An instrument-preserving or diamond-norm reconstruction of unknown states", "Lorentz structure, GR or the standard model"],
            "reference_scope": "Effect difference tensored with an arbitrary identity has the same operator norm. This bounds probabilities with a reference; it does not establish a complete quantum channel.",
        },
        "fixed_trace_hidden_fourth_coordinate": {
            "sample_count": len(points),
            "sample_antipodal_gap_min": min(hidden_sample_gaps),
            "sample_contrast_rank": int(np.linalg.matrix_rank(TETRA[1:] - TETRA[0])),
            "true_pole_antipodal_gap": opnorm(projection_effect(pole) - projection_effect(-pole)),
            "distance_of_missing_pole_to_sample_or_antipodes": float(np.sqrt(2)),
            "sample_gap_minus_cover_correction_at_necessary_radius": certified_gap(min(hidden_sample_gaps), 0.4, np.sqrt(2)),
        },
        "variable_trace_four_dimensional_counterexample": {
            "alpha": alpha,
            "legal_alpha_ceiling": float(1 / (2 * np.sqrt(2))),
            "global_antipodal_gap_lower": 2 * alpha,
            "global_antipodal_gap_upper": float(2 * np.sqrt(2) * alpha),
            "assumption_absent": "Trace and spectrum vary; this is not a whole-family unitary orbit.",
        },
        "finite_cover_certificate": {
            "m": cal["m"], "n": cal["n"], "sample_count": len(cal["points"]),
            "probability_settings": 2 * len(cal["points"]),
            "probability_settings_are_not_shot_counts": True,
            "analytic_chord_coverage_upper": float(cal["epsilon"]),
            "random_rounding_max_distance_for_crosscheck_only": float(max(distances)),
            "effect_operator_lipschitz": cal["lipschitz"],
            "readout_error_each": cal["readout_error_each"],
            "minimum_observed_probe_gap": float(cal["observed_gaps"].min()),
            "certified_global_gap_lower": cal["lower"],
            "exact_global_gap_of_this_known_calibration_model": 0.8,
            "coarse_net_certificate": coarse["lower"],
        },
        "strictly_weaker_than_injectivity": {
            "map": "S^1 angle phi -> 0.8*(cos(3phi),sin(3phi),0)",
            "covering_multiplicity": 3,
            "antipodal_gap": 0.8,
            "contrast_rank": 2,
            "scope": "This witnesses only the weakened upper-bound premise, not the full d=3 theorem.",
        },
        "primary_sources": [
            "https://arxiv.org/pdf/1008.1134",
            "https://arxiv.org/html/1203.0686",
            "https://arxiv.org/html/1206.0630v4",
        ],
    }


class Checks(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(382)

    def test_01_effect_positivity_and_gap_match_matrix_spectra(self):
        for dimension in (3, 4):
            for x in sphere_samples(self.rng, 40, dimension):
                e = projection_effect(x)
                self.assertGreaterEqual(np.linalg.eigvalsh(e).min(), -1e-14)
                self.assertLessEqual(np.linalg.eigvalsh(e).max(), 1 + 1e-14)
                self.assertAlmostEqual(np.trace(e).real, 1)
                self.assertAlmostEqual(opnorm(e - projection_effect(-x)),
                                       0.8 * np.linalg.norm(x[:3]))

    def test_02_antipodal_distinction_does_not_require_global_injection(self):
        def w(phi):
            return 0.8 * np.array([np.cos(3 * phi), np.sin(3 * phi), 0.])
        for phi in np.linspace(0, 2 * np.pi, 25):
            np.testing.assert_allclose(w(phi), w(phi + 2 * np.pi / 3), atol=5e-15)
            np.testing.assert_allclose(w(phi + np.pi), -w(phi), atol=5e-15)
            self.assertAlmostEqual(opnorm(effect(w(phi)) - effect(w(phi + np.pi))), 0.8)

    def test_03_four_dimensional_hidden_pair_despite_full_sample_contrasts(self):
        samples = np.column_stack((TETRA, np.zeros(4)))
        for x in samples:
            self.assertAlmostEqual(opnorm(projection_effect(x) - projection_effect(-x)), 0.8)
        self.assertEqual(np.linalg.matrix_rank(TETRA[1:] - TETRA[0]), 3)
        pole = np.array([0., 0., 0., 1.])
        np.testing.assert_array_equal(projection_effect(pole), projection_effect(-pole))

    def test_04_missing_poles_invalidate_an_asserted_fine_cover(self):
        samples = np.column_stack((TETRA, np.zeros(4)))
        samples = np.vstack((samples, -samples))
        pole = np.array([0., 0., 0., 1.])
        nearest = np.min(np.linalg.norm(samples - pole, axis=1))
        self.assertAlmostEqual(nearest, np.sqrt(2))
        self.assertGreater(certified_gap(0.8, 0.4, 0.1), 0)  # invalid coverage premise
        self.assertLessEqual(certified_gap(0.8, 0.4, nearest), 0)

    def test_05_variable_trace_allows_four_dimensional_antipodal_separation(self):
        alpha = 0.3
        special = np.vstack((np.eye(4), -np.eye(4),
                             [1 / np.sqrt(2), 1 / np.sqrt(2), 0, 0]))
        for x in np.vstack((special, sphere_samples(self.rng, 100, 4))):
            e = varying_trace_effect(x, alpha)
            self.assertGreaterEqual(np.linalg.eigvalsh(e).min(), 0)
            self.assertLessEqual(np.linalg.eigvalsh(e).max(), 1)
            np.testing.assert_allclose(varying_trace_effect(-x, alpha), I - e)
            gap = opnorm(e - varying_trace_effect(-x, alpha))
            expected = 2 * alpha * (abs(x[0]) + np.linalg.norm(x[1:]))
            self.assertAlmostEqual(gap, expected)
            self.assertGreaterEqual(gap, 2 * alpha - 1e-14)
        self.assertNotEqual(np.trace(varying_trace_effect(special[0])),
                            np.trace(varying_trace_effect(special[1])))

    def test_06_constructed_grid_contains_poles_and_antipodes(self):
        points = grid()
        self.assertEqual(len(points), 114)
        np.testing.assert_allclose(np.linalg.norm(points, axis=1), 1)
        for x in points:
            self.assertLess(np.min(np.linalg.norm(points + x, axis=1)), 2e-15)

    def test_07_rounding_gives_a_certified_cover_not_only_nearest_sample_fits(self):
        edge_cases = np.vstack((np.eye(3), -np.eye(3),
                                point(np.pi / 16, np.pi / 16),
                                point(np.pi - np.pi / 16, 2 * np.pi - np.pi / 16)))
        for x in np.vstack((edge_cases, sphere_samples(self.rng, 1200, 3))):
            q = rounded_grid_point(x)
            geodesic = np.arccos(np.clip(np.dot(x, q), -1, 1))
            self.assertLessEqual(geodesic, coverage_radius_bound() + 1e-14)
            self.assertLessEqual(np.linalg.norm(x - q), geodesic + 1e-14)
            self.assertLess(np.min(np.linalg.norm(grid() - q, axis=1)), 3e-15)

    def test_08_effect_error_deducts_two_per_pair_and_constant_is_attainable(self):
        e, f = projection_effect([0, 0, 1]), projection_effect([0, 0, -1])
        delta = 0.01
        e_hat, f_hat = e + delta * PAULI[2], f - delta * PAULI[2]
        self.assertAlmostEqual(opnorm(e_hat - e), delta)
        self.assertAlmostEqual(opnorm(f_hat - f), delta)
        self.assertAlmostEqual(opnorm(e_hat - f_hat) - opnorm(e - f), 2 * delta)
        for _ in range(50):
            a, b = sphere_samples(self.rng, 2, 3)
            estimate = opnorm(e + delta * sigma(a) - f - delta * sigma(b))
            self.assertLessEqual(abs(estimate - opnorm(e - f)), 2 * delta + 1e-14)

    def test_09_probability_probe_witnesses_the_operator_gap(self):
        for x in sphere_samples(self.rng, 60, 4):
            e, f = projection_effect(x), projection_effect(-x)
            state = best_probe(e - f)
            self.assertAlmostEqual(abs(probability(e, state) - probability(f, state)),
                                   opnorm(e - f))
            mixed = effect(0.5 * sphere_samples(self.rng, 1, 3)[0])
            self.assertLessEqual(abs(probability(e, mixed) - probability(f, mixed)),
                                 opnorm(e - f) + 1e-14)

    def test_10_noisy_covered_certificate_is_positive_and_globally_valid(self):
        cal = certificate_case()
        self.assertLessEqual(np.max(abs(cal["observed_gaps"] - cal["exact_gaps"])),
                             2 * cal["readout_error_each"] + 1e-14)
        self.assertGreater(cal["lower"], 0.47)
        self.assertLess(cal["lower"], 0.8)
        for x in sphere_samples(self.rng, 100, 3):
            self.assertLess(cal["lower"], opnorm(projection_effect(x) - projection_effect(-x)))

    def test_11_failing_a_sufficient_finite_certificate_does_not_show_a_collision(self):
        cal = certificate_case(2, 4)
        self.assertLess(cal["lower"], 0)
        np.testing.assert_allclose(cal["exact_gaps"], 0.8)

    def test_12_gap_lipschitz_constant_has_the_required_factor_two(self):
        # At the hidden pole of S^3 this ratio approaches 2K, not K.
        pole = np.array([0., 0., 0., 1.])
        for angle in (0.2, 0.05, 0.001):
            x = np.array([np.sin(angle), 0, 0, np.cos(angle)])
            h = opnorm(projection_effect(x) - projection_effect(-x))
            ratio = h / np.linalg.norm(x - pole)
            self.assertAlmostEqual(ratio, 0.8 * np.cos(angle / 2))
            self.assertGreater(ratio, 0.79)
            self.assertLessEqual(ratio, 0.8)

    def test_13_reference_cannot_increase_a_single_effect_gap(self):
        for x in sphere_samples(self.rng, 10, 4):
            delta = projection_effect(x) - projection_effect(-x)
            joint = np.kron(delta, np.eye(3))
            self.assertAlmostEqual(opnorm(joint), opnorm(delta))
            v = self.rng.normal(size=6) + 1j * self.rng.normal(size=6)
            v /= np.linalg.norm(v)
            self.assertLessEqual(abs(np.vdot(v, joint @ v)), opnorm(delta) + 1e-14)


if __name__ == "__main__":
    main(__name__, "antipodal_dimension_bound_audit", report)
