"""Round 386: strong-dilation interface and finite-scale reversal certificates.

The topology theorem is analytic, not inferred from this R^3 calibration model.
The model has base-dependent nonlinear charts: finite differences need not be
group differences. The old round-382 net is reused, not counted as a new result.
"""
import unittest

import numpy as np

from antipodal_dimension_bound_audit import grid, coverage_radius_bound
from growing_stream_audit import main


KAPPA = 0.4
BASE = np.array([0.7, -0.2, 0.1])
E2 = np.array([0., 1., 0.])
PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]],
                  [[1, 0], [0, -1]]], complex)


def coefficient(x):
    return KAPPA * np.sin(x[0])


def chart(x, y):
    v = np.asarray(y) - x
    return v + coefficient(x) * v[0] ** 2 * E2


def unchart(x, q):
    return x + q - coefficient(x) * q[0] ** 2 * E2


def delta(x, epsilon, y):
    return unchart(x, epsilon * chart(x, y))


def difference(x, epsilon, u, v):
    q = delta(x, epsilon, u)
    return delta(q, 1 / epsilon, delta(x, epsilon, v))


def limit_difference(x, u, v):
    return unchart(x, chart(x, v) - chart(x, u))


def approximate_sum(x, epsilon, u, v):
    return delta(x, 1 / epsilon, delta(delta(x, epsilon, u), epsilon, v))


def limit_sum(x, u, v):
    return unchart(x, chart(x, u) + chart(x, v))


def difference_bound(x, epsilon, u, v):
    a, b = chart(x, u), chart(x, v)
    step = np.linalg.norm(delta(x, epsilon, u) - x)
    return (step + abs(coefficient(x)) * epsilon * abs(b[0] ** 2 - a[0] ** 2)
            + KAPPA * epsilon * (1 + abs(a[0])) * (b[0] - a[0]) ** 2)


def raw_inverse_chart(epsilon, n, x=BASE):
    change = coefficient(x + epsilon * n) - coefficient(x)
    return (epsilon - 1) * (n + change * n[0] ** 2 * E2)


def reversal(epsilon, n, x=BASE):
    q = raw_inverse_chart(epsilon, n, x)
    return q / np.linalg.norm(q)


def reversal_bound(epsilon, raw_chart_error=0.):
    """Conditional chord bound, assuming the actually normalized vector is nonzero.

    A sufficient implementation guarantee for all errors of the given size is
    raw_chart_error < (1-epsilon)*(1-kappa*epsilon). A larger bound is only
    conditional; it does not certify that normalization can be performed.
    """
    # For epsilon < 1 and kappa*epsilon < 1, the raw point avoids the center.
    denominator = (1 - epsilon) * (1 - KAPPA * epsilon)
    if denominator <= 0:
        raise ValueError("The normalized inverse needs a nonzero annular margin.")
    return 2 * KAPPA * epsilon + 2 * raw_chart_error / denominator


def rotation(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])


def effect(f, scalar=0.5):
    return scalar * np.eye(2) + np.einsum('i,ijk->jk', f, PAULI)


def probe_probabilities(e):
    return np.array([[np.trace(e @ ((np.eye(2) + sign * p) / 2)).real
                      for sign in (1, -1)] for p in PAULI])


def reconstruct_traceless(probabilities):
    return (probabilities[:, 0] - probabilities[:, 1]) / 2


def gap_certificate(sample_gap, lipschitz, cover_radius, reversal_error, probability_error):
    return float(sample_gap - 2 * np.sqrt(3) * probability_error
                 - lipschitz * (2 * cover_radius + reversal_error))


def sufficient_shots(port_count, probability_error, failure_probability=0.01):
    # There are 12*port_count Bernoulli settings; Hoeffding and a union bound.
    return int(np.ceil(np.log(24 * port_count / failure_probability)
                       / (2 * probability_error ** 2)))


def certificate_case(epsilon, probability_error=None, noise_seed=386):
    visibility = epsilon ** 2
    tau = visibility / 100 if probability_error is None else probability_error
    points = grid()
    o = rotation(np.log(epsilon))  # No limiting orientation as epsilon -> 0.
    rng = np.random.default_rng(noise_seed)
    observed = []
    exact = []
    for n in points:
        r = reversal(epsilon, n)
        f, fr = visibility / 2 * (o @ n), visibility / 2 * (o @ r)
        pn, pr = probe_probabilities(effect(f)), probe_probabilities(effect(fr))
        pn += rng.uniform(-tau, tau, pn.shape)
        pr += rng.uniform(-tau, tau, pr.shape)
        observed.append(np.linalg.norm(reconstruct_traceless(pn) - reconstruct_traceless(pr)))
        exact.append(np.linalg.norm(f - fr))
    bound = gap_certificate(min(observed), visibility / 2, coverage_radius_bound(),
                            reversal_bound(epsilon), tau)
    shots = sufficient_shots(len(points), tau)
    return {"epsilon": epsilon, "visibility": visibility,
            "probability_error": tau, "ports": len(points),
            "reversal_chord_bound": reversal_bound(epsilon),
            "sample_observed_gap": float(min(observed)),
            "sample_true_approximate_reversal_gap": float(min(exact)),
            "certified_exact_inverse_gap": bound,
            "certificate_over_visibility": bound / visibility,
            "true_exact_inverse_gap": visibility,
            "sufficient_independent_shots_per_setting": shots,
            "sufficient_total_shots": 12 * len(points) * shots}


def high_frequency(k, theta):
    return np.array([np.cos(k * theta), np.sin(k * theta), 0.]) / 2


def report():
    u = unchart(BASE, np.array([0.6, -0.3, 0.4]))
    v = unchart(BASE, np.array([-0.2, 0.5, 0.1]))
    rows = []
    for epsilon in (0.2, 0.1, 0.05):
        actual = difference(BASE, epsilon, u, v)
        rows.append({"epsilon": epsilon,
                     "finite_difference_error": float(np.linalg.norm(actual - limit_difference(BASE, u, v))),
                     "analytic_error_bound": difference_bound(BASE, epsilon, u, v)})
    high_rows = [{"even_k": k, "reversal_angular_error": float(np.pi / k),
                  "actual_readout_gap": float(np.linalg.norm(high_frequency(k, 0.) - high_frequency(k, np.pi + np.pi / k))),
                  "exact_inverse_readout_gap": float(np.linalg.norm(high_frequency(k, 0.) - high_frequency(k, np.pi))),
                  "readout_modulus_at_reversal_error": 1.} for k in (20, 200, 2000)]
    return {
        "round": 386,
        "scope": {
            "analytic_geometric_input": "Complete locally compact separable position metric space; full nondegenerate strong dilatation axioms A0-A4 with scale group (0,infinity). These are extra inputs, not consequences of FUCP.",
            "topology": "On a compact physical neighborhood, continuity and nondegeneracy of the tangent metric give the same topology by compact-to-Hausdorff. Buliga yields a local positively graded Lie group. A small shell is a full boundary in the real position neighborhood, not automatically an equal-distance physical sphere.",
            "finite_scale_theorem": "A continuous traceless qubit readout on a full tangent shell transfers approximate-reversal separation to exact inverse separation when the former gap exceeds the readout modulus at the calibrated reversal error. One scale suffices; no readout limit or nonvanishing absolute visibility is required.",
            "full_dimension_three_requires": "The separate exact lower-bound contract of round 384 at that scale; approximate covariance or sampled full rank alone is insufficient.",
            "numerical_model": "The explicit base-dependent shear charts on a supplied R^3 are a nontrivial finite-operation calibration, not a derivation of its dimension or a demonstration on a presupposed curved manifold.",
            "resources": "Normalization, six calibrated probe preparations, complete shell coverage and finite error moduli are explicit. The numerical perturbations are synthetic bounded calibration errors. Shot counts are analytic sufficient budgets for a stated independent Bernoulli protocol; those shots were not simulated or executed, and the counts are not universal lower bounds.",
            "not_claimed": ["All positive graduations are Carnot stratifications", "The actual metric is Euclidean or shares tangent physical volume exponents", "A global displacement group or a smooth compatible atlas has been obtained", "Qualitative convergence supplies a usable finite convergence rate", "Small physical or angular endpoint error alone makes readout error small", "Derivation of SR or GR"],
        },
        "nonlinear_finite_operation": {"base": BASE.tolist(), "kappa": KAPPA,
                                        "difference_convergence": rows,
                                        "raw_inverse_norm_lower": "(1-epsilon)*(1-kappa*epsilon)",
                                        "normalized_reversal_chord_bound": "2*kappa*epsilon + 2*raw_error/((1-epsilon)*(1-kappa*epsilon))"},
        "relative_precision_certificates": [certificate_case(e) for e in (0.25, 0.125, 0.0625)],
        "fixed_precision_failure": certificate_case(0.001, probability_error=1e-5),
        "high_frequency_counterexample": high_rows,
        "finite_sampling": {"net_reused_from_round": 382,
                            "sphere_positive_calibration_ports": len(grid()),
                            "known_chord_cover_radius_upper": coverage_radius_bound(),
                            "unknown_dimension_coverage_warning": "This supplied S^2 net does not certify that an unknown physical boundary has no hidden directions."},
        "primary_sources": ["https://imar.ro/~mbuliga/buliga10.pdf",
                            "https://pi.math.cornell.edu/~hatcher/AT/ATch2.pdf"],
    }


class Checks(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(386)

    def test_01_actual_based_dilations_compose_exactly_with_inverses(self):
        for x, y in self.rng.normal(size=(25, 2, 3)):
            np.testing.assert_allclose(delta(x, 0.3, delta(x, 0.7, y)), delta(x, 0.21, y), atol=2e-15)
            np.testing.assert_allclose(delta(x, 1 / 0.3, delta(x, 0.3, y)), y, atol=3e-15)

    def test_02_rescaled_actual_metric_converges_with_an_explicit_uniform_compact_bound(self):
        for x, u, v in self.rng.uniform(-1, 1, size=(40, 3, 3)):
            a, b = chart(x, u), chart(x, v)
            epsilon = 0.07
            scaled = np.linalg.norm(delta(x, epsilon, u) - delta(x, epsilon, v)) / epsilon
            bound = abs(coefficient(x)) * epsilon * abs(a[0] ** 2 - b[0] ** 2)
            self.assertLessEqual(abs(scaled - np.linalg.norm(a - b)), bound + 6e-15)

    def test_03_finite_recentered_difference_converges_to_the_local_group_difference(self):
        for u, v in self.rng.uniform(-0.5, 0.5, size=(20, 2, 3)):
            u, v = unchart(BASE, u), unchart(BASE, v)
            for epsilon in (0.2, 0.05):
                error = np.linalg.norm(difference(BASE, epsilon, u, v) - limit_difference(BASE, u, v))
                self.assertLessEqual(error, difference_bound(BASE, epsilon, u, v) + 5e-14)
        n = np.array([0.6, -0.3, 0.4])
        u = unchart(BASE, n)
        errs = [np.linalg.norm(difference(BASE, e, u, BASE) - limit_difference(BASE, u, BASE))
                for e in (0.1, 0.01)]
        self.assertGreater(errs[0] / errs[1], 8.)

    def test_04_finite_sums_need_not_be_associative_although_the_limit_is(self):
        u, v, w = [unchart(BASE, q) for q in ([0.4, 0.2, 0.1], [-0.2, 0.3, 0.2], [0.1, -0.2, 0.3])]
        np.testing.assert_allclose(limit_sum(BASE, limit_sum(BASE, u, v), w),
                                   limit_sum(BASE, u, limit_sum(BASE, v, w)), atol=5e-16)
        defects = []
        for e in (0.2, 0.02):
            left = approximate_sum(BASE, e, approximate_sum(BASE, e, u, v), w)
            right = approximate_sum(BASE, e, u, approximate_sum(BASE, e, v, w))
            defects.append(np.linalg.norm(left - right))
        self.assertGreater(defects[0], 0.01)
        self.assertGreater(defects[0] / defects[1], 5.)

    def test_05_executed_raw_inverse_and_closed_form_agree_but_are_not_exact_involutions(self):
        e = 0.15
        for n in self.rng.normal(size=(25, 3)):
            n /= np.linalg.norm(n)
            u = unchart(BASE, n)
            np.testing.assert_allclose(chart(BASE, difference(BASE, e, u, BASE)),
                                       raw_inverse_chart(e, n), atol=2e-15)
        n = np.array([1., 0., 0.])
        twice = reversal(e, reversal(e, n))
        self.assertGreater(np.linalg.norm(twice - n), 1e-4)

    def test_06_normalization_and_physical_endpoint_error_have_different_scale_budgets(self):
        for n in self.rng.normal(size=(40, 3)):
            n /= np.linalg.norm(n)
            e = 0.12
            r = reversal(e, n)
            self.assertAlmostEqual(np.linalg.norm(r), 1.)
            self.assertLessEqual(np.linalg.norm(r + n), reversal_bound(e) + 1e-14)
            u, v = unchart(BASE, r), unchart(BASE, -n)
            physical = np.linalg.norm(delta(BASE, e, u) - delta(BASE, e, v))
            chord = np.linalg.norm(r + n)
            self.assertLessEqual(e * (1 - 2 * KAPPA * e) * chord, physical + 1e-15)
            self.assertLessEqual(physical, e * (1 + 2 * KAPPA * e) * chord + 1e-15)

    def test_07_six_probes_remove_an_unknown_direction_dependent_trace(self):
        for n in self.rng.normal(size=(25, 3)):
            n /= np.linalg.norm(n)
            f, scalar = 0.3 * n, 0.5 + 0.05 * n[2]
            e = effect(f, scalar)
            self.assertGreaterEqual(np.linalg.eigvalsh(e).min(), 0)
            self.assertLessEqual(np.linalg.eigvalsh(e).max(), 1)
            np.testing.assert_allclose(reconstruct_traceless(probe_probabilities(e)), f, atol=1e-16)

    def test_08_probability_error_vector_bound_is_sharp_and_enters_both_ports(self):
        tau = 0.004
        f = np.array([0.1, 0.05, -0.08])
        noise = np.array([[tau, -tau]] * 3)
        estimate = reconstruct_traceless(probe_probabilities(effect(f)) + noise)
        self.assertAlmostEqual(np.linalg.norm(estimate - f), np.sqrt(3) * tau)
        f2 = -f
        estimate2 = reconstruct_traceless(probe_probabilities(effect(f2)) - noise)
        self.assertLessEqual(abs(np.linalg.norm(estimate - estimate2) - np.linalg.norm(f - f2)), 2 * np.sqrt(3) * tau + 1e-15)

    def test_09_full_shell_certificate_stays_positive_as_visibility_vanishes_and_orientation_rotates(self):
        rows = [certificate_case(e) for e in (0.25, 0.125, 0.0625)]
        for row in rows:
            self.assertGreater(row['certificate_over_visibility'], 0.45)
            self.assertLess(row['certified_exact_inverse_gap'], row['true_exact_inverse_gap'])
        n = np.array([1., 0., 0.])
        self.assertGreater(np.linalg.norm(rotation(np.log(0.25)) @ n - rotation(np.log(0.0625)) @ n), 1.)

    def test_10_fixed_precision_can_fail_while_relative_precision_has_a_costed_certificate(self):
        self.assertLess(certificate_case(0.001, 1e-5)['certified_exact_inverse_gap'], 0)
        self.assertGreater(certificate_case(0.001)['certified_exact_inverse_gap'], 0)
        count = len(grid())
        tau = 0.0005
        n = sufficient_shots(count, tau)
        self.assertLessEqual(24 * count * np.exp(-2 * n * tau ** 2), 0.01)
        self.assertGreater(24 * count * np.exp(-2 * (n - 1) * tau ** 2), 0.01)
        self.assertAlmostEqual(sufficient_shots(count, tau / 2) / n, 4., places=6)

    def test_11_small_angular_error_is_not_small_readout_error_for_unbounded_frequency(self):
        for k in (20, 200, 2000):
            for t in (0., 0.123, -0.52):
                np.testing.assert_allclose(high_frequency(k, t + np.pi), high_frequency(k, t), atol=4e-13)
                self.assertAlmostEqual(np.linalg.norm(high_frequency(k, t + np.pi + np.pi / k) - high_frequency(k, t)), 1.)
            self.assertEqual(k % 2, 0)
        self.assertLess(np.pi / 2000, 0.002)

    def test_12_annular_nonzero_margin_is_needed_before_normalizing_a_noisy_inverse(self):
        epsilon = 1 - 1e-6
        n = np.array([0., 0., 1.])
        q = raw_inverse_chart(epsilon, n)
        noise = np.array([1 - epsilon, 0., 0.])
        actual = (q + noise) / np.linalg.norm(q + noise)
        self.assertLess(np.linalg.norm(noise), 2e-6)
        self.assertGreater(np.linalg.norm(actual + n), 0.7)
        self.assertGreater(reversal_bound(epsilon, np.linalg.norm(noise)), 1.)


if __name__ == '__main__':
    main(__name__, 'asymptotic_displacement_readout_audit', report)
