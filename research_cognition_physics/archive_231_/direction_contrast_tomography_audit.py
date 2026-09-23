"""Round 380: qubit direction contrasts, orbit topology and calibration bounds.

A local topological position boundary and its complete covariant effect labelling
are inputs. No position manifold, physical dimension or autonomous device is
obtained from finite calibration data alone.
"""
import itertools
import unittest

import numpy as np

from growing_stream_audit import main


I = np.eye(2, dtype=complex)
PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]],
                  [[1, 0], [0, -1]]], dtype=complex)
TETRA = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1],
                  [-1, -1, 1]], dtype=float) / np.sqrt(3)


def hermitian(v):
    return np.einsum("i,ijk->jk", np.asarray(v, dtype=float), PAULI)


def rho(r):
    return (I + hermitian(r)) / 2


def effect(w):
    return (I + hermitian(w)) / 2


def probabilities(vectors, states):
    return np.array([[np.trace(effect(w) @ rho(r)).real for r in states]
                     for w in vectors])


def latitude(angles, a=0.6, b=0.8):
    angles = np.asarray(angles)
    return np.column_stack((a * np.cos(angles), a * np.sin(angles),
                            np.full_like(angles, b)))


def contrast_matrix(vectors):
    vectors = np.asarray(vectors)
    return (vectors[1:] - vectors[0]) / 2


def double_difference(table):
    return table[1:, 1:] - table[1:, :1] - table[:1, 1:] + table[0, 0]


def singular_min(matrix):
    return float(np.linalg.svd(matrix, compute_uv=False)[-1])


def unitary(axis, angle):
    axis = np.asarray(axis, float)
    axis = axis / np.linalg.norm(axis)
    values, vectors = np.linalg.eigh(hermitian(axis))
    return (vectors * np.exp(-0.5j * angle * values)) @ vectors.conj().T


def rotation(axis, angle):
    n = np.asarray(axis, float)
    n = n / np.linalg.norm(n)
    cross = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]],
                      [-n[1], n[0], 0]])
    return (np.cos(angle) * np.eye(3) + (1 - np.cos(angle)) * np.outer(n, n)
            + np.sin(angle) * cross)


def tetrahedral_group():
    out = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((-1, 1), repeat=3):
            if np.prod(signs) != 1:
                continue
            r = np.diag(signs) @ np.eye(3, dtype=int)[list(perm)]
            if round(np.linalg.det(r)) == 1:
                out.append(r)
    return out


def latitude_reconstruct(readings, a=0.6, b=0.8):
    p0, p1, p2, p3 = readings
    return np.array([(p0 - p2) / a, (p1 - p3) / a,
                     (2 * np.mean(readings) - 1) / b])


def calibration_bound(eta, preparation_error, probability_error):
    if eta <= 0 or preparation_error < 0 or probability_error < 0:
        raise ValueError("Positive calibration visibility and nonnegative errors required.")
    return (6 * preparation_error + 12 * probability_error) / eta


def calibration_case():
    rng = np.random.default_rng(380)
    eta, prep_error, read_error = 0.75, 0.003, 0.0005
    nominal = np.vstack((np.zeros(3), eta * np.eye(3)))
    deviations = rng.normal(size=(4, 3))
    deviations *= prep_error / np.linalg.norm(deviations, axis=1)[:, None]
    actual = nominal + deviations
    table = probabilities(TETRA, actual)
    noise = rng.uniform(-read_error, read_error, size=(4, 4))
    measured = table + noise
    estimate = double_difference(measured) / eta
    actual_contrast = contrast_matrix(TETRA)
    delta = calibration_bound(eta, prep_error, read_error)
    target = np.array([0.31, -0.46, 0.54])
    output_error = 0.0007
    output_noise = rng.uniform(-output_error, output_error, size=4)
    readings = probabilities(TETRA, [target])[:, 0] + output_noise
    recovered = np.linalg.solve(estimate, readings[1:] - readings[0])
    recovery_bound = (delta + 2 * np.sqrt(3) * output_error) / singular_min(estimate)
    return {"eta": eta, "preparation_error": prep_error, "read_error": read_error,
            "actual_preparations": actual, "actual_table": table,
            "measured_table": measured, "estimate": estimate,
            "actual_contrast": actual_contrast, "delta": delta,
            "target": target, "output_error": output_error,
            "recovered": recovered, "recovery_bound": recovery_bound}


def report():
    angles = np.arange(4) * np.pi / 2
    lat = latitude(angles)
    target = np.array([0.2, -0.3, 0.4])
    readings = probabilities(lat, [target])[:, 0]
    ordinary = np.vstack((np.r_[1., np.zeros(3)],
                          np.column_stack((np.ones(4), lat)) / 2))
    states = np.array([[0.4, 0, 0.5], [0.4, 0, -0.5]])
    pair = probabilities(lat, states)
    cal = calibration_case()
    d = contrast_matrix(TETRA)
    return {
        "round": 380,
        "scope": (
            "A conditional bridge from a complete qubit effect interface to the "
            "dimension of an already supplied local topological position boundary. "
            "Continuous injective labelling, transitive physically implemented "
            "redirection covariance and contrast tomography are additional inputs. "
            "Finite calibration certifies selected qubit contrasts subject to error "
            "budgets; it cannot establish the global boundary or absence of hidden directions."),
        "latitude_counterexample": {
            "boundary": "S^1 (local position dimension d=2)",
            "a": 0.6, "b": 0.8,
            "normalized_state_tomography_rank": int(np.linalg.matrix_rank(lat)),
            "effects_plus_normalization_rank": int(np.linalg.matrix_rank(ordinary)),
            "contrast_rank": int(np.linalg.matrix_rank(contrast_matrix(lat))),
            "haar_mean_effect_bloch_vector": [0., 0., 0.8],
            "reconstruction_error": float(np.linalg.norm(
                latitude_reconstruct(readings) - target)),
            "same_contrast_state_bloch_vectors": states.tolist(),
            "contrast_disagreement": float(np.max(np.abs(
                (pair[1:, 0] - pair[0, 0]) - (pair[1:, 1] - pair[0, 1])))),
            "orientation_mean_probability_gap": float(abs(np.mean(pair[:, 0])
                                                          - np.mean(pair[:, 1]))),
        },
        "conditional_theorem": {
            "covariance_scope": "The same group action maps the entire effect family "
                "by E_(g.x)=U_g E_x U_g^dagger; pairwise unitary equivalence is insufficient.",
            "contrast_condition": "span_R{E_x-E_y}=Herm_0(2)",
            "equivalent_condition_within_one_unitary_orbit":
                "ordinary normalized-state tomography AND scalar Haar-mean effect",
            "topological_conclusion": "X homeomorphic to S^2; if X=S^(d-1), d=3",
            "d1_audit": "Two settings have contrast dimension at most one.",
            "without_sphere_topology": "The 12-element tetrahedral rotation group "
                "has a four-point, contrast-complete orbit with scalar mean.",
        },
        "four_setting_certificate": {
            "tetrahedral_vectors": TETRA.tolist(),
            "contrast_matrix": d.tolist(),
            "singular_values": np.linalg.svd(d, compute_uv=False).tolist(),
            "minimum_settings": 4,
            "exact_mean_effect": "I/2",
        },
        "finite_calibration": {
            "eta": cal["eta"],
            "preparation_bloch_error_bound": cal["preparation_error"],
            "each_probability_error_bound": cal["read_error"],
            "estimated_contrast_matrix": cal["estimate"].tolist(),
            "actual_matrix_error": float(np.linalg.norm(cal["estimate"] - d, 2)),
            "matrix_error_bound": cal["delta"],
            "measured_sigma_min": singular_min(cal["estimate"]),
            "certified_sigma_min_lower": singular_min(cal["estimate"]) - cal["delta"],
            "true_sigma_min": singular_min(d),
            "new_state_bloch_error": float(np.linalg.norm(cal["recovered"] - cal["target"])),
            "new_state_bloch_error_bound": cal["recovery_bound"],
            "calibration_setting_preparation_pairs": 16,
            "new_state_measurement_settings": 4,
        },
        "resource_scope": [
            "Complete qubit carrier and independently calibrated preparations are inputs.",
            "Repeated preparation is required; there is no single-copy state or direction readout.",
            "16 calibration probabilities and 4 new-state probabilities are not 20 shots.",
            "Error budgets must be justified by separate preparation and statistical evidence.",
            "Carrier tomography does not require all degrees of freedom of every particle to be directional.",
            "No derived spatial dynamics, metric scale, manifold existence or autonomous controller."],
        "sources": ["https://arxiv.org/html/1206.0630v4"],
    }


class Checks(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(380)

    def test_01_effect_validity_including_projective_latitude(self):
        vectors = np.vstack((latitude(np.linspace(0, 2 * np.pi, 23)), TETRA))
        for w in vectors:
            eig = np.linalg.eigvalsh(effect(w))
            np.testing.assert_allclose(eig, [0, 1], atol=4e-16)
            np.testing.assert_allclose(effect(w) @ effect(w), effect(w), atol=5e-16)

    def test_02_born_trace_and_bloch_probability_agree(self):
        states = self.rng.normal(size=(20, 3))
        states /= (1 + np.linalg.norm(states, axis=1))[:, None]
        vectors = np.vstack((TETRA, latitude(np.arange(7))))
        np.testing.assert_allclose(probabilities(vectors, states),
                                   (1 + vectors @ states.T) / 2, atol=3e-16)

    def test_03_latitude_transitive_qubit_covariance(self):
        for phi, angle in self.rng.normal(size=(16, 2)):
            u = unitary([0, 0, 1], angle)
            np.testing.assert_allclose(u @ effect(latitude([phi])[0]) @ u.conj().T,
                                       effect(latitude([phi + angle])[0]), atol=7e-16)
        self.assertGreater(np.linalg.norm(latitude([0])[0] - latitude([0.2])[0]), 0)

    def test_04_latitude_full_tomography_two_independent_inversions(self):
        vectors = latitude(np.arange(4) * np.pi / 2)
        self.assertEqual(np.linalg.matrix_rank(vectors), 3)
        for _ in range(15):
            r = self.rng.normal(size=3)
            r /= 1 + np.linalg.norm(r)
            p = probabilities(vectors, [r])[:, 0]
            fit = np.linalg.lstsq(vectors, 2 * p - 1, rcond=None)[0]
            np.testing.assert_allclose(fit, r, atol=6e-16)
            np.testing.assert_allclose(latitude_reconstruct(p), r, atol=6e-16)

    def test_05_latitude_same_direction_contrasts_hide_state_information(self):
        angles = np.linspace(0, 2 * np.pi, 97, endpoint=False)
        vectors = latitude(angles)
        states = np.array([[0.4, 0, 0.5], [0.4, 0, -0.5]])
        p = probabilities(vectors, states)
        np.testing.assert_allclose(p[:, 0] - p[0, 0], p[:, 1] - p[0, 1], atol=3e-16)
        self.assertAlmostEqual(np.mean(p[:, 0]) - np.mean(p[:, 1]), 0.4)
        opposite = probabilities(latitude(angles + np.pi), states)
        np.testing.assert_allclose(p - opposite,
                                   np.tile((0.6 * 0.4 * np.cos(angles))[:, None], (1, 2)),
                                   atol=4e-16)
        self.assertEqual(np.linalg.matrix_rank(contrast_matrix(vectors)), 2)

    def test_06_scalar_mean_without_tomography_is_insufficient(self):
        equator = latitude(np.arange(8) * np.pi / 4, a=1, b=0)
        np.testing.assert_allclose(np.mean([effect(w) for w in equator], axis=0), I / 2,
                                   atol=2e-16)
        self.assertEqual(np.linalg.matrix_rank(contrast_matrix(equator)), 2)
        np.testing.assert_allclose(probabilities(equator, [[0, 0, 1], [0, 0, -1]]), 0.5)
        two = np.array([[0, 0, 1], [0, 0, -1]])
        self.assertEqual(np.linalg.matrix_rank(contrast_matrix(two)), 1)

    def test_07_tetrahedral_contrast_spectrum_and_independent_matrix_recovery(self):
        d = contrast_matrix(TETRA)
        np.testing.assert_allclose(np.linalg.svd(d, compute_uv=False),
                                   np.array([2, 1, 1]) / np.sqrt(3), atol=5e-16)
        # Obtain the response using actual density-matrix probabilities.
        states = np.vstack((np.zeros(3), np.eye(3)))
        np.testing.assert_allclose(double_difference(probabilities(TETRA, states)), d,
                                   atol=3e-16)
        np.testing.assert_allclose(np.mean([effect(w) for w in TETRA], axis=0), I / 2)

    def test_08_three_settings_have_a_genuine_indistinguishable_contrast_pair(self):
        vectors = TETRA[:3]
        d = contrast_matrix(vectors)
        null = np.linalg.svd(d, full_matrices=True)[2][-1]
        p = probabilities(vectors, [0.8 * null, -0.8 * null])
        np.testing.assert_allclose(p[1:, 0] - p[0, 0], p[1:, 1] - p[0, 1], atol=3e-16)
        self.assertGreater(np.linalg.norm(rho(0.8 * null) - rho(-0.8 * null)), 1)

    def test_09_discrete_tetrahedral_orbit_preserves_the_topology_hypothesis(self):
        group = tetrahedral_group()
        self.assertEqual(len(group), 12)
        keys = {tuple(r.reshape(-1)) for r in group}
        for r in group:
            for s in group:
                self.assertIn(tuple((r @ s).reshape(-1)), keys)
        orbit = np.array([r @ TETRA[0] for r in group])
        self.assertEqual(len(np.unique(np.round(orbit, 12), axis=0)), 4)
        np.testing.assert_allclose(orbit.mean(axis=0), 0, atol=2e-16)
        self.assertEqual(np.linalg.matrix_rank(contrast_matrix(orbit)), 3)

    def test_10_contrast_certificate_is_independent_of_qubit_frame(self):
        for axis, angle in zip(self.rng.normal(size=(8, 3)), self.rng.normal(size=8)):
            r, u = rotation(axis, angle), unitary(axis, angle)
            rotated = TETRA @ r.T
            for old, new in zip(TETRA, rotated):
                np.testing.assert_allclose(u @ effect(old) @ u.conj().T, effect(new),
                                           atol=9e-16)
            np.testing.assert_allclose(np.linalg.svd(contrast_matrix(rotated), compute_uv=False),
                                       np.linalg.svd(contrast_matrix(TETRA), compute_uv=False),
                                       atol=9e-16)

    def test_11_double_difference_identity_with_nonideal_preparations(self):
        cal = calibration_case()
        states = cal["actual_preparations"]
        self.assertLess(np.max(np.linalg.norm(states, axis=1)), 1)
        r = (states[1:] - states[0]).T
        np.testing.assert_allclose(double_difference(cal["actual_table"]),
                                   cal["actual_contrast"] @ r, atol=3e-16)

    def test_12_readout_error_constant_is_attainable(self):
        eps = 0.0001
        errors = np.full((4, 4), eps)
        errors[0, 1:] = -eps
        errors[1:, 0] = -eps
        np.testing.assert_allclose(double_difference(errors), np.full((3, 3), 4 * eps))
        self.assertAlmostEqual(np.linalg.norm(double_difference(errors), 2), 12 * eps)

    def test_13_joint_preparation_readout_bound_and_positive_certificate(self):
        cal = calibration_case()
        actual = np.linalg.norm(cal["estimate"] - cal["actual_contrast"], 2)
        self.assertLessEqual(actual, cal["delta"])
        certified = singular_min(cal["estimate"]) - cal["delta"]
        self.assertGreater(certified, 0.5)
        self.assertLessEqual(certified, singular_min(cal["actual_contrast"]))

    def test_14_small_noise_fake_full_rank_does_not_pass_certificate(self):
        eta, eps = 0.75, 0.0005
        states = np.vstack((np.zeros(3), eta * np.eye(3)))
        vectors = latitude(np.arange(4) * np.pi / 2)
        table = probabilities(vectors, states)
        table[1, 3] += eps
        estimate = double_difference(table) / eta
        self.assertEqual(np.linalg.matrix_rank(estimate), 3)
        self.assertLessEqual(singular_min(estimate) - calibration_bound(eta, 0, eps), 0)

    def test_15_state_recovery_respects_certified_operator_bound(self):
        cal = calibration_case()
        dhat, delta = cal["estimate"], cal["delta"]
        eps = cal["output_error"]
        bound = (delta + 2 * np.sqrt(3) * eps) / singular_min(dhat)
        for _ in range(30):
            r = self.rng.normal(size=3)
            r /= 1 + np.linalg.norm(r)
            noise = self.rng.uniform(-eps, eps, size=4)
            p = probabilities(TETRA, [r])[:, 0] + noise
            estimate = np.linalg.solve(dhat, p[1:] - p[0])
            self.assertLessEqual(np.linalg.norm(estimate - r), bound)


if __name__ == "__main__":
    main(__name__, "direction_contrast_tomography_audit", report)
