"""Round 370: record fields versus operational radar dimension.

Inputs are flat d-dimensional space, a calibrated inertial clock, isotropic
speed c, immediate identifiable echoes and a specified directional response.
Nothing here selects spatial d, proves a microscopic ontology, or derives GR.
"""
import itertools
import unittest
import numpy as np
from growing_stream_audit import main


def norm_geometry(x):
    x = np.asarray(x, float)
    r = float(np.linalg.norm(x))
    if x.ndim != 1 or r <= 0:
        raise ValueError("The local radar chart excludes the observer itself.")
    u = x / r
    return r, u, np.eye(len(x)) - np.outer(u, u)


def record(time, x, c=1.):
    r, u, _ = norm_geometry(x)
    if c <= 0:
        raise ValueError("c must be positive.")
    return {"send": float(time-r/c), "receive": float(time+r/c),
            "port": u.tolist()}


def inverse_record(value, c=1.):
    time = (value["send"] + value["receive"]) / 2
    radius = c * (value["receive"] - value["send"]) / 2
    return time, radius * np.asarray(value["port"])


def event_data(z, c=1.):
    value = record(z[0], z[1:], c)
    return np.r_[value["send"], value["receive"], value["port"]]


def event_jacobian(x, c=1.):
    r, u, p = norm_geometry(x)
    result = np.zeros((len(x)+2, len(x)+1))
    result[:2, 0] = 1.
    result[0, 1:], result[1, 1:] = -u/c, u/c
    result[2:, 1:] = p/r
    return result


def numerical_jacobian(fun, x, h=2e-6):
    x = np.asarray(x, float)
    columns = []
    for j in range(len(x)):
        step = np.eye(len(x))[j]*h
        columns.append((fun(x+step)-fun(x-step))/(2*h))
    return np.column_stack(columns)


def tangent_rows(u):
    # An orthonormal basis of u-perp, with zero rows in d=1.
    return np.linalg.svd(np.asarray(u)[None, :], full_matrices=True)[2][1:]


def partial_port_jacobian(x, components, c=1., controlled_send=False):
    r, u, p = norm_geometry(x)
    a = tangent_rows(u)[:components]
    if controlled_send:
        return np.vstack((2*u[None, :]/c, a@p/r))
    out = np.zeros((2+components, len(x)+1))
    out[:2, 0] = 1.
    out[0, 1:], out[1, 1:] = -u/c, u/c
    out[2:, 1:] = a@p/r
    return out


def event_fisher(x, c=1., sigma_t=.2, sigma_u=.1):
    j = event_jacobian(x, c)
    weights = np.r_[np.repeat(1/sigma_t**2, 2),
                    np.repeat(1/sigma_u**2, len(x))]
    return j.T@(weights[:, None]*j)


def spatial_fisher(x, c=1., sigma_t=.2, sigma_u=.1,
                   controlled_send=False):
    r, u, p = norm_geometry(x)
    radial = (4. if controlled_send else 2.)/(c*sigma_t)**2
    return radial*np.outer(u, u) + p/(r*sigma_u)**2


def emission_nuisance_fisher(x, c=1., sigma_t=.2, sigma_u=.1):
    r, u, p = norm_geometry(x)
    j = np.zeros((len(x)+2, len(x)+1))
    j[:2, 0] = 1.  # unknown true emission time e, read twice with noise
    j[1, 1:] = 2*u/c
    j[2:, 1:] = p/r
    w = np.r_[np.repeat(1/sigma_t**2, 2),
              np.repeat(1/sigma_u**2, len(x))]
    return j.T@(w[:, None]*j)


def profile_emission(info):
    return info[1:, 1:] - np.outer(info[1:, 0], info[0, 1:])/info[0, 0]


def coarse_port(x):
    # Three fixed axis-pair sectors in d=3; ties excluded in local tests.
    return int(np.argmax(np.abs(x)))


def coarse_static_record(x, send=4., anchor=None, c=1.):
    anchor = np.zeros(len(x)) if anchor is None else np.asarray(anchor)
    relative = np.asarray(x)-anchor
    return {"send": send, "receive": send+2*np.linalg.norm(relative)/c,
            "port": coarse_port(relative)}


def range_jacobian(x, anchors):
    return np.array([norm_geometry(np.asarray(x)-a)[1] for a in anchors])


def pack_grid(coordinates, bits):
    # One integer stores a d-coordinate finite grid, without being a chart.
    return sum(((int(value) >> b) & 1) << (len(coordinates)*b+j)
               for j, value in enumerate(coordinates) for b in range(bits))


def unpack_grid(label, dimension, bits):
    return tuple(sum(((label >> (dimension*b+j)) & 1) << b
                     for b in range(bits)) for j in range(dimension))


def categorical_port(v, epsilon=.3):
    # EXTRA protocol: a random detected mark, not a chosen launch label.
    # v=(u_x,u_y) on the known upper-hemisphere direction chart.
    v = np.asarray(v, float)
    a = epsilon/3*np.array([[1., 0.], [0., 1.], [-1., -1.]])
    probabilities = np.ones(3)/3 + a@v
    return probabilities, a


def categorical_information(v, samples=1, epsilon=.3):
    probabilities, a = categorical_port(v, epsilon)
    return samples*a.T@(a/probabilities[:, None])


def categorical_estimator_covariance(v, samples=1, epsilon=.3):
    probabilities, _ = categorical_port(v, epsilon)
    p = probabilities[:2]
    return 9/(samples*epsilon**2)*(np.diag(p)-np.outer(p, p))


def dimension_cases():
    rows = []
    for d in range(1, 5):
        x = 2*np.arange(1., d+1)/np.linalg.norm(np.arange(1., d+1))
        value = record(4., x)
        info = event_fisher(x)
        rows.append({"spatial_dimension": d, "field_count": len(value),
                     "record": value, "direction_local_dimension": d-1,
                     "event_jacobian_rank": int(np.linalg.matrix_rank(event_jacobian(x))),
                     "time_only_event_rank": int(np.linalg.matrix_rank(event_jacobian(x)[:2])),
                     "controlled_send_spatial_rank": int(np.linalg.matrix_rank(
                         partial_port_jacobian(x, d-1, controlled_send=True))),
                     "event_fisher_eigenvalues": np.linalg.eigvalsh(info).tolist(),
                     "profiled_spatial_fisher_eigenvalues": np.linalg.eigvalsh(
                         profile_emission(emission_nuisance_fisher(x))).tolist()})
    return rows


def report():
    p, q = np.array([.8, .36, .48]), np.array([.8, .48, .36])
    second = np.array([0., .5, 0.])
    anchors = np.array([[0., 0., 0.], second, [0., 0., .5]])
    v = np.array([.2, .3])
    probability, a = categorical_port(v)
    binary_p = np.array([probability[0], probability[1:].sum()])
    binary_a = np.vstack((a[0], a[1:].sum(axis=0)))
    coarse_j = numerical_jacobian(
        lambda x: np.array([coarse_static_record(x)["receive"], coarse_port(x)]), p)
    return {
        "round": 370,
        "scope": {
            "baseline": "Frozen through 367; no dependence on 368 or 369.",
            "claim": "Three named record fields do not select spatial dimension. Directional response rank, not the number of labels, controls local spatial identifiability.",
            "inputs": ["flat space of specified dimension d", "known inertial observer and calibrated clock", "isotropic signal speed c", "immediate identifiable echoes", "specified port response and, where used, independent readout noise", "known extra observer baselines"],
            "not_proved": ["emergence of dimension three", "a physical identification of Bloch directions with propagation directions", "all geometrical facts are produced by consensus", "general relativity from F+U+C+P", "global identifiability from a local full-rank Fisher matrix"]},
        "dimensions": dimension_cases(),
        "general_rank": {"event": "2 + rank(D port on tangent sphere)",
                         "static_target_controlled_emission": "1 + rank(D port on tangent sphere)",
                         "full_direction": "d+1 event; d spatial"},
        "coarse_port_pair": {
            "point_p": p.tolist(), "point_q": q.tolist(),
            "first_record_p": coarse_static_record(p),
            "first_record_q": coarse_static_record(q),
            "first_record_local_spatial_rank": int(np.linalg.matrix_rank(coarse_j)),
            "second_anchor": second.tolist(),
            "second_record_p": coarse_static_record(p, anchor=second),
            "second_record_q": coarse_static_record(q, anchor=second),
            "second_roundtrip_gap": float(abs(coarse_static_record(p, anchor=second)["receive"]-coarse_static_record(q, anchor=second)["receive"])),
            "two_anchor_local_rank": int(np.linalg.matrix_rank(range_jacobian(p, anchors[:2]))),
            "three_anchor_local_rank": int(np.linalg.matrix_rank(range_jacobian(p, anchors))),
            "three_anchor_global_mirror_ambiguity": "p and (-p_x,p_y,p_z) have the same three ranges"},
        "finite_label_encoding": {"bits_per_coordinate": 3,
            "one_integer_codes_per_dimension": {str(d): 2**(3*d) for d in range(1, 5)},
            "interpretation": "Exact reversible finite-grid storage; no continuous scalar spatial coordinate is inferred."},
        "extra_categorical_detector": {
            "direction_chart": v.tolist(), "probabilities": probability.tolist(),
            "three_outcome_angular_fisher_rank": int(np.linalg.matrix_rank(categorical_information(v))),
            "merged_two_outcome_angular_fisher_rank": int(np.linalg.matrix_rank(binary_a.T@(binary_a/binary_p[:, None]))),
            "single_sample_estimator_covariance": categorical_estimator_covariance(v).tolist(),
            "10000_sample_estimator_covariance": categorical_estimator_covariance(v, 10000).tolist(),
            "controlled_port_label_fisher_rank": 0,
            "extra_resources": "Repeated probes of one fixed target, known directional likelihood, calibration and stored outcome counts; a chosen launch label alone has no target-dependent likelihood."},
        "direction_precision": {"sigma_time": .2, "sigma_direction": .1,
            "range_1_spatial_fisher_eigenvalues": np.linalg.eigvalsh(spatial_fisher(p)).tolist(),
            "range_10_spatial_fisher_eigenvalues": np.linalg.eigvalsh(spatial_fisher(10*p)).tolist(),
            "interpretation": "Transverse Fisher information falls as r^-2; full rank alone supplies no range-independent accuracy."},
        "sources": ["https://arxiv.org/html/0708.0170", "https://arxiv.org/html/gr-qc/0104077"]}


class Checks(unittest.TestCase):
    def test_01_same_three_fields_roundtrip_inverse_all_dimensions(self):
        for d in range(1, 5):
            for sign in (-1., 1.):
                x = sign*np.arange(1., d+1)
                value = record(3.7, x, c=1.4)
                t2, x2 = inverse_record(value, c=1.4)
                self.assertEqual(len(value), 3)
                self.assertAlmostEqual(t2, 3.7)
                np.testing.assert_allclose(x2, x, atol=2e-15)
                self.assertAlmostEqual(1.4**2*(3.7-value["send"])**2,
                                       float(x@x))
                self.assertAlmostEqual(1.4**2*(value["receive"]-3.7)**2,
                                       float(x@x))

    def test_02_analytic_radar_jacobian_matches_independent_differences(self):
        for d in range(1, 5):
            x = np.arange(.3, d+.3)
            actual = numerical_jacobian(lambda z: event_data(z, c=1.7), np.r_[2.3, x])
            np.testing.assert_allclose(actual, event_jacobian(x, c=1.7), atol=4e-10)

    def test_03_time_and_direction_have_different_ranks(self):
        for d in range(1, 5):
            j = event_jacobian(np.arange(1., d+1))
            self.assertEqual(np.linalg.matrix_rank(j), d+1)
            self.assertEqual(np.linalg.matrix_rank(j[:2]), 2)
            self.assertEqual(np.linalg.matrix_rank(j[2:], tol=1e-12), d-1)

    def test_04_partial_port_rank_counts_actual_tangent_responses(self):
        for d in range(1, 5):
            x = np.arange(1., d+1)
            for k in range(d):
                self.assertEqual(np.linalg.matrix_rank(partial_port_jacobian(x, k)), 2+k)
                self.assertEqual(np.linalg.matrix_rank(partial_port_jacobian(x, k, controlled_send=True)), 1+k)

    def test_05_gaussian_fisher_closed_form_and_time_nuisance(self):
        for d in range(1, 5):
            x = np.arange(1., d+1)
            info = event_fisher(x, c=1.7)
            expected = np.zeros_like(info)
            expected[0, 0] = 2/.2**2
            expected[1:, 1:] = spatial_fisher(x, c=1.7)
            np.testing.assert_allclose(info, expected, atol=2e-14)
            np.testing.assert_allclose(profile_emission(info), expected[1:, 1:], atol=2e-14)

    def test_06_controlled_emission_and_profiled_unknown_emission(self):
        for d in range(1, 5):
            x = np.arange(1., d+1)
            sigma_t, sigma_u, c = .2, .1, 1.7
            info = emission_nuisance_fisher(x, c)
            np.testing.assert_allclose(profile_emission(info), spatial_fisher(x, c), atol=2e-14)
            r, u, p = norm_geometry(x)
            known_j = np.vstack((2*u[None, :]/c/sigma_t, p/r/sigma_u))
            np.testing.assert_allclose(known_j.T@known_j, spatial_fisher(x, c, controlled_send=True), atol=2e-14)
            self.assertEqual(np.linalg.matrix_rank(known_j), d)

    def test_07_same_radius_and_same_coarse_port_leave_distinct_points(self):
        p, q = np.array([.8, .36, .48]), np.array([.8, .48, .36])
        self.assertGreater(np.linalg.norm(p-q), .1)
        self.assertEqual(coarse_static_record(p), coarse_static_record(q))
        j = numerical_jacobian(lambda x: np.array([coarse_static_record(x)["receive"], coarse_port(x)]), p)
        self.assertEqual(np.linalg.matrix_rank(j), 1)
        self.assertEqual(float(np.linalg.norm(j[1])), 0.)

    def test_08_extra_observer_resolves_pair_but_only_adds_one_local_rank(self):
        p, q = np.array([.8, .36, .48]), np.array([.8, .48, .36])
        anchors = np.array([[0., 0., 0.], [0., .5, 0.]])
        gap = abs(coarse_static_record(p, anchor=anchors[1])["receive"]-
                  coarse_static_record(q, anchor=anchors[1])["receive"])
        self.assertGreater(gap, .13)
        self.assertEqual(np.linalg.matrix_rank(range_jacobian(p, anchors)), 2)

    def test_09_repetition_does_not_add_rank_and_full_local_rank_is_not_global(self):
        p = np.array([.8, .36, .48])
        repeats = np.zeros((20, 3))
        self.assertEqual(np.linalg.matrix_rank(range_jacobian(p, repeats)), 1)
        anchors = np.array([[0., 0., 0.], [0., .5, 0.], [0., 0., .5]])
        self.assertEqual(np.linalg.matrix_rank(range_jacobian(p, anchors)), 3)
        mirror = p*np.array([-1., 1., 1.])
        np.testing.assert_allclose(np.linalg.norm(p-anchors, axis=1),
                                   np.linalg.norm(mirror-anchors, axis=1), atol=1e-15)

    def test_10_a_single_integer_can_encode_finite_grids_of_any_dimension(self):
        for d in range(1, 5):
            codes = set()
            for point in itertools.product(range(8), repeat=d):
                label = pack_grid(point, bits=3)
                self.assertEqual(unpack_grid(label, d, bits=3), point)
                codes.add(label)
            self.assertEqual(codes, set(range(2**(3*d))))

    def test_11_finite_random_outcomes_may_carry_continuous_direction_information(self):
        v = np.array([.2, .3])
        probabilities, a = categorical_port(v)
        self.assertAlmostEqual(float(probabilities.sum()), 1.)
        self.assertGreater(float(probabilities.min()), .2)
        np.testing.assert_allclose(numerical_jacobian(lambda w: categorical_port(w)[0], v), a, atol=2e-11)
        self.assertEqual(np.linalg.matrix_rank(categorical_information(v)), 2)
        binary_p = np.array([probabilities[0], probabilities[1:].sum()])
        binary_a = np.vstack((a[0], a[1:].sum(axis=0)))
        self.assertEqual(np.linalg.matrix_rank(binary_a.T@(binary_a/binary_p[:, None])), 1)
        # A launch distribution chosen independently of the target has zero derivative.
        independent = numerical_jacobian(lambda w: np.array([.2, .3, .5]), v)
        self.assertEqual(np.linalg.matrix_rank(independent), 0)

    def test_12_repeated_categorical_readout_has_exact_covariance_and_cost(self):
        v = np.array([.2, .3])
        probabilities, _ = categorical_port(v)
        np.testing.assert_allclose((3*probabilities[:2]-1)/.3, v, atol=1e-15)
        for n in (1, 100, 10000):
            cov = categorical_estimator_covariance(v, n)
            np.testing.assert_allclose(cov, np.linalg.inv(categorical_information(v, n)), atol=4e-14)
        np.testing.assert_allclose(categorical_estimator_covariance(v, 100),
                                   100*categorical_estimator_covariance(v, 10000), atol=1e-15)

    def test_13_joint_frame_change_preserves_ranges_in_every_dimension(self):
        for d in range(1, 5):
            x, anchor = np.arange(1., d+1), np.arange(1., d+1)*.1
            rotation = np.linalg.qr(np.eye(d)+.17*np.ones((d, d)))[0]
            shift = np.arange(.3, d+.3)
            transformed = rotation@x+shift - (rotation@anchor+shift)
            original = record(3., x-anchor)
            changed = record(3.+4., transformed)
            self.assertAlmostEqual(changed["receive"]-changed["send"],
                                   original["receive"]-original["send"])
            np.testing.assert_allclose(changed["port"], rotation@np.asarray(original["port"]), atol=4e-15)
            self.assertAlmostEqual(changed["send"]-original["send"], 4.)

    def test_14_full_rank_does_not_give_range_independent_transverse_accuracy(self):
        x = np.array([.8, .36, .48])
        _, u, _ = norm_geometry(x)
        rows = tangent_rows(u)
        near, far = spatial_fisher(x), spatial_fisher(10*x)
        np.testing.assert_allclose(rows@near@rows.T, 100*rows@far@rows.T, atol=2e-13)
        self.assertAlmostEqual(float(u@near@u), float(u@far@u))
        self.assertEqual(np.linalg.matrix_rank(far), 3)


if __name__ == "__main__":
    main(__name__, "radar_port_dimension_audit", report)
