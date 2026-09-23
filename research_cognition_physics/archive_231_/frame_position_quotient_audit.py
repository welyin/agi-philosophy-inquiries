"""Round 378: a specified frame bundle, position quotient and control labels.

SO(4), its unit S^3 position readout, the Levi-Civita lift and classical
control/calibration are inputs. This is not a selection of physical dimension.
"""
import unittest

import numpy as np

from growing_stream_audit import main


I4 = np.eye(4)
NORTH = I4[:, 3]


def plane(i, j, n=4):
    out = np.zeros((n, n))
    out[i, j], out[j, i] = 1.0, -1.0
    return out


def generator(u):
    out = np.zeros((4, 4))
    out[:3, 3] = u
    out[3, :3] = -np.asarray(u)
    return out


def horizontal_rotation(u):
    u = np.asarray(u, dtype=float)
    a = np.linalg.norm(u)
    if a < 1e-15:
        return I4.copy()
    k = generator(u / a)
    return I4 + np.sin(a) * k + (1.0 - np.cos(a)) * (k @ k)


def matrix_exponential_skew(k):
    vals, vecs = np.linalg.eigh(1j * k)
    out = (vecs * np.exp(-1j * vals)) @ vecs.conj().T
    if np.max(np.abs(out.imag)) > 1e-12:
        raise ArithmeticError("Unexpected imaginary roundoff.")
    return out.real


def fiber(r):
    out = I4.copy()
    out[:3, :3] = r
    return out


def random_so(rng, n):
    q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    q[:, 0] *= np.linalg.det(q)
    return q


def position(q):
    return q[:, 3]


def endpoint(q, controls):
    out = np.array(q, dtype=float)
    for u in controls:
        out = out @ horizontal_rotation(u)
    return out


def rk4_endpoint(q, controls, subdivisions=64):
    out = np.array(q, dtype=float)
    h = 1.0 / subdivisions
    for u in controls:
        k = generator(u)
        for _ in range(subdivisions):
            a = out @ k
            b = (out + 0.5 * h * a) @ k
            c = (out + 0.5 * h * b) @ k
            d = (out + h * c) @ k
            out += h * (a + 2 * b + 2 * c + d) / 6
    return out


def triangle(i, j, angle):
    e = np.eye(3)
    return [0.5 * np.pi * e[i], angle * e[j], -0.5 * np.pi * e[i]]


def cost(controls):
    return sum(float(np.linalg.norm(u)) for u in controls)


def position_effects(q):
    """Four separate binary effects, not four outcomes of one POVM."""
    return 0.5 * (1.0 + position(q))


def direction_effect(q, body_index, marker):
    return 0.5 * (1.0 + np.dot(marker, q[:, body_index]))


def lift(q, a):
    """Horizontal lift of the base field Y_a(x)=a-(a.x)x."""
    return q @ generator(q[:, :3].T @ a)


def base_field(x, a):
    return a - x * np.dot(x, a)


def lift_derivative(q, a, tangent):
    return (tangent @ generator(q[:, :3].T @ a)
            + q @ generator(tangent[:, :3].T @ a))


def lift_bracket(q, a, b):
    return (lift_derivative(q, b, lift(q, a))
            - lift_derivative(q, a, lift(q, b)))


def so3_factors(r):
    """R = product exp(phi K_ij); Givens elimination, no Euler-angle chart."""
    work = np.array(r, dtype=float)
    factors = []
    for i, j in ((0, 2), (0, 1), (1, 2)):
        angle = np.arctan2(work[j, i], work[i, i])
        work = matrix_exponential_skew(angle * plane(i, j, 3)) @ work
        factors.append((i, j, -float(angle)))
    if np.linalg.norm(work - np.eye(3)) > 2e-12:
        raise ArithmeticError("Givens elimination failed.")
    return factors


def synthesize(target):
    """At most one horizontal segment and three 3-segment position loops."""
    x = position(target)
    theta = np.arccos(np.clip(x[3], -1.0, 1.0))
    norm = np.linalg.norm(x[:3])
    if norm > 1e-13:
        u = theta * x[:3] / norm
    else:
        u = np.array([theta, 0.0, 0.0])
    qh = horizontal_rotation(u)
    r = (qh.T @ target)[:3, :3]
    controls = [u]
    for i, j, angle in so3_factors(r):
        controls.extend(triangle(i, j, -angle))
    return controls


def report():
    controls = triangle(0, 1, np.pi / 2)
    final = endpoint(I4, controls)
    numerical = rk4_endpoint(I4, controls)
    s = 0.4
    u = np.array([s, 0.0, 0.0])
    straight = endpoint(I4, [u])
    uncorrected = endpoint(final, [u])
    corrected = endpoint(final, [final[:3, :3].T @ u])
    rng = np.random.default_rng(378)
    target = random_so(rng, 4)
    recipe = synthesize(target)
    q = random_so(rng, 4)
    a, b = rng.normal(size=(2, 4))
    bracket = lift_bracket(q, a, b)
    x = position(q)
    projected = np.dot(x, a) * b - np.dot(x, b) * a
    return {
        "round": 378,
        "scope": (
            "Given SO(4) oriented frames over a unit S^3 with the specified "
            "Levi-Civita horizontal control, calibrated position effects and "
            "frame-dependent relabelling. Horizontal nonintegrability is compatible "
            "with a three-dimensional position quotient. Neither physical dimension "
            "nor an autonomous quantum dynamics is derived."),
        "dimensions": {"configuration": 6, "vertical_fiber": 3,
                       "horizontal_rank": 3, "position": 3, "bracket_rank": 6},
        "triangle": {
            "vertices": [NORTH.tolist(), I4[:, 0].tolist(), I4[:, 1].tolist(),
                         position(final).tolist()],
            "position_return_error": float(np.linalg.norm(position(final) - NORTH)),
            "final_frame": final.tolist(),
            "frame_rotation_radians": float(np.arccos(
                np.clip((np.trace(final[:3, :3]) - 1.0) / 2, -1.0, 1.0))),
            "length": cost(controls),
            "all_four_position_effect_error": float(np.max(np.abs(
                position_effects(final) - position_effects(I4)))),
            "orientation_effect_initial": direction_effect(I4, 0, I4[:, 0]),
            "orientation_effect_final": direction_effect(final, 0, I4[:, 0]),
            "rk4_full_frame_error": float(np.linalg.norm(numerical - final)),
        },
        "control_label_audit": {
            "subsequent_pulse_length": s,
            "uncorrected_endpoint_chord_gap": float(np.linalg.norm(
                position(straight) - position(uncorrected))),
            "uncorrected_marker_probability_gap": float(abs(
                position_effects(straight)[0] - position_effects(uncorrected)[0])),
            "compensated_endpoint_error": float(np.linalg.norm(
                position(straight) - position(corrected))),
            "relabelling": "Q -> Q diag(R,1), u -> R^T u",
            "same_untranslated_body_label_is_projectable": False,
        },
        "full_configuration_synthesis": {
            "pulses": len(recipe), "length": cost(recipe),
            "endpoint_residual": float(np.linalg.norm(endpoint(I4, recipe) - target)),
            "general_bound": "At most 10 prescribed segments; length <= 7*pi.",
        },
        "projectable_lift": {
            "bracket_projection_error": float(np.linalg.norm(bracket[:, 3] - projected)),
            "criterion": "D intersect V = 0 and d pi(D) = TM; D need not be involutive.",
        },
        "resource_scope": {
            "inputs": ["unit sphere geometry", "frame and marker calibration",
                       "signed horizontal control", "pulse clock and scheduling",
                       "reliable orientation-dependent label conversion"],
            "not_provided": ["autonomous controller", "energy budget",
                             "finite-noise position certification", "dimension selection"],
        },
        "sources": ["https://arxiv.org/pdf/1610.07359"],
    }


class Checks(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(378)

    def test_01_horizontal_rotation_independent_exponential(self):
        for u in [np.zeros(3), np.array([np.pi, 0, 0]),
                  *self.rng.normal(size=(12, 3))]:
            q = horizontal_rotation(u)
            np.testing.assert_allclose(q, matrix_exponential_skew(generator(u)),
                                       atol=2e-14)
            np.testing.assert_allclose(q.T @ q, I4, atol=2e-14)
            self.assertAlmostEqual(np.linalg.det(q), 1.0, places=13)

    def test_02_fibers_and_position_effect_completeness(self):
        for _ in range(12):
            q, r = random_so(self.rng, 4), random_so(self.rng, 3)
            other = q @ fiber(r)
            np.testing.assert_allclose(position(q), position(other), atol=2e-14)
            np.testing.assert_allclose(q.T @ other, fiber(r), atol=2e-14)
            np.testing.assert_allclose(2 * position_effects(q) - 1, position(q))

    def test_03_horizontal_vertical_and_quotient_ranks(self):
        q = random_so(self.rng, 4)
        horizontal = [q @ generator(e) for e in np.eye(3)]
        vertical = [q @ plane(i, j) for i, j in ((0, 1), (0, 2), (1, 2))]
        self.assertEqual(np.linalg.matrix_rank(
            np.column_stack([v.reshape(-1) for v in horizontal])), 3)
        self.assertEqual(np.linalg.matrix_rank(
            np.column_stack([v.reshape(-1) for v in horizontal + vertical])), 6)
        np.testing.assert_allclose(np.column_stack([v[:, 3] for v in vertical]), 0)
        projected = np.column_stack([v[:, 3] for v in horizontal])
        np.testing.assert_allclose(projected.T @ projected, np.eye(3), atol=2e-14)
        self.assertEqual(np.linalg.matrix_rank(0.5 * projected), 3)

    def test_04_brackets_and_vertical_curvature_sign(self):
        h = [generator(e) for e in np.eye(3)]
        commutators = []
        for i, j in ((0, 1), (0, 2), (1, 2)):
            comm = h[i] @ h[j] - h[j] @ h[i]
            np.testing.assert_allclose(comm, -plane(i, j), atol=0)
            commutators.append(comm)
        self.assertEqual(np.linalg.matrix_rank(
            np.column_stack([a.reshape(-1) for a in h + commutators])), 6)

    def test_05_geodesic_and_parallel_frame_equations(self):
        q, u = random_so(self.rng, 4), self.rng.normal(size=3)
        f, x = q[:, :3], position(q)
        tangent = q @ generator(u)
        np.testing.assert_allclose(tangent[:, 3], f @ u, atol=2e-14)
        np.testing.assert_allclose(tangent[:, :3], -np.outer(x, u), atol=2e-14)
        np.testing.assert_allclose((I4 - np.outer(x, x)) @ tangent[:, :3],
                                   0, atol=2e-14)
        np.testing.assert_allclose((tangent @ generator(u))[:, 3],
                                   -np.dot(u, u) * x, atol=2e-14)

    def test_06_exact_triangle_holonomy_all_coordinate_planes(self):
        for i, j in ((0, 1), (0, 2), (1, 2)):
            for angle in (-2.1, -0.4, 0.0, 0.7, np.pi / 2):
                end = endpoint(I4, triangle(i, j, angle))
                expected = matrix_exponential_skew(-angle * plane(i, j))
                np.testing.assert_allclose(end, expected, atol=2e-14)
                np.testing.assert_allclose(position(end), NORTH, atol=2e-14)

    def test_07_triangle_vertices_and_readout_separation(self):
        controls = triangle(0, 1, np.pi / 2)
        q1 = endpoint(I4, controls[:1])
        q2 = endpoint(I4, controls[:2])
        q3 = endpoint(I4, controls)
        np.testing.assert_allclose(position(q1), I4[:, 0], atol=2e-14)
        np.testing.assert_allclose(position(q2), I4[:, 1], atol=2e-14)
        np.testing.assert_allclose(position_effects(q3), position_effects(I4),
                                   atol=2e-14)
        self.assertAlmostEqual(direction_effect(I4, 0, I4[:, 0]), 1.0)
        self.assertAlmostEqual(direction_effect(q3, 0, I4[:, 0]), 0.5)
        self.assertAlmostEqual(cost(controls), 1.5 * np.pi)

    def test_08_independent_rk4_full_frame_convergence(self):
        controls = triangle(0, 1, np.pi / 2)
        truth = endpoint(I4, controls)
        e16 = np.linalg.norm(rk4_endpoint(I4, controls, 16) - truth)
        e32 = np.linalg.norm(rk4_endpoint(I4, controls, 32) - truth)
        e64 = np.linalg.norm(rk4_endpoint(I4, controls, 64) - truth)
        self.assertGreater(e16 / e32, 15)
        self.assertGreater(e32 / e64, 15)
        self.assertLess(e64, 2e-8)

    def test_09_fixed_body_label_fails_position_state_closure(self):
        loop = endpoint(I4, triangle(0, 1, np.pi / 2))
        for s in (0.1, 0.4, 0.9):
            u = np.array([s, 0, 0])
            p, q = position(endpoint(I4, [u])), position(endpoint(loop, [u]))
            self.assertAlmostEqual(np.linalg.norm(p - q), np.sqrt(2) * np.sin(s))
            self.assertAlmostEqual((p[0] - q[0]) / 2, np.sin(s) / 2)

    def test_10_relabelled_protocol_is_covariant_at_every_step(self):
        for _ in range(8):
            q, r = random_so(self.rng, 4), random_so(self.rng, 3)
            other = q @ fiber(r)
            for u in self.rng.normal(size=(6, 3)):
                q = endpoint(q, [u])
                other = endpoint(other, [r.T @ u])
                np.testing.assert_allclose(other, q @ fiber(r), atol=3e-14)
                np.testing.assert_allclose(position(other), position(q), atol=3e-14)

    def test_11_projectable_base_fields_and_lift_brackets(self):
        for _ in range(8):
            q, r = random_so(self.rng, 4), random_so(self.rng, 3)
            a, b = self.rng.normal(size=(2, 4))
            x, f = position(q), q[:, :3]
            np.testing.assert_allclose(lift(q, a)[:, 3], base_field(x, a), atol=2e-14)
            np.testing.assert_allclose(lift(q @ fiber(r), a), lift(q, a) @ fiber(r),
                                       atol=3e-14)
            bracket = lift_bracket(q, a, b)
            expected_base = np.dot(x, a) * b - np.dot(x, b) * a
            np.testing.assert_allclose(bracket[:, 3], expected_base, atol=3e-14)
            u, v = f.T @ a, f.T @ b
            np.testing.assert_allclose((q.T @ bracket)[:3, :3],
                                       -(np.outer(u, v) - np.outer(v, u)), atol=3e-14)

    def test_12_bracket_independent_curve_directional_difference(self):
        q = random_so(self.rng, 4)
        a, b = self.rng.normal(size=(2, 4))
        eps = 1e-6
        ua, ub = q[:, :3].T @ a, q[:, :3].T @ b
        db = (lift(endpoint(q, [eps * ua]), b)
              - lift(endpoint(q, [-eps * ua]), b)) / (2 * eps)
        da = (lift(endpoint(q, [eps * ub]), a)
              - lift(endpoint(q, [-eps * ub]), a)) / (2 * eps)
        np.testing.assert_allclose(db - da, lift_bracket(q, a, b), atol=2e-9)

    def test_13_full_configuration_finite_synthesis(self):
        targets = [I4, horizontal_rotation([np.pi, 0, 0]),
                   *[random_so(self.rng, 4) for _ in range(15)]]
        for target in targets:
            recipe = synthesize(target)
            np.testing.assert_allclose(endpoint(I4, recipe), target, atol=4e-13)
            self.assertLessEqual(len(recipe), 10)
            self.assertLessEqual(cost(recipe), 7 * np.pi + 1e-13)

    def test_14_same_position_protocol_from_arbitrary_initial_frames(self):
        for _ in range(8):
            q, r = random_so(self.rng, 4), random_so(self.rng, 3)
            protocol = triangle(0, 2, 0.7)
            end1 = endpoint(q, protocol)
            end2 = endpoint(q @ fiber(r), [r.T @ u for u in protocol])
            np.testing.assert_allclose(position(end1), position(q), atol=3e-14)
            np.testing.assert_allclose(position(end2), position(q), atol=3e-14)
            np.testing.assert_allclose(end2, end1 @ fiber(r), atol=3e-14)


if __name__ == "__main__":
    main(__name__, "frame_position_quotient_audit", report)
