"""Round 377: full horizontal isotropy does not identify all reachable directions."""
import unittest
import numpy as np
from growing_stream_audit import main


ZERO = np.zeros(6)
EYE = np.eye(3)


def product(g, h):
    g, h = np.asarray(g), np.asarray(h)
    return np.concatenate([g[:3]+h[:3],
                           g[3:]+h[3:]+np.cross(g[:3], h[:3])/2])


def inverse(g):
    return -np.asarray(g)


def flow(g, u, time=1.):
    return product(g, np.concatenate([time*np.asarray(u), np.zeros(3)]))


def horizontal(g, u):
    return np.concatenate([u, np.cross(np.asarray(g)[:3], u)/2])


def horizontal_frame(g):
    return np.column_stack([horizontal(g, e) for e in EYE])


def flow_jacobian(u, time=1.):
    out = np.eye(6)
    out[3:, :3] = np.column_stack([time*np.cross(e, u)/2 for e in EYE])
    return out


def endpoint(g, increments):
    out = np.array(g, dtype=float)
    for v in increments:
        out = flow(out, v)
    return out


def area_endpoint(g, increments):
    increments = np.asarray(increments, dtype=float)
    net = increments.sum(axis=0) if len(increments) else np.zeros(3)
    area = sum((np.cross(a, b)/2
                for j, a in enumerate(increments) for b in increments[j+1:]),
               start=np.zeros(3))
    return np.concatenate([g[:3]+net, g[3:]+np.cross(g[:3], net)/2+area])


def rk4_endpoint(g, increments, steps_per_pulse=7):
    # Integrate the vector field itself, without group multiplication.
    out = np.array(g, dtype=float)
    dt = 1./steps_per_pulse
    for u in increments:
        for _ in range(steps_per_pulse):
            k1 = horizontal(out, u)
            k2 = horizontal(out+dt*k1/2, u)
            k3 = horizontal(out+dt*k2/2, u)
            k4 = horizontal(out+dt*k3, u)
            out += dt*(k1+2*k2+2*k3+k4)/6
    return out


def path_cost(increments):
    return float(sum(np.linalg.norm(v) for v in increments))


def commutator_loop(a, b):
    return np.array([a, b, -np.asarray(a), -np.asarray(b)])


def synthesis(target):
    """A straight pulse and one four-pulse loop reach any (x,z) exactly."""
    target = np.asarray(target)
    x, z = target[:3], target[3:]
    nz = np.linalg.norm(z)
    if nz == 0:
        return np.array([x])
    n = z/nz
    seed = EYE[np.argmin(abs(n))]
    a0 = seed-np.dot(seed, n)*n
    a0 /= np.linalg.norm(a0)
    a = np.sqrt(nz)*a0
    b = np.sqrt(nz)*np.cross(n, a0)
    return np.vstack([x, commutator_loop(a, b)])


def orthogonal_action(g, rotation):
    """The central coordinates are pseudovectors under orientation reversal."""
    return np.concatenate([rotation @ g[:3],
                           np.linalg.det(rotation)*rotation @ g[3:]])


def between_points(q, source, target, rotation):
    return product(target, orthogonal_action(product(inverse(source), q), rotation))


def dilation(g, scale):
    return np.concatenate([scale*g[:3], scale**2*g[3:]])


def central_effect(g, covector=EYE[2]):
    return float((1+np.tanh(np.dot(covector, g[3:])))/2)


def numerical_jacobian(function, point, step=1e-6):
    basis = np.eye(len(point))
    return np.column_stack([(function(point+step*e)-function(point-step*e))/(2*step)
                            for e in basis])


def orthogonal(rng, sign=1):
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    q[:, 0] *= sign*np.linalg.det(q)
    return q


def bracket_from_vector_fields(g, a, b):
    da = numerical_jacobian(lambda p: horizontal(p, a), g)
    db = numerical_jacobian(lambda p: horizontal(p, b), g)
    return db @ horizontal(g, a)-da @ horizontal(g, b)


def report():
    g = np.array([.3, -.2, .4, .1, .5, -.3])
    a, b = EYE[0], EYE[1]
    positive = endpoint(ZERO, commutator_loop(a, b))
    negative = endpoint(ZERO, commutator_loop(b, a))
    target = np.array([.4, -.2, .6, .3, -.5, .8])
    pulses = synthesis(target)
    frame = horizontal_frame(g)
    brackets = np.column_stack([bracket_from_vector_fields(g, EYE[i], EYE[j])
                                 for i, j in [(1, 2), (2, 0), (0, 1)]])
    rng = np.random.default_rng(377)
    reflection = orthogonal(rng, -1)
    sample = rng.normal(size=(8, 3))
    expected = orthogonal_action(endpoint(g, sample), reflection)
    actual = endpoint(orthogonal_action(g, reflection), sample @ reflection.T)
    scales = [.25, .5, 1., 2., 3.]
    return {
        "round": 377,
        "scope": "A specified smooth free step-two rank-three Carnot control model. "
                 "Full O(3) horizontal isotropy coexists with six-dimensional reachability "
                 "and homogeneous volume exponent nine. Not a derivation of physical space, "
                 "a counterexample to the Muller-Masanes theorem, or autonomous quantum dynamics.",
        "dimensions": {"horizontal_rank": int(np.linalg.matrix_rank(frame)),
                       "bracket_span_rank": int(np.linalg.matrix_rank(np.column_stack([frame, brackets]))),
                       "manifold_dimension": 6, "homogeneous_volume_exponent": 9},
        "opposite_loops": {"positive_endpoint": positive.tolist(),
                           "negative_endpoint": negative.tolist(),
                           "same_horizontal_endpoint": bool(np.array_equal(positive[:3], negative[:3])),
                           "pulse_cost_each": path_cost(commutator_loop(a, b)),
                           "central_effect_positive": central_effect(positive),
                           "central_effect_negative": central_effect(negative),
                           "probability_gap": central_effect(positive)-central_effect(negative)},
        "constructive_reachability": {
            "target": target.tolist(), "endpoint": endpoint(ZERO, pulses).tolist(),
            "endpoint_residual": float(np.linalg.norm(endpoint(ZERO, pulses)-target)),
            "pulses": len(pulses), "cost": path_cost(pulses),
            "analytic_cost_upper_bound": float(np.linalg.norm(target[:3])+4*np.sqrt(np.linalg.norm(target[3:])))},
        "isotropy": {"group": "O(3): (x,z) -> (R x, det(R) R z)",
                     "reflection_endpoint_residual": float(np.linalg.norm(expected-actual)),
                     "proper_rotations_on_qubit": "Inherited SU(2)/SO(3) interface from round 371; no new qubit test.",
                     "reflection_is_not_assumed_a_qubit_unitary": True},
        "dilation": [{"scale": r, "volume_factor": float(np.linalg.det(np.diag([r]*3+[r*r]*3))),
                      "expected_r_to_9": r**9} for r in scales],
        "coordinate_readout_input": "The full endpoint is an observable configuration. "
                                   "With x-only effects and the listed controls, z can instead be quotiented out.",
        "quantum_interface": "Each prescribed pulse is a unit-Jacobian diffeomorphism, hence "
                             "its pullback is unitary on L2(G) and extends by identity to any reference. "
                             "This is an added infinite-dimensional representation, not FUCP's finite reconstruction.",
        "sources": ["https://arxiv.org/pdf/1610.07359",
                    "https://cvgmt.sns.it/media/doc/paper/3191/cut-free-v2.pdf"]}


class Checks(unittest.TestCase):
    def test_01_group_associativity_identity_and_inverse(self):
        rng = np.random.default_rng(1377)
        for _ in range(20):
            g, h, k = rng.normal(size=(3, 6))
            np.testing.assert_allclose(product(product(g, h), k), product(g, product(h, k)), atol=3e-15)
            np.testing.assert_array_equal(product(g, ZERO), g)
            np.testing.assert_allclose(product(g, inverse(g)), ZERO, atol=1e-15)

    def test_02_exact_flow_and_generator(self):
        g = np.array([.2, -.4, .7, .1, -.5, .3])
        u = np.array([.3, -.8, .4])
        np.testing.assert_allclose(flow(flow(g, u, .7), u, -.2), flow(g, u, .5), atol=1e-15)
        derivative = (flow(g, u, 1e-6)-flow(g, u, -1e-6))/2e-6
        np.testing.assert_allclose(derivative, horizontal(g, u), atol=5e-11)

    def test_03_independent_brackets_and_growth_rank(self):
        g = np.array([.2, -.4, .7, .1, -.5, .3])
        brackets = []
        for i, j in [(1, 2), (2, 0), (0, 1)]:
            bracket = bracket_from_vector_fields(g, EYE[i], EYE[j])
            np.testing.assert_allclose(bracket, np.concatenate([np.zeros(3), np.cross(EYE[i], EYE[j])]), atol=1e-10)
            brackets.append(bracket)
        self.assertEqual(np.linalg.matrix_rank(horizontal_frame(g)), 3)
        self.assertEqual(np.linalg.matrix_rank(np.column_stack([horizontal_frame(g)]+brackets)), 6)
        # Horizontal coefficients do not depend on any central coordinate.
        for j in range(3):
            self.assertLess(np.linalg.norm(numerical_jacobian(lambda p: horizontal(p, EYE[j]), g)[:, 3:]), 1e-12)

    def test_04_finite_commutator_loop_at_arbitrary_start(self):
        rng = np.random.default_rng(2377)
        for _ in range(20):
            g, a, b = rng.normal(size=6), rng.normal(size=3), rng.normal(size=3)
            expected = g+np.concatenate([np.zeros(3), np.cross(a, b)])
            np.testing.assert_allclose(endpoint(g, commutator_loop(a, b)), expected, atol=4e-15)

    def test_05_constructive_reachability_not_only_rank(self):
        rng = np.random.default_rng(3377)
        targets = [ZERO, np.array([1., -2., 3., 0, 0, 0])] + list(rng.normal(size=(20, 6)))
        for target in targets:
            pulses = synthesis(target)
            np.testing.assert_allclose(endpoint(ZERO, pulses), target, atol=3e-15)
            bound = np.linalg.norm(target[:3])+4*np.sqrt(np.linalg.norm(target[3:]))
            self.assertLessEqual(path_cost(pulses), bound+5e-15)

    def test_06_O3_automorphisms_preserve_flow_and_cost(self):
        rng = np.random.default_rng(4377)
        for sign in [-1, 1]:
            r = orthogonal(rng, sign)
            g, h = rng.normal(size=(2, 6))
            pulses = rng.normal(size=(9, 3))
            np.testing.assert_allclose(orthogonal_action(product(g, h), r),
                                       product(orthogonal_action(g, r), orthogonal_action(h, r)), atol=4e-15)
            np.testing.assert_allclose(orthogonal_action(endpoint(g, pulses), r),
                                       endpoint(orthogonal_action(g, r), pulses @ r.T), atol=7e-15)
            self.assertAlmostEqual(path_cost(pulses), path_cost(pulses @ r.T), places=13)

    def test_07_every_horizontal_linear_isometry_extends_between_points(self):
        rng = np.random.default_rng(5377)
        source, target = rng.normal(size=(2, 6))
        r = orthogonal(rng, -1)
        f = lambda q: between_points(q, source, target, r)
        np.testing.assert_allclose(f(source), target, atol=1e-15)
        derivative = numerical_jacobian(f, source)
        np.testing.assert_allclose(derivative @ horizontal_frame(source),
                                   horizontal_frame(target) @ r, atol=4e-10)

    def test_08_area_sum_and_ODE_independently_match_group_flow(self):
        rng = np.random.default_rng(6377)
        for _ in range(8):
            g = rng.normal(size=6)
            pulses = rng.normal(size=(12, 3))
            exact = endpoint(g, pulses)
            np.testing.assert_allclose(area_endpoint(g, pulses), exact, atol=7e-15)
            np.testing.assert_allclose(rk4_endpoint(g, pulses), exact, atol=3e-14)

    def test_09_central_readout_and_x_only_quotient_are_different_contracts(self):
        plus = endpoint(ZERO, commutator_loop(EYE[0], EYE[1]))
        minus = endpoint(ZERO, commutator_loop(EYE[1], EYE[0]))
        self.assertAlmostEqual(central_effect(plus)-central_effect(minus), np.tanh(1.))
        future = np.array([[.2, -.1, .8], [-.4, .6, -.2]])
        np.testing.assert_array_equal(endpoint(plus, future)[:3], endpoint(minus, future)[:3])
        r = np.diag([-1., 1., 1.])
        covector = np.linalg.det(r)*r @ EYE[2]
        self.assertAlmostEqual(central_effect(orthogonal_action(plus, r), covector), central_effect(plus))

    def test_10_reflection_requires_pseudovector_center(self):
        r = np.diag([-1., 1., 1.])
        g, h = np.r_[EYE[0], np.zeros(3)], np.r_[EYE[1], np.zeros(3)]
        wrong = lambda q: np.r_[r @ q[:3], r @ q[3:]]
        self.assertGreater(np.linalg.norm(wrong(product(g, h))-product(wrong(g), wrong(h))), .9)
        np.testing.assert_array_equal(orthogonal_action(product(g, h), r),
                                      product(orthogonal_action(g, r), orthogonal_action(h, r)))

    def test_11_dilation_preserves_products_and_scales_all_pulse_costs(self):
        rng = np.random.default_rng(7377)
        g, h = rng.normal(size=(2, 6))
        pulses = rng.normal(size=(7, 3))
        for scale in [.25, .5, 2., 3.]:
            np.testing.assert_allclose(dilation(product(g, h), scale),
                                       product(dilation(g, scale), dilation(h, scale)), atol=5e-15)
            np.testing.assert_allclose(dilation(endpoint(g, pulses), scale),
                                       endpoint(dilation(g, scale), scale*pulses), atol=2e-14)
            self.assertAlmostEqual(path_cost(scale*pulses), scale*path_cost(pulses), places=13)

    def test_12_Haar_measure_and_homogeneous_volume_jacobians(self):
        point = np.array([.3, -.7, .9, -.2, .1, .8])
        anchor = np.array([.4, .2, -.5, .9, -.6, .1])
        left_jac = numerical_jacobian(lambda q: product(anchor, q), point)
        self.assertAlmostEqual(np.linalg.det(left_jac), 1., places=8)
        flow_jac = numerical_jacobian(lambda q: flow(q, [.2, -.5, .8], .7), point)
        np.testing.assert_allclose(flow_jac, flow_jacobian([.2, -.5, .8], .7), atol=8e-11)
        self.assertAlmostEqual(np.linalg.det(flow_jac), 1., places=8)
        for scale in [.5, 2.]:
            jac = numerical_jacobian(lambda q: dilation(q, scale), point)
            self.assertAlmostEqual(np.linalg.det(jac)/scale**9, 1., places=8)

    def test_13_control_ball_bounds_and_nonzero_cost_of_central_motion(self):
        rng = np.random.default_rng(8377)
        for _ in range(20):
            pulses = rng.normal(size=(9, 3))
            final = endpoint(ZERO, pulses)
            length = path_cost(pulses)
            self.assertLessEqual(np.linalg.norm(final[:3]), length+1e-14)
            self.assertLessEqual(np.linalg.norm(final[3:]), length**2/4+1e-14)
        target = np.array([0, 0, 0, 0, 0, 1.])
        self.assertAlmostEqual(path_cost(synthesis(target)), 4.)
        self.assertEqual(2*np.sqrt(np.linalg.norm(target[3:])), 2.)


if __name__ == "__main__":
    main(__name__, "isotropic_direction_orbit_audit", report)
