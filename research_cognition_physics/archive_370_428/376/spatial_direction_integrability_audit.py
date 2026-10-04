"""Round 376: local three-direction rank versus controlled spatial reachability.

R4, its measure, smooth vector fields, reversible pulse access and position
readout are inputs. Pullback unitaries preserve arbitrary unknown references;
they are not derived from the formal Pauli symbol or an uncontrolled natural H.
"""
from functools import lru_cache
import unittest
import numpy as np
from growing_stream_audit import main


PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]],
                  [[1, 0], [0, -1]]], complex)


def fields(point, lam=0.):
    x, y, _, _ = np.asarray(point)
    return np.array([[1., 0., 0., lam*y],
                     [0., 1., 0., x], [0., 0., 1., 0.]]).T


def theta(point, lam=0.):
    x, y, _, _ = np.asarray(point)
    return np.array([-lam*y, -x, 0., 1.])


def bracket_by_differences(point, first=0, second=1, lam=0., step=1e-5):
    jac = [np.column_stack([(fields(point+step*v, lam)[:, a] -
                            fields(point-step*v, lam)[:, a])/(2*step)
                           for v in np.eye(4)]) for a in (first, second)]
    fs = fields(point, lam)
    return jac[1]@fs[:, first]-jac[0]@fs[:, second]


def flow(points, axis, duration, lam=0.):
    result = np.array(points, dtype=float, copy=True)
    if axis == 0:
        result[..., 3] += lam*duration*result[..., 1]
        result[..., 0] += duration
    elif axis == 1:
        result[..., 3] += duration*result[..., 0]
        result[..., 1] += duration
    elif axis == 2:
        result[..., 2] += duration
    else:
        raise ValueError("Only the three supplied channels are allowed.")
    return result


def flow_matrix(axis, duration, lam=0.):
    matrix = np.eye(5)
    matrix[axis, 4] = duration
    if axis == 0:
        matrix[3, 1] = lam*duration
    elif axis == 1:
        matrix[3, 0] = duration
    return matrix


def protocol(points, pulses, lam=0.):
    result = np.asarray(points, float)
    for axis, duration in pulses:
        result = flow(result, axis, duration, lam)
    return result


def protocol_matrix(pulses, lam=0.):
    matrix = np.eye(5)
    for axis, duration in pulses:
        matrix = flow_matrix(axis, duration, lam)@matrix
    return matrix


def reverse_protocol(pulses):
    return [(a, -s) for a, s in reversed(pulses)]


def loop(a, b):
    return [(0, a), (1, b), (0, -a), (1, -b)]


def cost(pulses):
    return float(sum(abs(s) for _, s in pulses))


def target_protocol(target, lam=0.):
    x, y, z, w = np.asarray(target)
    eta = 1-lam
    delta = w-x*y
    pulses = [(0, x), (1, y), (2, z)]
    if abs(delta) < 1e-15:
        return pulses
    if eta == 0:
        raise ValueError("The target is outside the integrable leaf w=x*y.")
    product = delta/eta
    a = np.sqrt(abs(product))
    return pulses+loop(a, np.sign(product)*a)


def rhs(point, control, lam=0.):
    return fields(point, lam)@np.asarray(control)


def rk4(initial, control, end, steps, lam=0.):
    q = np.asarray(initial, float).copy()
    dt = end/steps
    for n in range(steps):
        t = n*dt
        k1 = rhs(q, control(t), lam)
        k2 = rhs(q+dt*k1/2, control(t+dt/2), lam)
        k3 = rhs(q+dt*k2/2, control(t+dt/2), lam)
        k4 = rhs(q+dt*k3, control(t+dt), lam)
        q += dt*(k1+2*k2+2*k3+k4)/6
    return q


def integrated_protocol(initial, pulses, lam=0.):
    q = np.asarray(initial, float).copy()
    for axis, parameter in pulses:
        if parameter:
            direction = np.sign(parameter)*np.eye(3)[axis]
            q = rk4(q, lambda t: direction, abs(parameter), 5, lam)
    return q


def smooth_control(t):
    return np.array([np.cos(t), np.sin(2*t), .2+np.cos(3*t)])


def smooth_endpoint(t, lam=0.):
    x = np.sin(t)
    y = (1-np.cos(2*t))/2
    z = .2*t+np.sin(3*t)/3
    original_w = np.sin(t)/2-np.sin(3*t)/6
    return np.array([x, y, z, lam*x*y+(1-lam)*original_w])


def endpoint_map(parameters, lam=0., excursion=.2):
    x, y, z, b = parameters
    return protocol(np.zeros(4), [(0, x), (1, y), (2, z)]+loop(excursion, b), lam)


def endpoint_jacobian(parameters, lam=0., step=1e-5):
    return np.column_stack([(endpoint_map(parameters+step*v, lam) -
                             endpoint_map(parameters-step*v, lam))/(2*step)
                            for v in np.eye(4)])


def dilation(points, scale):
    return np.asarray(points)*np.array([scale, scale, scale, scale*scale])


def symbol(point, covector, lam=0.):
    return np.einsum("a,aij->ij", fields(point, lam).T@covector, PAULI)


def quantum_coefficients():
    # Six amplitudes are payload(2) times inactive reference(3).
    rng = np.random.default_rng(376)
    raw = rng.normal(size=(6, 2))+1j*rng.normal(size=(6, 2))
    orthogonal = np.linalg.qr(raw)[0]
    return orthogonal[:, 0].reshape(2, 3)/np.sqrt(2), \
        orthogonal[:, 1].reshape(2, 3)/np.sqrt(2)


C0, C1 = quantum_coefficients()


def wavefunction(points):
    q = np.asarray(points)
    envelope = np.exp(-np.sum(q*q, axis=-1)/2)/np.pi
    return envelope[..., None, None]*(C0+np.sqrt(2)*q[..., 0, None, None]*C1)


@lru_cache(maxsize=None)
def gaussian_quadrature(order):
    nodes, weights = np.polynomial.hermite.hermgauss(order)
    points = np.stack(np.meshgrid(*([nodes]*4), indexing="ij"), axis=-1).reshape(-1, 4)
    combined = np.prod(np.stack(np.meshgrid(*([weights]*4), indexing="ij"), axis=-1),
                       axis=-1).reshape(-1)/np.pi**2
    return points, combined


def quantum_statistics(pulses, lam=0., order=16):
    points, weights = gaussian_quadrature(order)
    # Divide the pulled-back amplitude by the reference Gaussian; the weights
    # integrate pi^-2 exp(-|q|^2) dq exactly for polynomial input moments.
    preimage = protocol(points, reverse_protocol(pulses), lam)
    ratio = np.exp((np.sum(points*points, axis=1)-np.sum(preimage*preimage, axis=1))/2)
    amplitudes = ratio[:, None, None]*(C0+np.sqrt(2)*preimage[:, 0, None, None]*C1)
    weighted = amplitudes.reshape(len(points), 6)*np.sqrt(weights)[:, None]
    payload_reference = weighted.T@weighted.conjugate()
    probabilities = np.sum(abs(amplitudes)**2, axis=(1, 2))
    mean = np.sum(weights[:, None]*probabilities[:, None]*points, axis=0)
    return payload_reference, mean


def transport_case():
    pulses = loop(.3, -.4)
    initial, initial_mean = quantum_statistics([], order=18)
    final, final_mean = quantum_statistics(pulses, order=18)
    coarse, _ = quantum_statistics(pulses, order=12)
    rng = np.random.default_rng(1376)
    points = rng.normal(size=(20, 4))
    seq = wavefunction(protocol(points, reverse_protocol(pulses)))
    direct = wavefunction(points-np.array([0., 0., 0., -.12]))
    return {
        "loop_parameters": [.3, -.4], "positive_elapsed_time_at_unit_axis_rate": cost(pulses),
        "net_spatial_translation": [0., 0., 0., -.12],
        "initial_norm": float(np.trace(initial).real), "final_norm": float(np.trace(final).real),
        "payload_reference_reduced_state_error": float(np.linalg.norm(final-initial, 2)),
        "centroid_displacement": (final_mean-initial_mean).tolist(),
        "quadrature_orders": [12, 18],
        "quadrature_reduced_state_difference": float(np.linalg.norm(final-coarse, 2)),
        "complete_wavefunction_sequence_translation_error": float(np.max(abs(seq-direct))),
        "scope": "The proof covers all L2 spatial states with arbitrary payload/reference; "
                 "quadrature only checks one spatially entangled six-amplitude example."
    }


def report():
    point = np.array([.3, -.2, .1, .4])
    parameters = np.array([.12, -.17, .21, .11])
    differential = []
    for lam in (0., .8, 1.):
        bracket = bracket_by_differences(point, lam=lam)
        jac = endpoint_jacobian(parameters, lam)
        differential.append({
            "lambda": lam, "instantaneous_rank": int(np.linalg.matrix_rank(fields(point, lam))),
            "bracket_X1_X2": bracket.tolist(),
            "rank_with_bracket": int(np.linalg.matrix_rank(
                np.column_stack((fields(point, lam), bracket)))),
            "finite_protocol_endpoint_rank": int(np.linalg.matrix_rank(jac, tol=1e-8)),
            "finite_protocol_endpoint_determinant": float(np.linalg.det(jac))
        })
    loops = []
    for lam in (0., .99, 1.):
        for a, b in ((.2, .3), (.1, .15), (.05, .075)):
            pulses = loop(a, b)
            exact = protocol(point, pulses, lam)
            numerical = integrated_protocol(point, pulses, lam)
            loops.append({
                "lambda": lam, "a": a, "b": b, "axis_time_cost": cost(pulses),
                "w_displacement": float(exact[3]-point[3]),
                "predicted_w_displacement": (1-lam)*a*b,
                "independent_RK_endpoint_error": float(np.linalg.norm(numerical-exact))
            })
    rk_errors = [float(np.linalg.norm(rk4(np.zeros(4), smooth_control, 1.3, n, .3) -
                                       smooth_endpoint(1.3, .3))) for n in (8, 16, 32)]
    return {
        "round": 376,
        "scope": {
            "baseline": "Frozen through 375; no dependence on other new rounds.",
            "claim": "Rank-three local transport directions do not force rank-three "
                     "controlled reachability. Involutivity gives three-dimensional "
                     "leaves; an explicit nonintegrable family gives exact fourth-direction loops.",
            "inputs": ["smooth R4 position manifold and Lebesgue measure",
                       "X1=partial_x+lambda*y*partial_w, X2=partial_y+x*partial_w, X3=partial_z",
                       "positive-duration access to both signs of three supplied channels",
                       "signed pulse amplitudes, scheduling, position readout and control cost",
                       "optional L2 pullback implementation with an unchanged payload"],
            "not_proved": ["spatial dimension is or is not three from cognition",
                           "FUCP selects a manifold, distribution or control fields",
                           "the formal Pauli symbol supplies the independent pulse controls",
                           "a single natural Hamiltonian generates the controlled orbit",
                           "finite-dimensional bounded-energy quantum realization",
                           "Lorentz propagation, metric dynamics or gravity"],
            "new_interface": "Cross-position integrability and full finite control "
                             "orbits, beyond graph degree, local Bloch rank and a prescribed R3 Weyl model."
        },
        "local_and_endpoint_ranks": differential,
        "finite_loops": loops,
        "smooth_ODE_refinement": {"steps": [8, 16, 32], "errors": rk_errors,
                                  "successive_error_ratios": [rk_errors[i]/rk_errors[i+1]
                                                              for i in range(2)]},
        "near_integrable_resource": {
            "fixed_w_resolution": 1e-4,
            "eta_values": [1., .01, .0001],
            "minimum_axis_time_for_that_vertical_displacement":
                [4*np.sqrt(1e-4/eta) for eta in (1., .01, .0001)],
            "formula": "closed xy loop: |delta_w| <= |1-lambda| L^2/16, attained by a square"
        },
        "dimensions_and_measure": {
            "instantaneous_direction_rank": 3,
            "lambda_zero_control_orbit_topological_dimension": 4,
            "lambda_zero_control_ball_volume_exponent": 5,
            "integrable_leaf_dimension": 3,
            "volume_bounds_lambda_zero": "R^5/6250 <= Lebesgue_volume(B_L1(R)) <= 16 R^5",
            "volume_scope": "Specified L1 control length and Lebesgue volume; not a universal physical or spectral dimension."
        },
        "unknown_quantum_payload": transport_case(),
        "sources": ["https://doi.org/10.1090/S0002-9947-1973-0321133-2",
                    "https://scispace.com/pdf/orbits-of-families-of-vector-fields-and-integrability-of-1eeq49bauq.pdf"]
    }


class Checks(unittest.TestCase):
    def test_01_rank_three_symbol_has_one_invisible_covector(self):
        for point in ([0., 0., 0., 0.], [.4, -.2, .7, .1]):
            for lam in (0., .4, 1.):
                f = fields(point, lam)
                self.assertEqual(np.linalg.matrix_rank(f), 3)
                np.testing.assert_allclose(f.T@theta(point, lam), 0, atol=1e-15)
                np.testing.assert_allclose(symbol(point, theta(point, lam), lam), 0, atol=1e-15)
                covector = np.array([.2, -.3, .4, .5])
                h = symbol(point, covector, lam)
                np.testing.assert_allclose(h@h,
                    np.linalg.norm(f.T@covector)**2*np.eye(2), atol=2e-15)

    def test_02_bracket_and_failure_of_involutivity(self):
        point = np.array([.4, -.2, .7, .1])
        for lam in (0., .4, 1.):
            bracket = bracket_by_differences(point, lam=lam)
            np.testing.assert_allclose(bracket, [0, 0, 0, 1-lam], atol=3e-12)
            f = fields(point, lam)
            residual = bracket-f@np.linalg.lstsq(f, bracket, rcond=None)[0]
            self.assertAlmostEqual(np.linalg.norm(residual),
                abs(1-lam)/np.linalg.norm(theta(point, lam)), places=11)

    def test_03_actual_flows_are_invertible_and_preserve_volume(self):
        points = np.random.default_rng(3).normal(size=(17, 4))
        for lam in (0., .4, 1.):
            for axis in range(3):
                np.testing.assert_allclose(flow(flow(points, axis, .3, lam), axis, -.3, lam),
                                           points, atol=3e-16)
                matrix = flow_matrix(axis, .3, lam)
                self.assertAlmostEqual(np.linalg.det(matrix[:4, :4]), 1.)
                expected = np.column_stack((points, np.ones(len(points))))@matrix.T
                np.testing.assert_allclose(flow(points, axis, .3, lam), expected[:, :4], atol=2e-16)

    def test_04_finite_loop_against_independent_ODE_and_exact_translation(self):
        for lam in (0., .8, 1.):
            for point in ([0., 0., 0., 0.], [.3, -.7, .1, .5]):
                for a, b in ((.3, .4), (-.2, .5), (.3, -.4)):
                    expected = np.array(point)+[0, 0, 0, (1-lam)*a*b]
                    pulses = loop(a, b)
                    np.testing.assert_allclose(protocol(point, pulses, lam), expected, atol=3e-16)
                    np.testing.assert_allclose(integrated_protocol(point, pulses, lam), expected, atol=8e-16)

    def test_05_closed_loop_area_cost_bound_and_optimal_square(self):
        rng = np.random.default_rng(5)
        for lam in (0., .8):
            for _ in range(30):
                pulses = [(int(rng.integers(2)), float(rng.uniform(-.3, .3))) for _ in range(8)]
                interim = protocol(np.zeros(4), pulses, lam)
                pulses += [(0, -interim[0]), (1, -interim[1])]
                end = protocol(np.zeros(4), pulses, lam)
                np.testing.assert_allclose(end[:3], 0, atol=1e-15)
                self.assertLessEqual(abs(end[3]), abs(1-lam)*cost(pulses)**2/16+1e-15)
            for side in (.1, .7):
                pulses = loop(side, side)
                self.assertAlmostEqual(protocol(np.zeros(4), pulses, lam)[3],
                                       (1-lam)*cost(pulses)**2/16, places=14)

    def test_06_seven_supplied_pulses_reach_arbitrary_R4_targets(self):
        targets = np.random.default_rng(6).uniform(-.4, .4, size=(20, 4))
        for lam in (0., .5, 1.2):
            for target in targets:
                pulses = target_protocol(target, lam)
                self.assertLessEqual(len(pulses), 7)
                np.testing.assert_allclose(integrated_protocol(np.zeros(4), pulses, lam),
                                           target, atol=2e-15)

    def test_07_endpoint_map_rank_distinguishes_three_and_four(self):
        parameters = np.array([.1, -.2, .3, .4])
        for lam in (0., .9, 1.):
            jac = endpoint_jacobian(parameters, lam)
            self.assertEqual(np.linalg.matrix_rank(jac, tol=1e-8), 3 if lam == 1 else 4)
            self.assertAlmostEqual(np.linalg.det(jac), .2*(1-lam), places=11)

    def test_08_nonconstant_integrable_comparison_preserves_leaf(self):
        rng = np.random.default_rng(8)
        for _ in range(20):
            point = rng.normal(size=4)
            pulses = [(int(rng.integers(3)), float(rng.uniform(-.4, .4))) for _ in range(12)]
            end = integrated_protocol(point, pulses, 1.)
            self.assertAlmostEqual(end[3]-end[0]*end[1], point[3]-point[0]*point[1], places=13)
        with self.assertRaises(ValueError):
            target_protocol([0., 0., 0., .1], 1.)

    def test_09_independent_smooth_ODE_refinement(self):
        errors = [np.linalg.norm(rk4(np.zeros(4), smooth_control, 1.3, n, .3) -
                                  smooth_endpoint(1.3, .3)) for n in (8, 16, 32)]
        for i in range(2):
            self.assertGreater(errors[i]/errors[i+1], 14)
            self.assertLess(errors[i]/errors[i+1], 18)

    def test_10_control_dilation_and_actual_inner_box_paths(self):
        rng = np.random.default_rng(10)
        for scale in (.1, .4, 2.):
            pulses = [(0, .2), (1, -.3), (2, .1)]+loop(.2, .15)
            scaled = [(axis, scale*s) for axis, s in pulses]
            np.testing.assert_allclose(protocol(np.zeros(4), scaled),
                dilation(protocol(np.zeros(4), pulses), scale), atol=2e-16)
            self.assertAlmostEqual(np.linalg.det(np.diag([scale]*3+[scale**2])), scale**5)
            for _ in range(12):
                target = rng.uniform(-1, 1, size=4)*np.array([scale/10]*3+[scale**2/100])
                route = target_protocol(target)
                self.assertLessEqual(cost(route), scale)
                np.testing.assert_allclose(protocol(np.zeros(4), route), target, atol=5e-16)

    def test_11_all_bounded_control_paths_obey_outer_box(self):
        rng = np.random.default_rng(11)
        for _ in range(40):
            pulses = [(int(rng.integers(3)), float(rng.uniform(-.2, .2))) for _ in range(10)]
            length = cost(pulses)
            end = integrated_protocol(np.zeros(4), pulses)
            self.assertTrue(np.all(abs(end[:3]) <= length+1e-15))
            self.assertLessEqual(abs(end[3]), length**2+1e-15)

    def test_12_near_involutive_budget_does_not_change_exact_orbit_rank(self):
        displacement = 1e-4
        for eta in (1., .01, .0001):
            side = np.sqrt(displacement/eta)
            pulses = loop(side, side)
            self.assertAlmostEqual(cost(pulses), 4*np.sqrt(displacement/eta))
            self.assertAlmostEqual(protocol(np.zeros(4), pulses, 1-eta)[3], displacement, places=13)
            self.assertEqual(np.linalg.matrix_rank(np.column_stack((
                fields(np.zeros(4), 1-eta), [0, 0, 0, eta]))), 4)

    def test_13_quantum_norm_and_unknown_reference_marginal_under_real_flows(self):
        original, _ = quantum_statistics([], order=18)
        self.assertAlmostEqual(np.trace(original).real, 1., places=13)
        for pulses in ([(0, .25)], [(1, -.3)], loop(.3, -.4)):
            final, _ = quantum_statistics(pulses, order=18)
            self.assertLess(np.linalg.norm(final-original, 2), 2e-10)
            self.assertAlmostEqual(np.trace(final).real, 1., places=10)

    def test_14_entangled_wavefunction_loop_equals_translation_and_has_no_postselection(self):
        row = transport_case()
        self.assertLess(row["complete_wavefunction_sequence_translation_error"], 2e-16)
        np.testing.assert_allclose(row["centroid_displacement"], [0, 0, 0, -.12], atol=2e-10)
        self.assertLess(row["quadrature_reduced_state_difference"], 2e-10)
        points = np.random.default_rng(14).normal(size=(10, 4))
        pulses = [(0, .2), (1, -.3), (2, .1)]+loop(.2, .15)
        recovered = protocol(protocol(points, pulses), reverse_protocol(pulses))
        np.testing.assert_allclose(wavefunction(recovered), wavefunction(points), atol=3e-16)
        matrix = protocol_matrix(pulses)
        homogeneous = np.column_stack((points, np.ones(len(points))))
        np.testing.assert_allclose(protocol(points, pulses), (homogeneous@matrix.T)[:, :4], atol=3e-16)


if __name__ == "__main__":
    main(__name__, "spatial_direction_integrability_audit", report)
