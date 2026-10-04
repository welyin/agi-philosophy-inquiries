"""Round 385: operational scaling, a spherical budget boundary and reversal.

The structure theorem is analytic and uses an exact extra displacement contract.
Heisenberg costs below are shortest horizontal path costs, not an arbitrary gauge.
A separate Jordan-flow case audits the given scaling versus a chosen grading.
"""
import unittest

import numpy as np

from growing_stream_audit import main


def multiply(g, h):
    x, y, z = np.asarray(g, float)
    a, b, c = np.asarray(h, float)
    return np.array([x + a, y + b, z + c + (x * b - y * a) / 2])


def matrix(g):
    x, y, z = g
    return np.array([[1., x, z + x * y / 2], [0., 1., y], [0., 0., 1.]])


def dilation(r, g):
    return np.array([r * g[0], r * g[1], r * r * g[2]])


def theta_minus_sin(t):
    if abs(t) < 0.001:
        return t ** 3 / 6 - t ** 5 / 120 + t ** 7 / 5040
    return t - np.sin(t)


def area_ratio(t):
    if abs(t) < 0.001:
        return t / 12 + t ** 3 / 360 + t ** 5 / 10080
    return theta_minus_sin(t) / (8 * np.sin(t / 2) ** 2)


def cc_parameters(g):
    """Return exact-formula CC length and signed turning angle, numerically.

    Our group has z += (x*y'-y*x')/2. The cited paper uses t=-4z.
    Monotonic inversion of |z|/rho^2=(theta-sin theta)/(4(1-cos theta))
    uses a bracket in (0,2*pi), with exact axis cases treated separately.
    """
    x, y, z = np.asarray(g, float)
    rho = np.hypot(x, y)
    if rho == 0:
        return np.sqrt(4 * np.pi * abs(z)), float(np.sign(z) * 2 * np.pi)
    if z == 0:
        return rho, 0.
    target = abs(z) / rho ** 2
    lo, hi = 0., 2 * np.pi
    for _ in range(65):
        mid = (lo + hi) / 2
        if area_ratio(mid) < target:
            lo = mid
        else:
            hi = mid
    theta = np.sign(z) * (lo + hi) / 2
    length = rho / np.sinc(theta / (2 * np.pi))
    return float(length), float(theta)


def cc_distance(g):
    return cc_parameters(g)[0]


def directional_cost(g, beta=0.3):
    if abs(beta) >= 1:
        raise ValueError("The horizontal directional cost must stay positive.")
    return cc_distance(g) + beta * g[0]


def two_way_budget(g, beta=0.3):
    return cc_distance(g) + abs(beta * g[0])


def arc_point(g, t):
    length, theta = cc_parameters(g)
    if length == 0:
        return np.zeros(3)
    if theta == 0:
        return t * np.asarray(g)
    alpha = (np.arctan2(g[1], g[0]) - theta / 2
             if np.hypot(g[0], g[1]) != 0 else 0.)
    xy = length * t * np.sinc(theta * t / (2 * np.pi)) * np.exp(1j * (alpha + theta * t / 2))
    z = length ** 2 * theta_minus_sin(theta * t) / (2 * theta ** 2)
    return np.array([xy.real, xy.imag, z])


def finite_pulses(g, count=160, beta=0.3):
    length, theta = cc_parameters(g)
    alpha = (np.arctan2(g[1], g[0]) - theta / 2
             if np.hypot(g[0], g[1]) != 0 else 0.)
    endpoint, cost = np.zeros(3), 0.
    for j in range(count):
        angle = alpha + theta * (j + 0.5) / count
        step = np.array([length * np.cos(angle) / count,
                         length * np.sin(angle) / count, 0.])
        endpoint = multiply(endpoint, step)
        cost += np.linalg.norm(step[:2]) + beta * step[0]
    return endpoint, float(cost)


def sphere_to_budget(v, beta=0.3):
    return dilation(1 / two_way_budget(v, beta), v)


def budget_to_sphere(g):
    xy2 = np.dot(g[:2], g[:2])
    scale2 = 2 / (xy2 + np.sqrt(xy2 ** 2 + 4 * g[2] ** 2))
    return dilation(np.sqrt(scale2), g)


def radial_redirection(g, rotation, beta=0.):
    if np.linalg.norm(g) == 0:
        return np.zeros(3)
    budget = two_way_budget(g, beta)
    direction = budget_to_sphere(g)
    return dilation(budget, sphere_to_budget(rotation @ direction, beta))


def heisenberg_automorphism(g, horizontal, central_shear):
    return np.r_[horizontal @ g[:2],
                 np.linalg.det(horizontal) * g[2] + np.dot(central_shear, g[:2])]


D = np.array([[1., 4.], [0., 1.]])


def lyapunov_matrix(d=D):
    n = len(d)
    basis = np.eye(n * n).reshape(n * n, n, n)
    linear = np.column_stack([(d.T @ b + b @ d).ravel() for b in basis])
    return np.linalg.solve(linear, np.eye(n).ravel()).reshape(n, n)


def jordan_flow(t, v):
    return np.exp(t) * np.array([v[0] + 4 * t * v[1], v[1]])


def jordan_radius(v):
    if np.linalg.norm(v) == 0:
        return 0.
    p = lyapunov_matrix()
    lo, hi = -50., 50.
    for _ in range(80):
        mid = (lo + hi) / 2
        w = jordan_flow(-mid, v)
        if w @ p @ w > 1:
            lo = mid
        else:
            hi = mid
    return float(np.exp((lo + hi) / 2))


def second_euclidean_crossing():
    v = np.array([1., -1.]) / np.sqrt(2)
    lo, hi = 0.25, 0.4
    for _ in range(70):
        mid = (lo + hi) / 2
        if np.dot(jordan_flow(mid, v), jordan_flow(mid, v)) < 1:
            lo = mid
        else:
            hi = mid
    return float((lo + hi) / 2)


def report():
    endpoint = np.array([0.8, -0.3, 0.2])
    exact_cost = directional_cost(endpoint)
    pulse_rows = []
    for count in (160, 320):
        q, cost = finite_pulses(endpoint, count)
        pulse_rows.append({"pulses": count,
                           "endpoint_error": float(np.linalg.norm(q - endpoint)),
                           "cost": cost,
                           "cost_minus_optimal_to_actual_endpoint": cost - directional_cost(q)})
    rotation = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])
    g, h = np.array([0., 0., 1.]), np.array([0., 0.6, 0.])
    lhs = radial_redirection(multiply(g, h), rotation)
    rhs = multiply(radial_redirection(g, rotation), radial_redirection(h, rotation))
    return {
        "round": 385,
        "scope": {
            "conditional_theorem": "A Hausdorff topological displacement group with a proper continuous positive homogeneous budget and a jointly continuous positive-real automorphism scaling has finite-dimensional simply connected nilpotent Lie structure. The max of forward and inverse budgets has a complete spherical unit boundary with inverse as a free involution.",
            "redundancies_removed": "Proper budget implies local compactness; homogeneity and properness imply contraction, hence the scaling paths imply connectedness. These are not separate axioms.",
            "extra_central_result": "If actual budget-preserving redirections are also displacement-group automorphisms and every unit-boundary orbit is dense, the nilpotent group is abelian. With the separately required qubit interfaces of rounds 383-384 this gives an R^3 displacement group conditionally.",
            "inputs": ["Positions are identified homeomorphically with a free transitive displacement torsor", "Exact group composition and common replay at each basepoint", "Exact positive-real scaling by group automorphisms", "Continuous proper budget with its declared operational meaning", "Direction/readout contracts are still separate from these assumptions"],
            "not_claimed": ["Derivation of those inputs from FUCP", "A finite sample proof of a structure theorem", "Every positive homogeneous budget is a metric or travel time", "Finite displacements in a curved spacetime obey this exact homogeneous contract", "A derivation of GR or a unique natural Hamiltonian"],
        },
        "actual_horizontal_protocol": {
            "model": "Heisenberg law with central increment (x*y'-y*x')/2",
            "running_cost": "sqrt(u_x^2+u_y^2)+beta*u_x, beta=0.3",
            "endpoint": endpoint.tolist(),
            "cc_distance": cc_distance(endpoint),
            "optimal_forward_cost": exact_cost,
            "optimal_inverse_cost": directional_cost(-endpoint),
            "two_way_max_budget": two_way_budget(endpoint),
            "pulse_crosschecks": pulse_rows,
            "vertical_unit_budget_coordinate": float(1 / (4 * np.pi)),
            "sphere_topology_does_not_determine_metric_volume_exponent": 4,
        },
        "given_scaling_need_not_be_a_pure_grading": {
            "generator": D.tolist(), "lyapunov_P": lyapunov_matrix().tolist(),
            "lyapunov_equation_residual": float(np.linalg.norm(D.T @ lyapunov_matrix() + lyapunov_matrix() @ D - np.eye(2))),
            "euclidean_squared_norm_initial_derivative": -2.,
            "second_euclidean_unit_sphere_crossing_t": second_euclidean_crossing(),
            "scope": "A separate additive R^2 Jordan-flow calibration, not the CC scaling of the Heisenberg protocol.",
        },
        "budget_preservation_is_not_composition_preservation": {
            "radial_SO3_redirection_product_defect": float(np.linalg.norm(lhs - rhs)),
            "moves_central_direction_to_noncentral_direction": True,
            "scope": "The exact radial conjugacy is a homeomorphism action preserving budgets. Calling it a displacement automorphism would add a false property.",
        },
        "circle_control_counterexample": {
            "unscaled_0_and_2pi_have_same_endpoint": True,
            "halved_control_endpoint_chord_difference": 2.,
            "scope": "Small continuous pulse amplitudes do not by themselves define a well-defined scaling on endpoint equivalence classes.",
        },
        "primary_sources": ["https://arxiv.org/pdf/0909.4565",
                            "https://arxiv.org/pdf/1512.04936",
                            "https://arxiv.org/pdf/1412.1797"],
    }


class Checks(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(385)

    def test_01_endpoint_group_law_has_an_independent_faithful_matrix_representation(self):
        for g, h in self.rng.normal(size=(30, 2, 3)):
            np.testing.assert_allclose(matrix(multiply(g, h)), matrix(g) @ matrix(h), atol=2e-15)
            np.testing.assert_allclose(matrix(-g), np.linalg.inv(matrix(g)), atol=2e-15)

    def test_02_path_return_and_inverse_displacement_from_same_basepoint_differ(self):
        p, g = np.array([0.2, 0.3, 0.4]), np.array([0.6, -0.2, 0.1])
        q = multiply(p, g)
        np.testing.assert_allclose(multiply(q, -g), p)
        opposite = multiply(p, -g)
        self.assertGreater(np.linalg.norm(opposite - p), 0.5)
        np.testing.assert_allclose(multiply(p, -multiply(-p, opposite)), q)

    def test_03_anisotropic_scaling_preserves_composition_but_naive_coordinate_scaling_fails(self):
        for g, h in self.rng.normal(size=(20, 2, 3)):
            np.testing.assert_allclose(dilation(0.4, multiply(g, h)),
                                       multiply(dilation(0.4, g), dilation(0.4, h)), atol=4e-16)
        g, h = np.array([1., 0., 0.]), np.array([0., 1., 0.])
        self.assertAlmostEqual(np.linalg.norm(0.4 * multiply(g, h) - multiply(0.4 * g, 0.4 * h)), 0.12)

    def test_04_cc_formula_and_arc_endpoints_match_including_axes(self):
        self.assertAlmostEqual(cc_distance([0., 0., 0.3]), np.sqrt(1.2 * np.pi))
        self.assertAlmostEqual(cc_distance([0.3, -0.4, 0.]), 0.5)
        for g in np.vstack((self.rng.normal(size=(30, 3)), [0, 0, 0.3], [0, 0, -0.2], [0.3, 0.4, 0])):
            np.testing.assert_allclose(arc_point(g, 1), g, atol=2e-14)
            np.testing.assert_allclose(arc_point(g, 0), 0)

    def test_05_actual_piecewise_controls_converge_and_pay_the_directional_cost(self):
        g = np.array([0.8, -0.3, 0.2])
        errors = []
        length = cc_distance(g)
        for count in (160, 320):
            q, cost = finite_pulses(g, count)
            errors.append(np.linalg.norm(q - g))
            self.assertAlmostEqual(cost - 0.3 * q[0], length)
            self.assertGreaterEqual(cost - directional_cost(q), -2e-14)
        self.assertGreater(errors[0] / errors[1], 3.99)
        self.assertLess(errors[1], 1e-5)

    def test_06_forward_and_inverse_costs_give_a_symmetric_homogeneous_subadditive_budget(self):
        for g, h in self.rng.normal(size=(25, 2, 3)):
            self.assertAlmostEqual(directional_cost(g) - directional_cost(-g), 0.6 * g[0])
            self.assertAlmostEqual(two_way_budget(g), max(directional_cost(g), directional_cost(-g)))
            self.assertAlmostEqual(two_way_budget(dilation(0.43, g)), 0.43 * two_way_budget(g))
            self.assertLessEqual(two_way_budget(multiply(g, h)), two_way_budget(g) + two_way_budget(h) + 1e-13)

    def test_07_complete_budget_boundary_is_reconstructed_bijectively_and_inverse_equivariantly(self):
        for v in self.rng.normal(size=(40, 3)):
            v /= np.linalg.norm(v)
            g = sphere_to_budget(v)
            self.assertAlmostEqual(two_way_budget(g), 1)
            np.testing.assert_allclose(budget_to_sphere(g), v, atol=2e-15)
            np.testing.assert_allclose(sphere_to_budget(-v), -g, atol=3e-16)
        self.assertAlmostEqual(sphere_to_budget([0., 0., 1.])[2], 1 / (4 * np.pi))

    def test_08_jordan_scaling_requires_an_adapted_lyapunov_norm(self):
        p = lyapunov_matrix()
        self.assertGreater(np.linalg.eigvalsh(p).min(), 0)
        np.testing.assert_allclose(D.T @ p + p @ D, np.eye(2))
        v = np.array([1., -1.]) / np.sqrt(2)
        self.assertAlmostEqual(v @ (D.T + D) @ v, -2)
        t = second_euclidean_crossing()
        self.assertGreater(t, 0.25)
        self.assertAlmostEqual(np.linalg.norm(jordan_flow(t, v)), 1)
        self.assertGreater(jordan_flow(t, v) @ p @ jordan_flow(t, v), v @ p @ v)

    def test_09_lyapunov_radial_coordinate_has_one_crossing_and_respects_the_given_flow(self):
        p = lyapunov_matrix()
        for v in self.rng.normal(size=(20, 2)):
            radius = jordan_radius(v)
            unit = jordan_flow(-np.log(radius), v)
            self.assertAlmostEqual(unit @ p @ unit, 1)
            for t in (-0.7, 0.4):
                self.assertAlmostEqual(jordan_radius(jordan_flow(t, v)), np.exp(t) * radius, places=12)

    def test_10_all_displayed_nonabelian_automorphisms_preserve_the_central_sector(self):
        for _ in range(20):
            a = np.eye(2) + 0.3 * self.rng.normal(size=(2, 2))
            shear = self.rng.normal(size=2)
            g, h = self.rng.normal(size=(2, 3))
            transform = lambda x: heisenberg_automorphism(x, a, shear)
            np.testing.assert_allclose(transform(multiply(g, h)), multiply(transform(g), transform(h)), atol=2e-15)
            np.testing.assert_allclose(transform([0., 0., 0.7])[:2], 0)

    def test_11_full_boundary_redirections_can_preserve_budgets_while_failing_group_composition(self):
        q = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])
        for g in self.rng.normal(size=(15, 3)):
            transformed = radial_redirection(g, q)
            self.assertAlmostEqual(two_way_budget(transformed, 0.), two_way_budget(g, 0.))
            np.testing.assert_allclose(radial_redirection(transformed, q.T), g, atol=2e-14)
        center, h = np.array([0., 0., 1.]), np.array([0., 0.6, 0.])
        defect = radial_redirection(multiply(center, h), q) - multiply(radial_redirection(center, q), radial_redirection(h, q))
        self.assertGreater(np.linalg.norm(defect), 0.5)

    def test_12_continuous_pulse_attenuation_need_not_descend_to_endpoint_classes(self):
        pulse_zero, pulse_loop = 0., 2 * np.pi
        self.assertLess(abs(np.exp(1j * pulse_zero) - np.exp(1j * pulse_loop)), 3e-16)
        self.assertAlmostEqual(abs(np.exp(0.5j * pulse_zero) - np.exp(0.5j * pulse_loop)), 2)


if __name__ == "__main__":
    main(__name__, "scalable_displacement_boundary_audit", report)
