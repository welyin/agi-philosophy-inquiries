"""Round 379: position identity is not a complete dynamical state.

The smooth configuration, canonical kinetic Hamiltonian, x-only positioning
interface, launch preparations and measurement noise are explicit inputs.
This is a restricted classical audit, not a model of all FUCP operations.
"""
import math
import unittest
import numpy as np
from growing_stream_audit import main


def cross_matrix(v):
    x, y, z = np.asarray(v, float)
    return np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])


def velocity(state):
    x, _, p, lam = np.asarray(state).reshape(4, 3)
    return p + .5*np.cross(lam, x)


def make_state(x, u, lam, z=None):
    x, u, lam = map(lambda v: np.asarray(v, float), (x, u, lam))
    if z is None:
        z = np.zeros(3)
    return np.concatenate((x, z, u-.5*np.cross(lam, x), lam))


def hamiltonian(state):
    u = velocity(state)
    return float(u@u/2)


def canonical_rhs(state):
    x, _, _, lam = np.asarray(state).reshape(4, 3)
    u = velocity(state)
    return np.concatenate((u, .5*np.cross(x, u), .5*np.cross(lam, u), np.zeros(3)))


def finite_gradient(function, point, step=1e-5):
    point = np.asarray(point, float)
    return np.array([(function(point+step*v)-function(point-step*v))/(2*step)
                     for v in np.eye(len(point))])


def canonical_matrix():
    return np.block([[np.zeros((6, 6)), np.eye(6)],
                     [-np.eye(6), np.zeros((6, 6))]])


def rk4(rhs, initial, t, steps):
    y = np.asarray(initial, float).copy()
    dt = t/steps
    for _ in range(steps):
        k1 = rhs(y)
        k2 = rhs(y+dt*k1/2)
        k3 = rhs(y+dt*k2/2)
        k4 = rhs(y+dt*k3)
        y += dt*(k1+2*k2+2*k3+k4)/6
    return y


def coefficients(theta):
    sinc = float(np.sinc(theta/np.pi))
    cosc = .5*float(np.sinc(theta/(2*np.pi)))**2
    if abs(theta) < 1e-3:
        sincc = 1/6-theta**2/120+theta**4/5040-theta**6/362880
    else:
        sincc = (theta-np.sin(theta))/theta**3
    return sinc, cosc, sincc


def propagation_matrices(lam, t):
    ell = cross_matrix(lam)
    sinc, cosc, sincc = coefficients(np.linalg.norm(lam)*t)
    rotation = np.eye(3)+t*sinc*ell+t*t*cosc*(ell@ell)
    integral = t*np.eye(3)+t*t*cosc*ell+t**3*sincc*(ell@ell)
    return rotation, integral


def exact_position_velocity(x, u, lam, t):
    rotation, integral = propagation_matrices(lam, t)
    return np.asarray(x)+integral@u, rotation@u


def full_reference_flow(initial, t, order=40):
    """Exact x/u/lambda formula plus independent Gauss integration of z."""
    x0, z0, _, lam = np.asarray(initial).reshape(4, 3)
    u0 = velocity(initial)
    x, u = exact_position_velocity(x0, u0, lam, t)
    nodes, weights = np.polynomial.legendre.leggauss(order)
    area = np.zeros(3)
    for node, weight in zip(nodes, weights):
        xs, us = exact_position_velocity(x0, u0, lam, t*(node+1)/2)
        area += weight*np.cross(xs, us)*t/4
    return make_state(x, u, lam, z0+area)


def reduced_rhs(xp, lam):
    x, p = np.asarray(xp).reshape(2, 3)
    u = p+.5*np.cross(lam, x)
    return np.concatenate((u, .5*np.cross(lam, u)))


def preparation_bounds(t, dx, du, dlambda, speed_bound):
    t = abs(t)
    return {"velocity": du+t*speed_bound*dlambda,
            "position": dx+t*du+.5*t*t*speed_bound*dlambda}


def position_tv(displacement, sigma):
    """TV of two R3 Gaussians with covariance sigma^2 I and displaced means."""
    return math.erf(float(np.linalg.norm(displacement))/(2*np.sqrt(2)*sigma))


def position_tv_bound(position_error, sigma):
    return min(1., position_error/(np.sqrt(2*np.pi)*sigma))


def launch_matrix(directions):
    # Acceleration is lambda cross u = -[u]_cross lambda.
    return np.vstack([-cross_matrix(u) for u in directions])


def acceleration_readout(lam, u, h, noise):
    samples = np.array([exact_position_velocity(np.zeros(3), u, lam, j*h)[0]
                        for j in range(3)])+noise
    return (samples[2]-2*samples[1]+samples[0])/h**2


def calibration_case():
    lam = np.array([.4, -.3, .7])
    directions = np.eye(3)[:2]
    h, epsilon, cap = .04, 2e-6, 1.
    rng = np.random.default_rng(1379)
    raw = rng.normal(size=(2, 3, 3))
    noise = epsilon*raw/np.linalg.norm(raw, axis=2, keepdims=True)
    observations = np.concatenate([acceleration_readout(lam, u, h, n)
                                    for u, n in zip(directions, noise)])
    matrix = launch_matrix(directions)
    fitted = np.linalg.lstsq(matrix, observations, rcond=None)[0]
    acceleration_bound = h*cap**2+4*epsilon/h**2
    parameter_bound = np.sqrt(2)*acceleration_bound
    t = 1.5
    x0, u0 = np.array([.2, -.1, .3]), np.array([.3, .4, .5])
    truth, _ = exact_position_velocity(x0, u0, lam, t)
    forecast, _ = exact_position_velocity(x0, u0, fitted, t)
    position_bound = .5*t*t*parameter_bound  # target speed is at most one
    return {
        "true_lambda": lam.tolist(), "fitted_lambda": fitted.tolist(),
        "known_lambda_norm_cap": cap, "positive_readout_times": [0., h, 2*h],
        "per_position_vector_error_bound": epsilon,
        "single_launch_acceleration_error_bound": acceleration_bound,
        "stacked_matrix_singular_values": np.linalg.svd(matrix, compute_uv=False).tolist(),
        "lambda_error": float(np.linalg.norm(fitted-lam)),
        "lambda_error_bound": parameter_bound,
        "prediction_time": t, "target_speed_bound": 1.,
        "position_prediction_error": float(np.linalg.norm(truth-forecast)),
        "position_prediction_error_bound": position_bound,
        "position_readout_sigma": .1,
        "actual_prediction_TV": position_tv(truth-forecast, .1),
        "prediction_TV_bound": position_tv_bound(position_bound, .1),
        "scope": "Repeated preparations share one conserved lambda; the two initial "
                 "unit velocities, clock, position error and lambda cap are supplied. "
                 "This is deterministic bounded-error calibration, not a finite-shot confidence proof."
    }


def counterexample_case():
    u0, x0 = np.array([1., 0., 0.]), np.zeros(3)
    t = 1.
    rows = []
    for b in (0., 1.):
        initial = make_state(x0, u0, [0., 0., b])
        reference = full_reference_flow(initial, t)
        numerical = rk4(canonical_rhs, initial, t, 128)
        rows.append({"lambda": initial[9:].tolist(), "energy": hamiltonian(initial),
                     "initial_full_configuration_velocity": canonical_rhs(initial)[:6].tolist(),
                     "position_at_t": reference[:3].tolist(),
                     "velocity_at_t": velocity(reference).tolist(),
                     "internal_coordinate_at_t": reference[3:6].tolist(),
                     "canonical_RK4_full_state_error": float(np.linalg.norm(numerical-reference))})
    delta = np.array(rows[1]["position_at_t"])-rows[0]["position_at_t"]
    return {"time": t, "initial_position": x0.tolist(), "initial_velocity": u0.tolist(),
            "cases": rows, "position_separation": float(np.linalg.norm(delta)),
            "position_readout_sigma": .1,
            "future_position_record_TV": position_tv(delta, .1),
            "initial_position_record_TV": 0.,
            "scope": "Different dynamical preparations at the same R3 position; "
                     "future distinguishability does not change their initial location identity."}


def report():
    initial = make_state([.2, -.3, .1], [.7, .2, -.4], [.3, -.5, .8], [.1, .4, -.2])
    reference = full_reference_flow(initial, 1.3)
    steps = [8, 16, 32]
    errors = [float(np.linalg.norm(rk4(canonical_rhs, initial, 1.3, n)-reference))
              for n in steps]
    x0, u0, lam0 = np.array([.2, -.1, .3]), np.array([.8, .3, -.2]), np.array([.2, -.4, .7])
    dx, du, dl = np.array([.0006, -.0008, 0.]), np.array([.002, 0., 0.]), np.array([0., -.01, 0.])
    deviations = []
    for t in np.linspace(0, 2., 41):
        x, u = exact_position_velocity(x0, u0, lam0, t)
        xp, up = exact_position_velocity(x0+dx, u0+du, lam0+dl, t)
        deviations.append((np.linalg.norm(xp-x), np.linalg.norm(up-u)))
    return {
        "round": 379,
        "scope": {
            "baseline": "Frozen through 377; independent of round 378.",
            "claim": "Position identity need not be a complete state for autonomous "
                     "trajectory prediction. A specified R3 position remains R3 when "
                     "a hidden conserved launch label changes later position records.",
            "model_inputs": ["smooth R3_x times R3_z configuration and canonical cotangent phase space",
                             "autonomous H=|p_x+lambda cross x/2|^2/2, mass unit one",
                             "position readouts depend only on x; z is internal configuration",
                             "repeated preparations, supplied velocity directions, clock and error budgets"],
            "not_proved": ["a full FUCP countermodel or a cognition-derived Hamiltonian",
                           "all internal preparations must have the same free trajectory",
                           "trajectory nonclosure adds spatial dimensions",
                           "natural or unique physical interpretation of the x-position interface",
                           "a normalizable quantum state with exact continuous lambda",
                           "three-dimensional space, relativity or gravity from cognition"],
            "new_interface": "Autonomous canonical flow, fixed-preparation prediction, "
                             "uniform launch-error bounds and two-launch calibration; "
                             "not a repetition of quotient closure or Landau spectra."
        },
        "same_location_different_dynamical_preparations": counterexample_case(),
        "independent_full_phase_space_ODE": {"steps": steps, "full_state_errors": errors,
            "refinement_ratios": [errors[i]/errors[i+1] for i in range(2)],
            "z_reference": "40-point Gauss integration of exact x,u; x,u,lambda have closed formulas"},
        "uniform_preparation_example": {
            "time_window": [0., 2.], "dx_bound": .001, "du_bound": .002,
            "dlambda_bound": .01, "speed_bound": 1.,
            "sampled_maximum_position_error": float(np.max(np.array(deviations)[:, 0])),
            "sampled_maximum_velocity_error": float(np.max(np.array(deviations)[:, 1])),
            "analytic_uniform_bounds": preparation_bounds(2., .001, .002, .01, 1.),
            "scope": "General bounds follow from orthogonal Duhamel evolution; samples only cross-check."},
        "calibration_and_prediction": calibration_case(),
        "sources": ["https://cvgmt.sns.it/media/doc/paper/3191/cut-free-v2.pdf"],
        "source_location": "Rizzi-Serres, section 2, Example 2 and section 2.1 equations (13)-(14)."
    }


class Checks(unittest.TestCase):
    def test_01_canonical_equations_from_independent_energy_gradient(self):
        point = make_state([.2, -.3, .1], [.7, .2, -.4], [.3, -.5, .8], [.1, .4, -.2])
        rhs = canonical_matrix()@finite_gradient(hamiltonian, point)
        np.testing.assert_allclose(rhs, canonical_rhs(point), atol=2e-11)

    def test_02_poisson_bracket_signs(self):
        point = make_state([.2, -.3, .1], [.7, .2, -.4], [.3, -.5, .8])
        jac = np.vstack([finite_gradient(lambda y: velocity(y)[i], point) for i in range(3)])
        # {u_i,u_j}=[lambda]_cross_ij = -epsilon_ijk lambda_k.
        np.testing.assert_allclose(jac@canonical_matrix()@jac.T,
                                   cross_matrix(point[9:]), atol=2e-11)
        np.testing.assert_allclose(jac@canonical_rhs(point),
                                   np.cross(point[9:], velocity(point)), atol=2e-11)

    def test_03_closed_rotation_against_matrix_series_and_group_law(self):
        for lam in (np.zeros(3), np.array([1e-10, 0., 0.]), np.array([.3, -.5, .8])):
            matrix = .7*cross_matrix(lam)
            term, series = np.eye(3), np.eye(3)
            for n in range(1, 45):
                term = term@matrix/n
                series += term
            rotation, _ = propagation_matrices(lam, .7)
            np.testing.assert_allclose(rotation, series, atol=4e-16)
            np.testing.assert_allclose(rotation.T@rotation, np.eye(3), atol=4e-16)
            ra, ja = propagation_matrices(lam, .2)
            rb, jb = propagation_matrices(lam, .5)
            _, integral = propagation_matrices(lam, .7)
            np.testing.assert_allclose(rotation, ra@rb, atol=4e-16)
            np.testing.assert_allclose(integral, ja+ra@jb, atol=4e-16)

    def test_04_full_canonical_ODE_and_independent_reference_flow(self):
        initial = make_state([.2, -.3, .1], [.7, .2, -.4], [.3, -.5, .8], [.1, .4, -.2])
        for t in (-.6, 0., 1.3):
            exact = full_reference_flow(initial, t)
            numerical = rk4(canonical_rhs, initial, t, 160)
            np.testing.assert_allclose(numerical, exact, atol=9e-11)
            np.testing.assert_allclose(exact, full_reference_flow(initial, t, order=24), atol=5e-16)
        circular = full_reference_flow(make_state([0, 0, 0], [1, 0, 0], [0, 0, 1]), 1.)
        np.testing.assert_allclose(circular[3:6], [0, 0, (1-np.sin(1.))/2], atol=3e-16)

    def test_05_ODE_refinement_and_autonomous_conserved_quantities(self):
        initial = make_state([.2, -.3, .1], [.7, .2, -.4], [.3, -.5, .8])
        exact = full_reference_flow(initial, 1.3)
        errors = [np.linalg.norm(rk4(canonical_rhs, initial, 1.3, n)-exact)
                  for n in (8, 16, 32)]
        for i in range(2):
            self.assertGreater(errors[i]/errors[i+1], 14.)
            self.assertLess(errors[i]/errors[i+1], 18.)
        for t in np.linspace(-2., 2., 9):
            final = full_reference_flow(initial, t)
            self.assertAlmostEqual(hamiltonian(final), hamiltonian(initial), places=14)
            self.assertAlmostEqual(velocity(final)@final[9:], velocity(initial)@initial[9:], places=14)

    def test_06_same_location_velocity_energy_different_future(self):
        row = counterexample_case()
        straight, curved = row["cases"]
        self.assertEqual(straight["energy"], curved["energy"])
        np.testing.assert_array_equal(straight["initial_full_configuration_velocity"],
                                      curved["initial_full_configuration_velocity"])
        np.testing.assert_allclose(straight["position_at_t"], [1, 0, 0], atol=1e-15)
        np.testing.assert_allclose(curved["position_at_t"], [np.sin(1), 1-np.cos(1), 0], atol=1e-15)
        self.assertGreater(row["position_separation"], .48)

    def test_07_fixed_label_R3_canonical_reduction_and_internal_z_independence(self):
        initial = make_state([.2, -.3, .1], [.7, .2, -.4], [.3, -.5, .8], [.1, .4, -.2])
        reduced = np.concatenate((initial[:3], initial[6:9]))
        reduced_final = rk4(lambda y: reduced_rhs(y, initial[9:]), reduced, 1.1, 128)
        full_final = rk4(canonical_rhs, initial, 1.1, 128)
        np.testing.assert_allclose(reduced_final, np.r_[full_final[:3], full_final[6:9]], atol=1e-15)
        shifted = initial.copy()
        shifted[3:6] += [100., -30., 7.]
        shifted_final = rk4(canonical_rhs, shifted, 1.1, 128)
        np.testing.assert_array_equal(shifted_final[:3], full_final[:3])
        np.testing.assert_array_equal(shifted_final[6:], full_final[6:])

    def test_08_energy_does_not_bound_label_and_one_launch_can_be_blind(self):
        for b in (0., 1., 1e4):
            initial = make_state([0, 0, 0], [1, 0, 0], [b, 0, 0])
            self.assertEqual(hamiltonian(initial), .5)
            x, u = exact_position_velocity(initial[:3], velocity(initial), initial[9:], .6)
            np.testing.assert_allclose(x, [.6, 0, 0], atol=1e-15)
            np.testing.assert_allclose(u, [1, 0, 0], atol=1e-15)
        self.assertEqual(np.linalg.matrix_rank(launch_matrix([[1., 0, 0]])), 2)

    def test_09_uniform_preparation_bound_without_label_magnitude_cap(self):
        rng = np.random.default_rng(379)
        for scale in (0., 1., 50.):
            for _ in range(7):
                x, u, lam = rng.normal(size=(3, 3))
                lam *= scale
                dx, du, dl = rng.normal(size=(3, 3))*np.array([.01, .02, .03])[:, None]
                speed = max(np.linalg.norm(u), np.linalg.norm(u+du))
                for t in (-1.2, .05, .8, 2.):
                    a, v = exact_position_velocity(x, u, lam, t)
                    b, w = exact_position_velocity(x+dx, u+du, lam+dl, t)
                    bound = preparation_bounds(t, np.linalg.norm(dx), np.linalg.norm(du), np.linalg.norm(dl), speed)
                    self.assertLessEqual(np.linalg.norm(b-a), bound["position"]+1e-13)
                    self.assertLessEqual(np.linalg.norm(w-v), bound["velocity"]+1e-13)

    def test_10_quadratic_preparation_position_bound_has_sharp_leading_constant(self):
        t, epsilon = .7, 1e-5
        straight, _ = exact_position_velocity(np.zeros(3), [1, 0, 0], np.zeros(3), t)
        curved, _ = exact_position_velocity(np.zeros(3), [1, 0, 0], [0, 0, epsilon], t)
        ratio = np.linalg.norm(curved-straight)/(.5*t*t*epsilon)
        self.assertGreater(ratio, .999999)
        self.assertLessEqual(ratio, 1.+1e-12)

    def test_11_current_position_identity_and_future_readout_are_different(self):
        row = counterexample_case()
        self.assertEqual(row["initial_position_record_TV"], 0.)
        self.assertGreater(row["future_position_record_TV"], .98)
        for displacement in (0., .001, .1, 1.):
            self.assertLessEqual(position_tv([displacement, 0, 0], .1),
                                 position_tv_bound(displacement, .1)+1e-15)

    def test_12_two_declared_launch_directions_identify_the_conserved_label(self):
        matrix = launch_matrix(np.eye(3)[:2])
        np.testing.assert_allclose(np.linalg.svd(matrix, compute_uv=False), [np.sqrt(2), 1, 1])
        lam = np.array([.4, -.3, .7])
        accelerations = np.concatenate([np.cross(lam, u) for u in np.eye(3)[:2]])
        np.testing.assert_allclose(np.linalg.lstsq(matrix, accelerations, rcond=None)[0], lam, atol=3e-16)

    def test_13_finite_positive_time_acceleration_and_readout_error_budget(self):
        lam = np.array([.4, -.3, .7])
        u, epsilon = np.array([1., 0, 0]), 2e-6
        noise = epsilon*np.array([[1., 0, 0], [-1., 0, 0], [1., 0, 0]])
        for h in (.02, .04, .08):
            estimate = acceleration_readout(lam, u, h, noise)
            exact = np.cross(lam, u)
            bound = h*np.linalg.norm(lam)**2*np.linalg.norm(u)+4*epsilon/h**2
            self.assertLessEqual(np.linalg.norm(estimate-exact), bound+1e-12)

    def test_14_calibration_uncertainty_propagates_to_finite_time_position_records(self):
        row = calibration_case()
        self.assertLess(row["lambda_error"], row["lambda_error_bound"])
        self.assertLess(row["position_prediction_error"], row["position_prediction_error_bound"])
        self.assertLess(row["actual_prediction_TV"], row["prediction_TV_bound"])
        self.assertLess(row["prediction_TV_bound"], 1.)


if __name__ == "__main__":
    main(__name__, "position_prediction_closure_audit", report)
