"""Round 304: conditional Jacobson closure and a self-consistent scalar FLRW branch."""
import unittest
import numpy as np
from conformal_volume_geometry_audit import ricci_from_connection, MINKOWSKI
from growing_stream_audit import main

INDICES = [(i, j) for i in range(4) for j in range(i, 4)]


def null_directions():
    axes = [np.r_[1., sign*np.eye(3)[i]] for i in range(3) for sign in (-1, 1)]
    mixed = [np.r_[1., (np.eye(3)[i]+np.eye(3)[j])/np.sqrt(2)]
             for i, j in ((0, 1), (0, 2), (1, 2))]
    return np.array(axes+mixed)


def null_measurement_matrix(directions):
    return np.array([[k[i]*k[j]*(1 if i == j else 2) for i, j in INDICES] for k in directions])


def pack(tensor):
    return np.array([tensor[i, j] for i, j in INDICES])


def unpack(vector):
    result = np.zeros((4, 4))
    for value, (i, j) in zip(vector, INDICES):
        result[i, j] = result[j, i] = value
    return result


def trace_free(tensor):
    return tensor-np.trace(MINKOWSKI@tensor)*MINKOWSKI/4


def horizon_budget(tkk, sigma, acceleration, hbar=1., area=0.7, epsilon=0.01):
    beta = 2*np.pi/(hbar*sigma)
    rkk = beta*tkk
    temperature = hbar*acceleration/(2*np.pi)
    # Leading local-equilibrium terms, not a finite-horizon exact integration.
    heat = acceleration*tkk*area*epsilon**2/2
    entropy = sigma*rkk*area*epsilon**2/2
    return heat, temperature*entropy, beta


def rhs(time, state, beta=1.):
    a, phi, velocity, hubble = state
    return np.array([a*hubble, velocity, -3*hubble*velocity, -beta*velocity**2/2])


def integrate_background(steps, right=4., beta=1., initial_hubble=1/3):
    times = np.linspace(1., right, steps+1)
    step = (right-1)/steps
    values = np.empty((steps+1, 4))
    values[0] = [1., 0., np.sqrt(2/(3*beta)), initial_hubble]
    for i, t in enumerate(times[:-1]):
        y = values[i]
        d1 = rhs(t, y, beta)
        d2 = rhs(t+step/2, y+step*d1/2, beta)
        d3 = rhs(t+step/2, y+step*d2/2, beta)
        d4 = rhs(t+step, y+step*d3, beta)
        values[i+1] = y+step*(d1+2*d2+2*d3+d4)/6
    return times, values


def exact_background(t, beta=1.):
    return np.array([t**(1/3), np.sqrt(2/(3*beta))*np.log(t),
                     np.sqrt(2/(3*beta))/t, 1/(3*t)])


def constraint(values, beta=1.):
    return values[:, 3]**2-beta*values[:, 2]**2/6


def conserved_test_field_counterexample(t=0.):
    h = 0.5
    velocity = np.exp(-3*h*t)
    rho = velocity**2/2
    return {'t': t, 'rho': float(rho), 'p': float(rho),
            'continuity_residual': float(-3*h*velocity**2+3*h*(rho+rho)),
            'null_closure_residual': float(-velocity**2)}


def report():
    matrix = null_measurement_matrix(null_directions())
    errors, constraints, momenta = [], [], []
    for steps in (32, 64, 128):
        times, states = integrate_background(steps)
        errors.append(float(np.max(abs(states[-1]-exact_background(times[-1])))))
        constraints.append(float(np.max(abs(constraint(states)))))
        momenta.append(float(np.max(abs(states[:, 0]**3*states[:, 2]-np.sqrt(2/3)))))
    rng = np.random.default_rng(304)
    raw = rng.normal(size=(4, 4))
    source = (raw+raw.T)/2
    restored = unpack(np.linalg.lstsq(matrix, matrix@pack(source), rcond=None)[0])
    shifted = []
    for initial_hubble in (1/3, 0.5):
        _, states = integrate_background(256, initial_hubble=initial_hubble)
        c = initial_hubble**2-1/9
        shifted.append({'initial_H': initial_hubble, 'Lambda': 3*c,
                        'final_a': float(states[-1, 0]),
                        'max_constant_drift': float(np.max(abs(constraint(states)-c)))})
    return {'round': 304,
            'scope': 'Known Jacobson argument under supplied area entropy, Unruh temperature, all-horizon equilibrium and conserved matter; classical scalar backreaction calibration, not first-principles cognition-to-gravity derivation.',
            'inputs': {'dimension': 4, 'entropy_density_sigma': float(2*np.pi), 'hbar': 1.,
                       'beta_8piG': 1., 'scalar_xi': 0., 'scalar_mass': 0.,
                       'spatial_curvature': 0., 'homogeneous_scalar': True},
            'null_tensor_interface': {'directions': 9, 'rank': int(np.linalg.matrix_rank(matrix)),
                                      'axis_only_rank': int(np.linalg.matrix_rank(matrix[:6])),
                                      'metric_kernel_error': float(np.max(abs(matrix@pack(MINKOWSKI)))),
                                      'trace_free_reconstruction_error': float(np.max(abs(restored-trace_free(source))))},
            'stiff_scalar_background': {'time_interval': [1., 4.], 'steps': [32, 64, 128],
                                         'endpoint_errors': errors,
                                         'error_ratios': [errors[i]/errors[i+1] for i in (0, 1)],
                                         'constraint_drifts': constraints, 'momentum_drifts': momenta,
                                         'exact_final_state': exact_background(4.).tolist(),
                                         'initial_R': -2/3, 'final_R': -1/24},
            'integration_constant_branches': shifted,
            'conserved_test_field_on_fixed_de_sitter': conserved_test_field_counterexample(),
            'remaining_inputs': ['area entropy law and coefficient', 'physical Unruh/KMS response',
                                 'local equilibrium for all null directions', 'state and cosmological constant selection']}


class BackreactionTests(unittest.TestCase):
    def test_01_null_measurements_rank_and_metric_kernel(self):
        directions = null_directions()
        self.assertLess(np.max(abs(np.einsum('ni,ij,nj->n', directions, MINKOWSKI, directions))), 1e-14)
        matrix = null_measurement_matrix(directions)
        self.assertEqual(np.linalg.matrix_rank(matrix), 9)
        self.assertLess(np.max(abs(matrix@pack(MINKOWSKI))), 1e-14)

    def test_02_tensor_reconstruction_up_to_trace(self):
        matrix = null_measurement_matrix(null_directions())
        rng = np.random.default_rng(304)
        for _ in range(10):
            raw = rng.normal(size=(4, 4))
            tensor = (raw+raw.T)/2
            restored = unpack(np.linalg.lstsq(matrix, matrix@pack(tensor), rcond=None)[0])
            np.testing.assert_allclose(restored, trace_free(tensor), atol=1e-13)

    def test_03_axis_directions_miss_spatial_shear(self):
        tensor = np.zeros((4, 4))
        tensor[1, 2] = tensor[2, 1] = 1
        matrix = null_measurement_matrix(null_directions())
        self.assertEqual(np.linalg.matrix_rank(matrix[:6]), 6)
        np.testing.assert_allclose(matrix[:6]@pack(tensor), 0)
        self.assertAlmostEqual((matrix@pack(tensor))[6], 1.)

    def test_04_cosmological_term_is_null_invisible(self):
        rng = np.random.default_rng(1304)
        spatial = rng.normal(size=(64, 3))
        spatial /= np.linalg.norm(spatial, axis=1)[:, None]
        matrix = null_measurement_matrix(np.column_stack((np.ones(64), spatial)))
        tensor = np.diag([2., 3., 4., 5.])
        for constant in (-3., 0., 8.):
            np.testing.assert_allclose(matrix@pack(tensor+constant*MINKOWSKI), matrix@pack(tensor), atol=1e-13)

    def test_05_clausius_coefficient_and_acceleration_cancellation(self):
        for sigma in (1., 2*np.pi, 10.):
            for acceleration in (0.1, 1., 20.):
                heat, tdS, beta = horizon_budget(0.3, sigma, acceleration)
                self.assertAlmostEqual(heat, tdS)
                self.assertAlmostEqual(beta*sigma, 2*np.pi)

    def test_06_geometric_null_projection_matches_scalar_source(self):
        for t in (1., 2., 4.):
            a, _, velocity, h = exact_background(t)
            dh = -1/(3*t*t)
            ricci, scalar = ricci_from_connection((1, 0, 0), (a, a*h, a*(dh+h*h)))
            for vector in null_directions():
                k = vector.copy()
                k[1:] /= a
                self.assertAlmostEqual(k@ricci@k, velocity**2)
            self.assertAlmostEqual(scalar, -velocity**2)

    def test_07_fourth_order_backreaction_convergence(self):
        data = report()['stiff_scalar_background']
        for ratio in data['error_ratios']:
            self.assertTrue(13. < ratio < 20., ratio)
        self.assertLess(data['endpoint_errors'][-1], 1e-7)

    def test_08_constraint_and_scalar_momentum_propagate(self):
        _, states = integrate_background(128)
        self.assertLess(np.max(abs(constraint(states))), 1e-8)
        self.assertLess(np.max(abs(states[:, 0]**3*states[:, 2]-np.sqrt(2/3))), 1e-7)
        # Analytic derivative of C=H^2-beta*v^2/6 vanishes even off C=0.
        for y in ([1, 0, 0.4, 0.8], [2, 1, -0.7, 0.3]):
            derivative = rhs(1., y)
            self.assertAlmostEqual(2*y[3]*derivative[3]-y[2]*derivative[2]/3, 0.)

    def test_09_different_lambda_and_vacuum_expansion_remain_allowed(self):
        data = report()['integration_constant_branches']
        self.assertGreater(abs(data[0]['final_a']-data[1]['final_a']), 0.2)
        for row in data:
            self.assertLess(row['max_constant_drift'], 1e-8)
        for h in (0.2, 0.5, 1.):
            np.testing.assert_allclose(rhs(1., [1, 0, 0, h])[2:], [0, 0])
            ricci, scalar = ricci_from_connection((1, 0, 0), (1, h, h*h))
            einstein = ricci-scalar*MINKOWSKI/2
            np.testing.assert_allclose(einstein+3*h*h*MINKOWSKI, 0, atol=1e-14)

    def test_10_conserved_matter_does_not_imply_backreaction(self):
        for t in (0., 0.5, 1.):
            data = conserved_test_field_counterexample(t)
            self.assertEqual(data['continuity_residual'], 0.)
            self.assertLess(data['null_closure_residual'], -0.01)


if __name__ == '__main__':
    main(__name__, 'horizon_backreaction_audit', report)
