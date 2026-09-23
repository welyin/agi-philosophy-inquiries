"""Round 322: massive test worldlines in the positive-F stealth background."""
from functools import lru_cache
import unittest
import numpy as np
from growing_stream_audit import main
from heat_kernel_area_audit import gauss

ETA = np.diag([-1., 1., 1., 1.])
NULL = np.array([1., -1., 0., 0.])


def conformal(u):
    u = np.asarray(u)
    if np.any(u <= 0):
        raise ValueError('Use only the connected u>0, F>0 branch.')
    f = 1-1/(1+u)**2
    return np.sqrt(f), 1/((1+u)**3*f)


def connection(gradient):
    eye = np.eye(4)
    return (np.einsum('ab,c->abc', eye, gradient)
            +np.einsum('ac,b->abc', eye, gradient)
            -np.einsum('bc,a->abc', ETA, ETA@gradient))


def primitive(u):
    return np.sqrt((1+u)**2-1)-np.arccos(1/(1+u))


def integrate(rhs, initial, end, steps):
    values = [np.asarray(initial, dtype=float)]
    h = end/steps
    for _ in range(steps):
        y = values[-1]
        a = rhs(y); b = rhs(y+h*a/2)
        c = rhs(y+h*b/2); d = rhs(y+h*c)
        values.append(y+h*(a+2*b+2*c+d)/6)
    return np.asarray(values)


def worldline(steps=128, transformed_mass=True, mass=2., velocity=.3, duration=1.2):
    gamma = 1/np.sqrt(1-velocity**2)
    j_velocity = gamma*np.array([1., velocity, 0., 0.])
    start = np.array([1., 0., 0., 0.])
    finish = start+duration*j_velocity
    rate = NULL@j_velocity
    end = (primitive(NULL@finish)-primitive(NULL@start))/rate
    omega0 = conformal(NULL@start)[0]
    initial = np.r_[start, j_velocity/omega0, 0., 0.]

    def rhs(y):
        omega, log_gradient = conformal(NULL@y[:4])
        w = y[4:8]
        gradient = log_gradient*NULL
        projector = ETA/omega**2+np.outer(w, w)
        force = projector@gradient if transformed_mass else np.zeros(4)
        acceleration = -np.einsum('abc,b,c->a', connection(gradient), w, w)+force
        return np.r_[w, acceleration, 1/omega, mass/omega]

    return integrate(rhs, initial, end, steps), finish, end


@lru_cache(maxsize=12)
def audit(steps=128, transformed_mass=True):
    y, target, end = worldline(steps, transformed_mass)
    omega = conformal(y[:, :4]@NULL)[0]
    norms = omega**2*np.einsum('na,ab,nb->n', y[:, 4:8], ETA, y[:, 4:8])
    path = y[:, 1]-.3*(y[:, 0]-1)
    return {'endpoint_error': float(np.max(abs(y[-1, :4]-target))),
            'straight_line_residual': float(np.max(abs(path))),
            'normalization_error': float(np.max(abs(norms+1))),
            'Jordan_elapsed_error': float(abs(y[-1, 8]-1.2)),
            'phase_error_against_Jordan': float(abs(y[-1, 9]-2.4)),
            'Einstein_proper_duration': float(end),
            'final_coordinate_velocity': float(y[-1, 5]/y[-1, 4]),
            'final_position': y[-1, :2].tolist()}


def report():
    return {'round': 322,
            'scope': 'Supplied four-dimensional classical scalar-tensor stealth solution, u=t-x>0, F=1-(1+u)^(-2). Massive structureless test probes are an added input; backreaction is neglected. The correct change of variables transforms both metric and mass. No universal matter coupling or quantum frame equivalence is derived.',
            'correct_mass_at_128_steps': audit(),
            'incorrect_constant_mass_control': audit(128, False),
            'convergence': [{'steps': n, 'endpoint_error': audit(n)['endpoint_error']} for n in (16, 32, 64, 128)],
            'next_interface': 'Compare actual dimensionless phase, radar and photon readouts with consistently transformed units.'}


class Checks(unittest.TestCase):
    def test_01_connection_from_metric_derivatives(self):
        gradient = conformal(.7)[1]*NULL
        f = conformal(.7)[0]**2
        g = f*ETA
        dg = 2*np.einsum('a,bc->abc', gradient, g)
        direct = np.zeros((4, 4, 4))
        for a in range(4):
            for b in range(4):
                for c in range(4):
                    direct[a, b, c] = sum(np.linalg.inv(g)[a, d]*(dg[b, d, c]+dg[c, d, b]-dg[d, b, c])/2 for d in range(4))
        np.testing.assert_allclose(direct, connection(gradient), atol=2e-15)

    def test_02_exact_parameter_against_quadrature(self):
        nodes, weights = gauss(48)
        rate = (1-.3)/np.sqrt(1-.3**2)
        numerical = .6*np.dot(weights, conformal(1+rate*.6*(nodes+1))[0])
        self.assertAlmostEqual(numerical, worldline(16)[2], places=13)

    def test_03_correct_orbit_recovers_Jordan_line(self):
        self.assertLess(audit()['endpoint_error'], 2e-10)
        self.assertLess(audit()['straight_line_residual'], 2e-11)

    def test_04_timelike_normalization_is_preserved(self):
        self.assertLess(audit()['normalization_error'], 2e-10)
        self.assertLess(audit(128, False)['normalization_error'], 2e-10)

    def test_05_RK4_convergence(self):
        errors = [audit(n)['endpoint_error'] for n in (16, 32, 64)]
        self.assertGreater(errors[0]/errors[1], 13)
        self.assertGreater(errors[1]/errors[2], 13)

    def test_06_old_proper_parameter_is_recovered(self):
        self.assertLess(audit()['Jordan_elapsed_error'], 2e-10)
        self.assertLess(audit()['phase_error_against_Jordan'], 4e-10)

    def test_07_constant_mass_changes_the_actual_path(self):
        self.assertGreater(audit(128, False)['straight_line_residual'], .01)
        self.assertGreater(abs(audit(128, False)['final_coordinate_velocity']-.3), .01)

    def test_08_overall_mass_does_not_change_free_fall(self):
        a = worldline(64, mass=1.)[0]
        b = worldline(64, mass=7.)[0]
        np.testing.assert_array_equal(a[:, :9], b[:, :9])
        np.testing.assert_allclose(b[:, 9], 7*a[:, 9], atol=2e-14)

    def test_09_transformed_force_is_orthogonal(self):
        w = np.array([1., .4, .1, 0.])/np.sqrt(1-.4**2-.1**2)
        omega, q = conformal(.8)
        w /= omega
        force = (ETA/omega**2+np.outer(w, w))@(q*NULL)
        self.assertLess(abs(w@(omega**2*ETA)@force), 2e-15)
        self.assertGreater(np.linalg.norm(force), .1)

    def test_10_positive_branch_guard(self):
        for value in (0., -.25):
            with self.assertRaises(ValueError):
                conformal(value)


if __name__ == '__main__':
    main(__name__, 'conformal_probe_motion_audit', report)
