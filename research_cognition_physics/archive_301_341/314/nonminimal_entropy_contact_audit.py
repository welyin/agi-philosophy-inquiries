"""Round 314: distinguish ordinary flat-space entropy and a conical/Wald contact term."""
import unittest
import numpy as np
from growing_stream_audit import main
from heat_kernel_area_audit import proper_moment, richardson_curvature, inverse_g_shift, gauss, e1_small


def terms(mass, cutoff, xi):
    j = proper_moment(mass, cutoff)
    phi_squared = j/(16*np.pi**2)
    statistical = j/(48*np.pi)
    contact = -2*np.pi*xi*phi_squared
    return {'xi': xi, 'regulated_phi_squared': phi_squared,
            'ordinary_entropy_per_area': statistical,
            'contact_entropy_per_area': contact,
            'geometric_entropy_per_area': statistical+contact,
            'inverse_G_shift_over_four': inverse_g_shift(mass, cutoff, xi)/4}


def coincident_propagator_momentum(mass, cutoff):
    nodes, weights = gauss(192)
    v = (nodes+1)*np.sqrt(40)/2
    x = (mass*cutoff)**2
    integral = np.dot(weights, v**3*np.exp(-v*v)/(v*v+x))*np.sqrt(40)/2
    return float(np.exp(-x)*integral/(8*np.pi**2*cutoff**2))


def regulated_action(matrix, cutoff):
    eigenvalues = np.linalg.eigvalsh(matrix)
    if eigenvalues.min() <= 0 or cutoff**2*eigenvalues.max() > 1:
        raise ValueError('Positive operator and small-argument E1 regime required.')
    return -0.5*sum(e1_small(cutoff**2*v) for v in eigenvalues)


def determinant_contact(xi=0.2, step=1e-3, position=2):
    n, cutoff = 6, 0.2
    laplacian = 2*np.eye(n)-np.roll(np.eye(n), 1, axis=0)-np.roll(np.eye(n), -1, axis=0)
    matrix = laplacian+np.diag(0.8+np.arange(n)/5)
    source = np.zeros((n, n)); source[position, position] = 4*np.pi
    action = lambda a: regulated_action(matrix+xi*(1-a)*source, cutoff)
    derivative = (action(1+step)-action(1-step))/(2*step)
    eigenvalues, modes = np.linalg.eigh(matrix)
    regulated_inverse = (modes*(np.exp(-cutoff**2*eigenvalues)/eigenvalues))@modes.T
    prediction = -0.5*xi*np.trace(regulated_inverse@source)
    return {'xi': xi, 'step': step, 'surface_position': position,
            'spectral_determinant_derivative': float(derivative),
            'local_contact_prediction': float(prediction),
            'derivative_error': float(abs(derivative-prediction))}


def report():
    return {'round': 314,
            'scope': 'One-loop nonminimal real scalar, linear conical response in the local area/EH sector. Contact is a geometric/Wald contribution, not a nonnegative von Neumann entropy.',
            'mass': 1., 'cutoff': 0.2,
            'entropy_and_coupling_terms': [terms(1., 0.2, xi) for xi in (-0.1, 0., 1/12, 1/6, 1/3)],
            'sphere_coefficient_checks': [{'xi': xi, 'spectral_R_coefficient': richardson_curvature(xi=xi),
                                           'expected': 1/6-xi} for xi in (-0.1, 0., 1/6, 1/3)],
            'determinant_contact_convergence': [determinant_contact(step=h) for h in (0.01, 0.005, 0.0025, 0.00125)],
            'local_contact_controls': [determinant_contact(step=1e-4, position=i) for i in (0, 4)],
            'flat_state_observation': 'At fixed flat metric R=0 and the same canonical algebra, mass, and state prescription, xi does not change the free scalar ground state. Curvature differentiation can still depend on xi.',
            'replica_scope': 'Only the first variation about a smooth cone is used. Singular conical potentials at finite deficit require a separate definition.'}


class Checks(unittest.TestCase):
    def test_01_spectral_curvature_response(self):
        for xi in (-0.1, 0., 1/6, 1/3):
            self.assertLess(abs(richardson_curvature(xi=xi)-(1/6-xi)), 2e-7)

    def test_02_conformal_coupling_does_not_zero_entanglement(self):
        r = terms(1., 0.2, 1/6)
        self.assertGreater(r['ordinary_entropy_per_area'], 0.1)
        self.assertLess(abs(r['geometric_entropy_per_area']), 1e-14)
        self.assertEqual(r['inverse_G_shift_over_four'], 0.)

    def test_03_negative_geometric_term(self):
        r = terms(1., 0.2, 1/3)
        self.assertGreater(r['ordinary_entropy_per_area'], 0.)
        self.assertLess(r['geometric_entropy_per_area'], 0.)
        self.assertAlmostEqual(r['geometric_entropy_per_area'], r['inverse_G_shift_over_four'])

    def test_04_positive_xi_determinant_derivative(self):
        r = determinant_contact(xi=0.2, step=1e-4)
        self.assertLess(r['local_contact_prediction'], 0.)
        self.assertLess(r['derivative_error'], 1e-7)

    def test_05_negative_xi_determinant_derivative(self):
        r = determinant_contact(xi=-0.1, step=1e-4)
        self.assertGreater(r['local_contact_prediction'], 0.)
        self.assertLess(r['derivative_error'], 1e-7)

    def test_06_independent_momentum_propagator(self):
        for mass in (0.3, 1., 3.):
            self.assertAlmostEqual(coincident_propagator_momentum(mass, 0.2), proper_moment(mass, 0.2)/(16*np.pi**2), places=10)

    def test_07_contact_derivative_converges_quadratically(self):
        values = [determinant_contact(step=h)['derivative_error'] for h in (0.01, 0.005, 0.0025)]
        self.assertTrue(all(3.9 < a/b < 4.1 for a, b in zip(values, values[1:])))

    def test_08_contact_samples_local_state(self):
        first = determinant_contact(step=1e-4, position=0)
        last = determinant_contact(step=1e-4, position=4)
        self.assertGreater(abs(first['local_contact_prediction']-last['local_contact_prediction']), 0.05)
        self.assertLess(max(first['derivative_error'], last['derivative_error']), 1e-7)


if __name__ == '__main__':
    main(__name__, 'nonminimal_entropy_contact_audit', report)
