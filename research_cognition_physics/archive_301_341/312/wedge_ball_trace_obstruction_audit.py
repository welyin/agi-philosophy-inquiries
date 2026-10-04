"""Round 312: exact stress-current distinction between boost and conformal ball flow."""
import unittest
import numpy as np
from growing_stream_audit import main

ETA = np.diag([-1., 1., 1., 1.])


def vector(point, region, radius=2.):
    t, *xyz = point
    xyz = np.array(xyz)
    if region == 'wedge':
        return np.array([xyz[0], t, 0., 0.])
    if region == 'ball':
        return np.r_[(radius**2-t*t-xyz@xyz)/(2*radius), -t*xyz/radius]
    raise ValueError('Choose wedge or ball.')


def jacobian(point, region, radius=2.):
    t, *xyz = point
    if region == 'wedge':
        result = np.zeros((4, 4)); result[0, 1] = result[1, 0] = 1.
        return result
    if region == 'ball':
        result = -t/radius*np.eye(4)
        result[0, 1:] = -np.array(xyz)/radius
        result[1:, 0] = -np.array(xyz)/radius
        return result
    raise ValueError('Choose wedge or ball.')


def derivative(function, point, axis, step=1e-4):
    shift = np.eye(4)[axis]*step
    return (-function(point+2*shift)+8*function(point+shift)-8*function(point-shift)+function(point-2*shift))/(12*step)


def scalar_stress(point, mass=0.8, amplitude=1.2):
    phase = mass*point[0]
    energy = 0.5*mass**2*amplitude**2
    pressure = -energy*np.cos(2*phase)
    return np.diag([energy, pressure, pressure, pressure])


def radiation_stress(point):
    return np.diag([1., 1/3, 1/3, 1/3])


def current(point, region, matter=scalar_stress, radius=2.):
    return matter(point)@ETA@vector(point, region, radius)


def direct_divergence(point, region, matter=scalar_stress, radius=2.):
    return float(sum(derivative(lambda p: current(p, region, matter, radius)[mu], point, mu) for mu in range(4)))


def ward_source(point, region, matter=scalar_stress, radius=2.):
    lowered = ETA@jacobian(point, region, radius)
    return float(0.5*np.sum(matter(point)*(lowered+lowered.T)))


def fixed_ball_balance(t, radius=2., mass=0.8, amplitude=1.2):
    energy = 0.5*mass**2*amplitude**2
    pressure = -energy*np.cos(2*mass*t)
    volume = 4*np.pi*radius**3/3
    charge_rate = -energy*volume*t/radius
    outward_flux = -3*pressure*volume*t/radius
    integrated_source = (energy-3*pressure)*volume*t/radius
    return {'time': t, 'charge_rate': charge_rate, 'outward_current_flux': outward_flux,
            'integrated_trace_source': integrated_source,
            'balance_error': abs(-charge_rate+outward_flux-integrated_source)}


def report():
    rows = []
    for t in (0., 0.2, 0.4, 0.8):
        point = np.array([t, 0.6, -0.2, 0.1])
        trace = float(np.sum(ETA*scalar_stress(point)))
        rows.append({'time': t, 'stress_trace': trace,
                     'wedge_current_divergence': direct_divergence(point, 'wedge'),
                     'ball_current_divergence': direct_divergence(point, 'ball'),
                     'ball_predicted_trace_source': -t*trace/2,
                     'traceless_radiation_ball_divergence': direct_divergence(point, 'ball', radiation_stress)})
    return {'round': 312,
            'scope': 'Analytic Ward identity with supplied Minkowski geometry and conserved matter. Classical scalar solutions are counterexamples/calibrations, not a quantum vacuum or a derivation of gravity.',
            'local_current_cases': rows,
            'finite_region_balance': [fixed_ball_balance(t) for t in (0.2, 0.4, 0.8)],
            'theorem_interface': {'wedge': 'BW requires an appropriate local relativistic vacuum theory; massive models can satisfy it.',
                                  'ball': 'The simple conformal stress charge needs additional structure. Conserved T alone leaves a trace source away from t=0.'},
            'research_decision': 'Prioritize the conditional local-horizon route for general matter. Keep ball corrections as a separate branch; next audit the common origin of area entropy and gravitational coupling.',
            'remaining_inputs': ['Verified Wightman or other sufficient BW assumptions, beyond locality and covariance alone', 'geometric scaling limit', 'area entropy coefficient and renormalization', 'local equilibrium condition']}


class Checks(unittest.TestCase):
    def test_01_vector_jacobian_independent(self):
        point = np.array([0.4, 0.6, -0.2, 0.1])
        for region in ('wedge', 'ball'):
            direct = np.column_stack([derivative(lambda p: vector(p, region), point, axis) for axis in range(4)])
            np.testing.assert_allclose(direct, jacobian(point, region), atol=1e-11)

    def test_02_boost_killing(self):
        lowered = ETA@jacobian(np.ones(4), 'wedge')
        np.testing.assert_allclose(lowered+lowered.T, 0, atol=1e-15)

    def test_03_ball_conformal_killing(self):
        point = np.array([0.4, 0.6, -0.2, 0.1]); radius = 2.
        lowered = ETA@jacobian(point, 'ball', radius)
        np.testing.assert_allclose(lowered+lowered.T, -2*point[0]/radius*ETA, atol=1e-15)

    def test_04_scalar_stress_is_conserved(self):
        point = np.array([0.4, 0.6, -0.2, 0.1])
        for nu in range(4):
            divergence = sum(derivative(lambda p: scalar_stress(p)[mu, nu], point, mu) for mu in range(4))
            self.assertLess(abs(divergence), 1e-11)

    def test_05_massive_current_ward_identity(self):
        for t in (0.2, 0.4, 0.8):
            point = np.array([t, 0.6, -0.2, 0.1])
            for region in ('wedge', 'ball'):
                self.assertAlmostEqual(direct_divergence(point, region), ward_source(point, region), places=10)

    def test_06_genuine_nonzero_ball_obstruction(self):
        point = np.array([0.4, 0.6, -0.2, 0.1])
        self.assertGreater(abs(direct_divergence(point, 'ball')), 0.2)
        self.assertLess(abs(direct_divergence(point, 'wedge')), 1e-11)

    def test_07_traceless_control(self):
        point = np.array([0.4, 0.6, -0.2, 0.1])
        self.assertLess(abs(direct_divergence(point, 'ball', radiation_stress)), 1e-11)

    def test_08_t0_can_hide_obstruction(self):
        point = np.array([0., 0.6, -0.2, 0.1])
        self.assertLess(abs(direct_divergence(point, 'ball')), 1e-11)
        self.assertGreater(abs(np.sum(ETA*scalar_stress(point))), 1)

    def test_09_integrated_balance(self):
        from numpy.polynomial.legendre import leggauss
        nodes, weights = leggauss(32)
        radius = 2.; radial = radius*(nodes+1)/2
        for t in (0.2, 0.4, 0.8):
            r = fixed_ball_balance(t)
            point = np.array([t, 0., 0., 0.])
            def charge(p):
                values = [-current(np.array([p[0], x, 0., 0.]), 'ball')[0] for x in radial]
                return 4*np.pi*np.dot(weights, radial**2*values)*radius/2
            rate = derivative(charge, point, 0)
            flux = 4*np.pi*radius**2*current(np.array([t, radius, 0., 0.]), 'ball')[1]
            source = 4*np.pi*np.dot(weights, radial**2*[ward_source(np.array([t, x, 0., 0.]), 'ball') for x in radial])*radius/2
            self.assertAlmostEqual(rate, r['charge_rate'], places=9)
            self.assertAlmostEqual(flux, r['outward_current_flux'], places=11)
            self.assertAlmostEqual(source, r['integrated_trace_source'], places=11)
            self.assertLess(r['balance_error'], 1e-12)
            self.assertGreater(abs(r['integrated_trace_source']), 1)


if __name__ == '__main__':
    main(__name__, 'wedge_ball_trace_obstruction_audit', report)
