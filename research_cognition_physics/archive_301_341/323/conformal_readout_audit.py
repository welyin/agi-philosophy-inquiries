"""Round 323: dimensionless test-matter readouts after a conformal frame change."""
import unittest
import numpy as np
from growing_stream_audit import main
from heat_kernel_area_audit import gauss
from conformal_probe_motion_audit import ETA, conformal, primitive


def quadrature(function, lower, upper, order=48):
    nodes, weights = gauss(order)
    x = (lower+upper)/2+(upper-lower)*nodes/2
    return (upper-lower)/2*np.dot(weights, function(x))


def radar(mass=2., separation=.4, emission=1.):
    arrival = emission+2*separation
    omega = lambda t: conformal(t)[0]
    proper = quadrature(omega, emission, arrival)
    phase = quadrature(lambda t: mass/omega(t)*omega(t), emission, arrival)
    return {'Jordan_elapsed': 2*separation, 'Einstein_proper_elapsed': float(proper),
            'Einstein_elapsed_primitive': float(primitive(arrival)-primitive(emission)),
            'roundtrip_phase': float(phase), 'dimensionless_distance': float(phase/2),
            'Jordan_dimensionless_distance': mass*separation,
            'incorrect_fixed_mass_distance': float(mass*proper/2)}


def photon(u, velocity=0., direction=-1., energy=3., mass=2.):
    omega = conformal(u)[0]
    covector = energy*np.array([-1., direction, 0., 0.])
    observer = np.array([1., velocity, 0., 0.])/np.sqrt(1-velocity**2)
    ej = -covector@observer
    ee = -covector@(observer/omega)
    return {'Jordan_energy': float(ej), 'Einstein_energy': float(ee),
            'Jordan_ratio': float(ej/mass), 'Einstein_ratio': float(ee/(mass/omega))}


def redshift():
    # Left-moving light: (t,x)=(1,.4) to (1.4,0), hence u=.6 to 1.4.
    emission, detection = photon(.6), photon(1.4)
    return {'raw_Einstein_energy_ratio': detection['Einstein_energy']/emission['Einstein_energy'],
            'dimensionless_frequency_ratio': detection['Einstein_ratio']/emission['Einstein_ratio'],
            'Jordan_frequency_ratio': detection['Jordan_ratio']/emission['Jordan_ratio']}


def eikonal_errors():
    momentum = np.array([-np.sqrt(2.**2+.7**2), .7, 0., 0.])
    sample = np.linspace(.1, 3., 41)
    omega = conformal(sample)[0]
    residual = (momentum@ETA@momentum)/omega**2+(2/omega)**2
    wrong = (momentum@ETA@momentum)/omega**2+2**2
    return float(np.max(abs(residual))), float(np.max(abs(wrong)))


def report():
    return {'round': 323,
            'scope': 'Classical positive-F frame redefinition with added massive phase counters and vacuum light rays. Phase counters use the point-particle action / leading WKB phase, not a derived atomic clock. Four-dimensional Maxwell density and massive eikonal are checked; quantum amplitudes, anomalies and microscopic clock ontology are outside scope.',
            'radar': radar(), 'photon_redshift': redshift(),
            'massive_eikonal_residual': eikonal_errors()[0],
            'incorrect_constant_mass_eikonal_residual': eikonal_errors()[1],
            'moving_detector': photon(1.4, velocity=.35),
            'next_interface': 'Test when multiple species admit a common matter metric; mass ratios are invariant under a common frame change.'}


class Checks(unittest.TestCase):
    def test_01_worldline_phase_is_frame_invariant(self):
        self.assertAlmostEqual(radar()['roundtrip_phase'], 1.6, places=13)

    def test_02_metric_proper_elapsed_is_not_the_phase(self):
        r = radar()
        self.assertAlmostEqual(r['Einstein_proper_elapsed'], r['Einstein_elapsed_primitive'], places=13)
        self.assertGreater(abs(r['Einstein_proper_elapsed']-r['Jordan_elapsed']), .05)

    def test_03_radar_distance_in_matter_units_agrees(self):
        for mass, separation in ((1., .2), (2., .4), (5., .7)):
            r = radar(mass, separation)
            self.assertAlmostEqual(r['dimensionless_distance'], r['Jordan_dimensionless_distance'], places=13)

    def test_04_fixed_units_control_is_distinguishable(self):
        self.assertGreater(abs(radar()['incorrect_fixed_mass_distance']-.8), .05)
        self.assertGreater(abs(redshift()['raw_Einstein_energy_ratio']-1), .05)

    def test_05_normalized_photon_reading_is_invariant(self):
        self.assertAlmostEqual(redshift()['dimensionless_frequency_ratio'], 1., places=13)
        for u in (.2, .6, 1.4, 3.):
            p = photon(u)
            self.assertAlmostEqual(p['Jordan_ratio'], p['Einstein_ratio'], places=13)

    def test_06_Doppler_effect_survives_frame_change(self):
        for velocity in (-.6, 0., .35):
            for direction in (-1., 1.):
                p = photon(.8, velocity, direction)
                expected = 1.5*(1-direction*velocity)/np.sqrt(1-velocity**2)
                self.assertAlmostEqual(p['Einstein_ratio'], expected, places=13)

    def test_07_massive_Hamilton_Jacobi_equation(self):
        correct, wrong = eikonal_errors()
        self.assertLess(correct, 3e-14)
        self.assertGreater(wrong, 10.)

    def test_08_vacuum_Maxwell_density_in_four_dimensions(self):
        field = np.array([[0., 1., .3, .2], [-1., 0., .5, -.4], [-.3, -.5, 0., .7], [-.2, .4, -.7, 0.]])
        original = np.einsum('ac,bd,ab,cd', ETA, ETA, field, field)
        for u in (.1, .5, 2.):
            omega = conformal(u)[0]
            metric = omega**2*ETA
            inverse = np.linalg.inv(metric)
            density = np.sqrt(-np.linalg.det(metric))*np.einsum('ac,bd,ab,cd', inverse, inverse, field, field)
            self.assertAlmostEqual(density, original, places=13)

    def test_09_radar_paths_are_null_in_both_frames(self):
        for t, x, direction in ((1.1, .1, 1.), (1.6, .2, -1.)):
            tangent = np.array([1., direction, 0., 0.])
            metric = conformal(t-x)[0]**2*ETA
            self.assertLess(abs(tangent@metric@tangent), 1e-15)


if __name__ == '__main__':
    main(__name__, 'conformal_readout_audit', report)
