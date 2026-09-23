"""Round 327: a scalar force survives universal spin-2 coupling; its range is explicit."""
import unittest
import numpy as np
from growing_stream_audit import main

KAPPA = 1.
G = KAPPA/(8*np.pi)
SOURCE_MASS = .05
SOURCE_ALPHA = .6
ALPHAS = (.3, .8)


def scalar_profile(r, mu, alpha=SOURCE_ALPHA, mass=SOURCE_MASS):
    return -alpha*mass*np.exp(-mu*np.asarray(r))/(4*np.pi*np.asarray(r))


def potential(r, mu, alpha):
    # Potential per unit target inertial mass, to linear order in the weak fields.
    return -G*SOURCE_MASS/np.asarray(r)+alpha*scalar_profile(r, mu)


def acceleration(r, mu, alpha):
    r = np.asarray(r)
    shape = (1+mu*r)*np.exp(-mu*r)
    relative = alpha*SOURCE_ALPHA/(4*np.pi*G)
    return -G*SOURCE_MASS/r**2*(1+relative*shape)


def eotvos(r, mu, alphas=ALPHAS):
    a, b = [acceleration(r, mu, value) for value in alphas]
    return 2*np.abs(a-b)/np.abs(a+b)


def radial_equation_error(mu=.7, step=.001):
    r = np.linspace(1., 3., 65)
    f = scalar_profile(r, mu)
    d1 = (scalar_profile(r-2*step, mu)-8*scalar_profile(r-step, mu)
          +8*scalar_profile(r+step, mu)-scalar_profile(r+2*step, mu))/(12*step)
    d2 = (-scalar_profile(r+2*step, mu)+16*scalar_profile(r+step, mu)-30*f
          +16*scalar_profile(r-step, mu)-scalar_profile(r-2*step, mu))/(12*step**2)
    return float(np.max(abs(d2+2*d1/r-mu**2*f)))


def report():
    return {'round': 327,
            'scope': 'Leading weak-field static limit of Einstein gravity plus a canonically normalized scalar with quadratic mass and species-dependent m_i(phi). The spin-2 Newton term is universal; scalar charges are independent inputs. Point-source exterior and structureless targets only, no binding-energy composition model, nonlinear screening, radiative stability or experimental fit.',
            'normalization': {'kappa_g': KAPPA, 'G': G, 'source_mass': SOURCE_MASS,
                              'source_alpha': SOURCE_ALPHA, 'target_alphas': list(ALPHAS)},
            'range_scan_at_r_1': [{'mu_r': x, 'Eotvos_parameter': float(eotvos(1., x))} for x in (0., 1., 5., 10., 20.)],
            'massless_accelerations_at_r_1': [float(acceleration(1., 0., a)) for a in ALPHAS],
            'equal_scalar_charge_control': float(eotvos(1., 0., (.3, .3))),
            'radial_Yukawa_equation_error': radial_equation_error(),
            'maximum_Newton_potential_in_test_domain': G*SOURCE_MASS,
            'maximum_target_log_mass_shift_in_test_domain': float(max(ALPHAS)*abs(scalar_profile(1., 0.))),
            'next_interface': 'Test the long-distance Einstein limit by integrating out the massive scalar with explicit error bounds; record the scalar mass gap and coupling alignment as selection inputs, not conclusions of conservation.'}


class Checks(unittest.TestCase):
    def test_01_Yukawa_equation_outside_source(self):
        self.assertLess(radial_equation_error(), 1e-10)

    def test_02_point_source_flux_normalization(self):
        r = 1e-5; mu = .7
        derivative = SOURCE_ALPHA*SOURCE_MASS*np.exp(-mu*r)*(1+mu*r)/(4*np.pi*r*r)
        self.assertAlmostEqual(4*np.pi*r*r*derivative, SOURCE_ALPHA*SOURCE_MASS, places=11)

    def test_03_force_is_gradient_of_the_same_potential(self):
        h = 1e-4
        for alpha in ALPHAS:
            for r in (1., 2., 3.):
                derivative = (potential(r-2*h,.7,alpha)-8*potential(r-h,.7,alpha)
                              +8*potential(r+h,.7,alpha)-potential(r+2*h,.7,alpha))/(12*h)
                self.assertAlmostEqual(-derivative, acceleration(r,.7,alpha), places=12)

    def test_04_universal_Newton_term_does_not_remove_scalar_difference(self):
        self.assertAlmostEqual(float(eotvos(1., 0.)), 1.2/3.32, places=14)
        self.assertGreater(abs(acceleration(1.,0.,.3)-acceleration(1.,0.,.8)), .001)

    def test_05_equal_scalar_charges_restore_this_test(self):
        self.assertEqual(float(eotvos(1., 0., (.3, .3))), 0.)
        self.assertNotEqual(float(acceleration(1., 0., .3)), -G*SOURCE_MASS)

    def test_06_large_mass_suppresses_long_distance_violation(self):
        x = np.linspace(0., 20., 101)
        values = eotvos(1., x)
        self.assertTrue(np.all(np.diff(values)<0))
        self.assertLess(float(values[-1]), 3e-8)

    def test_07_massless_and_short_distance_limits(self):
        self.assertAlmostEqual(float(eotvos(1., 1e-6)), float(eotvos(1., 0.)), places=11)
        self.assertAlmostEqual(float(eotvos(2., .5)), float(eotvos(1., 1.)), places=14)

    def test_08_scalar_force_obeys_pair_reciprocity(self):
        a, b, mass_a, mass_b, r, mu = .3, .8, .05, .07, 2., .4
        pair = -a*b*mass_a*mass_b*np.exp(-mu*r)/(4*np.pi*r)
        source_a = mass_b*b*scalar_profile(r,mu,alpha=a,mass=mass_a)
        source_b = mass_a*a*scalar_profile(r,mu,alpha=b,mass=mass_b)
        self.assertAlmostEqual(pair, source_a, places=15)
        self.assertAlmostEqual(pair, source_b, places=15)

    def test_09_weak_field_regime_is_explicit(self):
        self.assertLess(G*SOURCE_MASS, .003)
        self.assertLess(max(ALPHAS)*abs(scalar_profile(1.,0.)), .003)


if __name__ == '__main__':
    main(__name__, 'scalar_fifth_force_audit', report)
