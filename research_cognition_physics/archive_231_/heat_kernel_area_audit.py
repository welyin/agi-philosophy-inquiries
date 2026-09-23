"""Round 313: two independent geometries calibrate the same scalar one-loop coefficient."""
from functools import lru_cache
import math
import unittest
import numpy as np
from numpy.polynomial.legendre import leggauss
from growing_stream_audit import main


@lru_cache(maxsize=8)
def gauss(order):
    return leggauss(order)


def proper_moment(mass, cutoff, power=2, upper=None, order=192):
    """Integral ds s**(-power) exp(-m**2*s), cutoff**2 <= s <= upper**2."""
    if mass < 0 or cutoff <= 0 or power <= 1 or (upper is not None and upper < cutoff):
        raise ValueError('Nonnegative mass, positive cutoff, power>1, ordered endpoints required.')
    if upper == cutoff:
        return 0.
    if mass == 0:
        return float((cutoff**(2-2*power)-(0 if upper is None else upper**(2-2*power)))/(power-1))
    if mass*cutoff > 4:
        raise ValueError('This quadrature is calibrated only for mass*cutoff <= 4.')
    nodes, weights = gauss(order)
    limit = 45/(power-1) if upper is None else 2*np.log(upper/cutoff)
    t = (nodes+1)*limit/2
    integrand = np.exp(-(power-1)*t-(mass*cutoff)**2*np.exp(t))
    return float(cutoff**(2-2*power)*np.dot(weights, integrand)*limit/2)


def e1_small(x):
    if not 0 < x <= 1:
        raise ValueError('Convergent series check restricted to 0<x<=1.')
    term, total = -x, -x
    for k in range(2, 100):
        term *= -x/k
        add = term/k
        total += add
        if abs(add) < 1e-18:
            break
    return float(-np.euler_gamma-np.log(x)-total)


def cone_surface_coefficient(alpha):
    if alpha <= 0:
        raise ValueError('Positive cone opening ratio required.')
    return (1/alpha-alpha)/12


def cone_images_numeric(n, order=96):
    if n < 1 or int(n) != n:
        raise ValueError('Orbifold image check uses integer n>=1.')
    nodes, weights = gauss(order)
    radial = (nodes+1)*np.sqrt(40)/2
    integral = float(np.dot(weights, radial*np.exp(-radial**2))*np.sqrt(40)/2)
    return float(sum(integral/(2*np.sin(np.pi*k/n)**2) for k in range(1, n))/n)


def sphere_normalized_heat(time, radius=2., xi=0., maximum=None):
    if time <= 0 or radius <= 0:
        raise ValueError('Positive heat parameter and radius required.')
    if maximum is None:
        maximum = math.ceil(np.sqrt(45*radius**2/time))+8
    ell = np.arange(maximum+1)
    return float(time/radius**2*np.sum((2*ell+1)*np.exp(-time*ell*(ell+1)/radius**2))
                 *np.exp(-xi*2*time/radius**2))


def curvature_coefficient(time, radius=2., xi=0.):
    return (sphere_normalized_heat(time, radius, xi)-1)/(time*2/radius**2)


def richardson_curvature(time=0.01, xi=0.):
    return 2*curvature_coefficient(time/2, xi=xi)-curvature_coefficient(time, xi=xi)


def replica_action(alpha, mass, cutoff, area=1.):
    return -area*cone_surface_coefficient(alpha)*proper_moment(mass, cutoff)/(8*np.pi)


def replica_entropy_numeric(mass, cutoff, area=1., step=1e-3):
    action = lambda a: replica_action(a, mass, cutoff, area)
    derivative = (-action(1+2*step)+8*action(1+step)-8*action(1-step)+action(1-2*step))/(12*step)
    return float(derivative-action(1))


def inverse_g_shift(mass, cutoff, xi=0.):
    # W_R = -J*(1/6-xi)/(32*pi**2) int R = -delta(1/G)/(16*pi) int R.
    return (1/6-xi)*proper_moment(mass, cutoff)/(2*np.pi)


def report():
    rows = []
    for mass in (0., 0.3, 1., 3., 10.):
        j = proper_moment(mass, 0.2)
        ent = replica_entropy_numeric(mass, 0.2)
        q = inverse_g_shift(mass, 0.2)
        rows.append({'mass': mass, 'cutoff': 0.2, 'proper_integral_J': j,
                     'replica_entropy_per_area': ent, 'inverse_G_shift': q,
                     'matching_error': abs(ent-q/4)})
    return {'round': 313,
            'scope': 'Known one-loop free minimal scalar, supplied four-dimensional geometry and common proper-time regulator. Local Einstein-Hilbert/planar area sector only; no absolute G or gravitational dynamics derived.',
            'sphere_spectral_calibration': [{'heat_parameter': s, 'coefficient_R': curvature_coefficient(s)} for s in (0.08, 0.04, 0.02, 0.01)],
            'sphere_extrapolated_coefficient': richardson_curvature(),
            'cone_image_calibration': [{'orbifold_order': n, 'images_numeric': cone_images_numeric(n),
                                        'known_cone_coefficient': cone_surface_coefficient(1/n)} for n in (2, 3, 4, 7)],
            'massive_area_EH_matching': rows,
            'proper_parameter_warning': 'Heat-kernel proper time is a spectral integration parameter, not a microscopic clock model.',
            'IR_and_curvature_scope': 'Planar cone area term is explicit; the smooth-background calculation retains the local heat-kernel R term. Higher curvature and nonlocal finite terms are not included.'}


class Checks(unittest.TestCase):
    def test_01_massless_integral_independent(self):
        nodes, weights = gauss(128)
        t = (nodes+1)*22.5
        numerical = np.dot(weights, np.exp(-t))*22.5/0.2**2
        self.assertAlmostEqual(numerical, proper_moment(0, 0.2), places=10)

    def test_02_exponential_integral_identity(self):
        for mass in (0.3, 1., 3.):
            cutoff = 0.2; x = (mass*cutoff)**2
            independent = np.exp(-x)/cutoff**2-mass**2*e1_small(x)
            self.assertAlmostEqual(independent, proper_moment(mass, cutoff), places=9)

    def test_03_quadrature_resolution(self):
        for mass in (0.3, 1., 3., 10.):
            self.assertAlmostEqual(proper_moment(mass, 0.2, order=128), proper_moment(mass, 0.2), places=9)

    def test_04_orbifold_images(self):
        for n in (2, 3, 4, 7, 11):
            self.assertAlmostEqual(cone_images_numeric(n), cone_surface_coefficient(1/n), places=12)

    def test_05_sphere_heat_coefficient(self):
        self.assertLess(abs(richardson_curvature()-1/6), 3e-8)

    def test_06_sphere_spectral_tail(self):
        for s in (0.01, 0.04):
            self.assertAlmostEqual(sphere_normalized_heat(s), sphere_normalized_heat(s, maximum=400), places=13)

    def test_07_independent_entropy_and_curvature(self):
        spectral_q = richardson_curvature()*proper_moment(1., 0.2)/(2*np.pi)
        self.assertLess(abs(replica_entropy_numeric(1., 0.2)-spectral_q/4), 3e-8)

    def test_08_mass_decoupling(self):
        values = [proper_moment(m, 0.2) for m in (0., 0.3, 1., 3., 10.)]
        self.assertTrue(all(a > b for a, b in zip(values, values[1:])))
        self.assertLess(values[-1]/values[0], 0.004)

    def test_09_dimensional_scaling(self):
        for scale in (0.5, 2., 3.):
            self.assertAlmostEqual(proper_moment(0.7/scale, 0.2*scale)*scale**2, proper_moment(0.7, 0.2), places=11)


if __name__ == '__main__':
    main(__name__, 'heat_kernel_area_audit', report)
