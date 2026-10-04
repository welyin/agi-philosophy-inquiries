"""Round 306: fixed-volume small-ball geometry and conditional CFT entropy response."""
import unittest
import numpy as np
from growing_stream_audit import main

NODES, WEIGHTS = np.polynomial.legendre.leggauss(48)


def radial_scale(r, curvature):
    r = np.asarray(r)
    if curvature > 0:
        root = np.sqrt(curvature)
        return np.sin(root*r)/root
    if curvature < 0:
        root = np.sqrt(-curvature)
        return np.sinh(root*r)/root
    return r


def integrate_radially(function, radius):
    r = radius*(NODES+1)/2
    return float(radius/2*np.dot(WEIGHTS, function(r)))


def volume(radius, curvature):
    return 4*np.pi*integrate_radially(lambda r: radial_scale(r,curvature)**2, radius)


def area(radius, curvature):
    return float(4*np.pi*radial_scale(radius, curvature)**2)


def equal_volume_radius(ell, curvature):
    if ell <= 0 or abs(curvature)*ell*ell > 0.2:
        raise ValueError('This calibration is restricted to small positive-radius balls.')
    target = 4*np.pi*ell**3/3
    low, high = ell/2, 2*ell
    if not (volume(low,curvature) < target < volume(high,curvature)):
        raise ValueError('Root not bracketed.')
    for _ in range(70):
        middle = (low+high)/2
        if volume(middle,curvature) < target:
            low = middle
        else:
            high = middle
    return (low+high)/2


def ball_row(ell, curvature):
    radius = equal_volume_radius(ell, curvature)
    flat_area = 4*np.pi*ell*ell
    fixed_volume_delta = area(radius,curvature)-flat_area
    fixed_radius_delta = area(ell,curvature)-flat_area
    g00 = 3*curvature  # Time-symmetric slice: spatial R=6K=2G_uu.
    w = 4*np.pi*ell**4/15
    return {'ell': ell, 'curvature_K': curvature, 'fixed_volume_radius': radius,
            'relative_volume_error': abs(volume(radius,curvature)/(4*np.pi*ell**3/3)-1),
            'fixed_volume_area_change': fixed_volume_delta,
            'fixed_radius_area_change': fixed_radius_delta,
            'fixed_volume_response_ratio': float(-fixed_volume_delta/(w*g00)),
            'fixed_radius_response_ratio': float(-fixed_radius_delta/(w*g00)),
            'fixed_volume_G00_estimate': float(-fixed_volume_delta/w),
            'true_G00': g00}


def modular_energy(ell, density0=1., density2=0., hbar=1.):
    # Supplied CFT vacuum-ball modular kernel; not inferred from a finite state.
    integral = integrate_radially(lambda r: 4*np.pi*r*r*(ell*ell-r*r)/(2*ell)*
                                 (density0+density2*r*r),ell)
    return 2*np.pi*integral/hbar


def report():
    geometry = [ball_row(ell,0.3) for ell in (0.4,0.2,0.1,0.05)]
    gradient = []
    for ell in (0.4,0.2,0.1,0.05):
        uniform = modular_energy(ell)
        nonuniform = modular_energy(ell,density2=0.7)
        gradient.append({'ell':ell,'uniform_delta_K':uniform,
                         'relative_gradient_correction':(nonuniform-uniform)/uniform})
    # Compare derivatives at K=0, holding volume fixed, not finite area with an exact first law.
    derivative_rows = []
    ell = 0.2
    w = 4*np.pi*ell**4/15
    for step in (0.02,0.01,0.005):
        plus = ball_row(ell,step)['fixed_volume_area_change']
        minus = ball_row(ell,-step)['fixed_volume_area_change']
        derivative = (plus-minus)/(2*step)
        derivative_rows.append({'curvature_step':step,'area_derivative':derivative,
                                'linear_closure_relative_error':abs(derivative/(-3*w)-1)})
    return {'round':306,
            'scope':'Known Jacobson small-ball/CFT interface with supplied smooth 3+1 geometry, conformal vacuum kernel, area entropy and fixed-volume stationarity; no microscopic area-law derivation or extension to generic matter.',
            'geometry_rows':geometry,
            'geometry_error_halving_ratios':[(geometry[i]['fixed_volume_response_ratio']-1)/
                                              (geometry[i+1]['fixed_volume_response_ratio']-1) for i in range(3)],
            'fixed_volume_flat_derivative':derivative_rows,
            'inhomogeneous_matter_rows':gradient,
            'conditional_coupling':{'entropy_density_s':float(2*np.pi),'hbar':1.,
                                     'expected_beta_8piG':1.,
                                     'fixed_radius_wrong_beta':0.6,
                                     'fixed_radius_area_coefficient_factor':5/3},
            'distinguished_assumptions':['CFT vacuum modular Hamiltonian is geometric',
                                         'UV entropy variation equals s times area variation',
                                         'total entropy stationary at fixed volume',
                                         'small curvature and slowly varying matter',
                                         'first-order variations around a specified reference vacuum']}


class SmallBallTests(unittest.TestCase):
    def test_01_flat_geometry_and_area_volume_derivative(self):
        for ell in (0.1,0.4):
            self.assertAlmostEqual(volume(ell,0),4*np.pi*ell**3/3)
            self.assertAlmostEqual(area(ell,0),4*np.pi*ell*ell)
            step=1e-5
            for k in (-0.3,0.3):
                finite=(volume(ell+step,k)-volume(ell-step,k))/(2*step)
                self.assertLess(abs(finite-area(ell,k)),1e-8)

    def test_02_volume_preservation_positive_and_negative_curvature(self):
        for curvature in (-0.3,0.3):
            for ell in (0.1,0.4):
                row=ball_row(ell,curvature)
                self.assertLess(row['relative_volume_error'],1e-13)
                self.assertEqual(np.sign(row['fixed_volume_radius']-ell),np.sign(curvature))

    def test_03_small_ball_area_coefficient_and_second_order_scale_error(self):
        rows=report()['geometry_rows']
        self.assertLess(abs(rows[-1]['fixed_volume_response_ratio']-1),1e-3)
        for ratio in report()['geometry_error_halving_ratios']:
            self.assertTrue(3.8 < ratio < 4.2,ratio)

    def test_04_fixed_radius_produces_wrong_coupling(self):
        for curvature in (-0.3,0.3):
            row=ball_row(0.05,curvature)
            self.assertLess(abs(row['fixed_radius_response_ratio']-5/3),1e-3)
            self.assertGreater(abs(row['fixed_radius_response_ratio']-1),0.6)

    def test_05_modular_energy_constant_density_integral(self):
        for ell in (0.1,0.4):
            for density in (-0.2,1.,3.):
                expected=2*np.pi*(4*np.pi*ell**4/15)*density
                self.assertAlmostEqual(modular_energy(ell,density),expected)

    def test_06_nonuniform_density_remainder(self):
        for ell in (0.1,0.4):
            actual=modular_energy(ell,1.,0.7)
            exact=2*np.pi*4*np.pi*(ell**4/15+0.7*ell**6/35)
            self.assertAlmostEqual(actual,exact)
        rows=report()['inhomogeneous_matter_rows']
        for i in range(3):
            self.assertAlmostEqual(rows[i]['relative_gradient_correction']/rows[i+1]['relative_gradient_correction'],4.)

    def test_07_first_variation_closes_with_given_entropy_density(self):
        for row in report()['fixed_volume_flat_derivative']:
            self.assertLess(row['linear_closure_relative_error'],1e-6)
        ell=0.2
        w=4*np.pi*ell**4/15
        for entropy_density in (1.,2*np.pi,10.):
            beta=2*np.pi/entropy_density
            curvature_derivative=beta/3
            area_derivative=-3*w*curvature_derivative
            self.assertAlmostEqual(entropy_density*area_derivative+modular_energy(ell),0.)

    def test_08_finite_curvature_requires_corrections(self):
        row=ball_row(0.4,0.3)
        self.assertGreater(abs(row['fixed_volume_response_ratio']-1),1e-3)
        # Defining sA alone does not enforce stationarity for an arbitrary matter perturbation.
        s=2*np.pi
        mismatch=s*row['fixed_volume_area_change']+modular_energy(0.4,density0=0.)
        self.assertGreater(abs(mismatch),1e-3)

    def test_09_calibration_domain_and_scale_power(self):
        with self.assertRaises(ValueError):
            equal_volume_radius(2.,1.)
        for ell in (0.1,0.2):
            self.assertAlmostEqual(modular_energy(2*ell)/modular_energy(ell),16.)


if __name__ == '__main__':
    main(__name__,'fixed_volume_entanglement_audit',report)
