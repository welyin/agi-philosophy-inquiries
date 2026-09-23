"""Round 318: canonical null boost, improved stress and Wald cut charge in one convention."""
import math
import unittest
import numpy as np
from growing_stream_audit import main
from heat_kernel_area_audit import gauss


def profile(lam, center=0.7, width=0.6, amplitude=1.):
    if width <= 0:
        raise ValueError('Positive profile width required.')
    d = np.asarray(lam)-center
    phi = amplitude*np.exp(-d*d/(2*width**2))
    derivative = -d*phi/width**2
    square = phi**2
    square_prime = -2*d*square/width**2
    square_second = (4*d*d/width**4-2/width**2)*square
    return phi, derivative, square, square_prime, square_second


def canonical_integral_exact(cut, end, center=0.7, width=0.6, amplitude=1.):
    def primitive(lam):
        y = (lam-center)/width
        i3 = -(y*y+1)*np.exp(-y*y)/2
        i2 = np.sqrt(np.pi)*math.erf(y)/4-y*np.exp(-y*y)/2
        return i3+(center-cut)*i2/width
    return float(2*np.pi*amplitude**2*(primitive(end)-primitive(cut)))


def account(xi, cut=0., end=6., center=0.7, width=0.6, amplitude=1., order=128):
    if end <= cut:
        raise ValueError('Ordered null cuts required.')
    nodes, weights = gauss(order)
    lam = cut+(nodes+1)*(end-cut)/2
    _, derivative, square, square_prime, square_second = profile(lam, center, width, amplitude)
    measure = weights*(end-cut)/2*(lam-cut)*2*np.pi
    canonical = float(np.dot(measure, derivative**2))
    hilbert = float(np.dot(measure, derivative**2-xi*square_second))
    f0 = float(profile(cut, center, width, amplitude)[2])
    _, _, f1, df1, _ = profile(end, center, width, amplitude)
    endpoint = float((end-cut)*df1-f1)
    contact = -2*np.pi*xi*f0
    completion = hilbert-contact+2*np.pi*xi*endpoint
    return {'xi': xi, 'cut': cut, 'end': end, 'canonical_boost': canonical,
            'hilbert_boost': hilbert, 'wald_cut_charge': float(contact),
            'far_endpoint_B': endpoint, 'completed_boost': float(completion),
            'matching_error': abs(completion-canonical),
            'error_if_far_endpoint_is_omitted': abs(hilbert-contact-canonical),
            'error_if_cut_charge_is_omitted': abs(hilbert+2*np.pi*xi*endpoint-canonical)}


def report():
    rows = [account(xi) for xi in (-0.1, 0., 1/6, 0.4)]
    return {'round': 318,
            'scope': 'Flat null profile of a free massless real scalar, per unit transverse area. Analytic integration by parts plus independent quadrature checks; not a simulation of quantum entropy or a derivation of backreaction.',
            'long_interval_accounts': rows,
            'finite_endpoint_control': account(1/6, cut=0., end=1.2),
            'negative_improved_boost_control': account(0.8, cut=0.7),
            'canonical_integral_independent': canonical_integral_exact(0., 6.),
            'conditional_gravity_interface': 'Only if linearized geometric response is separately established, delta(A/4G)=-K_H at a cut relative to a settled future; adding the Wald cut charge then gives -K_can when the far endpoint vanishes.',
            'next_interface': 'Nonlinear Wald entropy focusing with a field-dependent curvature coefficient: test the positivity condition and boundary data explicitly; do not infer full quantum GSL from a classical balance.'}


class Checks(unittest.TestCase):
    def test_01_profile_second_derivative(self):
        h = 2e-4
        for lam in (0., 0.7, 1.2):
            f = lambda x: profile(x)[2]
            numeric = (-f(lam+2*h)+16*f(lam+h)-30*f(lam)+16*f(lam-h)-f(lam-2*h))/(12*h*h)
            self.assertLess(abs(numeric-profile(lam)[4]), 2e-8)

    def test_02_canonical_integral_independent_primitive(self):
        for cut, end in ((0., 6.), (0.25, 1.2), (0.7, 6.)):
            self.assertAlmostEqual(account(0., cut, end)['canonical_boost'], canonical_integral_exact(cut, end), places=11)

    def test_03_improvement_integrates_to_both_endpoints(self):
        for xi in (-0.1, 0., 1/6, 0.4):
            for cut, end in ((0., 6.), (0.25, 1.2), (0.7, 6.)):
                self.assertLess(account(xi, cut, end)['matching_error'], 1e-12)

    def test_04_wald_cut_cancels_long_interval_improvement(self):
        for row in report()['long_interval_accounts']:
            self.assertLess(abs(row['far_endpoint_B']), 1e-29)
            self.assertAlmostEqual(row['hilbert_boost']-row['wald_cut_charge'], row['canonical_boost'], places=11)

    def test_05_finite_far_boundary_cannot_be_dropped(self):
        row = report()['finite_endpoint_control']
        self.assertGreater(row['error_if_far_endpoint_is_omitted'], 1.)
        self.assertLess(row['matching_error'], 1e-12)

    def test_06_completed_generator_is_xi_independent(self):
        rows = report()['long_interval_accounts']
        values = [r['completed_boost'] for r in rows]
        self.assertLess(max(values)-min(values), 1e-12)
        self.assertGreater(max(r['hilbert_boost'] for r in rows)-min(r['hilbert_boost'] for r in rows), 0.1)

    def test_07_bulk_hilbert_boost_alone_need_not_be_positive(self):
        row = report()['negative_improved_boost_control']
        self.assertLess(row['hilbert_boost'], -1.)
        self.assertAlmostEqual(row['canonical_boost'], np.pi, places=11)
        self.assertLess(row['matching_error'], 1e-12)

    def test_08_affine_relabeling_does_not_change_the_charge(self):
        first = account(0.2, 0.2, 1.2)
        scale = 3.
        second = account(0.2, 0.2*scale, 1.2*scale, center=0.7*scale, width=0.6*scale)
        for key in ('canonical_boost', 'hilbert_boost', 'wald_cut_charge', 'far_endpoint_B'):
            self.assertAlmostEqual(first[key], second[key], places=11)

    def test_09_field_amplitude_scales_all_charges_quadratically(self):
        first, second = account(0.2), account(0.2, amplitude=0.3)
        for key in ('canonical_boost', 'hilbert_boost', 'wald_cut_charge'):
            self.assertAlmostEqual(second[key], 0.09*first[key], places=11)


if __name__ == '__main__':
    main(__name__, 'improved_null_energy_audit', report)
