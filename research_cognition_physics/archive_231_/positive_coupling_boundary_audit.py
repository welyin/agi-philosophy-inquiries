"""Round 321: a regular stealth crossing versus the degeneracy of the positive-F description."""
import unittest
import numpy as np
from growing_stream_audit import main


ETA = np.diag([-1., 1., 1., 1.])
NULL = np.array([1., 1., 0., 0.])/np.sqrt(2)


def stealth(lam, offset=0.):
    u = 1+np.asarray(lam)
    if np.any(u <= 0):
        raise ValueError('The scalar pole u=0 is outside this regular patch.')
    a = np.sqrt(6.)
    phi, dp, ddp = a/u+offset, -a/u**2, 2*a/u**3
    # Algebraically reduced expressions keep the undeformed F=0 crossing exact.
    f = 1-1/u**2-(2*a*offset/u+offset**2)/6
    df = -phi*dp/3
    ddf = -(dp*dp+phi*ddp)/3
    numerator = -2*a*offset/(3*u**3)
    return phi, dp, ddp, f, df, ddf, numerator


def flat_stress(lam):
    phi, dp, ddp, _, _, _, _ = stealth(lam)
    gradient = dp*NULL
    square_hessian = 2*(dp*dp+phi*ddp)*np.outer(NULL, NULL)
    squared_norm = float(gradient@ETA@gradient)
    box_square = float(np.einsum('ab,ab', ETA, square_hessian))
    return np.outer(gradient, gradient)-ETA*squared_norm/2+(ETA*box_square-square_hessian)/6


def shifted_crossing(offset):
    if not 0 < offset < np.sqrt(6.):
        raise ValueError('Use a small positive constant deformation.')
    return np.sqrt(6.)/(np.sqrt(6.)-offset)-1


def crossing_sensitivity(offset=0.01):
    cut = shifted_crossing(offset)
    at = stealth(cut, offset)
    limiting_residue = at[6]/at[4]
    rows = []
    for distance in (1e-2, 1e-3, 1e-4):
        phi, dp, ddp, f, df, ddf, numerator = stealth(cut+distance, offset)
        ricci = numerator/f
        rows.append({'distance_from_crossing': distance, 'F': float(f),
                     'required_Rkk': float(ricci), 'distance_times_Rkk': float(distance*ricci)})
    return {'offset': offset, 'shifted_crossing': float(cut),
            'nonzero_numerator_at_F_zero': float(at[6]),
            'limiting_pole_residue': float(limiting_residue), 'approach': rows}


def report():
    rows = []
    for lam in (-0.25, 0., 0.25):
        phi, dp, ddp, f, df, ddf, numerator = stealth(lam)
        rows.append({'lambda': lam, 'phi': float(phi), 'F': float(f), 'Wald_weight_rate': float(df),
                     'flat_stress_maximum': float(np.max(abs(flat_stress(lam)))),
                     'Einstein_metric_determinant': float(np.linalg.det(f*ETA))})
    # Local negative-F data: theta=shear=0, phi=5, phi'=1, xi=0.1, kappa=1.
    f = -1.5; z = -1/f
    theta_prime = -1/f-z*z
    return {'round': 321,
            'scope': 'Known massless conformal scalar stealth branch on flat spacetime, xi=1/6 and kappa=1; exact classical counterexample only. Regular Jordan background crossing does not establish perturbative stability, a positive graviton kinetic term or a well-posed crossing prescription.',
            'regular_stealth_crossing': rows,
            'negative_F_focusing_counterexample': {'F': f, 'theta': 0., 'shear_squared': 0.,
                                                   'phi': 5., 'phi_prime': 1., 'xi': 0.1,
                                                   'generalized_expansion_derivative': float(theta_prime)},
            'deformed_crossing': crossing_sensitivity(),
            'near_crossing_logarithmic_expansion': [{'lambda': v, 'Theta': float(stealth(v)[4]/stealth(v)[3])} for v in (0.1, 0.01, 0.001)],
            'next_interface': 'Operational meaning of conformal frames: transform massive probe matter as well as gravity, and distinguish a field redefinition from changing the matter coupling.'}


class Checks(unittest.TestCase):
    def test_01_full_improved_stress_vanishes(self):
        for lam in (-0.4, -0.25, 0., 0.25, 1.):
            self.assertLess(float(np.max(abs(flat_stress(lam)))), 2e-13)

    def test_02_massless_scalar_equation_is_satisfied(self):
        for lam in (-0.25, 0., 0.25):
            hessian = stealth(lam)[2]*np.outer(NULL, NULL)
            self.assertAlmostEqual(float(np.einsum('ab,ab', ETA, hessian)), 0., places=13)

    def test_03_background_fields_regular_across_F_zero(self):
        self.assertLess(stealth(-0.25)[3], 0.)
        self.assertEqual(stealth(0.)[3], 0.)
        self.assertGreater(stealth(0.25)[3], 0.)
        for lam in np.linspace(-0.25, 0.25, 21):
            self.assertTrue(all(np.isfinite(x) for x in stealth(lam)))

    def test_04_Wald_weight_regular_while_logarithm_fails(self):
        h = 1e-4
        slope = (stealth(-2*h)[3]-8*stealth(-h)[3]+8*stealth(h)[3]-stealth(2*h)[3])/(12*h)
        self.assertAlmostEqual(slope, stealth(0.)[4], places=10)
        self.assertAlmostEqual(slope, 2., places=10)

    def test_05_conformal_metric_is_degenerate_at_the_crossing(self):
        self.assertEqual(np.linalg.det(stealth(0.)[3]*ETA), 0.)
        with self.assertRaises(np.linalg.LinAlgError):
            np.linalg.inv(stealth(0.)[3]*ETA)
        self.assertEqual(np.linalg.det(ETA), -1.)

    def test_06_negative_F_can_reverse_the_focusing_sign(self):
        self.assertGreater(report()['negative_F_focusing_counterexample']['generalized_expansion_derivative'], 0.2)

    def test_07_generic_deformation_breaks_zero_over_zero_cancellation(self):
        data = crossing_sensitivity()
        self.assertLess(abs(stealth(data['shifted_crossing'], 0.01)[3]), 1e-14)
        self.assertGreater(abs(data['nonzero_numerator_at_F_zero']), 0.01)
        self.assertEqual(stealth(0.)[6], 0.)

    def test_08_required_curvature_has_a_resolved_simple_pole(self):
        data = crossing_sensitivity(); rows = data['approach']
        self.assertTrue(all(abs(b['required_Rkk']) > 9*abs(a['required_Rkk']) for a, b in zip(rows, rows[1:])))
        self.assertLess(abs(rows[-1]['distance_times_Rkk']/data['limiting_pole_residue']-1), 0.001)

    def test_09_logarithmic_expansion_can_diverge_with_flat_geometry(self):
        for lam in (0.01, 0.001):
            f, df = stealth(lam)[3:5]
            self.assertLess(abs(lam*df/f-1), 0.02)
        self.assertLess(float(np.max(abs(flat_stress(0.001)))), 1e-13)


if __name__ == '__main__':
    main(__name__, 'positive_coupling_boundary_audit', report)
