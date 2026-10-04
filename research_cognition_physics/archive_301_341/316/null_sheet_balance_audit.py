"""Round 316: exact null-sheet area balance before imposing a matter field equation."""
import unittest
import numpy as np
from growing_stream_audit import main
from heat_kernel_area_audit import gauss


def transverse(lam, tide, end=1.):
    """Positive Rosen scale with a(end)=1, a'(end)=0 and a''=-tide*a."""
    distance = end-np.asarray(lam)
    if np.any(distance < 0):
        raise ValueError('This patch is restricted to lambda <= end.')
    if tide > 0:
        z = np.sqrt(tide)*distance
        if np.any(z >= np.pi/2):
            raise ValueError('A caustic would invalidate this positive transverse patch.')
        return np.cos(z), np.sqrt(tide)*np.sin(z)
    if tide < 0:
        z = np.sqrt(-tide)*distance
        return np.cosh(z), -np.sqrt(-tide)*np.sinh(z)
    return np.ones_like(distance), np.zeros_like(distance)


def optics(lam, tx, ty, end=1.):
    a, da = transverse(lam, tx, end)
    b, db = transverse(lam, ty, end)
    px, py = da/a, db/b
    return {'area': a*b, 'theta': px+py, 'shear_squared': (px-py)**2/2,
            'ricci_kk': tx+ty}


def connection(lam, tx, ty, end=1.):
    a, da = transverse(lam, tx, end)
    b, db = transverse(lam, ty, end)
    metric = np.zeros((4, 4)); metric[0, 1] = metric[1, 0] = -1.
    metric[2, 2] = a*a; metric[3, 3] = b*b
    derivative = np.zeros((4, 4, 4))
    derivative[0, 2, 2] = 2*a*da; derivative[0, 3, 3] = 2*b*db
    inverse = np.linalg.inv(metric)
    gamma = np.zeros((4, 4, 4))
    for r in range(4):
        for m in range(4):
            for n in range(4):
                gamma[r, m, n] = sum(inverse[r, s]*(derivative[m, n, s]+derivative[n, m, s]-derivative[s, m, n])/2 for s in range(4))
    return gamma


def curvature_from_connection(lam, tx, ty, end=1., step=1e-4):
    gamma = connection(lam, tx, ty, end)
    derivative = (connection(lam-2*step, tx, ty, end)-8*connection(lam-step, tx, ty, end)
                  +8*connection(lam+step, tx, ty, end)-connection(lam+2*step, tx, ty, end))/(12*step)
    # Only lambda derivatives are nonzero. R_00 is evaluated from the connection definition.
    return float(derivative[0, 0, 0]-sum(derivative[r, 0, r] for r in range(4))
                 +sum(gamma[r, r, s]*gamma[s, 0, 0]-gamma[r, 0, s]*gamma[s, 0, r]
                      for r in range(4) for s in range(4)))


def balance(tx, ty, cut=1., end=1., order=96):
    if not 0 < cut <= end:
        raise ValueError('Require 0 < cut <= final equilibrium slice.')
    nodes, weights = gauss(order)
    lam = (nodes+1)*cut/2
    data = optics(lam, tx, ty, end)
    measure = weights*cut/2*lam*data['area']
    ricci = float(np.sum(measure)*(tx+ty))
    shear = float(np.dot(measure, data['shear_squared']))
    expansion = float(-np.dot(measure, data['theta']**2)/2)
    initial, final = optics(0., tx, ty, end), optics(cut, tx, ty, end)
    boundary = float(cut*final['area']*final['theta'])
    delta = float(final['area']-initial['area'])
    return {'tx': tx, 'ty': ty, 'cut': cut, 'end': end,
            'area_change': delta, 'ricci_integral': ricci, 'shear_integral': shear,
            'expansion_integral': expansion, 'endpoint_term': boundary,
            'balance_error': abs(delta-ricci-shear-expansion-boundary),
            'error_if_only_ricci_is_kept': abs(delta-ricci)}


def report():
    ricci = []
    vacuum = []
    for amplitude in (0.4, 0.2, 0.1, 0.05):
        r = balance(amplitude/2, amplitude/2)
        r['amplitude'] = amplitude
        r['first_order_area'] = amplitude/2
        r['first_order_error'] = abs(r['area_change']-amplitude/2)
        ricci.append(r)
        w = balance(amplitude, -amplitude)
        w['amplitude'] = amplitude
        w['area_over_amplitude_squared'] = w['area_change']/amplitude**2
        vacuum.append(w)
    return {'round': 316,
            'scope': 'Supplied four-dimensional Rosen null-sheet metric with constant transverse tidal coefficients; purely geometric identities, no matter Einstein equation, quantum entropy or global event-horizon construction.',
            'ricci_perturbations': ricci, 'ricci_flat_shear_perturbations': vacuum,
            'finite_slice_control': balance(0.3, 0.1, cut=0.6),
            'independent_connection_Rkk': curvature_from_connection(0.3, 0.3, 0.1),
            'vacuum_quadratic_limit': 1/6,
            'interpretation': 'First-order balance requires controlled endpoint and nonlinear optical terms. Shear, area corrections and quantum relative entropy are different objects; no universal nonnegative entropy production is assigned to their sum.'}


class Checks(unittest.TestCase):
    def test_01_independent_metric_curvature(self):
        for tx, ty in ((0.3, 0.1), (0.2, -0.2), (-0.1, -0.3)):
            self.assertAlmostEqual(curvature_from_connection(0.3, tx, ty), tx+ty, places=10)

    def test_02_raychaudhuri_from_scale_derivatives(self):
        lam = 0.3; h = 1e-5
        for tx, ty in ((0.3, 0.1), (0.2, -0.2)):
            d = optics(lam, tx, ty)
            slope = (optics(lam+h, tx, ty)['theta']-optics(lam-h, tx, ty)['theta'])/(2*h)
            self.assertLess(abs(slope+d['theta']**2/2+d['shear_squared']+tx+ty), 1e-10)

    def test_03_exact_finite_area_balance(self):
        for tx, ty in ((0.3, 0.1), (0.2, -0.2), (-0.1, -0.3)):
            for cut in (0.4, 0.8, 1.):
                self.assertLess(balance(tx, ty, cut)['balance_error'], 1e-13)

    def test_04_endpoint_cannot_be_omitted_at_an_arbitrary_cut(self):
        full, part = balance(0.3, 0.1), balance(0.3, 0.1, 0.6)
        self.assertEqual(full['endpoint_term'], 0.)
        self.assertGreater(abs(part['endpoint_term']), 0.05)
        self.assertLess(part['balance_error'], 1e-13)

    def test_05_first_order_ricci_limit_has_quadratic_error(self):
        rows = report()['ricci_perturbations']
        for first, second in zip(rows, rows[1:]):
            self.assertTrue(3.9 < first['first_order_error']/second['first_order_error'] < 4.1)
        self.assertLess(abs(rows[-1]['area_change']/0.05-0.5), 0.005)

    def test_06_ricci_flat_shear_changes_area(self):
        rows = report()['ricci_flat_shear_perturbations']
        for row in rows:
            self.assertEqual(row['ricci_integral'], 0.)
            self.assertGreater(row['area_change'], 0.)
        self.assertLess(abs(rows[-1]['area_over_amplitude_squared']-1/6), 2e-6)

    def test_07_shear_omission_is_a_detectable_failure(self):
        row = balance(0.4, -0.4)
        self.assertGreater(row['shear_integral'], 0.02)
        self.assertGreater(row['error_if_only_ricci_is_kept'], 0.02)
        self.assertLess(row['balance_error'], 1e-13)

    def test_08_expansion_is_a_separate_signed_correction(self):
        row = balance(0.3, 0.3)
        self.assertEqual(row['shear_integral'], 0.)
        self.assertLess(row['expansion_integral'], -0.005)
        self.assertAlmostEqual(row['area_change']-row['ricci_integral'], row['expansion_integral'])

    def test_09_affine_rescaling_preserves_balance(self):
        initial = balance(0.3, 0.1, 0.6)
        scale = 3.
        rescaled = balance(0.3/scale**2, 0.1/scale**2, 0.6*scale, scale)
        for key in ('area_change', 'ricci_integral', 'shear_integral', 'expansion_integral', 'endpoint_term'):
            self.assertAlmostEqual(initial[key], rescaled[key], places=12)

    def test_10_caustic_patch_is_rejected(self):
        with self.assertRaises(ValueError):
            balance(3., 0.)


if __name__ == '__main__':
    main(__name__, 'null_sheet_balance_audit', report)
