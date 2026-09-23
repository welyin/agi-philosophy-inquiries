"""Round 320: classical conformal frame change including null affine normalization."""
from functools import lru_cache
import unittest
import numpy as np
from growing_stream_audit import main
from heat_kernel_area_audit import gauss
from wald_focusing_audit import solve, derivative


def kinetic(phi, xi=0.8, kappa=1.):
    f = 1-kappa*xi*np.asarray(phi)**2
    if np.any(f <= 0) or kappa <= 0:
        raise ValueError('This real Einstein-frame map requires F>0 and kappa>0.')
    f_phi = -2*kappa*xi*np.asarray(phi)
    return 1/f+3/(2*kappa)*(f_phi/f)**2


def canonical_scalar(phi, xi=0.8, kappa=1., order=48):
    nodes, weights = gauss(order)
    values = np.asarray(phi)[..., None]*(nodes+1)/2
    return np.asarray(phi)/2*np.sum(weights*np.sqrt(kinetic(values, xi, kappa)), axis=-1)


def metric_connection(data):
    n = len(data['lambda'])
    g = np.zeros((n, 4, 4)); dg = np.zeros((n, 4, 4, 4))
    f, df = data['F'], data['dF']
    g[:, 0, 1] = g[:, 1, 0] = -f
    dg[:, 0, 0, 1] = dg[:, 0, 1, 0] = -df
    for coordinate, sign in ((2, 1), (3, -1)):
        g[:, coordinate, coordinate] = f*data['area']*np.exp(sign*2*data['beta'])
        dg[:, 0, coordinate, coordinate] = g[:, coordinate, coordinate]*(df/f+data['theta']+sign*2*data['dbeta'])
    term = dg+dg.swapaxes(1, 2)-dg.transpose(0, 2, 3, 1)
    gamma = np.einsum('nrs,nabs->nrab', np.linalg.inv(g), term)/2
    return g, gamma


@lru_cache(maxsize=12)
def audit(steps=512, xi=0.8):
    d = solve(steps, xi=xi); f = d['F']; g, gamma = metric_connection(d)
    h = 1/steps
    matrix = gamma[:, :, 0, :]
    derivative_part = derivative(gamma[:, 0, 0, 0]-np.trace(matrix, axis1=1, axis2=2), h)
    product = (np.einsum('nrrs,ns->n', gamma, gamma[:, :, 0, 0])
               -np.einsum('nrs,nsr->n', matrix, matrix))[2:-2]
    ricci_affine = (derivative_part+product)/f[2:-2]**2
    kin = kinetic(d['phi'], xi)
    energy_affine = kin*d['dphi']**2/f**2
    expansion = d['Theta']/f
    residual = derivative(expansion, h)/f[2:-2]+(expansion**2/2+2*d['dbeta']**2/f**2+energy_affine)[2:-2]
    affine_acceleration = -d['dF']/f**3+gamma[:, 0, 0, 0]/f**2
    scalar = canonical_scalar(d['phi'], xi)
    return {'area_identity_error': float(np.max(abs(np.sqrt(g[:, 2, 2]*g[:, 3, 3])-d['W']))),
            'affine_acceleration_error': float(np.max(abs(affine_acceleration))),
            'wrong_affine_acceleration': float(np.max(abs(gamma[:, 0, 0, 0]))),
            'curvature_equation_error': float(np.max(abs(ricci_affine-energy_affine[2:-2]))),
            'focusing_equation_error': float(np.max(abs(residual))),
            'canonical_scalar_derivative_error': float(np.max(abs(derivative(scalar, h)-(np.sqrt(kin)*d['dphi'])[2:-2]))),
            'missing_kinetic_correction': float(np.max(abs((kin-1/f)*d['dphi']**2/f**2))),
            'minimum_kinetic_coefficient': float(np.min(kin)),
            'maximum_kinetic_coefficient': float(np.max(kin)),
            'minimum_Einstein_expansion': float(np.min(expansion))}


def report():
    return {'round': 320,
            'scope': 'Positive-F classical field redefinition g_E=F*g_J in the same supplied scalar-tensor plane-wave solution. Area, affine generators and scalar kinetic coefficient checked jointly; no quantum frame equivalence or full global horizon is established.',
            'checks_at_512_steps': audit(),
            'independent_connection_convergence': [{'steps': n, 'Rkk_equation_error': audit(n)['curvature_equation_error']} for n in (128, 256, 512)],
            'minimal_coupling_control': audit(256, xi=0.),
            'next_interface': 'Audit the F=0 boundary without confusing failure of this field redefinition with a mandatory curvature singularity.'}


class Checks(unittest.TestCase):
    def test_01_Wald_weighted_area_is_Einstein_area(self):
        self.assertLess(audit()['area_identity_error'], 1e-13)

    def test_02_affine_generator_must_be_rescaled(self):
        self.assertLess(audit()['affine_acceleration_error'], 1e-13)
        self.assertGreater(audit()['wrong_affine_acceleration'], 1.)

    def test_03_field_equation_from_Einstein_metric_connection(self):
        self.assertLess(audit()['curvature_equation_error'], 2e-5)

    def test_04_connection_difference_converges(self):
        errors = [audit(n)['curvature_equation_error'] for n in (128, 256, 512)]
        self.assertTrue(all(12 < a/b < 20 for a, b in zip(errors, errors[1:])))

    def test_05_affine_Raychaudhuri_matches_canonical_scalar(self):
        self.assertLess(audit()['focusing_equation_error'], 2e-5)
        self.assertGreaterEqual(audit()['minimum_Einstein_expansion'], -1e-12)

    def test_06_scalar_redefinition_is_independently_integrated(self):
        self.assertLess(audit()['canonical_scalar_derivative_error'], 2e-6)

    def test_07_omitting_induced_kinetic_energy_fails(self):
        self.assertGreater(audit()['missing_kinetic_correction'], 1.)
        self.assertGreater(audit()['minimum_kinetic_coefficient'], 0.)

    def test_08_minimal_coupling_map_is_identity(self):
        row = audit(256, xi=0.)
        self.assertEqual(row['wrong_affine_acceleration'], 0.)
        self.assertEqual(row['missing_kinetic_correction'], 0.)
        np.testing.assert_allclose(canonical_scalar(np.array([-0.5, 0., 0.5]), xi=0.), [-0.5, 0., 0.5], atol=1e-13)

    def test_09_nonpositive_F_not_an_invertible_positive_frame(self):
        with self.assertRaises(ValueError):
            kinetic(1., xi=1.)
        with self.assertRaises(ValueError):
            kinetic(2., xi=1.)


if __name__ == '__main__':
    main(__name__, 'einstein_frame_optics_audit', report)
