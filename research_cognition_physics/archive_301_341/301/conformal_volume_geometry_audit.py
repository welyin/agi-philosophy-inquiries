"""Round 301: continuum conformal-volume calibration, not emergent geometry."""
import unittest
import numpy as np
from growing_stream_audit import main

MINKOWSKI = np.diag([-1., 1., 1., 1.])


def reconstruct_metric(base, volume_density):
    """A specified conformal representative plus a physical volume density."""
    base = np.asarray(base, dtype=float)
    if base.shape != (4, 4) or not np.allclose(base, base.T):
        raise ValueError('A symmetric four-dimensional representative is required.')
    if np.count_nonzero(np.linalg.eigvalsh(base) < 0) != 1 or np.linalg.det(base) >= 0:
        raise ValueError('Lorentz signature required.')
    if not np.isfinite(volume_density) or volume_density <= 0:
        raise ValueError('Positive finite physical volume density required.')
    return (volume_density / np.sqrt(-np.linalg.det(base)))**0.5 * base


def de_sitter_jet(eta, hubble=0.5):
    if eta >= 0 or hubble <= 0:
        raise ValueError('Expanding conformal patch requires eta < 0, H > 0.')
    return (-1/(hubble*eta), 1/(hubble*eta**2), -2/(hubble*eta**3))


def ricci_from_connection(lapse, scale):
    """Compute R_ab from Gamma and its derivative; jets are (f,f',f'')."""
    n, dn, ddn = lapse
    a, da, dda = scale
    gamma, derivative = np.zeros((4, 4, 4)), np.zeros((4, 4, 4))
    gamma[0, 0, 0] = dn/n
    derivative[0, 0, 0] = ddn/n-(dn/n)**2
    for i in range(1, 4):
        gamma[0, i, i] = a*da/n**2
        derivative[0, i, i] = (da**2+a*dda)/n**2-2*a*da*dn/n**3
        gamma[i, 0, i] = gamma[i, i, 0] = da/a
        derivative[i, 0, i] = derivative[i, i, 0] = dda/a-(da/a)**2
    ricci = np.zeros((4, 4))
    for mu in range(4):
        for nu in range(4):
            ricci[mu, nu] = derivative[0, mu, nu]
            if nu == 0:
                ricci[mu, nu] -= sum(derivative[lam, mu, lam] for lam in range(4))
            for lam in range(4):
                for sig in range(4):
                    ricci[mu, nu] += (gamma[lam, mu, nu]*gamma[sig, lam, sig]
                                      - gamma[sig, mu, lam]*gamma[lam, nu, sig])
    inverse_diag = np.array([-1/n**2, 1/a**2, 1/a**2, 1/a**2])
    return ricci, float(np.dot(inverse_diag, np.diag(ricci)))


def integrate(function, left, right):
    nodes, weights = np.polynomial.legendre.leggauss(48)
    points = (right-left)*nodes/2+(right+left)/2
    return float((right-left)/2*np.dot(weights, function(points)))


def causal_sample():
    rng = np.random.default_rng(301)
    points = rng.uniform(-0.4, 0.4, size=(256, 4))
    points[:, 0] = rng.uniform(-2., -1., 256)
    delta = points[None, :, :]-points[:, None, :]
    norm = np.einsum('...i,ij,...j->...', delta, MINKOWSKI, delta)
    flat = (delta[:, :, 0] > 0) & (norm < 0)
    a = -1/(0.5*points[:, 0])
    curved = (delta[:, :, 0] > 0) & (a[:, None]**2*norm < 0)
    # This uses the analytic conformal causal criterion, not geodesic intervals.
    return int(flat.sum()), int(np.count_nonzero(flat != curved))


def report():
    h = 0.5
    grid = np.linspace(-2, -1, 17)
    curvature, reconstruction, covariance = [], [], []
    for eta in grid:
        jet = de_sitter_jet(eta, h)
        a = jet[0]
        _, scalar = ricci_from_connection(jet, jet)
        curvature.append(abs(scalar-12*h*h))
        reconstruction.append(np.max(abs(reconstruct_metric(MINKOWSKI, a**4)-a*a*MINKOWSKI)))
        jac = np.diag([-eta, 1, 1, 1])
        base_z = jac.T @ MINKOWSKI @ jac
        w_z = a**4*(-eta)
        covariance.append(np.max(abs(reconstruct_metric(base_z, w_z)-jac.T@(a*a*MINKOWSKI)@jac)))
    pairs, disagreements = causal_sample()
    return {
        'round': 301,
        'scope': 'Given 3+1 Lorentz manifold, conformal class and exact physical volume; no derivation of dimension, volume law, emergent manifold or Einstein dynamics.',
        'inputs': {'spacetime_dimension': 4, 'H': h, 'eta_interval': [-2, -1],
                   'comoving_spatial_volume': 1, 'physical_density_known_for_count_to_volume': True},
        'causal_order': {'points': 256, 'timelike_ordered_pairs': pairs, 'disagreements': disagreements},
        'geometry': {'flat_R': 0., 'de_sitter_R': 12*h*h,
                     'max_connection_scalar_error': float(max(curvature)),
                     'max_metric_reconstruction_error': float(max(reconstruction)),
                     'max_coordinate_covariance_error': float(max(covariance)),
                     'flat_proper_time': 1.,
                     'de_sitter_proper_time': integrate(lambda t: -1/(h*t), -2, -1),
                     'flat_slab_volume': 1.,
                     'de_sitter_slab_volume': integrate(lambda t: 1/(h*t)**4, -2, -1),
                     'slab_volume_in_z': integrate(lambda z: np.exp(3*z)/h**4, -np.log(2), 0)},
        'count_degeneracy': {'same_poisson_intensity': 'rho_curved=1, rho_flat=a^4',
                             'max_intensity_difference': 0.,
                             'meaning': 'Identical intensity measure and order imply identical finite Poisson order laws; unknown density can hide curvature.'},
        'scale_ambiguity': {'metric_multiplier': 9., 'compensating_density_multiplier': 1/81,
                            'curvature_after_rescaling': 12*h*h/9},
    }


class GeometryTests(unittest.TestCase):
    def test_01_reconstruct_and_validate(self):
        for eta in np.linspace(-4, -0.5, 13):
            a = de_sitter_jet(eta)[0]
            np.testing.assert_allclose(reconstruct_metric(MINKOWSKI, a**4), a*a*MINKOWSKI)
        for volume in (0, -1, np.nan):
            with self.assertRaises(ValueError):
                reconstruct_metric(MINKOWSKI, volume)
        with self.assertRaises(ValueError):
            reconstruct_metric(np.eye(4), 1)

    def test_02_causal_order_and_null_directions(self):
        pairs, difference = causal_sample()
        self.assertGreater(pairs, 0)
        self.assertEqual(difference, 0)
        v = np.array([1., 0.6, 0.8, 0.])
        for eta in (-2, -1):
            a = de_sitter_jet(eta)[0]
            self.assertAlmostEqual(v @ (a*a*MINKOWSKI) @ v, 0)

    def test_03_connection_curvature_known_geometries(self):
        _, flat_r = ricci_from_connection((1, 0, 0), (1, 0, 0))
        self.assertEqual(flat_r, 0)
        for h in (0.2, 0.5, 1.1):
            for eta in (-3, -1.2, -0.7):
                jet = de_sitter_jet(eta, h)
                ricci, scalar = ricci_from_connection(jet, jet)
                np.testing.assert_allclose(ricci, 3*h*h*jet[0]**2*MINKOWSKI, atol=1e-13)
                self.assertAlmostEqual(scalar, 12*h*h)

    def test_04_coordinate_covariance_of_volume_reconstruction(self):
        for z in np.linspace(-np.log(2), 0, 7):
            jac = np.diag([np.exp(-z), 1, 1, 1])
            a = np.exp(z)/0.5
            expected = np.diag([-4, a*a, a*a, a*a])
            actual = reconstruct_metric(jac.T @ MINKOWSKI @ jac, a**4*np.exp(-z))
            np.testing.assert_allclose(actual, expected)

    def test_05_scalar_invariant_in_nonlinear_coordinates(self):
        for z in (-0.7, -0.3, 0.2):
            a, n = np.exp(z)/0.5, 2.
            _, curved_r = ricci_from_connection((n, 0, 0), (a, a, a))
            self.assertAlmostEqual(curved_r, 3.)
            n0 = np.exp(-z)
            _, flat_r = ricci_from_connection((n0, -n0, n0), (1, 0, 0))
            self.assertAlmostEqual(flat_r, 0.)

    def test_06_integrated_clocks_and_volumes(self):
        g = report()['geometry']
        self.assertAlmostEqual(g['de_sitter_proper_time'], 2*np.log(2))
        self.assertAlmostEqual(g['de_sitter_slab_volume'], 14/3)
        self.assertAlmostEqual(g['slab_volume_in_z'], 14/3)

    def test_07_density_curvature_degeneracy(self):
        for left, right in ((-3., -2.), (-2., -1.), (-1., -0.5)):
            curved_mean = integrate(lambda t: (-1/(0.5*t))**4, left, right)
            flat_variable_density_mean = integrate(lambda t: 16/t**4, left, right)
            self.assertEqual(curved_mean, flat_variable_density_mean)
        self.assertNotEqual(report()['geometry']['de_sitter_R'], 0.)

    def test_08_global_scale_and_density_ambiguity(self):
        jet = de_sitter_jet(-1.7)
        c = 3.
        _, r1 = ricci_from_connection(jet, jet)
        scaled = tuple(c*x for x in jet)
        _, r2 = ricci_from_connection(scaled, scaled)
        self.assertAlmostEqual(r2, r1/c**2)
        self.assertAlmostEqual((c*jet[0])**4/c**4, jet[0]**4)

    def test_09_density_is_not_a_coordinate_scalar(self):
        z = -np.log(2)
        a = np.exp(z)/0.5
        wrong = reconstruct_metric(MINKOWSKI, a**4*np.exp(-z))
        actual = np.diag([-4., a*a, a*a, a*a])
        self.assertGreater(np.max(abs(wrong-actual)), 1.)


if __name__ == '__main__':
    main(__name__, 'conformal_volume_geometry_audit', report)
