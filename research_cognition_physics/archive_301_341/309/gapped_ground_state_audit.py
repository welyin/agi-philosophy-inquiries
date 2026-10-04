"""Round 309: dimerized ground states, independently controlled gap and quadrature."""
from functools import lru_cache
import unittest
import numpy as np
from growing_stream_audit import main
from ground_state_modular_audit import ground_projector, interval_correlation, slater_state, reduced_prefix, gaussian_density, modular_matrix


def dispersion(q, delta):
    return np.sqrt(np.cos(q)**2+delta**2*np.sin(q)**2)


def correlation_length_sites(delta):
    d = abs(delta)
    if not 0 < d < 1:
        raise ValueError('Use 0 < |delta| < 1 for a finite nonzero correlation length.')
    return 2/np.log((1+d)/(1-d))


def physical_matrix(n, delta, antiperiodic=False):
    if n % 2 or n < 4 or abs(delta) > 1:
        raise ValueError('Use an even chain of at least four sites and |delta| <= 1.')
    i = np.arange(1, n)
    bonds = -(1+np.where(i % 2 == 0, delta, -delta))/2
    h = np.diag(bonds, 1)+np.diag(bonds, -1)
    if antiperiodic:
        h[0, -1] = h[-1, 0] = (1+delta)/2
    return h


@lru_cache(maxsize=96)
def correlation(n, delta, quadrature=4096):
    if n % 2 or n < 2 or abs(delta) > 1 or quadrature < n:
        raise ValueError('Even region, |delta| <= 1, and enough momentum samples required.')
    if delta == 0:
        c = interval_correlation(n)
    else:
        q = 2*np.pi*(np.arange(quadrature)+0.5)/quadrature
        v, w = (1-delta)/2, (1+delta)/2
        f = (v+w*np.exp(-1j*q))/(2*abs(v+w*np.exp(-1j*q)))
        r = np.arange(-n//2, n//2+1)
        coefficients = (np.exp(1j*r[:, None]*q)@f/quadrature).real
        c = np.eye(n)*0.5
        for i in range(0, n, 2):
            for j in range(1, n, 2):
                c[i, j] = c[j, i] = coefficients[i//2-j//2+n//2]
    c.setflags(write=False)
    return c


def entropy_from_c(c):
    z = np.linalg.eigvalsh(c)
    if z.min() < -1e-12 or z.max() > 1+1e-12:
        raise ValueError('Correlation is not a positive contraction.')
    # Endpoint limits of the entropy are zero; this never constructs log(C).
    resolved = z[(z > 0) & (z < 1)]
    return float(-np.sum(resolved*np.log(resolved)+(1-resolved)*np.log1p(-resolved)))


def report():
    rows = []
    for delta in (0.02, 0.1, 0.3, -0.3):
        c = correlation(32, delta)
        rows.append({'delta': delta, 'band_gap': 2*abs(delta),
                     'correlation_length_in_sites': correlation_length_sites(delta),
                     'correlation_length_in_two_site_cells': correlation_length_sites(delta)/2,
                     'entropy_N32': entropy_from_c(c),
                     'quadrature_2048_4096_error_fro': float(np.linalg.norm(c-correlation(32, delta, 2048)))})
    reference = correlation(16, 0.01)
    quadrature = [{'momentum_points': m, 'error_fro': float(np.linalg.norm(correlation(16, 0.01, m)-reference))}
                  for m in (128, 256, 512, 1024)]
    c, occupied, _ = ground_projector(physical_matrix(8, 0.2))
    direct = reduced_prefix(slater_state(occupied), 3)
    rebuilt = gaussian_density(modular_matrix(c[:3, :3]))
    return {'round': 309,
            'scope': 'Explicit 1D dimerized free-fermion branch; mass/gap and cell geometry are inputs.',
            'gapped_cases': rows, 'quadrature_convergence': quadrature,
            'independent_fock_partial_trace_error': float(np.linalg.norm(direct-rebuilt)),
            'fully_dimerized_entropy': {'cut_strong_bonds_delta_plus_1': entropy_from_c(correlation(8, 1.)),
                                       'cut_weak_bonds_delta_minus_1': entropy_from_c(correlation(8, -1.))},
            'cut_convention': 'Site 1--2 bond is (1-delta)/2. Even interval cuts strong bonds for positive delta.',
            'endpoint_policy': 'Entropy uses its continuous endpoint value; unresolved full modular logarithms remain forbidden.'}


class Checks(unittest.TestCase):
    def test_01_gap_and_mass_sign(self):
        for d in (0.02, 0.1, 0.3):
            self.assertAlmostEqual(dispersion(np.pi/2, d), d)
            np.testing.assert_allclose(dispersion(np.arange(7)/7, d), dispersion(np.arange(7)/7, -d))

    def test_02_bloch_projector_against_finite_h(self):
        c, _, _ = ground_projector(physical_matrix(64, 0.2, True))
        np.testing.assert_allclose(c[:16, :16], correlation(16, 0.2, 32), atol=3e-14)

    def test_03_partial_trace_of_actual_ground_state(self):
        self.assertLess(report()['independent_fock_partial_trace_error'], 1e-12)

    def test_04_quadrature_convergence(self):
        errors = [r['error_fro'] for r in report()['quadrature_convergence']]
        self.assertTrue(all(a > b for a, b in zip(errors, errors[1:])))
        self.assertLess(errors[-1], 1e-9)

    def test_05_positive_contraction_and_parity(self):
        for d in (0.02, 0.1, -0.3):
            c = correlation(32, d); z = np.linalg.eigvalsh(c)
            self.assertGreater(z.min(), -1e-13); self.assertLess(z.max(), 1+1e-13)
            parity = np.diag((-1.)**np.arange(32))
            np.testing.assert_allclose(parity@c@parity, np.eye(32)-c, atol=1e-14)

    def test_06_critical_reference(self):
        np.testing.assert_allclose(correlation(16, 0), interval_correlation(16), atol=1e-15)
        self.assertLess(np.linalg.norm(correlation(16, 0.002)-correlation(16, 0)), 0.04)

    def test_07_identical_gap_different_cut_entropy(self):
        self.assertGreater(entropy_from_c(correlation(32, 0.3))-entropy_from_c(correlation(32, -0.3)), 0.9)

    def test_08_exact_dimer_limits(self):
        self.assertAlmostEqual(entropy_from_c(correlation(8, 1)), 2*np.log(2), places=11)
        self.assertLess(entropy_from_c(correlation(8, -1)), 1e-11)

    def test_09_complex_momentum_length(self):
        for d in (0.02, 0.1, 0.3):
            q = np.pi/2+1j/correlation_length_sites(d)
            self.assertLess(abs(np.cos(q)**2+d*d*np.sin(q)**2), 1e-14)


if __name__ == '__main__':
    main(__name__, 'gapped_ground_state_audit', report)
