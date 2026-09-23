"""Round 311: half-chain modular state and a controlled massive continuum scaling."""
import unittest
import numpy as np
from growing_stream_audit import main
from gapped_ground_state_audit import correlation, correlation_length_sites
from massive_interval_geometry_audit import candidate, covariance_of_generator, ctm_coefficient, elliptic_k


def boundary_error(length, delta, prefix=8, scale=1.):
    k = scale*candidate(length, delta, 'halfline')
    predicted = covariance_of_generator(k)[:prefix, :prefix]
    return float(np.linalg.norm(predicted-correlation(prefix, delta)))


def scaled_dispersion_squared(spacing, mass, momentum):
    if abs(mass*spacing) >= 1 or abs(momentum*spacing) > 1:
        raise ValueError('This controlled expansion uses |m a| < 1 and |p a| <= 1.')
    return float(np.sin(spacing*momentum)**2/spacing**2+mass**2*np.cos(spacing*momentum)**2)


def continuum_case(spacing, mass=0.8, momentum=0.7):
    delta = mass*spacing
    exact = scaled_dispersion_squared(spacing, mass, momentum)
    target = mass**2+momentum**2
    return {'spacing': spacing, 'mass': mass, 'momentum': momentum, 'lattice_dimerization': delta,
            'energy_squared': exact, 'relativistic_target': target,
            'energy_squared_error': target-exact,
            'analytic_error_upper': spacing**2*(momentum**4/3+mass**2*momentum**2),
            'halfspace_coefficient_relative_error': ctm_coefficient(delta)/(2*np.pi)-1,
            'physical_correlation_length': spacing*correlation_length_sites(delta)}


def report():
    rows = []
    for d in (0.1, 0.3, -0.3):
        for n in (16, 32, 64, 128):
            rows.append({'delta': d, 'halfchain_truncation': n, 'observed_prefix_sites': 8,
                         'prefix_covariance_error_fro': boundary_error(n, d)})
    return {'round': 311,
            'scope': 'Reuse exact CTM half-chain result in an integrable 1D model; verify local state convergence and a separate low-momentum continuum expansion. Not a derivation of universal Lorentz symmetry.',
            'boundary_state_checks': rows,
            'normalization_control_delta_0_3_N128': [{'scale': scale, 'error': boundary_error(128, 0.3, scale=scale)} for scale in (0.5, 1., 2.)],
            'massive_continuum_scaling': [continuum_case(a) for a in (0.2, 0.1, 0.05, 0.025)],
            'known_ctm_identity': 'gamma(delta)=4 K(|delta|)=4 K(sqrt(1-r^2))/(1+|delta|), r=(1-|delta|)/(1+|delta|).',
            'limits_kept_separate': ['Fixed lattice gap, half-chain truncation to infinity, fixed observed prefix.',
                                    'Spacing to zero with delta=m*a and fixed physical momentum; neither alone proves the full QFT scaling limit.']}


class Checks(unittest.TestCase):
    def test_01_strong_cut_boundary_state(self):
        self.assertLess(boundary_error(128, 0.3), 2e-12)

    def test_02_weak_cut_boundary_state(self):
        self.assertLess(boundary_error(128, -0.3), 2e-12)

    def test_03_longer_correlation_length_needs_larger_truncation(self):
        values = [boundary_error(n, 0.1) for n in (16, 32, 64, 128)]
        self.assertTrue(all(a > b for a, b in zip(values, values[1:])))
        self.assertLess(values[-1], 1e-8)

    def test_04_wrong_normalization_detected(self):
        for scale in (0.5, 2.):
            self.assertGreater(boundary_error(128, 0.3, scale=scale), 1e-5)

    def test_05_landen_identity(self):
        for d in (0., 0.02, 0.1, 0.3):
            self.assertAlmostEqual(ctm_coefficient(d), 4*elliptic_k(d), places=12)

    def test_06_relativistic_dispersion_error_bound(self):
        for a in (0.2, 0.1, 0.05, 0.025):
            for m in (0.2, 0.8, 1.2):
                r = continuum_case(a, m)
                self.assertGreater(r['energy_squared_error'], 0)
                self.assertLessEqual(r['energy_squared_error'], r['analytic_error_upper']+1e-14)

    def test_07_quadratic_scaling(self):
        rows = [continuum_case(a) for a in (0.2, 0.1, 0.05, 0.025)]
        ratios = [a['energy_squared_error']/b['energy_squared_error'] for a, b in zip(rows, rows[1:])]
        self.assertTrue(all(3.9 < r < 4.1 for r in ratios))

    def test_08_boost_coefficient_small_mass_expansion(self):
        for d in (0.02, 0.04, 0.08):
            measured = (ctm_coefficient(d)/(2*np.pi)-1)/d**2
            self.assertLess(abs(measured-0.25), 0.001)

    def test_09_physical_correlation_length(self):
        rows = [continuum_case(a) for a in (0.2, 0.1, 0.05, 0.025)]
        errors = [abs(r['physical_correlation_length']-1/0.8) for r in rows]
        self.assertTrue(all(a > b for a, b in zip(errors, errors[1:])))
        self.assertLess(errors[-1], 0.0002)


if __name__ == '__main__':
    main(__name__, 'massive_halfspace_bridge_audit', report)
