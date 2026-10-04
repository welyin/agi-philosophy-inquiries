"""Round 315: common proper-time split, contact accounting and unfixed gravitational normalization."""
import unittest
import numpy as np
from growing_stream_audit import main
from heat_kernel_area_audit import proper_moment

# (mass, nonminimal curvature coupling, number of real scalar species)
SPECIES = ((0., 0., 2), (0.7, 0.1, 1), (1.3, 1/6, 1), (2., 1/3, 1))


def shift(species, cutoff, upper=None):
    return float(sum(n*(1-6*xi)*proper_moment(m, cutoff, upper=upper)/(12*np.pi) for m, xi, n in species))


def scale_account(scale, cutoff=0.1, bare_inverse_g=2., species=SPECIES):
    running = bare_inverse_g+shift(species, cutoff, upper=scale)
    residual_statistical = sum(n*proper_moment(m, scale)/(48*np.pi) for m, xi, n in species)
    residual_contact = sum(-n*xi*proper_moment(m, scale)/(8*np.pi) for m, xi, n in species)
    total = running/4+residual_statistical+residual_contact
    target = (bare_inverse_g+shift(species, cutoff))/4
    return {'split_length': scale, 'running_inverse_G': float(running),
            'residual_statistical_area_coefficient': float(residual_statistical),
            'residual_contact_area_coefficient': float(residual_contact),
            'generalized_area_coefficient': float(total),
            'all_scales_target': float(target), 'split_error': float(abs(total-target)),
            'coefficient_if_contact_is_omitted': float(running/4+residual_statistical)}


def running_beta(scale, species=SPECIES):
    return float(sum(n*(1-6*xi)*np.exp(-(m*scale)**2)/(6*np.pi*scale**2) for m, xi, n in species))


def report():
    rows = [scale_account(ell) for ell in (0.1, 0.15, 0.25, 0.5, 1., 2.)]
    renormalization = []
    for eps in (0.2, 0.1, 0.05, 0.025):
        bare = 8.-shift(SPECIES, eps)
        account = scale_account(0.5, eps, bare)
        renormalization.append({'UV_cutoff': eps, 'adjusted_bare_inverse_G': bare,
                                'fixed_renormalized_inverse_G': 8.,
                                'generalized_area_coefficient': account['generalized_area_coefficient']})
    families = [{'bare_inverse_G': bare, 'renormalized_inverse_G': bare+shift(SPECIES, 0.1),
                 'split_error': scale_account(0.5, bare_inverse_g=bare)['split_error']} for bare in (0., 2., 5.)]
    return {'round': 315,
            'scope': 'One-loop local area/EH sector with a common proper-time split, free real scalars and stated contact terms. Not an exact Hilbert-space partial trace, full interacting RG flow, or a prediction of the measured G.',
            'species': [{'mass': m, 'xi': xi, 'multiplicity': n} for m, xi, n in SPECIES],
            'scale_accounts': rows, 'UV_renormalization_controls': renormalization,
            'unfixed_bare_coupling_family': families,
            'physical_scope': 'Curvature-squared, nonlocal, state-dependent finite entropy and outer-boundary terms remain outside this planar area-sector audit.',
            'next_interface': 'Local equilibrium and entropy production: determine which first-order balance is justified and what second-order terms remain.'}


class Checks(unittest.TestCase):
    def test_01_split_integral_composes(self):
        for m in (0., 0.7, 2.):
            direct = proper_moment(m, 0.1, upper=0.5)
            parts = proper_moment(m, 0.1, upper=0.25)+proper_moment(m, 0.25, upper=0.5)
            self.assertAlmostEqual(direct, parts, places=10)

    def test_02_generalized_coefficient_is_scale_independent(self):
        for ell in (0.1, 0.15, 0.25, 0.5, 1., 2.):
            self.assertLess(scale_account(ell)['split_error'], 1e-10)

    def test_03_dropping_contact_produces_drift(self):
        first, last = scale_account(0.1), scale_account(1.)
        self.assertGreater(abs(first['coefficient_if_contact_is_omitted']-last['coefficient_if_contact_is_omitted']), 2.)

    def test_04_differential_running_independent(self):
        step = 1e-4
        for ell in (0.2, 0.5, 1.):
            direct = (scale_account(ell*np.exp(step))['running_inverse_G']-scale_account(ell*np.exp(-step))['running_inverse_G'])/(2*step)
            self.assertLess(abs(direct-running_beta(ell)), 1e-6)

    def test_05_regulator_mismatch_detected(self):
        # Only minimal massless fields: using a different split in the entropy breaks cancellation.
        ell = 0.5; species = ((0., 0., 1),)
        correct = scale_account(ell, species=species)
        mismatched = correct['running_inverse_G']/4+proper_moment(0., ell*2)/(48*np.pi)
        self.assertGreater(abs(mismatched-correct['all_scales_target']), 0.01)

    def test_06_UV_counterterm_keeps_renormalized_value(self):
        for row in report()['UV_renormalization_controls']:
            self.assertAlmostEqual(row['generalized_area_coefficient'], 2., places=10)

    def test_07_matching_does_not_select_bare_coupling(self):
        rows = report()['unfixed_bare_coupling_family']
        self.assertAlmostEqual(rows[-1]['renormalized_inverse_G']-rows[0]['renormalized_inverse_G'], 5.)
        self.assertLess(max(r['split_error'] for r in rows), 1e-10)

    def test_08_threshold_running_decouples_massive_species(self):
        light, heavy = ((0., 0., 1),), ((2., 0., 1),)
        self.assertAlmostEqual(running_beta(1., heavy)/running_beta(1., light), np.exp(-4), places=13)

    def test_09_induced_scale_is_still_an_input(self):
        # Bare q=0 is an added boundary condition; in a massless minimal model G=12*pi*eps**2/N.
        for n in (1, 3, 7):
            for eps in (0.1, 0.2):
                calculated_g = 1/shift(((0., 0., n),), eps)
                self.assertAlmostEqual(calculated_g, 12*np.pi*eps**2/n, places=12)


if __name__ == '__main__':
    main(__name__, 'generalized_entropy_scale_audit', report)
