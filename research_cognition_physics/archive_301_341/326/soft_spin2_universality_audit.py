"""Round 326: conditional leading soft spin-2 Ward constraint on connected species."""
import unittest
import numpy as np
from growing_stream_audit import main
from conformal_probe_motion_audit import ETA

SIGNS = np.array([-1., -1., 1., 1.])


def scattering(masses=(1., 1.7), momentum=.8, angle=.9):
    energy = np.sqrt(np.asarray(masses)**2+momentum**2)
    before = np.array([0., 0., momentum])
    after = momentum*np.array([np.sin(angle), 0., np.cos(angle)])
    return np.array([np.r_[energy[0], before], np.r_[energy[1], -before],
                     np.r_[energy[0], after], np.r_[energy[1], -after]])


def soft_tensor(p, q, couplings=(1., 1.)):
    coupling = np.tile(couplings, 2)
    denominator = p@ETA@q
    if np.min(abs(denominator)) < 1e-14:
        raise ValueError('Soft denominators must be nonzero.')
    return np.einsum('i,ia,ib->ab', SIGNS*coupling/denominator, p, p)


def ward_vector(p, couplings):
    return (SIGNS*np.tile(couplings, 2))@p


def channel_constraints(species_count, edges):
    rows = []
    for left, right in edges:
        p = scattering((1.+.2*left, 1.+.2*right))
        block = np.zeros((4, species_count))
        block[:, left] = p[2]-p[0]
        block[:, right] = p[3]-p[1]
        rows.extend(block)
    return np.asarray(rows)


def report():
    p = scattering(); q = np.array([1., .6, .8, 0.])*.01
    universal = q@ETA@soft_tensor(p, q)
    differing = q@ETA@soft_tensor(p, q, (1., 1.4))
    graph = channel_constraints(4, [(0, 1), (1, 2), (2, 3)])
    split = channel_constraints(4, [(0, 1), (2, 3)])
    return {'round': 326,
            'scope': 'Application of Weinberg leading soft graviton consistency in a supplied Lorentz-invariant asymptotic scattering setting with a massless helicity-2 mode, external-leg pole factorization and pure-gauge decoupling. These are additional inputs, not consequences of cognition or of the finite FLRW solution. Only leading pole consistency and species connectivity are tested; no complete S matrix or unique gravitational action is derived.',
            'universal_Ward_norm': float(np.linalg.norm(universal)),
            'nonuniversal_Ward_vector': differing.tolist(),
            'nonuniversal_Ward_norm': float(np.linalg.norm(differing)),
            'connected_channel_matrix_rank': int(np.linalg.matrix_rank(graph)),
            'connected_nullity': int(4-np.linalg.matrix_rank(graph)),
            'disconnected_audited_channel_nullity': int(4-np.linalg.matrix_rank(split)),
            'forward_scattering_Ward_norm': float(np.linalg.norm(ward_vector(scattering(angle=0.), (1., 1.4)))),
            'next_interface': 'Universal spin-2 coupling still allows independent scalar charges. Compute their weak-field observable and range dependence in the same covariant model class.'}


class Checks(unittest.TestCase):
    def test_01_physical_hard_kinematics(self):
        p = scattering()
        np.testing.assert_allclose(SIGNS@p, 0., atol=2e-15)
        norm = np.einsum('ia,ab,ib->i', p, ETA, p)
        np.testing.assert_allclose(norm, -np.tile(np.array([1., 1.7])**2, 2), atol=2e-15)

    def test_02_pole_contraction_is_Ward_vector(self):
        p = scattering()
        for direction in ([.6, .8, 0.], [0., 0., 1.], [1., 0., 0.]):
            q = .01*np.r_[1., direction]
            np.testing.assert_allclose(q@ETA@soft_tensor(p, q, (1., 1.4)), ward_vector(p, (1., 1.4)), atol=2e-14)

    def test_03_universal_coupling_decouples_pure_gauge(self):
        p = scattering(); q = .01*np.array([1., .6, .8, 0.])
        zeta = np.array([0., .8, -.6, 0.])
        self.assertLess(abs(q@ETA@zeta), 1e-16)
        polarization = np.outer(ETA@q, ETA@zeta)+np.outer(ETA@zeta, ETA@q)
        self.assertLess(abs(np.sum(soft_tensor(p, q)*polarization)), 2e-14)

    def test_04_nonuniversal_coupling_has_gauge_response(self):
        p = scattering(); q = .01*np.array([1., .6, .8, 0.])
        zeta = np.array([0., .8, -.6, 0.])
        self.assertLess(abs(q@ETA@zeta), 1e-16)
        gauge = np.outer(ETA@q, ETA@zeta)+np.outer(ETA@zeta, ETA@q)
        observed = np.sum(soft_tensor(p, q, (1., 1.4))*gauge)
        expected = 2*ward_vector(p, (1., 1.4))@ETA@zeta
        self.assertAlmostEqual(observed, expected, places=13)
        self.assertGreater(abs(observed), .1)

    def test_05_connected_species_leave_one_coupling(self):
        matrix = channel_constraints(4, [(0, 1), (1, 2), (2, 3)])
        self.assertEqual(np.linalg.matrix_rank(matrix), 3)
        np.testing.assert_allclose(matrix@np.ones(4), 0., atol=2e-15)
        null = np.linalg.svd(matrix)[2][-1]
        np.testing.assert_allclose(null/null[0], np.ones(4), atol=2e-14)

    def test_06_missing_channels_leave_independent_constants(self):
        matrix = channel_constraints(4, [(0, 1), (2, 3)])
        self.assertEqual(np.linalg.matrix_rank(matrix), 2)
        np.testing.assert_allclose(matrix@np.array([1., 1., 2., 2.]), 0., atol=2e-15)

    def test_07_forward_channel_cannot_identify_couplings(self):
        np.testing.assert_allclose(ward_vector(scattering(angle=0.), (1., 7.)), 0., atol=2e-15)
        self.assertGreater(np.linalg.norm(ward_vector(scattering(), (1., 7.))), 1.)

    def test_08_soft_scaling_does_not_repair_nonuniversality(self):
        p = scattering(); q = .01*np.array([1., .6, .8, 0.])
        s = soft_tensor(p, q, (1., 1.4))
        np.testing.assert_allclose(soft_tensor(p, q/10, (1., 1.4)), 10*s, atol=2e-12)
        np.testing.assert_allclose((q/10)@ETA@soft_tensor(p, q/10, (1., 1.4)), q@ETA@s, atol=2e-14)

    def test_09_Ward_vector_transforms_as_Lorentz_vector(self):
        speed = .35; gamma = 1/np.sqrt(1-speed**2)
        boost = np.eye(4); boost[0, 0] = boost[1, 1] = gamma
        boost[0, 1] = boost[1, 0] = -gamma*speed
        np.testing.assert_allclose(boost.T@ETA@boost, ETA, atol=3e-16)
        p = scattering()
        np.testing.assert_allclose(ward_vector(p@boost.T, (1., 1.4)), boost@ward_vector(p, (1., 1.4)), atol=3e-15)


if __name__ == '__main__':
    main(__name__, 'soft_spin2_universality_audit', report)
