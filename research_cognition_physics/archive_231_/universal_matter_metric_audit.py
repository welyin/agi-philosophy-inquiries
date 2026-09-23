"""Round 324: a shared point-particle metric and invariant nonuniversal couplings."""
import unittest
import numpy as np
from growing_stream_audit import main
from conformal_probe_motion_audit import ETA, NULL, connection


def velocity(v):
    v = np.asarray(v, dtype=float)
    return np.r_[1., v]/np.sqrt(1-v@v)


def projector(u, metric=ETA):
    return np.linalg.inv(metric)+np.outer(u, u)


def acceleration(alpha, gradient, u, metric=ETA):
    return -alpha*projector(u, metric)@gradient


def transformed_derivative(u, a, factor, frame_gradient):
    # Coordinate derivative of Ubar=U/B with dsbar=B ds.
    return (a-u*(frame_gradient@u))/factor**2


def example():
    phi, slope, frame_alpha = .4, .2, .5
    couplings = np.array([.3, .8])
    gradient = np.array([0., slope, 0., 0.])
    u = velocity([0., 0., 0.])
    factor = np.exp(frame_alpha*phi)
    base = np.array([acceleration(alpha, gradient, u)[1] for alpha in couplings])
    changed = np.array([acceleration(alpha-frame_alpha, gradient, u/factor, factor**2*ETA)[1] for alpha in couplings])
    return {'couplings': couplings.tolist(), 'base_spatial_accelerations': base.tolist(),
            'changed_frame_accelerations': changed.tolist(),
            'difference_rescaled_to_base': float(factor**2*(changed[0]-changed[1])),
            'base_difference': float(base[0]-base[1]),
            'mass_ratio_at_phi_0': 1., 'mass_ratio_at_phi_1': float(np.exp(couplings[0]-couplings[1])),
            'frame_factor': float(factor)}


def projection_ranks():
    rest = velocity([0., 0., 0.]); moving = velocity([.4, 0., 0.])
    p0, p1 = projector(rest), projector(moving)
    hidden = np.array([1., 0., 0., 0.])
    return {'single_velocity_rank': int(np.linalg.matrix_rank(p0)),
            'two_velocity_rank': int(np.linalg.matrix_rank(np.vstack([p0, p1]))),
            'hidden_gradient_rest_acceleration': float(np.linalg.norm(p0@hidden)),
            'hidden_gradient_moving_acceleration': float(np.linalg.norm(p1@hidden))}


def stealth_probe_example():
    phi = np.sqrt(6)/2  # Same u=1 point in the frozen round-321 background.
    gradient = -np.sqrt(6)/4*NULL
    u = velocity([0., 0., 0.])
    values = [float(acceleration(alpha, gradient, u)[1]) for alpha in (.3, .8)]
    return {'phi': float(phi), 'F': float(1-phi**2/6),
            'spatial_accelerations': values, 'difference': values[0]-values[1]}


def report():
    return {'round': 324,
            'scope': 'Positive smooth position-dependent masses of structureless test particles on a connected Lorentzian domain of dimension >=2. A common conformally related matter metric exists iff all mass ratios are spatially and temporally constant, equivalently universal free fall for every timelike initial velocity. This is a restricted classical criterion, not a derivation of the full Einstein equivalence principle or a backreacted multispecies solution.',
            'nonuniversal_example': example(), 'velocity_identifiability': projection_ranks(),
            'same_stealth_background_test_probes': stealth_probe_example(),
            'next_interface': 'Does closing the matter-scalar energy account enforce universal coupling? Test covariant exchange identities with independent species couplings before adding a new equivalence axiom.'}


class Checks(unittest.TestCase):
    def test_01_projector_kernel_is_velocity_covector(self):
        u = velocity([.3, -.2, .1]); p = projector(u)
        self.assertLess(np.linalg.norm(p@(ETA@u)), 4e-16)
        self.assertEqual(np.linalg.matrix_rank(p), 3)

    def test_02_two_initial_velocities_determine_gradient(self):
        r = projection_ranks()
        self.assertEqual((r['single_velocity_rank'], r['two_velocity_rank']), (3, 4))
        self.assertEqual(r['hidden_gradient_rest_acceleration'], 0.)
        self.assertGreater(r['hidden_gradient_moving_acceleration'], .4)

    def test_03_common_mass_factor_becomes_constant_masses(self):
        phi = np.linspace(-1., 1., 31)
        common = np.exp(.3*phi+.2*phi**2)
        constants = np.array([1., 2., 7.])
        masses = constants[:, None]*common
        np.testing.assert_allclose(masses/common, np.broadcast_to(constants[:, None], masses.shape), atol=2e-15)
        np.testing.assert_allclose(masses[1]/masses[0], 2., atol=2e-15)

    def test_04_common_metric_geodesic_equation(self):
        u = velocity([.2, .1, -.3])
        gradient = np.array([.1, -.2, .05, .3])
        factor, alpha = 1.3, .7
        a = acceleration(alpha, gradient, u)
        du = transformed_derivative(u, a, factor, alpha*gradient)
        connection_term = np.einsum('abc,b,c->a', connection(alpha*gradient), u/factor, u/factor)
        self.assertLess(np.linalg.norm(du+connection_term), 4e-16)

    def test_05_different_couplings_have_different_free_fall(self):
        e = example()
        np.testing.assert_allclose(e['base_spatial_accelerations'], [-.06, -.16], atol=1e-15)
        self.assertAlmostEqual(e['base_difference'], .1, places=14)

    def test_06_frame_change_does_not_remove_difference(self):
        e = example()
        self.assertAlmostEqual(e['difference_rescaled_to_base'], e['base_difference'], places=14)
        self.assertGreater(abs(e['changed_frame_accelerations'][0]-e['changed_frame_accelerations'][1]), .05)

    def test_07_mass_ratio_is_invariant(self):
        phi = np.linspace(-1., 2., 37)
        m1, m2, factor = np.exp(.3*phi), np.exp(.8*phi), np.exp(.5*phi+.1*phi**2)
        np.testing.assert_allclose((m1/factor)/(m2/factor), m1/m2, atol=3e-15)
        self.assertGreater(np.ptp(m1/m2), 1.)

    def test_08_force_transformation_for_general_initial_data(self):
        u = velocity([.4, .1, -.2])
        gradient = np.array([.3, -.1, .2, .4]); factor = 1.7
        for alpha in (.3, .8):
            a = acceleration(alpha, gradient, u)
            frame_gradient = .5*gradient
            du = transformed_derivative(u, a, factor, frame_gradient)
            abar = du+np.einsum('abc,b,c->a', connection(frame_gradient), u/factor, u/factor)
            expected = acceleration(alpha-.5, gradient, u/factor, factor**2*ETA)
            np.testing.assert_allclose(abar, expected, atol=3e-16)

    def test_09_constant_scalar_background_can_hide_couplings(self):
        u = velocity([.3, 0., 0.])
        np.testing.assert_array_equal(acceleration(.3, np.zeros(4), u), acceleration(.8, np.zeros(4), u))
        self.assertNotEqual(.3, .8)


    def test_10_existing_background_does_not_select_probe_universality(self):
        r = stealth_probe_example()
        self.assertAlmostEqual(r['F'], .75, places=14)
        self.assertAlmostEqual(r['difference'], np.sqrt(6)/8, places=14)


if __name__ == '__main__':
    main(__name__, 'universal_matter_metric_audit', report)
