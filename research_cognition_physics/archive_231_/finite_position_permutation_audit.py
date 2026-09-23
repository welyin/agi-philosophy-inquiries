"""Round 388: a robust four-position obstruction to a circular boundary.

The complete position shell, one actual homeomorphism, and connected-region or
continuous-path calibration are inputs. Finite endpoint matching is not a
certificate for arbitrary metrics, program labels, or uncontrolled path gaps.
"""
import itertools
import unittest

import numpy as np

from direction_contrast_tomography_audit import (
    TETRA, probabilities, rotation, tetrahedral_group,
)
from growing_stream_audit import main


SIGMA = np.array([1, 2, 0, 3], dtype=int)
AXIS = np.array([0.23, -0.37, 0.59])
AXIS /= np.linalg.norm(AXIS)
ANGLE = 0.14


def tetrahedral_seed():
    """Reuse round 380's matrices, without recounting group closure as new work."""
    return next(r.astype(float) for r in tetrahedral_group()
                if np.max(abs(TETRA @ r.T - TETRA[SIGMA])) < 1e-12)


def dihedral_permutations(order):
    order = tuple(order)
    n = len(order)
    out = set()
    for sign in (1, -1):
        for shift in range(n):
            p = [0] * n
            for k, label in enumerate(order):
                p[label] = order[(shift + sign * k) % n]
            out.add(tuple(p))
    return out


def geodesic(u, v):
    u, v = np.asarray(u), np.asarray(v)
    return np.arccos(np.clip(np.sum(u * v, axis=-1), -1., 1.))


def circle_distance(u, v):
    return np.abs((np.asarray(u) - np.asarray(v) + 0.5) % 1. - 0.5)


def metric_margin(center_distances, matching_distances, eta_s, eta_e):
    """A numerical margin only; connected metric balls are a separate input."""
    if eta_s < 0 or eta_e < 0:
        raise ValueError('Error bounds must be nonnegative.')
    return float(np.min(center_distances) - eta_s
                 - 2 * (np.max(matching_distances) + eta_e))


def path_margin(readings, reading_error, modulus_at_mesh):
    """Bound pairwise path separation in one common position-only readout."""
    readings = np.asarray(readings)
    gaps = [np.min(np.linalg.norm(readings[i, :, None, :]
                                 - readings[j, None, :, :], axis=-1))
            for i in range(len(readings)) for j in range(i)]
    return float(min(gaps) - 2 * (reading_error + modulus_at_mesh))


def sphere_case():
    seed = tetrahedral_seed()
    actual = rotation(AXIS, ANGLE) @ seed
    images = TETRA @ actual.T
    center_distances = np.array([geodesic(TETRA[i], TETRA[j])
                                 for i in range(4) for j in range(i)])
    matching = geodesic(images, TETRA[SIGMA])
    eta_s, eta_e = 0.01, 0.005
    center_estimates = center_distances + eta_s * np.linspace(-0.8, 0.8, 6)
    matching_estimates = matching + eta_e * np.array([0.8, -0.6, 0.4, -0.2])
    margin = metric_margin(center_estimates, matching_estimates, eta_s, eta_e)
    s_lower = float(np.min(center_estimates) - eta_s)
    e_upper = float(np.max(matching_estimates) + eta_e)
    return {'seed': seed, 'actual': actual, 'images': images,
            'center_distances': center_distances, 'matching': matching,
            'eta_s': eta_s, 'eta_e': eta_e, 'center_estimates': center_estimates,
            'matching_estimates': matching_estimates, 's_lower': s_lower,
            'e_upper': e_upper, 'margin': margin,
            'ball_radius': (e_upper + s_lower / 2) / 2,
            'cube_defect': float(np.linalg.norm(np.linalg.matrix_power(actual, 3)
                                                - np.eye(3), ord=2))}


def sphere_paths(times):
    """Path j goes from x_j to a(x_sigma_inverse(j))."""
    return np.stack([TETRA @ rotation(AXIS, ANGLE * t).T for t in times], axis=1)


def path_case():
    intervals, visibility, error = 8, 0.7, 0.002
    grid = np.linspace(0., 1., intervals + 1)
    paths = sphere_paths(grid)
    means = (1 + visibility * paths) / 2
    rng = np.random.default_rng(388)
    noise = rng.normal(size=means.shape)
    noise *= error / np.linalg.norm(noise, axis=-1, keepdims=True)
    readings = means + noise
    mesh = 1 / (2 * intervals)
    modulus = visibility * ANGLE * mesh / 2
    return {'intervals': intervals, 'grid': grid, 'paths': paths,
            'visibility': visibility, 'means': means, 'readings': readings,
            'error': error, 'mesh': mesh, 'modulus': modulus,
            'margin': path_margin(readings, error, modulus)}


def sharp_lift(t, kappa):
    """A decreasing degree-minus-one lift, for 0 < kappa < 1/8."""
    if not 0 < kappa < 0.125:
        raise ValueError('Require 0 < kappa < 1/8.')
    t = np.asarray(t, dtype=float)
    integer = np.floor(t)
    fraction = t - integer
    ordinates = [0.375 + kappa, 0.375 - kappa, 0., -0.25,
                 -0.625 + kappa]
    return np.interp(fraction, np.linspace(0., 1., 5), ordinates) - integer


def sharp_case(kappa):
    centers = np.arange(4) / 4
    images = sharp_lift(centers, kappa) % 1.
    errors = circle_distance(images, centers[SIGMA])
    slopes = np.array([-8 * kappa, -1.5 + 4 * kappa,
                       -1., -1.5 + 4 * kappa])
    return {'kappa': kappa, 's_min': 0.25, 'e_max': float(max(errors)),
            'margin': float(0.25 - 2 * max(errors)),
            'slopes': slopes, 'Lipschitz': float(max(abs(slopes))),
            'inverse_Lipschitz': float(1 / min(abs(slopes)))}


def circle_paired_readout(t):
    """Continuous R^4 readout, equal on four disconnected two-point sets."""
    centers = np.arange(4) / 4
    shifted = (centers + 0.07) % 1
    knots = np.concatenate((centers, shifted))
    labels = np.concatenate((np.arange(4), SIGMA))
    order = np.argsort(knots)
    knots, labels = knots[order], labels[order]
    values = np.eye(4)[labels]
    knots = np.r_[knots, 1.]
    values = np.vstack((values, values[0]))
    t = np.asarray(t) % 1.
    return np.stack([np.interp(t, knots, values[:, j]) for j in range(4)], axis=-1)


def compatible_nonlength_metric(u, v, weight=0.1):
    return (weight * circle_distance(u, v)
            + np.linalg.norm(circle_paired_readout(u) - circle_paired_readout(v), axis=-1))


def metric_counterexample():
    centers = np.arange(4) / 4
    images = (centers + 0.07) % 1
    center_distances = np.array([compatible_nonlength_metric(centers[i], centers[j])
                                 for i in range(4) for j in range(i)])
    matching = compatible_nonlength_metric(images, centers[SIGMA])
    return {'centers': centers, 'images': images,
            's_min': float(min(center_distances)), 'e_max': float(max(matching)),
            'naive_margin': metric_margin(center_distances, matching, 0., 0.),
            'radius': 0.1,
            'included_pair_distances_from_zero': compatible_nonlength_metric(0., [0., 0.57]),
            'excluded_separators_distances_from_zero': compatible_nonlength_metric(0., [0.07, 0.75])}


def report():
    s, p, c = sphere_case(), path_case(), metric_counterexample()
    return {
        'round': 388,
        'scope': {
            'theorem': 'A homeomorphism of one full shell S^(n-1), four distinct positions, and four pairwise disjoint connected sets pairing sigma=(123)(4) imply n>=3.',
            'required_inputs': ['A complete actual position shell and its sphere identity',
                                'One actual global position homeomorphism',
                                'Connected regions or actual continuous connecting paths',
                                'Calibrated error bounds and, for path sampling, a known common position-readout continuity modulus'],
            'metric_certificate': 's_hat-eta_s-2*(e_hat+eta_e)>0 is sufficient only with separately certified connected balls, for example a geodesic metric.',
            'path_certificate': 'The same position-only readout is used on all paths, with matched internal/controller preparation. Finite samples need a uniform parameter coverage modulus.',
            'reused': 'Round-380 tetrahedral_group and TETRA supply the seed rotation and points; their old group closure is not counted as new science.',
            'not_claimed': ['The whole map has order three', 'The whole A4 presentation holds',
                            'The map sends entire connected regions into each other',
                            'Arbitrary compatible metric balls are connected',
                            'Finite calibration alone proves actual position identity or a global homeomorphism',
                            'FUCP alone supplies the geometry or redirection',
                            'This lower bound alone selects exactly three dimensions'],
        },
        'circle_order_obstruction': {
            'sigma_zero_based': SIGMA.tolist(), 'all_label_orders_checked': 24,
            'compatible_dihedral_orders_for_sigma': sum(tuple(SIGMA) in dihedral_permutations(o)
                                                       for o in itertools.permutations(range(4))),
            'dimension_one': 'Four disjoint nonempty sets cannot fit in S^0.'},
        'perturbed_actual_rotation': {
            'angle': ANGLE, 'axis': AXIS.tolist(), 'seed': s['seed'].tolist(),
            'actual': s['actual'].tolist(), 'a_cubed_operator_defect': s['cube_defect'],
            's_min': float(min(s['center_distances'])), 'e_max': float(max(s['matching'])),
            'eta_s': s['eta_s'], 'eta_e': s['eta_e'], 's_lower': s['s_lower'],
            'e_upper': s['e_upper'], 'certified_margin': s['margin'],
            'connected_ball_radius': s['ball_radius']},
        'continuous_path_readout_certificate': {
            'paths': 4, 'samples_per_path': len(p['grid']), 'binary_settings_per_sample': 3,
            'mean_probabilities_to_estimate': 4 * len(p['grid']) * 3,
            'position_effect_visibility': p['visibility'],
            'vector_readout_error_bound': p['error'], 'parameter_cover_radius': p['mesh'],
            'readout_modulus_at_cover_radius': p['modulus'],
            'certified_all_times_cross_path_readout_separation': p['margin'],
            'samples_are_not_measurement_shot_count': True},
        'sharpness_from_circle_homeomorphisms': [
            {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in sharp_case(z).items()}
            for z in (0.025, 0.005, 0.001)],
        'insufficient_contracts': {
            'only_three_positions': 'A circle rotation by 1/3 realizes an exact three-cycle.',
            'continuous_but_not_homeomorphism': 't -> 2t mod 1 cycles 1/7,2/7,4/7 and fixes 0.',
            'compatible_metric_without_connected_balls': {
                k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in c.items()},
            'uncontrolled_path_sampling': 'At times 0 and 1, path 0 is read as zero; its triangular excursion reaches the constant path 1 at time 1/2. The certified modulus rejects the false endpoint-only certificate.'},
        'source': 'https://arxiv.org/pdf/math/0607481v3',
    }


class Checks(unittest.TestCase):
    def test_01_three_cycle_with_fixed_point_fails_every_circular_order(self):
        for order in itertools.permutations(range(4)):
            permitted = dihedral_permutations(order)
            self.assertEqual(len(permitted), 8)
            self.assertNotIn(tuple(SIGMA), permitted)

    def test_02_actual_perturbation_breaks_order_three_but_preserves_margin(self):
        s = sphere_case()
        np.testing.assert_allclose(s['actual'].T @ s['actual'], np.eye(3), atol=1e-14)
        self.assertGreater(s['cube_defect'], 0.05)
        self.assertLessEqual(max(s['matching']), ANGLE + 1e-14)
        self.assertGreater(s['margin'], 1.5)
        self.assertLess(s['e_upper'], s['ball_radius'])
        self.assertLess(s['ball_radius'], s['s_lower'] / 2)

    def test_03_distance_error_certificate_has_the_declared_worst_case_bound(self):
        s = sphere_case()
        exact = min(s['center_distances']) - 2 * max(s['matching'])
        for es, ee in itertools.product((-1., 1.), repeat=2):
            bound = metric_margin(s['center_distances'] + es * s['eta_s'],
                                  s['matching'] + ee * s['eta_e'], s['eta_s'], s['eta_e'])
            self.assertLessEqual(bound, exact + 1e-14)
        self.assertLess(metric_margin([0.3], [0.12], 0.03, 0.03), 0.)

    def test_04_actual_paths_connect_the_required_endpoints_and_read_one_function(self):
        p, s = path_case(), sphere_case()
        np.testing.assert_allclose(p['paths'][:, 0], TETRA, atol=1e-14)
        np.testing.assert_allclose(p['paths'][SIGMA, -1], s['images'], atol=1e-14)
        by_born = probabilities(p['paths'].reshape(-1, 3), p['visibility'] * np.eye(3))
        np.testing.assert_allclose(by_born.reshape(p['means'].shape), p['means'], atol=1e-14)

    def test_05_known_path_modulus_yields_positive_all_times_separation(self):
        p = path_case()
        self.assertGreater(p['margin'], 0.48)
        self.assertLessEqual(np.max(np.linalg.norm(p['readings'] - p['means'], axis=-1)),
                             p['error'] + 1e-14)
        fine_t = np.linspace(0., 1., 129)
        fine = (1 + p['visibility'] * sphere_paths(fine_t)) / 2
        for i in range(4):
            increments = np.linalg.norm(fine[i, 1:] - fine[i, :-1], axis=-1)
            self.assertLessEqual(max(increments), p['visibility'] * ANGLE / 256 + 1e-14)
        fine_min = min(np.min(np.linalg.norm(fine[i, :, None] - fine[j, None, :], axis=-1))
                       for i in range(4) for j in range(i))
        self.assertGreaterEqual(fine_min, p['margin'] - 1e-14)

    def test_06_nondistinguishing_readout_does_not_certify_paths(self):
        self.assertLess(path_margin(np.ones((4, 3, 2)), 0.01, 0.01), 0.)

    def test_07_factor_two_is_sharp_without_a_uniform_inverse_lipschitz_bound(self):
        grid = np.linspace(-2., 3., 2001)
        for kappa in (0.025, 0.005, 0.001):
            c = sharp_case(kappa)
            values = sharp_lift(grid, kappa)
            self.assertTrue(np.all(np.diff(values) < 0))
            np.testing.assert_allclose(sharp_lift(grid + 1, kappa), values - 1, atol=2e-15)
            self.assertAlmostEqual(c['e_max'], 0.125 + kappa)
            self.assertAlmostEqual(c['margin'], -2 * kappa)
            self.assertGreater(c['inverse_Lipschitz'], 1.)

    def test_08_three_positions_alone_are_realized_on_the_circle(self):
        centers = np.arange(3) / 3
        np.testing.assert_allclose(circle_distance((centers + 1 / 3) % 1, centers[[1, 2, 0]]), 0, atol=1e-15)
        self.assertIn((1, 2, 0), dihedral_permutations(range(3)))

    def test_09_noninvertible_circle_map_has_the_exact_forbidden_pattern(self):
        centers = np.array([1 / 7, 2 / 7, 4 / 7, 0.])
        np.testing.assert_allclose(circle_distance(2 * centers, centers[SIGMA]), 0, atol=1e-15)
        self.assertEqual(float(circle_distance(2 * 0., 2 * 0.5)), 0.)

    def test_10_compatible_nonlength_metric_gives_a_false_ball_certificate(self):
        c = metric_counterexample()
        self.assertGreater(c['naive_margin'], 1.)
        self.assertTrue(np.all(c['included_pair_distances_from_zero'] < c['radius']))
        self.assertTrue(np.all(c['excluded_separators_distances_from_zero'] > c['radius']))
        eight = np.r_[c['centers'], c['images']]
        self.assertEqual(len(np.unique(eight)), 8)
        rng = np.random.default_rng(1388)
        triples = rng.random((200, 3))
        self.assertTrue(np.all(compatible_nonlength_metric(triples[:, 0], triples[:, 2])
                               <= compatible_nonlength_metric(triples[:, 0], triples[:, 1])
                               + compatible_nonlength_metric(triples[:, 1], triples[:, 2]) + 1e-14))

    def test_11_uncontrolled_path_excursion_can_fool_endpoint_samples(self):
        readings = np.repeat(np.arange(4.)[:, None, None], 2, axis=1)
        self.assertEqual(path_margin(readings, 0., 0.), 1.)
        self.assertEqual(path_margin(readings, 0., 1.), -1.)
        triangle = lambda t: 1 - abs(2 * t - 1)
        self.assertEqual(triangle(0.), 0.)
        self.assertEqual(triangle(1.), 0.)
        self.assertEqual(triangle(0.5), 1.)


if __name__ == '__main__':
    main(__name__, 'finite_position_permutation_audit', report)
