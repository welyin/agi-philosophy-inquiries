"""Round 237: calibrate existing Myrheim-Meyer/chain statistics with an exact false positive.

Fixed N uniform points in null-coordinate square = Poisson sprinkling conditioned
on N in a given 1+1 Minkowski diamond. Geometry is a calibration input.
"""
import argparse
from functools import lru_cache
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np


def relation_from_points(points):
    points = np.asarray(points)
    return np.all(points[:, None, :] < points[None, :, :], axis=2)


def permutation_order(permutation):
    p = np.asarray(permutation)
    return np.triu(p[:, None] < p[None, :], 1)


def chain_counts(relation, maximum=4):
    a = relation.astype(np.int64)
    previous = np.ones(len(a), dtype=np.int64)
    counts = []
    for _ in range(maximum):
        counts.append(int(previous.sum()))
        previous = a.T@previous
    return counts


def ordering_fraction_model(d):
    return math.exp(math.lgamma(d+1)+math.lgamma(d/2)-math.log(2)-math.lgamma(3*d/2))


def dimension_from_fraction(r):
    if not 0 < r <= 1:
        raise ValueError('A finite dimension estimate requires 0 < r <= 1.')
    lower, upper = 1., 2.
    while ordering_fraction_model(upper) > r:
        upper *= 2
    for _ in range(64):
        mid = (lower+upper)/2
        if ordering_fraction_model(mid) > r:
            lower = mid
        else:
            upper = mid
    return (lower+upper)/2


def matched_two_layer(n):
    if n < 12 or n % 4:
        raise ValueError('Use N >= 12 divisible by four.')
    a = np.zeros((n, n), dtype=bool)
    a[:n//2, n//2:] = True
    for i in range(n//4):
        a[i, n//2+i] = False
    return a


def height(relation):
    # Any finite poset can be ordered by past size; predecessors have smaller size.
    order = np.argsort(relation.sum(axis=0), kind='stable')
    h = np.ones(len(relation), dtype=int)
    for j in order:
        predecessors = np.flatnonzero(relation[:, j])
        if len(predecessors):
            h[j] = 1+h[predecessors].max()
    return int(h.max()) if len(h) else 0


def metrics(relation):
    n = len(relation)
    c = chain_counts(relation)
    r = c[1]/math.comb(n, 2)
    return {'N': n, 'chain_counts_1_to_4': c, 'r': r,
            'three_chain_fraction': c[2]/math.comb(n, 3),
            'MM_dimension': dimension_from_fraction(r), 'height': height(relation)}


def describe(values):
    values = np.array(values, dtype=float)
    sd = float(values.std(ddof=1))
    return {'mean': float(values.mean()), 'sample_sd': sd,
            'mean_standard_error': sd/math.sqrt(len(values)),
            'min': float(values.min()), 'max': float(values.max())}


@lru_cache(maxsize=1)
def calibration():
    rng = np.random.default_rng(237001)
    results = []
    for n in [64, 128, 256]:
        samples = [metrics(relation_from_points(rng.random((n, 2)))) for _ in range(64)]
        r_theory_se = math.sqrt((2*n+5)/(18*n*(n-1))/len(samples))
        row = {'N': n, 'replicates': len(samples), 'seed_stream': 237001,
               'r_theoretical_mean': .5, 'r_theoretical_mean_SE': r_theory_se,
               'three_chain_theoretical_mean': 1/6}
        for key in ['r', 'three_chain_fraction', 'MM_dimension', 'height']:
            row[key] = describe([s[key] for s in samples])
        row['MM_dimension_of_mean_r'] = dimension_from_fraction(row['r']['mean'])
        row['r_mean_standardized_error'] = (row['r']['mean']-.5)/r_theory_se
        row['three_chain_mean_standardized_error'] = (
            row['three_chain_fraction']['mean']-1/6)/row['three_chain_fraction']['mean_standard_error']
        results.append(row)
    return results


@lru_cache(maxsize=1)
def exact_permutation_calibration():
    n = 6
    counts = np.array([chain_counts(permutation_order(p)) for p in itertools.permutations(range(n))])
    r = counts[:, 1]/math.comb(n, 2)
    return {'N': n, 'permutations': math.factorial(n),
            'mean_chain_counts': counts.mean(axis=0).tolist(),
            'exact_expected_chain_counts': [math.comb(n, k)/math.factorial(k) for k in range(1, 5)],
            'r_population_variance': float(r.var()),
            'r_analytic_variance': (2*n+5)/(18*n*(n-1))}


@lru_cache(maxsize=1)
def obstruction_certificate():
    s = matched_two_layer(12)[np.ix_([0, 1, 2, 6, 7, 8], [0, 1, 2, 6, 7, 8])]
    extensions = []
    for p in itertools.permutations(range(6)):
        positions = np.argsort(p)
        order = positions[:, None] < positions[None, :]
        if np.all(order[s]):
            extensions.append(order)
    intersections = sum(np.array_equal(a & b, s) for a in extensions for b in extensions)
    return {'induced_standard_example_S3': s.astype(int).tolist(),
            'linear_extensions': len(extensions), 'two_extension_pairs_checked': len(extensions)**2,
            'pairs_realizing_the_poset': intersections}


def scale_certificate():
    points = np.random.default_rng(237002).random((48, 2))
    return {'N': len(points), 'rescaling_factor': 3.,
            'same_order_after_rescaling': bool(np.array_equal(
                relation_from_points(points), relation_from_points(3*points))),
            'given_diamond_volumes': [1., 9.],
            'given_endpoint_proper_times': [math.sqrt(2), 3*math.sqrt(2)],
            'meaning': 'Coordinates and metric are supplied; order alone does not fix their scale.'}


def report():
    return {'round': 237, 'date': '2026-09-22',
            'hypothesis': 'Does matching MM dimension two certify a 1+1 Minkowski causal order?',
            'fixed_N_Minkowski_calibration': calibration(),
            'exact_small_ensemble': exact_permutation_calibration(),
            'exact_matched_negative_controls': [metrics(matched_two_layer(n)) for n in [64, 128, 256]],
            'two_dimensional_order_obstruction': obstruction_certificate(),
            'scale_control': scale_certificate(),
            'conclusion': 'The negative controls have r=1/2 exactly and MM dimension two but contain S3, obstructing a 1+1 order embedding. Three-chains reject them.',
            'extra_inputs': ['Finite event identities', 'Given 1+1 Minkowski diamond and uniform conditional sampling for calibration',
                             'No physical volume density derived'],
            'not_claimed': ['A new dimension estimator', 'Quantum gravity or spacetime emergence',
                            'Three-chain matching proves manifoldlikeness'],
            'runtime': {'python': platform.python_version(), 'numpy': np.__version__}}


class DimensionTests(unittest.TestCase):
    def test_known_dimension_values_and_inverse(self):
        for d, r in [(1, 1), (2, .5), (3, 8/35), (4, .1)]:
            self.assertAlmostEqual(ordering_fraction_model(d), r)
            self.assertAlmostEqual(dimension_from_fraction(r), d)
        with self.assertRaises(ValueError):
            dimension_from_fraction(0)

    def test_chain_counter_against_brute_force(self):
        rng = np.random.default_rng(237003)
        for _ in range(12):
            r = permutation_order(rng.permutation(7))
            brute = [sum(all(r[i, j] for i, j in itertools.combinations(c, 2))
                         for c in itertools.combinations(range(7), k)) for k in range(1, 5)]
            self.assertEqual(chain_counts(r), brute)

    def test_exact_small_ensemble_and_correlated_pair_variance(self):
        c = exact_permutation_calibration()
        np.testing.assert_allclose(c['mean_chain_counts'], c['exact_expected_chain_counts'], atol=1e-14)
        self.assertAlmostEqual(c['r_population_variance'], c['r_analytic_variance'])

    def test_negative_controls_match_dimension_but_have_no_three_chains(self):
        for n in [64, 128, 256]:
            r = matched_two_layer(n)
            c = metrics(r)
            self.assertEqual(c['r'], .5)
            self.assertAlmostEqual(c['MM_dimension'], 2.)
            self.assertEqual(c['height'], 2)
            self.assertEqual(c['three_chain_fraction'], 0.)
            self.assertFalse(np.any((r.astype(int)@r.astype(int)>0) & ~r))

    def test_induced_S3_has_no_two_linear_extension_realizer(self):
        c = obstruction_certificate()
        self.assertGreater(c['linear_extensions'], 0)
        self.assertEqual(c['pairs_realizing_the_poset'], 0)

    def test_seeded_calibration_matches_analytic_means(self):
        for c in calibration():
            self.assertLess(abs(c['r_mean_standardized_error']), 6)
            self.assertLess(abs(c['three_chain_mean_standardized_error']), 6)
            self.assertGreater(c['height']['min'], 2)

    def test_metric_scale_is_not_fixed_by_order(self):
        c = scale_certificate()
        self.assertTrue(c['same_order_after_rescaling'])
        self.assertNotEqual(*c['given_diamond_volumes'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(DimensionTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = report()
    data['checks'] = {'run': checks.testsRun, 'failures': 0, 'errors': 0}
    if args.write_results:
        target = Path(__file__).with_name('causal_dimension_audit_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8')) != json.loads(json.dumps(data)):
            raise RuntimeError('Existing result differs; inspect before replacing.')
        target.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(data, ensure_ascii=False, indent=2))
