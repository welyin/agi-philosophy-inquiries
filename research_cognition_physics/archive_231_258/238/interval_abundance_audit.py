"""Round 238: existing interval-abundance diagnostics, finite-N calibration and thinning.

Reuse round237's seeded fixed-cardinality ensemble. The conditional finite-N
formula is derived explicitly in the note rather than substituted for a Poisson
ensemble formula. Finite profile collisions do not refute asymptotic rigidity.
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
import causal_dimension_audit as dim


def profile(relation):
    n = len(relation)
    if n < 2:
        return np.zeros(0, dtype=np.int64)
    counts = relation.astype(np.int64)@relation.astype(np.int64)
    return np.bincount(counts[relation], minlength=n-1)


def expected_profile(n):
    if n < 2:
        return np.zeros(0)
    harmonic = np.r_[0., np.cumsum(1/np.arange(1, n+1, dtype=float))]
    m = np.arange(n-1)
    return n*(harmonic[n-1]-harmonic[m]) + (m+1)*(harmonic[n]-harmonic[m+1]) - 2*n+2*(m+1)


def thinning_transform(abundance, q):
    transformed = np.zeros(len(abundance), dtype=float)
    for n, count in enumerate(abundance):
        for m in range(n+1):
            transformed[m] += count*q*q*math.comb(n, m)*q**m*(1-q)**(n-m)
    return transformed


@lru_cache(maxsize=1)
def calibration():
    # Exactly the same original point draws as round237, not a new ensemble.
    rng = np.random.default_rng(237001)
    output = []
    for n in [64, 128, 256]:
        samples = np.array([profile(dim.relation_from_points(rng.random((n, 2)))) for _ in range(64)])
        means = samples.mean(axis=0)
        sem = samples.std(axis=0, ddof=1)/math.sqrt(len(samples))
        theory = expected_profile(n)
        low_z = (means[:8]-theory[:8])/sem[:8]
        neg = profile(dim.matched_two_layer(n))
        output.append({'N': n, 'replicates': 64, 'same_point_ensemble_as_round237': True,
                       'mean_profile': means.tolist(), 'mean_standard_errors': sem.tolist(),
                       'conditional_theory': theory.tolist(),
                       'low_m_0_to_7_standardized_mean_errors': low_z.tolist(),
                       'mean_C2_from_profile': float(means.sum()),
                       'mean_C3_from_profile': float(np.arange(n-1)@means),
                       'mean_links': float(means[0]),
                       'theoretical_links': float(theory[0]),
                       'negative_control_profile': neg.tolist(),
                       'negative_to_expected_links_ratio': float(neg[0]/theory[0])})
    return output


@lru_cache(maxsize=1)
def exact_ensemble_certificate():
    n = 6
    samples = [profile(dim.permutation_order(p)) for p in itertools.permutations(range(n))]
    return {'N': n, 'permutations': math.factorial(n),
            'enumerated_mean_profile': np.mean(samples, axis=0).tolist(),
            'conditional_formula': expected_profile(n).tolist(),
            'maximum_error': float(np.max(np.abs(np.mean(samples, axis=0)-expected_profile(n))))}


def enumerate_thinning(relation, q):
    n = len(relation)
    mean = np.zeros(n-1)
    total = 0.
    for mask in range(1 << n):
        keep = [i for i in range(n) if (mask >> i) & 1]
        probability = q**len(keep)*(1-q)**(n-len(keep))
        total += probability
        sub = relation[np.ix_(keep, keep)]
        abundance = profile(sub)
        mean[:len(abundance)] += probability*abundance
    return mean, total


@lru_cache(maxsize=1)
def thinning_certificate():
    r = dim.permutation_order([0, 2, 1, 3, 5, 4, 7, 6])
    p = profile(r)
    cases = []
    for q in [.25, .5, .75, 1.]:
        exact, total = enumerate_thinning(r, q)
        predicted = thinning_transform(p, q)
        cases.append({'keep_probability': q, 'all_subsets': 256,
                      'total_probability_error': float(abs(total-1)),
                      'profile_transform_max_error': float(np.max(np.abs(exact-predicted))),
                      'C2_scaling_error': float(abs(exact.sum()-q*q*p.sum())),
                      'C3_scaling_error': float(abs(np.arange(7)@exact-q**3*(np.arange(7)@p)))})
    return cases


def conditional_thinning_mixture(n, q):
    mixed = np.zeros(n-1)
    for k in range(2, n+1):
        p = math.comb(n, k)*q**k*(1-q)**(n-k)
        mixed[:k-1] += p*expected_profile(k)
    return mixed


def collision_certificate():
    permutations = [[0, 1, 3, 5, 4, 2], [0, 1, 4, 5, 2, 3]]
    graphs = [dim.permutation_order(p) for p in permutations]
    degrees = [sorted(zip(r.sum(0).tolist(), r.sum(1).tolist())) for r in graphs]
    dual = sorted((b, a) for a, b in degrees[0])
    return {'permutations': permutations,
            'relations': [r.astype(int).tolist() for r in graphs],
            'profiles': [profile(r).tolist() for r in graphs],
            'sorted_past_future_counts': degrees,
            'profiles_equal': bool(np.array_equal(profile(graphs[0]), profile(graphs[1]))),
            'different_degree_multisets': degrees[0] != degrees[1],
            'not_just_order_duals_by_degree_test': dual != degrees[1],
            'scope': 'Finite invariant is not injective. Both examples are two-dimensional orders; no asymptotic or faithful-embedding conjecture is refuted.'}


def report():
    return {'round': 238, 'date': '2026-09-22',
            'hypothesis': 'Does the existing interval profile resolve round237 false positives and determine a unique finite geometry?',
            'interval_profile_calibration': calibration(),
            'exact_small_ensemble': exact_ensemble_certificate(),
            'exact_thinning_checks': thinning_certificate(),
            'finite_profile_collision': collision_certificate(),
            'conclusion': 'Interval abundance rejects the matched two-layer controls and has an exact thinning law, but equal finite profiles do not identify the whole order.',
            'extra_inputs': ['Same supplied conditional Minkowski ensemble as round237',
                             'Independent per-event thinning if used'],
            'not_claimed': ['New general geometry reconstruction theorem',
                            'Proof or disproof of asymptotic CST rigidity',
                            'Thinning identities alone prove manifoldlikeness'],
            'runtime': {'python': platform.python_version(), 'numpy': np.__version__}}


class IntervalTests(unittest.TestCase):
    def test_profile_counts_and_chain_moments_by_independent_enumeration(self):
        r = dim.permutation_order([2, 0, 4, 1, 6, 3, 7, 5])
        brute = np.zeros(7, dtype=int)
        for i in range(8):
            for j in range(8):
                if r[i, j]:
                    interior = sum(bool(r[i, k] and r[k, j]) for k in range(8))
                    brute[interior] += 1
        np.testing.assert_array_equal(profile(r), brute)
        c = dim.chain_counts(r)
        self.assertEqual(int(brute.sum()), c[1])
        self.assertEqual(int(np.arange(7)@brute), c[2])

    def test_conditional_formula_exact_permutations_and_moments(self):
        self.assertLess(exact_ensemble_certificate()['maximum_error'], 1e-13)
        for n in [2, 3, 6, 64, 128, 256]:
            theory = expected_profile(n)
            self.assertGreater(theory.min(), -1e-10)
            self.assertAlmostEqual(theory.sum(), math.comb(n, 2)/2, delta=1e-8)
            self.assertAlmostEqual(np.arange(n-1)@theory, math.comb(n, 3)/6, delta=1e-7)

    def test_seeded_low_interval_profile_and_round237_consistency(self):
        previous = dim.calibration()
        for row, old in zip(calibration(), previous):
            self.assertLess(max(abs(z) for z in row['low_m_0_to_7_standardized_mean_errors']), 6)
            self.assertAlmostEqual(row['mean_C2_from_profile']/math.comb(row['N'], 2), old['r']['mean'])
            self.assertAlmostEqual(row['mean_C3_from_profile']/math.comb(row['N'], 3), old['three_chain_fraction']['mean'])

    def test_negative_controls_are_all_links(self):
        for n in [64, 128, 256]:
            p = profile(dim.matched_two_layer(n))
            self.assertEqual(p[0], n*(n-1)//4)
            self.assertEqual(int(p[1:].sum()), 0)
            self.assertGreater(p[0]/expected_profile(n)[0], 4)

    def test_thinning_transform_against_all_subsets(self):
        for row in thinning_certificate():
            for key in ['total_probability_error', 'profile_transform_max_error', 'C2_scaling_error', 'C3_scaling_error']:
                self.assertLess(row[key], 1e-12)

    def test_thinning_uses_binomial_mixture_not_mean_N_substitution(self):
        for n in [8, 16, 32]:
            for q in [.25, .5, .75]:
                np.testing.assert_allclose(thinning_transform(expected_profile(n), q),
                                           conditional_thinning_mixture(n, q), atol=2e-11, rtol=2e-11)

    def test_equal_profiles_do_not_identify_order_even_up_to_duality(self):
        c = collision_certificate()
        self.assertTrue(c['profiles_equal'])
        self.assertTrue(c['different_degree_multisets'])
        self.assertTrue(c['not_just_order_duals_by_degree_test'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(IntervalTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = report()
    data['checks'] = {'run': checks.testsRun, 'failures': 0, 'errors': 0}
    if args.write_results:
        target = Path(__file__).with_name('interval_abundance_audit_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8')) != json.loads(json.dumps(data)):
            raise RuntimeError('Existing result differs; inspect before replacing.')
        target.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(data, ensure_ascii=False, indent=2))
