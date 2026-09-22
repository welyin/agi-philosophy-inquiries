"""Round 239: finite 2-order sampling versus order-invariant, prefix-consistent growth.

All finite probability comparisons use Fraction. The fixed-N law is the
pushforward of a uniform permutation, not uniform over unlabeled posets.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
import causal_dimension_audit as dim
import interval_abundance_audit as intervals


def key(r):
    return tuple(tuple(bool(x) for x in row) for row in r)


def matrix(k):
    return np.array(k, dtype=bool).reshape(len(k), len(k))


@lru_cache(maxsize=None)
def extensions(k):
    n = len(k)
    past = [{i for i in range(n) if k[i][j]} for j in range(n)]
    found = []
    def visit(prefix, remaining):
        if not remaining:
            found.append(tuple(prefix))
        for j in sorted(remaining):
            if not past[j] & remaining:
                visit(prefix+[j], remaining-{j})
    visit([], set(range(n)))
    return tuple(found)


@lru_cache(maxsize=None)
def natural_relabelings(k):
    r = matrix(k)
    return tuple(sorted({key(r[np.ix_(order, order)]) for order in extensions(k)}))


def canonical(k):
    return min(natural_relabelings(k))


def insert_rank(p, rank):
    return tuple(x+(x >= rank) for x in p)+(rank,)


def insertion_permutations(n):
    out = [()]
    for size in range(n):
        out = [insert_rank(p, rank) for p in out for rank in range(size+1)]
    return out


@lru_cache(maxsize=None)
def target_law(n):
    out = defaultdict(F)
    for p in itertools.permutations(range(n)):
        out[canonical(key(dim.permutation_order(p)))] += F(1, math.factorial(n))
    return dict(out)


@lru_cache(maxsize=None)
def symmetrized_law(n):
    out = defaultdict(F)
    for k, probability in target_law(n).items():
        # Each distinct naturally labeled matrix has the same number of
        # linear extensions of a fixed representative (its automorphism count).
        labels = natural_relabelings(k)
        for label in labels:
            out[label] += probability/len(labels)
    return dict(out)


def unlabeled(law):
    out = defaultdict(F)
    for k, probability in law.items():
        out[canonical(k)] += probability
    return dict(out)


def prefix_law(law, size):
    out = defaultdict(F)
    for k, probability in law.items():
        cut = key(matrix(k)[:size, :size])
        out[canonical(cut)] += probability
    return dict(out)


def random_delete_law(n):
    out = defaultdict(F)
    for k, probability in target_law(n).items():
        r = matrix(k)
        for omitted in range(n):
            keep = [i for i in range(n) if i != omitted]
            out[canonical(key(r[np.ix_(keep, keep)]))] += probability/n
    return dict(out)


def total_variation(a, b):
    return sum((abs(a.get(k, F(0))-b.get(k, F(0))) for k in a.keys()|b.keys()), F(0))/2


def shape3(k):
    r = matrix(k)
    c = int(r.sum())
    if c == 0:
        return 'antichain'
    if c == 1:
        return 'chain_plus_isolate'
    if c == 3:
        return 'chain'
    return 'split' if r.sum(1).max() == 2 else 'join'


def covariance_witness():
    matrices = []
    distribution = defaultdict(F)
    for p in itertools.permutations(range(3)):
        distribution[key(dim.permutation_order(p))] += F(1, 6)
    for edge in [(0, 1), (0, 2), (1, 2)]:
        r = np.zeros((3, 3), dtype=bool)
        r[edge] = True
        k = key(r)
        matrices.append({'edge': list(edge), 'probability': str(distribution[k]),
                         'symmetrized_probability': str(symmetrized_law(3)[k])})
    return matrices


def quantum_sampler(weights):
    # Instrument I_r(rho)=w_r rho preserves any unknown data and reference.
    weights = np.array(weights, dtype=float)
    bell = np.array([1, 0, 0, 1], dtype=complex)/math.sqrt(2)
    rho = np.outer(bell, bell.conj())
    kraus = [math.sqrt(w)*np.eye(2) for w in weights]
    blocks = []
    for k in kraus:
        joint = np.kron(k, np.eye(2))
        blocks.append(joint@rho@joint.conj().T)
    return {'weights': weights.tolist(),
            'completeness_error': float(np.max(np.abs(sum(k.conj().T@k for k in kraus)-np.eye(2)))),
            'unknown_state_and_reference_error': float(np.max(np.abs(sum(blocks)-rho))),
            'outcome_probability_error': float(np.max(np.abs(np.array([np.trace(b).real for b in blocks])-weights)))}


@lru_cache(maxsize=1)
def report():
    chain2 = key(dim.permutation_order([0, 1]))
    q3 = target_law(3)
    prefix = prefix_law(symmetrized_law(3), 2)
    return {'round': 239, 'date': '2026-09-22',
            'question': 'Can fixed-N uniform 2-order laws be exact marginals of one internally temporal, order-invariant growth?',
            'fixed_N_law': 'Uniform permutation pushforward; not uniform unlabeled posets',
            'N3_target_shapes': {shape3(k):str(v) for k,v in q3.items()},
            'insertion_covariance_witness': covariance_witness(),
            'prefix_obstruction': {
                'target_N2_chain_probability': str(target_law(2)[chain2]),
                'N3_covariant_prefix_chain_probability': str(prefix[canonical(chain2)]),
                'exact_gap': str(target_law(2)[chain2]-prefix[canonical(chain2)]),
                'sum_of_TV_errors_lower_bound': '1/18',
                'max_of_TV_errors_lower_bound': '1/36'},
            'prefix_and_random_deletion_checks': [
                {'N': n, 'prefix_TV':str(total_variation(prefix_law(symmetrized_law(n), n-1), target_law(n-1))),
                 'random_deletion_TV':str(total_variation(random_delete_law(n), target_law(n-1)))}
                for n in range(3, 6)],
            'quantum_sampler_certificates': [quantum_sampler(w) for w in [[1/3]*3, [.2, .3, .5]]],
            'extra_inputs': ['Two independent total orders and their uniform measure',
                             'Order-invariance and internal temporality, if requiring physical sequential growth'],
            'conclusion': 'Exact target marginals, order-invariance and prefix consistency already conflict between N=2 and N=3.',
            'not_claimed': ['All two-dimensional spacetime growth is impossible',
                            'Quantum amplitudes or histories must obey this classical probability model',
                            'The old cognition contract forces order-invariance'],
            'runtime': {'python': platform.python_version(), 'numpy': np.__version__}}


class GrowthTests(unittest.TestCase):
    def test_insertion_is_bijection_and_preserves_old_order(self):
        for n in range(1, 7):
            ps = insertion_permutations(n)
            self.assertEqual(len(ps), math.factorial(n))
            self.assertEqual(set(ps), set(itertools.permutations(range(n))))
        for p in itertools.permutations(range(4)):
            old = dim.permutation_order(p)
            for rank in range(5):
                new = dim.permutation_order(insert_rank(p, rank))
                np.testing.assert_array_equal(new[:4, :4], old)
                self.assertFalse(new[4].any())

    def test_insertion_matches_previous_fixed_N_interval_calibration(self):
        ps = insertion_permutations(5)
        avg = np.mean([intervals.profile(dim.permutation_order(p)) for p in ps], axis=0)
        np.testing.assert_allclose(avg, intervals.expected_profile(5), atol=1e-13)

    def test_exact_covariance_witness(self):
        self.assertEqual([x['probability'] for x in covariance_witness()], ['1/6', '0', '1/6'])
        self.assertEqual([x['symmetrized_probability'] for x in covariance_witness()], ['1/9']*3)

    def test_symmetrization_preserves_shapes_and_equal_label_weights(self):
        for n in range(2, 6):
            law = symmetrized_law(n)
            self.assertEqual(sum(law.values()), 1)
            self.assertEqual(unlabeled(law), target_law(n))
            for k in target_law(n):
                self.assertEqual(len({law[x] for x in natural_relabelings(k)}), 1)

    def test_prefix_obstruction_is_exact(self):
        q = report()['prefix_obstruction']
        self.assertEqual(q['target_N2_chain_probability'], '1/2')
        self.assertEqual(q['N3_covariant_prefix_chain_probability'], '4/9')
        self.assertEqual(q['exact_gap'], '1/18')

    def test_random_deletion_consistency_is_not_birth_prefix_consistency(self):
        for n in range(3, 6):
            self.assertEqual(random_delete_law(n), target_law(n-1))
        self.assertNotEqual(prefix_law(symmetrized_law(3), 2), target_law(2))

    def test_instrument_does_not_select_record_weights_or_damage_unknown_state(self):
        for row in report()['quantum_sampler_certificates']:
            for name in ['completeness_error','unknown_state_and_reference_error','outcome_probability_error']:
                self.assertLess(row[name], 1e-14)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(GrowthTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = dict(report(), checks={'run':checks.testsRun,'failures':0,'errors':0})
    if args.write_results:
        target = Path(__file__).with_name('random_order_growth_audit_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8')) != data:
            raise RuntimeError('Existing result differs; inspect before replacing.')
        target.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n',encoding='utf-8')
    print(json.dumps(data, ensure_ascii=False, indent=2))
