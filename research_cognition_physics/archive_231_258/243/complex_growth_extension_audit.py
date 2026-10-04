"""Round 243: a cited complex sequential growth family and its extension.

Reuse the finite labelled-poset implementation from round 240. Single-coupling
models t0=1, t_k=tau have a parent-independent absolute transition sum.
Infinite extension conclusions use analysis and the cited scalar-measure
theorem, not extrapolation of these finite enumerations.
"""
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
import csg_support_audit as csg
import random_order_growth_audit as growth


def coefficient(size, maxima, k):
    return math.comb(size-maxima,k-maxima) if maxima <= k <= size else 0


def transition(r, subset, k, tau=1j):
    size = len(r)
    numerator = (1 if not subset else 0)+tau*coefficient(
        len(subset),csg.maximal_count(r,subset),k)
    denominator = 1+tau*math.comb(size,k) if size >= k else 1
    if denominator == 0:
        raise ValueError('Singular coupling: lambda(n,0)=0')
    return numerator/denominator


@lru_cache(maxsize=None)
def level(n,k,tau=1j):
    if n == 0:
        return {():1+0j}
    out = {}
    for parent,a in level(n-1,k,tau).items():
        r = growth.matrix(parent)
        for subset in csg.ideals(r):
            child = growth.key(csg.append_event(r,subset))
            out[child] = a*transition(r,subset,k,tau)
    return out


def absolute_transition_sum(n,k,tau=1j):
    r = math.comb(n,k) if n >= k else 0
    return (1+abs(tau)*r)/abs(1+tau*r)


def log_variation(n,k):
    # tau=i. Stable at large n: log((1+r)/sqrt(1+r^2)).
    return math.fsum(math.log1p(1/math.comb(j,k))-.5*math.log1p(1/math.comb(j,k)**2)
                     for j in range(k,n))


def variation(n,k):
    return math.exp(log_variation(n,k))


def extension_bounds(n=10000):
    if n < 2:
        raise ValueError('Tail bound needs N>=2')
    low = variation(n,2)
    return {'N':n,'S_N':low,'S_infinity_lower':low,
            'S_infinity_upper':low*math.exp(2/(n-1)),
            'log_tail_bound':2/(n-1),
            'basis':'For k=2,tau=i: log q_n <= 1/binom(n,2); telescoping tail =2/(N-1).'}


def parent_key(child):
    return tuple(tuple(row[:-1]) for row in child[:-1])


def born_prefix_gap():
    a3,a4 = level(3,2),level(4,2)
    p3 = {h:abs(a)**2 for h,a in a3.items()}
    norm3 = math.fsum(p3.values())
    p3 = {h:x/norm3 for h,x in p3.items()}
    marg = {h:0. for h in a3}
    norm4 = math.fsum(abs(a)**2 for a in a4.values())
    for h,a in a4.items():
        marg[parent_key(h)] += abs(a)**2/norm4
    support = [h for h,a in a3.items() if abs(a)>0]
    return {'N3_probabilities':[p3[h] for h in support],
            'N4_marginal_on_N3':[marg[h] for h in support],
            'TV':sum(abs(marg[h]-p3[h]) for h in p3)/2}


def complex_pair(z):
    return [float(z.real),float(z.imag)]


@lru_cache(maxsize=1)
def report():
    rows = []
    for k in [1,2]:
        for n in range(1,6):
            law = level(n,k)
            values = np.array(list(law.values()))
            rows.append({'k':k,'N':n,'labelled_orders':len(law),
                         'nonzero_cylinders':int(np.count_nonzero(values)),
                         'amplitude_sum_error':float(abs(values.sum()-1)),
                         'S_N_enumerated':float(np.abs(values).sum()),
                         'S_N_product':variation(n,k),
                         'diagonal_sum':float(np.vdot(values,values).real)})
    anti = np.zeros((2,2),dtype=bool)
    groups = {}
    for subset in csg.ideals(anti):
        h = growth.key(csg.append_event(anti,subset))
        groups.setdefault(growth.canonical(h),[]).append(transition(anti,subset,1))
    return {'round':243,'date':'2026-09-22',
            'literature':{'family_and_extension':'https://arxiv.org/html/2003.11311',
                          'vector_measure_interface':'https://arxiv.org/html/1007.2725'},
            'family':'t0=1, t_k=i, all other t_j=0; ordinary labelled causal-set growth tree',
            'finite_levels':rows,
            'finite_born_renormalization_counterexample':born_prefix_gap(),
            'multiplicity_counterexample':{
                'parent':'two-event antichain; k=1,t1=i',
                'sum_all_labelled_children':complex_pair(sum(sum(v) for v in groups.values())),
                'sum_one_child_per_isomorphism_class':complex_pair(sum(v[0] for v in groups.values()))},
            'infinite_classification':{
                'k1':'S_N=N/prod_{n=1}^{N-1}sqrt(1+1/n^2); N/e <= S_N <= N; no scalar countably additive extension.',
                'k2':'Exact product has finite limit; cited bounded-variation theorem gives unique complex measure extension.',
                'k2_certificate':extension_bounds(),
                'numerically_proved_infinite_extension':False},
            'extra_inputs':['Specified complex CSG transition family and natural labelled event tree',
                            'Scalar rank-one decoherence functional',
                            'Full cylinder-generated sigma-algebra as the extension target'],
            'not_claimed':['This scalar history model represents arbitrary quantum matter or generates spacetime',
                           'A quantum measure on every event is an ordinary outcome probability',
                           'Isomorphic labelled cylinder events may be deleted from the sum rule'],
            'runtime':{'python':platform.python_version(),'numpy':np.__version__}}


class ExtensionTests(unittest.TestCase):
    def test_transition_sum_and_absolute_sum_for_all_small_parents(self):
        for n in range(5):
            for h in csg.all_natural_orders(n):
                r = growth.matrix(h)
                for k,tau in [(1,1j),(2,1j),(2,1.),(2,.4+.3j)]:
                    weights = [transition(r,s,k,tau) for s in csg.ideals(r)]
                    self.assertLess(abs(sum(weights)-1),1e-13)
                    self.assertAlmostEqual(sum(abs(w) for w in weights),absolute_transition_sum(n,k,tau))

    def test_finite_amplitudes_and_cross_step_biadditivity(self):
        for k in [1,2]:
            for n in range(1,5):
                parents,children = level(n,k),level(n+1,k)
                sums = {h:0j for h in parents}
                for child,a in children.items():
                    sums[parent_key(child)] += a
                a = np.array(list(parents.values()))
                b = np.array([sums[h] for h in parents])
                np.testing.assert_allclose(a,b,atol=1e-13)
                np.testing.assert_allclose(np.outer(a.conj(),a),np.outer(b.conj(),b),atol=1e-13)
                self.assertLess(abs(sum(parents.values())-1),1e-13)

    def test_rank_one_strong_positivity_and_normalization(self):
        for k in [1,2]:
            a = np.array(list(level(4,k).values()))
            d = np.outer(a.conj(),a)
            self.assertGreater(np.linalg.eigvalsh(d).min(),-1e-13)
            self.assertEqual(np.linalg.matrix_rank(d,tol=1e-12),1)
            self.assertLess(abs(d.sum()-1),1e-13)

    def test_same_amplitude_for_all_isomorphic_N5_cylinders(self):
        for k in [1,2]:
            representatives = {}
            for h,a in level(5,k).items():
                canonical = growth.canonical(h)
                if canonical in representatives:
                    self.assertLess(abs(a-representatives[canonical]),1e-13)
                representatives[canonical] = a

    def test_exact_variation_product_matches_independent_enumeration(self):
        for row in report()['finite_levels']:
            self.assertAlmostEqual(row['S_N_enumerated'],row['S_N_product'])
        self.assertEqual([len(level(n,2)) for n in range(6)],[1,1,2,7,40,357])

    def test_label_multiplicity_cannot_be_removed_in_this_sample_space(self):
        row = report()['multiplicity_counterexample']
        np.testing.assert_allclose(row['sum_all_labelled_children'],[1,0],atol=1e-14)
        np.testing.assert_allclose(row['sum_one_child_per_isomorphism_class'],[.6,-.2],atol=1e-14)

    def test_square_modulus_renormalization_breaks_prefix_consistency(self):
        row = born_prefix_gap()
        np.testing.assert_allclose(row['N3_probabilities'],[.5,.5],atol=1e-14)
        np.testing.assert_allclose(row['N4_marginal_on_N3'],[.4,.6],atol=1e-14)
        self.assertAlmostEqual(row['TV'],.1)

    def test_explicit_tail_and_divergence_bounds(self):
        for n in [2,10,100,10000]:
            self.assertGreaterEqual(variation(n,1),n/math.e)
            self.assertLessEqual(variation(n,1),n+1e-10)
        bound = extension_bounds(1000)
        further = variation(10000,2)
        self.assertLessEqual(bound['S_N'],further)
        self.assertLessEqual(further,bound['S_infinity_upper'])

    def test_singular_denominators_are_rejected(self):
        with self.assertRaises(ValueError):
            transition(np.zeros((2,2),dtype=bool),(),2,-1.)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ExtensionTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = dict(report(),checks={'run':checks.testsRun,'failures':0,'errors':0})
    if args.write_results:
        target = Path(__file__).with_name('complex_growth_extension_audit_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8')) != data:
            raise RuntimeError('Existing result differs; inspect before replacing.')
        target.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(data,ensure_ascii=False,indent=2))
