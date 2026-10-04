"""Round 240: exact two-order support in the regular Rideout-Sorkin CSG family.

Use the established CSG transition formula. A padded standard example supplies
an obstruction for every positive t_k with k>=2. Remaining support is forests.
This is classical model analysis, not a derivation of physical or quantum growth.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
import random_order_growth_audit as growth
import causal_dimension_audit as dim


def lam(size, maxima, couplings):
    return sum((F(math.comb(size-maxima, k-maxima))*couplings[k]
                for k in range(maxima, size+1)), F(0))


def is_ideal(r, subset):
    return all(not r[i,j] or i in subset for j in subset for i in range(len(r)))


def ideals(r):
    for mask in range(1 << len(r)):
        subset = tuple(i for i in range(len(r)) if (mask >> i) & 1)
        if is_ideal(r, subset):
            yield subset


def maximal_count(r, subset):
    return sum(not any(r[i,j] for j in subset) for i in subset)


def transition(r, subset, couplings):
    return lam(len(subset), maximal_count(r, subset), couplings)/lam(len(r), 0, couplings)


def append_event(r, subset):
    new = np.zeros((len(r)+1, len(r)+1), dtype=bool)
    new[:-1, :-1] = r
    new[list(subset), -1] = True
    return new


@lru_cache(maxsize=None)
def all_natural_orders(n):
    if n == 0:
        return ((),)
    out = []
    for k in all_natural_orders(n-1):
        r = growth.matrix(k)
        out.extend(growth.key(append_event(r,s)) for s in ideals(r))
    return tuple(out)


def path_probability(r, couplings):
    out = F(1)
    for n in range(len(r)):
        subset = tuple(np.flatnonzero(r[:n,n]))
        out *= transition(r[:n,:n], subset, couplings)
    return out


def invariant_probability(r, couplings):
    numerator, denominator = F(1), F(1)
    for j in range(len(r)):
        subset = tuple(np.flatnonzero(r[:,j]))
        numerator *= lam(len(subset), maximal_count(r,subset), couplings)
        denominator *= lam(j,0,couplings)
    return numerator/denominator


def padded_obstruction(k):
    # k+1 minimal points: a0,a1,a2 plus k-2 shared padding points.
    n = k+4
    r = np.zeros((n,n), dtype=bool)
    for j in range(3):
        for i in range(k+1):
            if i != j:
                r[i,k+1+j] = True
    return r


def parents(r):
    ps = []
    for j in range(len(r)):
        past = tuple(np.flatnonzero(r[:,j]))
        maximal = [i for i in past if not any(r[i,h] for h in past)]
        if len(maximal)>1:
            raise ValueError('Not a rooted forest order')
        ps.append(maximal[0] if maximal else -1)
    return ps


def forest_realizer(r):
    ps = parents(r)
    children = {i:[] for i in range(-1,len(r))}
    for j,p in enumerate(ps):
        children[p].append(j)
    def traversal(reverse):
        out = []
        def visit(i):
            if i != -1:
                out.append(i)
            for j in sorted(children[i], reverse=reverse):
                visit(j)
        visit(-1)
        return out
    orders = [traversal(False),traversal(True)]
    ranks = np.array([[order.index(i) for i in range(len(r))] for order in orders]).T
    return dim.relation_from_points(ranks)


def forest_expected_relations(n, a=F(1)):
    value = F(0)
    for size in range(n):
        value += a*(value+size)/(1+size*a)
    return value


@lru_cache(maxsize=1)
def report():
    cases = []
    for k in range(2,7):
        r = padded_obstruction(k)
        t = [F(1),F(1)]+[F(0)]*(k+2)
        t[k] = F(2)
        probability = path_probability(r,t)
        denominator = math.prod(lam(n,0,t) for n in range(k+4))
        cases.append({'positive_coupling_index':k,'events':k+4,
                      't_k':'2','labeled_history_probability':str(probability),
                      'formula_probability':str(t[k]**3/denominator),
                      'induced_S3_vertices':[0,1,2,k+1,k+2,k+3]})
    forest = [F(1),F(1)]+[F(0)]*5
    all5 = all_natural_orders(5)
    supported5 = [k for k in all5 if path_probability(growth.matrix(k),forest)>0]
    means = []
    for n in [16,64,256,1024]:
        exact = forest_expected_relations(n)
        means.append({'N':n,'expected_C2':float(exact),
                      'expected_r':float(exact/F(math.comb(n,2))),
                      'target_flat_diamond_expected_r':.5})
    return {'round':240,'date':'2026-09-22',
            'family':'Regular CSG: t0=1, finite nonnegative t_k, lambda transition probabilities',
            'question':'Can regular CSG retain exact two-order support while allowing events to join independent past branches?',
            'conditional_theorem':'All finite outputs have order dimension <=2 with probability one iff t_k=0 for every k>=2.',
            'padded_standard_example_witnesses':cases,
            'forest_support_check':{'all_naturally_labeled_N5_orders':len(all5),
                                    'positive_weight_forest_orders':len(supported5),
                                    'all_positive_outputs_have_two_order_realizer':True},
            'forest_expectations_a1':means,
            'small_join_probability':{'t0':'1','t1':'1','t2':'1',
                                      'CSG_three_event_join':'1/8','target_flat_diamond_join':'1/6'},
            'extra_inputs':['Classical sequential growth formula and its causality/covariance conditions',
                            'Exact two-order support at every finite stage, stronger than a continuum limit'],
            'conclusion':'Exact two-order support leaves only rooted forests in this family; their relation fraction vanishes in probability for any fixed finite t1.',
            'not_claimed':['CSG cannot have any continuum approximation',
                           'A nonzero finite forbidden-order probability rules out coarse-grained spacetime',
                           'Cognitive composability forces classical Bell causality'],
            'runtime':{'python':platform.python_version(),'numpy':np.__version__}}


class CSGTests(unittest.TestCase):
    def test_existing_formula_normalizes_every_small_parent(self):
        families = [[F(1)]*7,[F(1),F(1)]+[F(0)]*5,
                    list(map(F,[1,2,0,3,1,0,0]))]
        self.assertEqual([len(all_natural_orders(n)) for n in range(6)],[1,1,2,7,40,357])
        for n in range(6):
            for k in all_natural_orders(n):
                r = growth.matrix(k)
                for t in families:
                    weights = [transition(r,s,t) for s in ideals(r)]
                    self.assertEqual(sum(weights),1)
                    self.assertGreaterEqual(min(weights),0)

    def test_covariance_by_all_N5_natural_labels(self):
        t = list(map(F,[1,2,3,1,0,0]))
        groups = {}
        for k in all_natural_orders(5):
            r = growth.matrix(k)
            probability = path_probability(r,t)
            self.assertEqual(probability,invariant_probability(r,t))
            c = growth.canonical(k)
            if c in groups:
                self.assertEqual(probability,groups[c])
            groups[c] = probability

    def test_bell_ratio_ignores_explicit_spectator(self):
        r = growth.matrix(((False,True,False),(False,False,False),(False,False,False)))
        small = r[:2,:2]
        t = list(map(F,[1,2,3,4]))
        self.assertEqual(transition(r,(0,1),t)/transition(r,(0,),t),
                         transition(small,(0,1),t)/transition(small,(0,),t))

    def test_padded_obstructions_have_positive_weight_and_induced_S3(self):
        s3 = padded_obstruction(2)
        for row in report()['padded_standard_example_witnesses']:
            k = row['positive_coupling_index']
            self.assertGreater(F(row['labeled_history_probability']),0)
            self.assertEqual(row['labeled_history_probability'],row['formula_probability'])
            r = padded_obstruction(k)
            ids = row['induced_S3_vertices']
            np.testing.assert_array_equal(r[np.ix_(ids,ids)],s3)
            for j in range(k+1,k+4):
                self.assertEqual(int(r[:,j].sum()),k)

    def test_remaining_support_is_forest_with_constructive_two_orders(self):
        t = [F(1),F(1)]+[F(0)]*5
        for k in all_natural_orders(5):
            r = growth.matrix(k)
            if path_probability(r,t)>0:
                np.testing.assert_array_equal(forest_realizer(r),r)
                self.assertTrue(all(maximal_count(r,tuple(np.flatnonzero(r[:,j])))<=1 for j in range(5)))

    def test_forest_expectation_against_full_distribution_and_closed_form(self):
        for a in [F(0),F(1,2),F(1),F(2)]:
            t = [F(1),a]+[F(0)]*5
            for n in range(2,6):
                exact = sum((path_probability(growth.matrix(k),t)*sum(map(sum,k))
                             for k in all_natural_orders(n)),F(0))
                self.assertEqual(exact,forest_expected_relations(n,a))
        for n in [2,3,16,64]:
            h = sum((F(1,i) for i in range(1,n+1)),F(0))
            self.assertEqual(forest_expected_relations(n),(n+1)*h-2*n)

    def test_join_has_positive_weight_only_with_t2_at_three_events(self):
        r = np.zeros((3,3),dtype=bool)
        r[0,2]=r[1,2]=True
        for a,b in [(F(1),F(0)),(F(1),F(1)),(F(2),F(3))]:
            probability = path_probability(r,[F(1),a,b])
            self.assertEqual(probability,b/((1+a)*(1+2*a+b)))
        self.assertEqual(path_probability(r,[F(1)]*3),F(1,8))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CSGTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = dict(report(),checks={'run':checks.testsRun,'failures':0,'errors':0})
    if args.write_results:
        target = Path(__file__).with_name('csg_support_audit_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8')) != data:
            raise RuntimeError('Existing result differs; inspect before replacing.')
        target.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(data,ensure_ascii=False,indent=2))
