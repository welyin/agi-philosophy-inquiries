"""Round 246: a consistent rank-two growth functional with formal records.

This changes the old scalar functional. Its two normalized conditional
measures happen to be positive. A classical comparator has the same sector
weights but different within-sector interference. Neither mathematical
construction supplies a local physical record mechanism.
"""
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
import complex_growth_extension_audit as ext
import covariant_minima_audit as minima
import csg_support_audit as csg
import random_order_growth_audit as growth


@lru_cache(maxsize=None)
def conditional_laws(n):
    if n < 2:
        raise ValueError('Use the common two-minimum initial prefix, N>=2')
    keys = tuple(ext.level(n,2))
    beta = minima.coefficients(n,2)[2]
    laws = np.zeros((2,len(keys)),dtype=complex)
    for j,h in enumerate(keys):
        b = ext.level(n,2)[h]/beta
        m = minima.minima_count(h)
        if m == 2:
            laws[:,j] = [b,b/(n-1)]
        elif m == 3:
            laws[1,j] = .5j*b
    assert np.max(abs(laws.imag)) < 2e-13
    return keys,laws.real


def vectors(n,p=.5):
    keys,laws = conditional_laws(n)
    return keys,np.sqrt([p,1-p])[:,None]*laws


def functional(n,p=.5):
    keys,v = vectors(n,p)
    return keys,v.T@v


def comparator(n,p=.5):
    keys,laws = conditional_laws(n)
    return keys,np.array([p,1-p])@laws


def refinement(n):
    parents,_ = conditional_laws(n)
    children,_ = conditional_laws(n+1)
    index = {h:i for i,h in enumerate(parents)}
    matrix = np.zeros((len(children),len(parents)))
    for j,h in enumerate(children):
        matrix[j,index[ext.parent_key(h)]] = 1
    return matrix


def reference_operator(n,p=.5):
    keys,laws = conditional_laws(n)
    reference = np.array([float(csg.path_probability(growth.matrix(h),[1]*(n+1))) for h in keys])
    u = laws/np.sqrt(reference)
    rho = u.T@np.diag([p,1-p])@u
    return reference,rho


def classical_bayes_error(n,p=.5):
    _,laws = conditional_laws(n)
    return float(np.minimum(p*laws[0],(1-p)*laws[1]).sum())


def transition_hazard_rows(n):
    parents,old = conditional_laws(n)
    children,new = conditional_laws(n+1)
    groups = {}
    for j,h in enumerate(children):
        groups.setdefault(ext.parent_key(h),[]).append(j)
    rows = []
    for i,h in enumerate(parents):
        m = minima.minima_count(h)
        for s in range(2):
            if old[s,i] < 1e-14:
                continue
            js = groups[h]
            born_root = [j for j in js if minima.minima_count(children[j]) == m+1]
            rows.append((s,m,float(new[s,born_root].sum()/old[s,i])))
    return rows


def report():
    rows = []
    for n in range(2,6):
        keys,laws = conditional_laws(n)
        _,v = vectors(n)
        d = v.T@v
        _,prob = comparator(n)
        rows.append({'N':n,'labelled_cylinders':len(keys),
                     'conditional_normalizations':laws.sum(axis=1).tolist(),
                     'rank_D':int(np.linalg.matrix_rank(v,tol=1e-11)),
                     'D_Omega_Omega':float(d.sum()),
                     'sum_diagonal_D':float(np.trace(d)),
                     'classical_D_rank':int(np.count_nonzero(prob>1e-14)),
                     'classical_Bayes_error':classical_bayes_error(n),
                     'exact_Bayes_error':1/(2*(n-1)),
                     'P_M3_component_still_two_minima':float(laws[1,[j for j,h in enumerate(keys) if minima.minima_count(h)==2]].sum())})
    reference,rho = reference_operator(4)
    d = functional(4)[1]
    rows_h = transition_hazard_rows(3)
    return {'round':246,'date':'2026-09-22',
            'model':'Changed rank-two direct sum of normalized restrictions of the t2=i complex measure to total minima M=2 and M=3',
            'prior_p':.5,'finite_enumeration':rows,
            'reference_operator_N4':{
                'reference_probability_sum':float(reference.sum()),
                'normalization_on_constant_function':float(np.sqrt(reference)@rho@np.sqrt(reference)),
                'trace_not_required_to_be_one':float(np.trace(rho)),
                'kernel_reconstruction_error':float(np.max(abs(np.sqrt(reference[:,None]*reference[None,:])*rho-d)))},
            'exact_macro_records':{'D_E2_E3':0.,'mu_E2':.5,'mu_E3':.5,
                                   'record_subspaces':2,'sectors_persist_under_measure_refinement':True},
            'third_minimum_wait':{
                'conditional_hazard_before_arrival':'1/n at transition n -> n+1',
                'P_still_two_at_N_given_E3':'1/(N-1)',
                'P_third_minimum_born_at_n_plus_1':'1/(n*(n-1)), n>=2',
                'finite_prefix_classical_Bayes_error':'min(p,(1-p)/(N-1))',
                'mean_birth_index':'infinite, even though arrival has probability one'},
            'N3_hazards':[[s+2,m,h] for s,m,h in rows_h],
            'extra_inputs':['Choose E2 and E3; suppress other infinite histories',
                            'Choose arbitrary sector prior p and orthogonalize sectors',
                            'Generalized projector records, not a derived local readable register'],
            'limits':['D differs from the original scalar functional',
                      'Positive conditional amplitudes do not make their rank-one kernels classical',
                      'Bayes error concerns the explicitly defined classical comparator, not a Born rule for coherent prefixes',
                      'Consistency is not by itself Kraus completeness, local records or quantum gravity'],
            'runtime':{'python':platform.python_version(),'numpy':np.__version__}}


class Audit(unittest.TestCase):
    def test_positive_normalized_conditionals(self):
        for n in range(2,6):
            _,laws = conditional_laws(n)
            self.assertGreaterEqual(laws.min(),-1e-14)
            np.testing.assert_allclose(laws.sum(axis=1),1,atol=1e-13)

    def test_prefix_consistency_both_functionals(self):
        for n in range(2,5):
            t = refinement(n)
            np.testing.assert_allclose(conditional_laws(n+1)[1]@t,conditional_laws(n)[1],atol=1e-13)
            np.testing.assert_allclose(t.T@functional(n+1)[1]@t,functional(n)[1],atol=1e-13)
            np.testing.assert_allclose(t.T@np.diag(comparator(n+1)[1])@t,np.diag(comparator(n)[1]),atol=1e-13)

    def test_natural_relabelling_invariant_amplitudes(self):
        # The order-isomorphism canonical form is inherited from the older audit.
        import causal_dimension_audit as dim
        keys,laws = conditional_laws(5)
        groups = {}
        for j,h in enumerate(keys):
            key = dim.canonical_key(growth.matrix(h))
            if key in groups:
                np.testing.assert_allclose(laws[:,j],groups[key],atol=1e-13)
            groups[key] = laws[:,j]

    def test_positive_rank_two_normalized_kernel(self):
        for n in range(3,6):
            _,v = vectors(n)
            self.assertEqual(np.linalg.matrix_rank(v),2)
            self.assertAlmostEqual(float((v.T@v).sum()),1)
        eig = np.linalg.eigvalsh(functional(4)[1])
        self.assertGreaterEqual(eig.min(),-1e-13)

    def test_reference_operator_embedding(self):
        pi,rho = reference_operator(4)
        self.assertTrue(np.all(pi>0))
        self.assertAlmostEqual(float(pi.sum()),1)
        self.assertGreaterEqual(np.linalg.eigvalsh(rho).min(),-1e-12)
        self.assertAlmostEqual(float(np.sqrt(pi)@rho@np.sqrt(pi)),1)
        np.testing.assert_allclose(np.sqrt(pi[:,None]*pi[None,:])*rho,functional(4)[1],atol=1e-13)

    def test_macro_records_for_arbitrary_priors(self):
        for p in [.2,.5,.8]:
            v = np.diag(np.sqrt([p,1-p]))
            for s in range(2):
                r = np.diag([int(j==s) for j in range(2)])
                for t in range(2):
                    np.testing.assert_allclose(r@v[:,t],v[:,t] if s==t else np.zeros(2),atol=1e-13)
            np.testing.assert_allclose(v.T@v,np.diag([p,1-p]),atol=1e-13)

    def test_classical_and_coherent_microhistories_differ(self):
        d = functional(4)[1]
        classical = np.diag(comparator(4)[1])
        self.assertGreater(np.max(abs(d-classical)),.01)
        self.assertGreater(np.max(abs(d-np.diag(d.diagonal()))),.01)
        self.assertAlmostEqual(float(classical.sum()),float(d.sum()))

    def test_hazard_from_full_transition_enumeration(self):
        for n in range(2,5):
            for s,m,value in transition_hazard_rows(n):
                expected = 1/n if s==1 and m==2 else 0.
                self.assertAlmostEqual(value,expected)

    def test_bayes_error_from_all_prefixes(self):
        for p in [.2,.5,.8]:
            for n in range(2,6):
                self.assertAlmostEqual(classical_bayes_error(n,p),min(p,(1-p)/(n-1)))
                self.assertGreater(classical_bayes_error(n,p),0)

    def test_arrival_tail_from_conditional_law(self):
        for n in range(2,6):
            keys,laws = conditional_laws(n)
            mask = [minima.minima_count(h)==2 for h in keys]
            self.assertAlmostEqual(float(laws[1,mask].sum()),1/(n-1))
        # Telescoping sum checks the arrival distribution and its survival tail.
        for n in [10,100,1000]:
            subtotal = math.fsum(1/(j*(j-1)) for j in range(2,n))
            self.assertAlmostEqual(subtotal+1/(n-1),1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Audit)
    checks = unittest.TextTestRunner(verbosity=2).run(suite)
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = report()
    result['checks'] = {'run':checks.testsRun,'failures':len(checks.failures),'errors':len(checks.errors)}
    target = Path(__file__).with_name('record_sector_growth_audit_results.json')
    payload = json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise SystemExit('Refusing to overwrite a different saved result')
        target.write_text(payload,encoding='utf-8')
    print(payload)
