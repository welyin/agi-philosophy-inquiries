"""Round 244: first covariant observable in the extendible t2=i model.

Compute the amplitude generating function for the total number of minimal
events. Coefficients and their squared moduli are not an ordinary probability
distribution. Tail certificates refer to the existing scalar complex measure.
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
import random_order_growth_audit as growth


def root_transition(n):
    return 1/(1+1j*math.comb(n,2))


def coefficients(n, max_minima=None):
    if n < 2:
        raise ValueError('The audited branch starts from two minimal events')
    size = n if max_minima is None else min(n,max_minima)
    a = np.zeros(size+1,dtype=complex)
    a[2] = 1
    for j in range(2,n):
        q = root_transition(j)
        shifted = np.zeros_like(a)
        shifted[1:] = a[:-1]
        a = (1-q)*a+q*shifted
    return a


def minima_count(h):
    r = growth.matrix(h)
    return int(np.count_nonzero(r.sum(axis=0)==0))


def enumerated_coefficients(n):
    out = np.zeros(n+1,dtype=complex)
    for h,a in ext.level(n,2).items():
        out[minima_count(h)] += a
    return out


def exactly_two_quantum_measure(n):
    return math.exp(-math.fsum(math.log1p(1/math.comb(j,2)**2) for j in range(2,n)))


def two_minima_bounds(n=1000):
    if n < 3:
        raise ValueError('Integral tail certificate needs N>=3')
    upper = exactly_two_quantum_measure(n)
    log_tail = 4/(3*(n-2)**3)
    return {'N':n,'quantum_measure_lower':upper*math.exp(-log_tail),
            'quantum_measure_upper':upper,'log_tail_bound':log_tail,
            'ordinary_outcome_probability':False}


def coefficient_tail_bound(n):
    return ext.variation(n,2)*math.expm1(4/(n-1))


@lru_cache(maxsize=1)
def report():
    cutoff = 50000
    a = coefficients(cutoff,8)
    epsilon = coefficient_tail_bound(cutoff)
    rows = []
    for m in range(2,9):
        x = a[m]
        rows.append({'minimal_events':m,'amplitude_real':float(x.real),'amplitude_imag':float(x.imag),
                     'quantum_measure_estimate':float(abs(x)**2),
                     'amplitude_absolute_tail_bound':epsilon,
                     'quantum_measure_absolute_tail_bound':2*abs(x)*epsilon+epsilon**2})
    pair_interference = float(2*(a[2].conjugate()*a[4]).real)
    union = float(abs(a[2]+a[4])**2)
    individual = float(abs(a[2])**2+abs(a[4])**2)
    return {'round':244,'date':'2026-09-22',
            'model':'Extendible scalar complex sequential growth: t0=1,t2=i',
            'observable':'Total number M of minimal elements in the infinite causet; invariant under natural relabelling',
            'finite_enumeration':[{'N':n,'max_coefficient_difference':float(np.max(np.abs(coefficients(n)-enumerated_coefficients(n))))}
                                  for n in range(2,6)],
            'N_cutoff':cutoff,'low_coefficients':rows,
            'exactly_two_minima':two_minima_bounds(),
            'exactly_three_minima':{
                'amplitude_ratio_to_M2':'-2i, exactly in the infinite extension',
                'quantum_measure_lower':4*two_minima_bounds()['quantum_measure_lower'],
                'quantum_measure_upper':4*two_minima_bounds()['quantum_measure_upper'],
                'ordinary_outcome_probability':False},
            'covariant_nonadditivity':{
                'mu_M2_union_M4_estimate':union,
                'mu_M2_plus_mu_M4_estimate':individual,
                'interference_estimate':pair_interference,
                'interference_tail_error_bound':2*((abs(a[2])+abs(a[4]))*epsilon+epsilon**2)},
            'analytic_support':{
                'M1_amplitude':'0: the N=2 chain has zero cylinder amplitude',
                'M_infinite_total_variation':'0: bounded by 2*S_infinity/(N-1) for arbitrarily large N',
                'M_ge_2_and_finite':'full complex measure; not an assertion of classical almost-sure sampling'},
            'extra_inputs':['The particular t2=i coupling branch selected as a worked example',
                            'The infinite scalar complex-measure extension established in round 243'],
            'not_claimed':['Squared coefficients are mutually exclusive outcome probabilities',
                           'This toy model selects two physical observers, dimensions, or a realistic cosmology',
                           'Nonzero quantum measure supplies a complete interpretation or a measurement instrument'],
            'runtime':{'python':platform.python_version(),'numpy':np.__version__}}


class MinimaTests(unittest.TestCase):
    def test_generating_function_matches_labelled_causet_enumeration(self):
        for n in range(2,6):
            np.testing.assert_allclose(coefficients(n),enumerated_coefficients(n),atol=1e-13)

    def test_complete_finite_coefficient_sum_is_one(self):
        for n in [2,3,5,16,64]:
            self.assertLess(abs(coefficients(n).sum()-1),1e-13)

    def test_truncation_does_not_change_low_coefficients(self):
        np.testing.assert_allclose(coefficients(80,8),coefficients(80)[:9],atol=1e-14)

    def test_two_minima_product_and_phase_are_separated(self):
        for n in [3,5,64]:
            self.assertAlmostEqual(abs(coefficients(n,2)[2])**2,exactly_two_quantum_measure(n))
        self.assertGreater(abs(coefficients(64,2)[2].imag),.1)

    def test_real_infinite_tail_certificate_contains_further_products(self):
        b = two_minima_bounds(100)
        further = exactly_two_quantum_measure(10000)
        self.assertLessEqual(b['quantum_measure_lower'],further)
        self.assertLessEqual(further,b['quantum_measure_upper'])

    def test_l1_amplitude_tail_bound_on_two_complete_finite_levels(self):
        small,large = coefficients(64),coefficients(256)
        padded = np.pad(small,(0,len(large)-len(small)))
        self.assertLess(np.abs(large-padded).sum(),coefficient_tail_bound(64))

    def test_covariant_partition_has_nonzero_interference(self):
        row = report()['covariant_nonadditivity']
        self.assertGreater(abs(row['interference_estimate']),.1)
        self.assertAlmostEqual(row['mu_M2_union_M4_estimate']-row['mu_M2_plus_mu_M4_estimate'],
                               row['interference_estimate'])
        self.assertGreater(abs(row['interference_estimate']),row['interference_tail_error_bound'])

    def test_adjacent_amplitudes_are_in_quadrature_and_mu_M3_exceeds_one(self):
        for n in [3,20,1000]:
            a = coefficients(n,3)
            np.testing.assert_allclose(a[3]/a[2],-1j*(2-2/(n-1)),atol=1e-12)
            self.assertAlmostEqual((a[2].conjugate()*a[3]).real,0)
        self.assertGreater(report()['exactly_three_minima']['quantum_measure_lower'],1)

    def test_no_single_minimum_and_summable_new_minimum_amplitudes(self):
        self.assertEqual(coefficients(20)[1],0)
        for start in [3,20,100]:
            tail = math.fsum(abs(root_transition(n)) for n in range(start,10000))
            self.assertLessEqual(tail,2/(start-1))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(MinimaTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = dict(report(),checks={'run':checks.testsRun,'failures':0,'errors':0})
    if args.write_results:
        target = Path(__file__).with_name('covariant_minima_audit_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8')) != data:
            raise RuntimeError('Existing result differs; inspect before replacing.')
        target.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(data,ensure_ascii=False,indent=2))
