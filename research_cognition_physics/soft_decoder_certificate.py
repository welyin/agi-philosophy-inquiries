"""Round 163: a strict global improvement over the certified old recovery.

The new upper bound covers every real effective density matrix, hence all
promised independent sources and external extensions. The lower witness is
a legal independent source. No optimality of the new decoder is claimed.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import unittest

import numpy as np

from certified_intervals import Interval as I, SCALE
from global_history_recovery import LOWER as OLD_LOWER
from joint_history_certificate import positive_ldl, stereographic_bloch, exact_product_sigma
from joint_relation_records import interval_fields
from noisy_classical_history import majority_error
from soft_history_decoder import (PARAMETER, rational_circle, branch_errors, branch_weights,
    supporting_coefficients, gram_from_sigma, schur_root)


DENOMINATOR = 1 << 40
UPPER = F(450678542, 10**9)
LOWER = F(449863269918, 10**12)
UPPER_INTEGERS = (76059903675, 1032253769, 17240626388, 4926202340, 9688559780,
    3500760390, 60383776, 76509393374, 861889912, 13170559195, 3590251484,
    9630015688, 12614911677, 1091759124, 101603734349, 13442104640, 3366859916,
    1238527508, 3666003468, 9950899858, 13442902795, 101803666422, 372745888,
    4376538885, 9623968928, 3590251484, 13183908478, 874869195, 76509393374,
    39229249, 3500558779, 9696374875, 4924414134, 16792693909, 1005880917, 75910473652)
LOWER_INTEGERS = (76025749428, 595884450, 17171745611, 4951692352, 9884391182,
    3638185577, -582415709, 76509393374, 210349530, 13517282767, 3764973933,
    9873776245, 12721251393, 623052926, 101652287020, 13817396382, 3468018154,
    1239797830, 3678830794, 10184932526, 13823610231, 101809552038, 228959486,
    4388331711, 9857807892, 3764973933, 13613033997, 279342519, 76509393374,
    -602312637, 3636615983, 9912088033, 4937770702, 16819405399, 551475459, 75908208931)
STEREOGRAPHIC = (-1574528, 1457762, -1451869, -1200266)


def coefficients(integers):
    if len(integers) != 36 or any(type(n) is not int for n in integers):
        raise ValueError('Use 36 exact integer coefficient numerators.')
    return tuple(F(n, DENOMINATOR) for n in integers)


def residual_blocks(integers, parameter=PARAMETER):
    a = coefficients(integers)
    residual = np.eye(8, dtype=object)
    for weight, w in zip(a, branch_errors(parameter)):
        residual = residual-weight*w
    return 1-sum(a), residual


def partial_reference(matrix):
    return np.trace(matrix.reshape(2, 4, 2, 4), axis1=0, axis2=2)


def block_expectation(sigma, scalar, negative):
    partial = partial_reference(negative)
    return scalar/2+sum(sigma[j, i]*partial[i, j] for i in range(4) for j in range(4))/4


def dual_data(integers=UPPER_INTEGERS, parameter=PARAMETER):
    a = coefficients(integers)
    p = branch_weights(majority_error(1))
    d = 1-sum(I.exact(x*x)/weight for x, weight in zip(a, p))
    if d.lo <= 0:
        raise ValueError('The positive completion denominator failed.')
    b, residual = residual_blocks(integers, parameter)
    numerator = b*b*np.eye(4, dtype=object)+partial_reference(residual.T @ residual)/2
    return d, numerator


def certify_upper(upper=UPPER):
    d, numerator = dual_data()
    matrix = [[2*upper*d*int(i == j)-numerator[i, j] for j in range(4)] for i in range(4)]
    pivots = positive_ldl(matrix)
    return {'upper_exact': str(upper), 'positive_denominator': interval_fields(d),
            'positive_ldl_pivots': [interval_fields(x) for x in pivots],
            'dual_numerator_exact': [[str(x) for x in row] for row in numerator]}


def witness_states():
    p = tuple(F(n, 10**6) for n in STEREOGRAPHIC)
    return stereographic_bloch(*p[:2]), stereographic_bloch(*p[2:])


def rayleigh_lower(sigma, integers=LOWER_INTEGERS):
    b, residual = residual_blocks(integers)
    norm = block_expectation(sigma, b*b, residual.T @ residual)
    if norm <= 0:
        raise ValueError('The external test vector has no positive norm.')
    target = block_expectation(sigma, b, residual.T)
    overlaps = [block_expectation(sigma, b, residual.T @ w) for w in branch_errors()]
    p = branch_weights(majority_error(1))
    return (target*target-sum(weight*x*x for weight, x in zip(p, overlaps)))/norm


def noiseless_interval():
    c, s = rational_circle()
    return (s*s+(I.exact(s**4+4*(1-c)**2)).sqrt())/6


@lru_cache(None)
def certificate():
    upper = certify_upper()
    a, b = witness_states()
    if sum(x*x for x in a) != 1 or sum(x*x for x in b) != 1:
        raise AssertionError('The lower witness is not exactly a pair of pure states.')
    lower = rayleigh_lower(exact_product_sigma(a, b))
    if F(lower.lo, SCALE) <= LOWER or UPPER >= OLD_LOWER:
        raise AssertionError('The claimed lower bound or strict improvement failed.')
    return {'round': 163, 'parameter_exact': str(PARAMETER),
            'new_strict_lower_exact': str(LOWER), 'new_strict_upper_exact': str(UPPER),
            'old_certified_strict_lower_exact': str(OLD_LOWER),
            'strict_worst_error_improvement_at_least_exact': str(OLD_LOWER-UPPER),
            'upper_certificate': upper, 'legal_witness_lower': interval_fields(lower),
            'coefficient_denominator': DENOMINATOR, 'upper_coefficient_numerators': UPPER_INTEGERS,
            'lower_coefficient_numerators': LOWER_INTEGERS,
            'witness_stereographic_numerators_over_one_million': STEREOGRAPHIC,
            'legal_witness_bloch_vectors_exact': [[str(x) for x in row] for row in (a, b)],
            'noiseless_error': interval_fields(noiseless_interval()),
            'all_independent_sources_and_external_extensions_covered': True,
            'new_upper_uses_the_larger_all_real_effective_states_domain': True,
            'new_actual_worst_error_exactly_solved': False,
            'decoder_parameter_or_general_recovery_globally_optimal': False,
            'three_reads_same_reports_and_one_coherent_history_bit': True,
            'certificate_requires_scipy': False}


class SoftDecoderCertificateTests(unittest.TestCase):
    def test_uniform_upper_has_strict_interval_positive_pivots(self):
        result = certify_upper()
        self.assertEqual(len(result['positive_ldl_pivots']), 4)
        self.assertLess(UPPER, OLD_LOWER)

    def test_legal_lower_witness_and_full_gram_agree(self):
        a, b = witness_states()
        self.assertEqual(sum(x*x for x in a), 1)
        self.assertEqual(sum(x*x for x in b), 1)
        sigma = exact_product_sigma(a, b)
        lower = rayleigh_lower(sigma)
        self.assertGreater(F(lower.lo, SCALE), LOWER)
        self.assertAlmostEqual(lower.floats()[0], schur_root(gram_from_sigma(np.array(sigma, dtype=float)))[0], places=12)

    def test_global_support_on_unseen_relaxed_real_states(self):
        d, numerator = dual_data()
        h = np.array(numerator, dtype=float)/(2*np.mean(d.floats()))
        rng = np.random.default_rng(163)
        for rank in (1, 2, 4):
            c = rng.normal(size=(4, rank)); sigma = c @ c.T/np.sum(c*c)
            actual = schur_root(gram_from_sigma(sigma))[0]
            self.assertLessEqual(actual, np.trace(h @ sigma)+2e-14)
            self.assertLess(np.trace(h @ sigma), float(UPPER))

    def test_support_expectation_matches_direct_unweighted_gram(self):
        from soft_history_decoder import gram_kernel
        sigma = np.eye(4)/4
        gram = np.einsum('abij,ji->ab', gram_kernel(), sigma)
        v = np.r_[-np.array(coefficients(LOWER_INTEGERS), dtype=float), 1.]
        overlap = gram @ v
        b, residual = residual_blocks(LOWER_INTEGERS)
        self.assertAlmostEqual(float(block_expectation(sigma, b*b, residual.T @ residual)), v @ overlap)
        for j, w in enumerate(branch_errors()):
            self.assertAlmostEqual(float(block_expectation(sigma, b, residual.T @ w)), overlap[j], places=13)

    def test_corrupt_or_overstrong_certificates_are_rejected(self):
        with self.assertRaises(ValueError): coefficients([0.]*36)
        with self.assertRaises(ValueError): dual_data((2*DENOMINATOR,)+(0,)*35)
        with self.assertRaises(ValueError): certify_upper(F(45, 100))

    def test_improvement_compares_new_upper_to_old_legal_lower(self):
        result = certificate()
        self.assertGreater(F(result['strict_worst_error_improvement_at_least_exact']), F(2311, 10**6))
        self.assertFalse(result['new_actual_worst_error_exactly_solved'])
        self.assertFalse(result['decoder_parameter_or_general_recovery_globally_optimal'])

    def test_noiseless_price_is_strictly_positive_and_certified(self):
        from soft_history_decoder import noiseless_error
        value = noiseless_interval()
        self.assertGreater(F(value.lo, SCALE), F(8, 1000))
        self.assertAlmostEqual(value.floats()[0], noiseless_error(), places=14)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SoftDecoderCertificateTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = certificate()
    report['automated_checks'] = {'run': checks.testsRun, 'failures': len(checks.failures), 'errors': len(checks.errors)}
    if args.write_results:
        Path(__file__).with_name('soft_decoder_certificate_results.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: report[k] for k in ('round', 'new_strict_lower_exact', 'new_strict_upper_exact',
        'old_certified_strict_lower_exact', 'strict_worst_error_improvement_at_least_exact', 'automated_checks')}, indent=2))


if __name__ == '__main__':
    main()
