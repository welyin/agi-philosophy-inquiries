"""Round 159: rigorous bounds for the fixed three-read recovery channel.

Saved rational witnesses are checked with 100-bit outward-rounded intervals.
Candidate discovery used the optional SciPy runtime; verification uses only
NumPy and exact integer/fraction arithmetic, with no optimizer tolerance.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import unittest

import numpy as np

from certified_intervals import Interval as I, SCALE
from joint_history_channel import (BIT_ERROR, bloch_state, error_from_sigma,
    gram_kernel_numerators, grouped_errors, grouped_weights, source_sigma)
from joint_relation_records import interval_fields
from noisy_classical_history import majority_error


DENOMINATOR = 1 << 40
UPPER_INTEGERS = (510771285238, 1019182055, 14145828109, 8696543109,
    19624097634, 7069396106, 1731454411, 23595442411, 7304821694,
    19652464021, 10453895973, 27159397963, 16946974460, 6224283919,
    1223397392, 3623284635, 1067567048, 13788584129)
LOWER_INTEGERS = (510771285238, -158111310, 14065056208, 8725032098,
    19869616575, 7388920099, 516363092, 24005576664, 7575672460,
    19953740967, 10621918534, 27958203165, 17330031520, 6220598806,
    242948163, 3655565519, -123196600, 13804823982)
STEREOGRAPHIC_INTEGERS = (-2123005, 1003654, -2100803, -855037)
RELAXED_VECTOR = (10008187, 40733053, 46174438, -78157379)
LOWER = F(452989542, 10**9)
UPPER = F(453436664, 10**9)
RELAXED_LOWER = F(453436662, 10**9)


def coefficients(integers):
    if len(integers) != 18 or any(type(n) is not int for n in integers):
        raise ValueError('Use 18 exact integer coefficient numerators.')
    return tuple(F(n, DENOMINATOR) for n in integers)


def stereographic_bloch(a, b):
    a, b = F(a), F(b)
    norm = 1+a*a+b*b
    return (2*a/norm, 2*b/norm, (1-a*a-b*b)/norm)


def promised_bloch_vectors():
    q = tuple(F(n, 10**6) for n in STEREOGRAPHIC_INTEGERS)
    return stereographic_bloch(*q[:2]), stereographic_bloch(*q[2:])


def exact_product_sigma(a, b):
    """Re(rho(a) tensor conjugate(rho(b))), without floating complex numbers."""
    ra = np.array([[1+a[2], a[0]], [a[0], 1-a[2]]], dtype=object)/2
    rb = np.array([[1+b[2], b[0]], [b[0], 1-b[2]]], dtype=object)/2
    j = np.array([[0, -1], [1, 0]], dtype=object)
    return np.kron(ra, rb) + a[1]*b[1]*np.kron(j, j)/4


def relaxed_sigma():
    v = np.array([F(n) for n in RELAXED_VECTOR], dtype=object)
    return np.outer(v, v)/sum(x*x for x in v)


def dual_data(integers=UPPER_INTEGERS):
    """B=diag(b I,D), d=1-sum a_j^2/p_j, N=b^2 I+Tr_R(D^T D)/2."""
    a = coefficients(integers)
    weights = grouped_weights(majority_error(1))
    if any(p.lo <= 0 for p in weights):
        raise ValueError('This certificate requires strictly positive group weights.')
    d = 1-sum(I.exact(x*x)/p for x, p in zip(a, weights))
    if d.lo <= 0:
        raise ValueError('The positive completion denominator is not certified.')
    b = 1-sum(a)
    residual = np.eye(8, dtype=object)
    for x, w in zip(a, grouped_errors()[0]):
        residual = residual-x*w
    squared = residual.T @ residual
    partial = np.trace(squared.reshape(2, 4, 2, 4), axis1=0, axis2=2)
    numerator = b*b*np.eye(4, dtype=object) + partial/2
    return d, numerator


def positive_ldl(matrix):
    """Interval LDL^T: positive pivot lower endpoints prove positive definiteness."""
    size = len(matrix)
    lower = [[I.exact(int(i == j)) for j in range(size)] for i in range(size)]
    pivots = []
    for i in range(size):
        pivot = matrix[i][i]-sum(lower[i][k]**2*pivots[k] for k in range(i))
        if pivot.lo <= 0:
            raise ValueError('A strictly positive LDL pivot could not be certified.')
        pivots.append(pivot)
        for j in range(i+1, size):
            lower[j][i] = (matrix[j][i]-sum(lower[j][k]*lower[i][k]*pivots[k]
                                          for k in range(i)))/pivot
    return pivots


def certify_upper(upper=UPPER, integers=UPPER_INTEGERS):
    d, numerator = dual_data(integers)
    residual = [[2*upper*d*int(i == j)-numerator[i, j] for j in range(4)]
                for i in range(4)]
    pivots = positive_ldl(residual)
    return {'upper_exact': str(upper), 'completion_denominator': interval_fields(d),
            'positive_ldl_pivots': [interval_fields(x) for x in pivots],
            'dual_numerator_exact': [[str(x) for x in row] for row in numerator]}


def rayleigh_lower(sigma, integers):
    """A legal external test vector B|psi>; its Rayleigh quotient is a lower bound."""
    kernel = gram_kernel_numerators()
    gram = np.array([[sum(F(int(kernel[a, b, i, j]))*sigma[j, i]
                         for i in range(4) for j in range(4))/8
                      for b in range(19)] for a in range(19)], dtype=object)
    vector = np.array(tuple(-a for a in coefficients(integers))+(F(1),), dtype=object)
    overlap = gram @ vector
    norm = vector @ overlap
    if norm <= 0:
        raise ValueError('The external test vector must have positive norm.')
    weights = grouped_weights(majority_error(1))
    numerator = overlap[-1]**2-sum(p*x*x for p, x in zip(weights, overlap[:-1]))
    return numerator/norm


@lru_cache(None)
def certificate():
    upper = certify_upper()
    a, b = promised_bloch_vectors()
    if sum(x*x for x in a) != 1 or sum(x*x for x in b) != 1:
        raise AssertionError('The saved independent source is not exactly pure.')
    witness = rayleigh_lower(exact_product_sigma(a, b), LOWER_INTEGERS)
    relaxed = rayleigh_lower(relaxed_sigma(), UPPER_INTEGERS)
    if F(witness.lo, SCALE) <= LOWER or F(relaxed.lo, SCALE) <= RELAXED_LOWER:
        raise AssertionError('A claimed strict lower endpoint did not verify.')
    # Every Re(product) is the average of two product density operators, so
    # its partial transpose is PSD. The saved relaxed pure state fails this.
    v = tuple(F(n) for n in RELAXED_VECTOR)
    partial_transpose_determinant = -(v[0]*v[3]-v[1]*v[2])**4/sum(x*x for x in v)**4
    if partial_transpose_determinant >= 0:
        raise AssertionError('The relaxed witness must be outside the source promise.')
    return {'round': 159, 'readout_contrast_exact': '1/2',
        'fixed_channel': 'round 147 once per label bit; 110->010, 111->011',
        'promised_worst_error_strict_lower_exact': str(LOWER),
        'promised_worst_error_strict_upper_exact': str(UPPER),
        'certified_bracket_width_exact': str(UPPER-LOWER),
        'upper_certificate': upper, 'promised_witness_rayleigh_lower': interval_fields(witness),
        'promised_witness_bloch_vectors_exact': [[str(x) for x in row] for row in (a, b)],
        'saved_denominator': DENOMINATOR, 'upper_coefficient_numerators': UPPER_INTEGERS,
        'lower_coefficient_numerators': LOWER_INTEGERS,
        'stereographic_numerators_over_one_million': STEREOGRAPHIC_INTEGERS,
        'relaxed_worst_error_strict_lower_exact': str(RELAXED_LOWER),
        'relaxed_witness_rayleigh_lower': interval_fields(relaxed),
        'relaxed_vector_before_normalization': RELAXED_VECTOR,
        'relaxed_witness_partial_transpose_determinant_exact': str(partial_transpose_determinant),
        'relaxed_witness_is_a_legal_independent_source': False,
        'original_source_promise_minimax_exactly_solved': False,
        'any_external_extension_included': True,
        'all_readout_outcomes_accepted': True, 'coherent_history_rebits': 1,
        'noisy_reads_and_fresh_pointers': 3, 'raw_record_bits': 3, 'decoded_message_bits': 3,
        'new_native_gates': 413, 'gates_including_initial_join': 579,
        'gates_and_coherent_storage_ideal': True,
        'source_transport_calibration_classical_compute_and_microscopic_bath_costs_included': False,
        'source_purification_pointers_and_baths_remain_in_the_total_system': True,
        'candidate_search_only': 'SciPy BFGS, floating Gram eigenvalues; not the proof',
        'certificate_requires_scipy': False}


class JointHistoryCertificateTests(unittest.TestCase):
    def test_upper_is_strictly_certified_by_four_interval_pivots(self):
        result = certify_upper()
        self.assertEqual(len(result['positive_ldl_pivots']), 4)
        self.assertLess(UPPER, F(454, 1000))

    def test_exact_pure_product_witness_and_complete_external_distance(self):
        a, b = promised_bloch_vectors()
        self.assertEqual(sum(x*x for x in a), 1)
        self.assertEqual(sum(x*x for x in b), 1)
        sigma = exact_product_sigma(a, b)
        numeric = source_sigma(bloch_state(tuple(map(float, a))), bloch_state(tuple(map(float, b))))
        np.testing.assert_allclose(np.array(sigma, dtype=float), numeric, atol=1e-16)
        value = rayleigh_lower(sigma, LOWER_INTEGERS)
        self.assertGreater(F(value.lo, SCALE), LOWER)
        self.assertAlmostEqual(value.floats()[0], error_from_sigma(numeric), places=12)

    def test_completion_inequality_on_unrestricted_complex_vectors(self):
        # This is the proof's state-independent matrix completion, tested
        # independently of encoding and of the reduced Gram implementation.
        rng = np.random.default_rng(159)
        weights = np.array(grouped_weights(BIT_ERROR), dtype=float)
        a = np.array(coefficients(UPPER_INTEGERS), dtype=float)
        d = 1-sum(a*a/weights)
        outputs = rng.normal(size=(7, 18)) + 1j*rng.normal(size=(7, 18))
        target = rng.normal(size=7) + 1j*rng.normal(size=7)
        residual = target-outputs @ a
        completed = (outputs*weights) @ outputs.conj().T-np.outer(target, target.conj())
        completed += np.outer(residual, residual.conj())/d
        self.assertGreater(np.linalg.eigvalsh(completed).min(), -1e-13)

    def test_relaxed_optimum_bracket_and_illegal_lower_witness(self):
        sigma = relaxed_sigma()
        lower = rayleigh_lower(sigma, UPPER_INTEGERS)
        self.assertGreater(F(lower.lo, SCALE), RELAXED_LOWER)
        self.assertLess(UPPER-RELAXED_LOWER, F(3, 10**9))
        pt = np.array(sigma, dtype=float).reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
        self.assertLess(np.linalg.det(pt), -.001)

    def test_uniform_dual_upper_for_unseen_real_density_matrices(self):
        d, numerator = dual_data()
        h = np.array(numerator, dtype=float)/(2*np.mean(d.floats()))
        rng = np.random.default_rng(259)
        for rank in (1, 2, 4):
            for _ in range(4):
                c = rng.normal(size=(4, rank))
                sigma = c @ c.T/np.sum(c*c)
                self.assertLessEqual(error_from_sigma(sigma), np.trace(h @ sigma)+2e-14)
                self.assertLess(np.trace(h @ sigma), float(UPPER))

    def test_invalid_or_overstrong_certificates_are_rejected(self):
        with self.assertRaises(ValueError):
            coefficients([0.] * 18)
        with self.assertRaises(ValueError):
            dual_data((2*DENOMINATOR,)+(0,)*17)
        with self.assertRaises(ValueError):
            certify_upper(F(453, 1000))
        with self.assertRaises(ValueError):
            positive_ldl([[I.exact(1), I.exact(2)], [I.exact(2), I.exact(1)]])

    def test_existing_resource_ledger_and_strict_bracket(self):
        result = certificate()
        self.assertEqual(result['new_native_gates'], 368+3*4+33)
        self.assertEqual(result['gates_including_initial_join'], 579)
        self.assertLess(UPPER-LOWER, F(448, 10**6))
        self.assertFalse(result['original_source_promise_minimax_exactly_solved'])
        self.assertFalse(result['relaxed_witness_is_a_legal_independent_source'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(JointHistoryCertificateTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = certificate()
    report['automated_checks'] = {'run': checks.testsRun, 'failures': len(checks.failures),
                                 'errors': len(checks.errors)}
    if args.write_results:
        Path(__file__).with_name('joint_history_certificate_results.json').write_text(
            json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: report[k] for k in ('round', 'promised_worst_error_strict_lower_exact',
        'promised_worst_error_strict_upper_exact', 'certified_bracket_width_exact',
        'relaxed_worst_error_strict_lower_exact', 'new_native_gates', 'automated_checks')}, indent=2))


if __name__ == '__main__':
    main()
