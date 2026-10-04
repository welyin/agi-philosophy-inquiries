"""Round 423: sharp finite source dimension for exact Bernoulli records.

The analytic lower bound uses positivity and successive endpoint limits.
These checks verify the finite instruments, constructions, and scope examples.
--check writes nothing; default exclusively creates the result JSON.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)
OBS = {}


def norm(a):
    return float(np.linalg.norm(a, 2))


def binomial(n, p):
    return np.array([math.comb(n, k) * p**k * (1-p)**(n-k) for k in range(n+1)])


def words(n):
    return list(itertools.product((0, 1), repeat=n))


def instrument(n, remaining):
    # The clipping defines a CPTP instrument even on unreachable k > remaining.
    q = np.minimum(np.arange(n+1) / remaining, 1.)
    a0 = np.diag(np.sqrt(1-q)).astype(complex)
    a1 = np.zeros((n+1, n+1), complex)
    for k in range(1, n+1):
        a1[k-1, k] = math.sqrt(q[k])
    return a0, a1


def sequential_branches(n):
    current = {(): np.eye(n+1, dtype=complex)}
    for remaining in range(n, 0, -1):
        a = instrument(n, remaining)
        current = {prefix+(b,): a[b] @ old
                   for prefix, old in current.items() for b in (0, 1)}
    return current


def coherent_decoder(n):
    branch = sequential_branches(n)
    return np.stack([branch[w] for w in words(n)], axis=1).reshape((n+1)*2**n, n+1)


def dicke_decoder(n):
    out = np.zeros(((n+1)*2**n, n+1), complex)
    for j, w in enumerate(words(n)):
        k = sum(w)
        out[j, k] = 1 / math.sqrt(math.comb(n, k))
    return out


def generalized_bernstein(c, p):
    m = len(c)-1
    values = [Fraction(ck)*p**k*(1-p)**(m-k) for k, ck in enumerate(c)]
    total = sum(values)
    return [v/total for v in values]


class Audit(unittest.TestCase):
    def test_01_literal_sequential_instruments_are_complete(self):
        worst = 0.
        for n in range(1, 8):
            for remaining in range(1, n+1):
                a0, a1 = instrument(n, remaining)
                worst = max(worst, norm(a0.conj().T @ a0 + a1.conj().T @ a1 - np.eye(n+1)))
                pointer_isometry = np.stack([a0, a1], axis=1).reshape(2*(n+1), n+1)
                self.assertLess(norm(pointer_isometry.conj().T @ pointer_isometry - np.eye(n+1)), 1e-12)
        self.assertLess(worst, 1e-12)
        OBS['instrument_completeness_error'] = worst

    def test_02_complete_record_distribution_from_count_source(self):
        maximum = 0.
        rows = []
        for n in (1, 2, 3, 5, 6):
            branches = sequential_branches(n)
            for p in (0., .001, .17, .5, .83, .999, 1.):
                rho = np.diag(binomial(n, p))
                actual = np.array([np.trace(k @ rho @ k.conj().T).real for k in branches.values()])
                expected = np.array([p**sum(w) * (1-p)**(n-sum(w)) for w in branches])
                maximum = max(maximum, float(np.max(np.abs(actual-expected))))
                self.assertLess(abs(float(actual.sum())-1), 1e-12)
            rows.append(dict(reads=n, source_dimension=n+1, output_words=2**n))
        self.assertLess(maximum, 1e-12)
        OBS['whole_record_probability_error'] = maximum
        OBS['exact_achieving_resources'] = rows

    def test_03_independent_rational_sampling_without_replacement(self):
        checked = 0
        for n in (2, 3, 5, 7):
            for p in (Fraction(1, 7), Fraction(2, 5), Fraction(5, 6)):
                for w in words(n):
                    k = sum(w)
                    probability = Fraction(math.comb(n, k))*p**k*(1-p)**(n-k)
                    remaining_ones = k
                    for i, b in enumerate(w):
                        r = n-i
                        probability *= Fraction(remaining_ones if b else r-remaining_ones, r)
                        remaining_ones -= b
                    self.assertEqual(probability, p**k*(1-p)**(n-k))
                    checked += 1
        OBS['exact_rational_full_words_checked'] = checked

    def test_04_full_unknown_input_and_reference_isometry(self):
        rng = np.random.default_rng(423)
        worst = 0.
        reference_error = 0.
        for n in (2, 3, 6):
            v, target = coherent_decoder(n), dicke_decoder(n)
            worst = max(worst, norm(v-target), norm(v.conj().T @ v-np.eye(n+1)))
            vector = rng.normal(size=(n+1, 3)) + 1j*rng.normal(size=(n+1, 3))
            vector /= np.linalg.norm(vector)
            evolved = v @ vector
            reference_error = max(reference_error, norm(evolved.conj().T @ evolved-vector.conj().T @ vector))
            self.assertLess(np.linalg.norm(evolved.reshape(n+1, 2**n, 3)[1:]), 1e-12)
            # A known pure source with binomial amplitudes decodes to N product qubits.
            p = .31
            expected = np.array([math.sqrt(p**sum(w)*(1-p)**(n-sum(w))) for w in words(n)])
            self.assertLess(np.linalg.norm((v @ np.sqrt(binomial(n, p)))[:2**n]-expected), 1e-12)
        self.assertLess(worst, 1e-12)
        self.assertLess(reference_error, 1e-12)
        OBS['dicke_isometry_and_decoder_error'] = worst
        OBS['arbitrary_input_reference_error'] = reference_error

    def test_05_endpoint_filter_induction_with_exact_coefficients(self):
        count = 0
        for coefficients in ([1, 4, 6, 4, 1], [2, 7, 3, 11, 5, 13, 17], [1, 8, 2]):
            c = list(coefficients)
            while len(c) > 1:
                for p in (Fraction(1, 100), Fraction(2, 7), Fraction(4, 5)):
                    distribution = generalized_bernstein(c, p)
                    filtered = [v/(1-distribution[0]) for v in distribution[1:]]
                    self.assertEqual(filtered, generalized_bernstein(c[1:], p))
                    count += 1
                c = c[1:]
        # Nonprojective filter: validate the inverse on supp B, not on its kernel.
        effects = [np.diag([1., .2, .4, .3]), np.diag([0., .5, .1, .2]),
                   np.diag([0., .3, .5, .5])]
        b = effects[1]+effects[2]
        sqrt_b = np.diag(np.sqrt(np.diag(b)))
        inverse = np.diag(1/np.sqrt(np.diag(b)[1:]))
        filtered_effects = [inverse @ e[1:, 1:] @ inverse for e in effects[1:]]
        vector = np.array([1, 2j, -1+.4j, .7], complex)
        rho = np.outer(vector, vector.conj()) / np.vdot(vector, vector).real
        success = np.trace(b @ rho).real
        conditional = (sqrt_b @ rho @ sqrt_b)[1:, 1:]/success
        errors = [abs(np.trace(conditional @ f)-np.trace(rho @ e)/success)
                  for f, e in zip(filtered_effects, effects[1:])]
        self.assertLess(norm(sum(filtered_effects)-np.eye(3)), 1e-12)
        self.assertLess(max(errors), 1e-12)
        OBS['exact_endpoint_filter_identities_checked'] = count
        OBS['nonprojective_filter_error'] = float(max(errors))

    def test_06_restricted_bias_qubit_counterexample(self):
        vectors = np.array([[math.cos(2*math.pi*k/3), math.sin(2*math.pi*k/3), 0.] for k in range(3)])
        matrices = [v[0]*X+v[1]*Y+v[2]*Z for v in vectors]
        effects = [(I2+m)/3 for m in matrices]
        worst = 0.
        largest_norm_squared = 0.
        for p in np.linspace(.25, .75, 101):
            q = binomial(2, float(p))
            r = 2*q @ vectors
            rho = (I2+r[0]*X+r[1]*Y+r[2]*Z)/2
            actual = np.array([np.trace(rho @ e).real for e in effects])
            worst = max(worst, float(np.max(np.abs(actual-q))))
            squared = float(r @ r)
            z = float(p)-.5
            self.assertAlmostEqual(squared, .25+6*z*z+36*z**4)
            self.assertGreater(float(np.linalg.eigvalsh(rho).min()), 0.)
            largest_norm_squared = max(largest_norm_squared, squared)
            # count 1 is randomly ordered; the full four-record law is iid.
            full = actual[[0, 1, 1, 2]] * np.array([1, .5, .5, 1])
            self.assertTrue(np.allclose(full, [(1-p)**2, p*(1-p), p*(1-p), p*p]))
        self.assertLess(worst, 1e-12)
        self.assertAlmostEqual(largest_norm_squared, 49/64)
        at_zero = (I2+2*matrices[0])/2
        self.assertAlmostEqual(float(np.linalg.eigvalsh(at_zero).min()), -.5)
        OBS['restricted_bias_probability_error'] = worst
        OBS['restricted_bias_max_bloch_norm_squared'] = largest_norm_squared
        OBS['same_encoding_at_zero_minimum_eigenvalue'] = float(np.linalg.eigvalsh(at_zero).min())

    def test_07_achieving_preparation_is_not_a_free_single_coin_channel(self):
        # An affine preparation channel acting on diag(1-p,p) must preserve this midpoint.
        rows = []
        for n in (2, 3, 6):
            middle = binomial(n, .5)
            endpoint_mixture = (binomial(n, 0.)+binomial(n, 1.))/2
            defect = float(np.sum(np.abs(middle-endpoint_mixture))/2)
            self.assertAlmostEqual(defect, 1-2**(1-n))
            rows.append(dict(reads=n, preparation_affinity_defect=defect))
        OBS['source_preparation_must_be_supplied'] = rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise SystemExit(1)
    saved = dict(round=423, baseline_round=422, status='verified_conditional_result',
                 tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
                 python=platform.python_version(), numpy=np.__version__, observations=OBS,
                 scope=dict(exact_all_bias_iid_records_are_extra_task=True,
                            arbitrary_consumption_of_source_allowed=True,
                            encoding_continuity_required=False,
                            full_parameter_dependent_resources_counted=True,
                            sharp_source_dimension_for_N_reads='N+1',
                            blank_record_storage_counted_separately=True,
                            noisy_optimal_dimension_proved=False,
                            restricted_bias_same_bound_claimed=False,
                            source_preparation_autonomously_generated=False,
                            original_round413_algorithm_invalidated=False,
                            three_dimensional_space_derived=False,
                            full_cognitive_countermodel_completed=False,
                            full_GR_goal_completed=False, phase_closure_triggered=False))
    if not args.check:
        with Path(__file__).with_name('finite_sampling_source_audit_results.json').open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(saved, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print(json.dumps(saved, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
