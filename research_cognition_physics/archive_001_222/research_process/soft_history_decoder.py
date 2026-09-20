"""Round 162: continuously soften four conditional X/Z history decoders.

Only the report-dependent recovery changes. The original six true branches,
three noisy reads, invalid-code fallback and coherent history remain fixed.
The scalar parameter is tan(theta/2), giving exact rational orthogonal maps.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import unittest

import numpy as np

from centered_history_recovery import purification
from classical_joining_history import independent_encoding, record_unitary
from coarsened_history_bound import pure_extension_metrics
from flagged_subject_joining import joining_kraus, random_state
from joint_history_channel import I2, X, Y, Z, BIT_ERROR, source_sigma, bloch_state
from minimal_coherent_joining_history import conditional_unitary
from noisy_coherent_history_recovery import confusion_matrix, recovery_kraus as old_kraus


PARAMETER = F(3, 40)


def rational_circle(parameter=PARAMETER):
    if type(parameter) not in (int, F) or not 0 <= parameter <= 1:
        raise ValueError('Use an exact rational tan(theta/2) in [0,1].')
    t = F(parameter)
    return (1-t*t)/(1+t*t), 2*t/(1+t*t)


@lru_cache(None)
def soft_unitaries(parameter=PARAMETER):
    c, s = rational_circle(parameter)
    return (np.kron(Z, np.kron(c*Z+s*X, I2)),
            np.kron(Z, np.eye(4, dtype=int)),
            np.kron(Z, np.kron(c*X+s*Z, I2)),
            np.kron(I2, np.kron(I2, c*Z+s*X)),
            np.eye(8, dtype=int),
            np.kron(I2, np.kron(I2, c*X+s*Z)))


@lru_cache(None)
def branch_errors(parameter=PARAMETER):
    decoders = soft_unitaries(parameter)
    true = soft_unitaries(F(0))
    return tuple(decoders[g].T @ true[r] for r in range(6) for g in range(6))


def branch_weights(error=BIT_ERROR):
    matrix = confusion_matrix(error)
    return tuple(matrix[r][g]/6 for r in range(6) for g in range(6))


def decoder(report, parameter=PARAMETER):
    if type(report) is not int or not 0 <= report < 6:
        raise ValueError('Use a reported label in 0..5.')
    ks = joining_kraus((2, 2))
    w = np.array(soft_unitaries(parameter)[report], dtype=float)
    return np.vstack((ks[(0, 0)], w @ ks[(0, 1)])).T


def recovery_kraus(parameter=PARAMETER, error=BIT_ERROR, report=None):
    if report is not None and (type(report) is not int or not 0 <= report < 6):
        raise ValueError('Use a reported label in 0..5.')
    matrix = confusion_matrix(error)
    return tuple(np.sqrt(float(matrix[r][g])/6)*decoder(g, parameter) @ conditional_unitary(r+1)
                 for r in range(6) for g in (range(6) if report is None else (report,)))


@lru_cache(None)
def gram_kernel(parameter=PARAMETER):
    ws = np.array(branch_errors(parameter)+(np.eye(8),), dtype=float)
    products = np.einsum('aki,bkj->abij', ws, ws).reshape(37, 37, 2, 4, 2, 4)
    partial = np.trace(products, axis1=2, axis2=4)
    return np.eye(4)[None, None]/2+(partial+partial.transpose(0, 1, 3, 2))/8


@lru_cache(None)
def gram_basis(parameter=PARAMETER):
    scale = np.r_[np.sqrt(np.array(branch_weights(), dtype=float)), 1]
    kernel = gram_kernel(parameter)*scale[:, None, None, None]*scale[None, :, None, None]
    pauli = (I2, X, Y, Z)
    return np.array([[np.einsum('abij,ji->ab', kernel, np.kron(a, b.conj()).real/4)
                      for b in pauli] for a in pauli])


def gram_from_sigma(sigma, parameter=PARAMETER):
    scale = np.r_[np.sqrt(np.array(branch_weights(), dtype=float)), 1]
    return np.einsum('abij,ji->ab', gram_kernel(parameter), sigma)*scale[:, None]*scale[None, :]


def schur_root(gram):
    m, g = gram[:-1, :-1], gram[:-1, -1]
    t = 1.
    for _ in range(35):
        c = np.linalg.solve(m+t*np.eye(len(g)), g)
        d = 1-c @ c
        if d <= 0:
            raise ArithmeticError('The positive-root denominator failed.')
        step = (1-t-g @ c)/d
        t += step
        if abs(step) < 2e-15:
            break
    else:
        raise ArithmeticError('The positive Schur root did not converge.')
    c = np.linalg.solve(m+t*np.eye(len(g)), g)
    return t, c, 1-c @ c


def error_jet(point, parameter=PARAMETER):
    point = np.asarray(point, dtype=float)
    if point.shape != (6,):
        raise ValueError('Use six Bloch coordinates.')
    u, v = np.r_[1., point[:3]], np.r_[1., point[3:]]
    basis = gram_basis(parameter)
    gram = np.einsum('i,j,ijab->ab', u, v, basis)
    value, c, d = schur_root(gram)
    z = np.r_[-c, 1.]
    first = np.concatenate((np.einsum('j,ijab->iab', v, basis[1:]),
                            np.einsum('i,ijab->jab', u, basis[:, 1:])))
    return {'value': float(value), 'gradient': np.einsum('a,iab,b->i', z, first, z)/d,
            'coefficients': c*np.sqrt(np.array(branch_weights(), dtype=float))}


def supporting_coefficients(sigma, parameter=PARAMETER):
    value, c, d = schur_root(gram_from_sigma(sigma, parameter))
    return value, c*np.sqrt(np.array(branch_weights(), dtype=float))


def noiseless_error(parameter=PARAMETER):
    c, s = rational_circle(parameter)
    return (float(s*s)+np.sqrt(float(s**4+4*(1-c)**2)))/6


class SoftHistoryDecoderTests(unittest.TestCase):
    def test_exact_rational_reflections_and_original_endpoint(self):
        self.assertEqual(rational_circle(), (F(1591, 1609), F(240, 1609)))
        for parameter in (F(0), PARAMETER, F(1)):
            for w in soft_unitaries(parameter):
                np.testing.assert_array_equal(w.T @ w, np.eye(8, dtype=int))
        for r, w in enumerate(soft_unitaries(F(0))):
            np.testing.assert_array_equal(w, record_unitary(r+1))

    def test_channel_is_trace_preserving_and_zero_parameter_recovers_old_map(self):
        ks = recovery_kraus()
        np.testing.assert_allclose(sum(k.T @ k for k in ks), np.eye(16), atol=2e-15)
        np.testing.assert_allclose(recovery_kraus(F(0)), old_kraus(BIT_ERROR), atol=3e-16)

    def test_gram_computes_complete_external_error_for_unknown_sources(self):
        rng = np.random.default_rng(162)
        for _ in range(6):
            states = (random_state(rng, 2), random_state(rng, 2))
            gram = gram_from_sigma(source_sigma(*states))
            actual = pure_extension_metrics(recovery_kraus(), purification(independent_encoding(states)))['trace_error']
            self.assertAlmostEqual(schur_root(gram)[0], actual, places=13)

    def test_analytic_source_gradient_matches_finite_differences(self):
        point = np.array([.2, -.1, .3, -.2, .4, .1]); eps = 1e-5
        finite = [(error_jet(point+eps*v)['value']-error_jet(point-eps*v)['value'])/(2*eps) for v in np.eye(6)]
        np.testing.assert_allclose(error_jet(point)['gradient'], finite, atol=3e-10)

    def test_noiseless_limit_has_a_nonzero_uniform_price(self):
        rng = np.random.default_rng(262)
        for _ in range(5):
            omega = independent_encoding((random_state(rng, 2), random_state(rng, 2)))
            value = pure_extension_metrics(recovery_kraus(error=F(0)), purification(omega))['trace_error']
            self.assertAlmostEqual(value, noiseless_error(), places=13)
        self.assertGreater(noiseless_error(), .008)
        self.assertEqual(noiseless_error(F(0)), 0.)

    def test_report_probabilities_remain_state_independent(self):
        matrix = confusion_matrix(BIT_ERROR)
        for g in range(6):
            probability = sum(matrix[r][g] for r in range(6))/6
            ks = recovery_kraus(report=g)
            np.testing.assert_allclose(sum(k.T @ k for k in ks), probability*np.eye(16), atol=3e-16)

    def test_guards_and_simultaneous_conjugation_symmetry(self):
        for bad in (-1, F(3, 2), True, .075):
            with self.assertRaises(ValueError): rational_circle(bad)
        with self.assertRaises(ValueError): decoder(6)
        x = np.array([.2, .3, -.4, -.1, -.2, .5]); reflected = x.copy(); reflected[[1, 4]] *= -1
        self.assertAlmostEqual(error_jet(x)['value'], error_jet(reflected)['value'], places=14)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SoftHistoryDecoderTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    c, s = rational_circle()
    report = {'round': 162, 'tan_half_angle_exact': str(PARAMETER),
              'cos_theta_exact': str(c), 'sin_theta_exact': str(s),
              'angle_radians_diagnostic': float(2*np.arctan(float(PARAMETER))),
              'modified_reports_zero_based': [0, 2, 3, 5],
              'readout_and_invalid_code_fallback_unchanged': True,
              'noiseless_error_formula': '(s^2+sqrt(s^4+4(1-c)^2))/6',
              'noiseless_error_diagnostic': noiseless_error(),
              'parameter_or_decoder_family_globally_optimal_claimed': False,
              'automated_checks': {'run': checks.testsRun, 'failures': len(checks.failures), 'errors': len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name('soft_history_decoder_results.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
