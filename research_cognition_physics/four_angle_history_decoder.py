"""Round 165: independent rational angles for four reported X/Z decoders.

Readout, six true branches, all raw outcomes and coherent capacity are fixed.
Floating derivatives discover candidates; round 166 verifies a global cover.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import unittest

import numpy as np

from centered_history_recovery import purification
from classical_joining_history import independent_encoding
from coarsened_history_bound import pure_extension_metrics
from compiled_soft_history_recovery import soft_decoder_circuit, report_probabilities
from compiled_subject_joining import circuit_matrix
from compiled_three_subject_joining import gate_counts
from flagged_subject_joining import joining_kraus, random_state
from joint_history_channel import I2, X, Y, Z, BIT_ERROR, source_sigma
from minimal_coherent_joining_history import conditional_unitary
from noisy_coherent_history_recovery import confusion_matrix
from soft_history_decoder import (PARAMETER as COMMON_PARAMETER, rational_circle,
    soft_unitaries, branch_weights, schur_root)


MODIFIED = (0, 2, 3, 5)
PARAMETERS = tuple(F(n, 10**6) for n in (87728, 84862, 84638, 25145))


def validate_parameters(parameters):
    if not isinstance(parameters, tuple) or len(parameters) != 4:
        raise ValueError('Use four exact half-angle tangents, in a tuple.')
    for value in parameters:
        rational_circle(value)
    return parameters


@lru_cache(None)
def unitaries(parameters=PARAMETERS):
    validate_parameters(parameters)
    cs = [rational_circle(t) for t in parameters]
    (a, b), (c, d), (e, f), (g, h) = cs
    return (np.kron(Z, np.kron(a*Z+b*X, I2)), np.kron(Z, np.eye(4, dtype=int)),
            np.kron(Z, np.kron(c*X+d*Z, I2)), np.kron(I2, np.kron(I2, e*Z+f*X)),
            np.eye(8, dtype=int), np.kron(I2, np.kron(I2, g*X+h*Z)))


@lru_cache(None)
def branch_errors(parameters=PARAMETERS):
    repair, true = unitaries(parameters), soft_unitaries(F(0))
    return tuple(repair[g].T @ true[r] for r in range(6) for g in range(6))


def decoder(report, parameters=PARAMETERS):
    if type(report) is not int or not 0 <= report < 6:
        raise ValueError('Use a report in 0..5.')
    ks = joining_kraus((2, 2))
    return np.vstack((ks[(0, 0)], np.array(unitaries(parameters)[report], dtype=float) @ ks[(0, 1)])).T


def recovery_kraus(parameters=PARAMETERS, error=BIT_ERROR):
    matrix = confusion_matrix(error)
    return tuple(np.sqrt(float(matrix[r][g])/6)*decoder(g, parameters) @ conditional_unitary(r+1)
                 for r in range(6) for g in range(6))


@lru_cache(None)
def gram_kernel(parameters=PARAMETERS):
    ws = np.array(branch_errors(parameters)+(np.eye(8),), dtype=float)
    products = np.einsum('aki,bkj->abij', ws, ws).reshape(37, 37, 2, 4, 2, 4)
    partial = np.trace(products, axis1=2, axis2=4)
    return np.eye(4)[None, None]/2+(partial+partial.transpose(0, 1, 3, 2))/8


@lru_cache(None)
def gram_basis(parameters=PARAMETERS):
    scale = np.r_[np.sqrt(np.array(branch_weights(), dtype=float)), 1]
    kernel = gram_kernel(parameters)*scale[:, None, None, None]*scale[None, :, None, None]
    return np.array([[np.einsum('abij,ji->ab', kernel, np.kron(a, b.conj()).real/4)
                      for b in (I2, X, Y, Z)] for a in (I2, X, Y, Z)])


def gram_from_sigma(sigma, parameters=PARAMETERS):
    scale = np.r_[np.sqrt(np.array(branch_weights(), dtype=float)), 1]
    return np.einsum('abij,ji->ab', gram_kernel(parameters), sigma)*scale[:, None]*scale[None, :]


def error_jet(point, parameters=PARAMETERS):
    point = np.asarray(point, dtype=float)
    if point.shape != (6,):
        raise ValueError('Use two three-coordinate Bloch vectors.')
    u, v = np.r_[1., point[:3]], np.r_[1., point[3:]]
    basis = gram_basis(parameters)
    gram = np.einsum('i,j,ijab->ab', u, v, basis)
    value, c, d = schur_root(gram)
    z = np.r_[-c, 1.]
    first = np.concatenate((np.einsum('j,ijab->iab', v, basis[1:]),
                            np.einsum('i,ijab->jab', u, basis[:, 1:])))
    return {'value':float(value), 'gradient':np.einsum('a,iab,b->i', z, first, z)/d,
            'coefficients':c*np.sqrt(np.array(branch_weights(), dtype=float))}


def decoder_circuit(report, parameters=PARAMETERS, angle_errors=(0., 0.)):
    validate_parameters(parameters)
    if type(report) is not int or not 0 <= report < 6:
        raise ValueError('Use a report in 0..5.')
    t = parameters[MODIFIED.index(report)] if report in MODIFIED else F(0)
    return soft_decoder_circuit(report, t, angle_errors)


class FourAngleHistoryDecoderTests(unittest.TestCase):
    def test_exact_orthogonality_and_common_angle_embedding(self):
        for parameters in (PARAMETERS, (F(0),)*4, (F(1), F(0), F(1, 2), F(1, 4))):
            for w in unitaries(parameters):
                np.testing.assert_array_equal(w.T @ w, np.eye(8, dtype=int))
        np.testing.assert_array_equal(unitaries((COMMON_PARAMETER,)*4), soft_unitaries())

    def test_complete_external_distance_matches_gram(self):
        rng = np.random.default_rng(165)
        for _ in range(5):
            states = (random_state(rng, 2), random_state(rng, 2))
            actual = pure_extension_metrics(recovery_kraus(), purification(independent_encoding(states)))['trace_error']
            self.assertAlmostEqual(schur_root(gram_from_sigma(source_sigma(*states)))[0], actual, places=13)
        ks = recovery_kraus()
        np.testing.assert_allclose(sum(k.T @ k for k in ks), np.eye(16), atol=3e-15)

    def test_source_gradient(self):
        x = np.array([.1, .2, -.3, -.2, .1, .4]); eps = 1e-5
        finite = [(error_jet(x+eps*v)['value']-error_jet(x-eps*v)['value'])/(2*eps) for v in np.eye(6)]
        np.testing.assert_allclose(error_jet(x)['gradient'], finite, atol=3e-10)

    def test_analytic_angle_gradient_and_search_kernel(self):
        from search_four_angle_history import kernel, angle_gradient
        angles = 2*np.arctan(np.array(PARAMETERS, dtype=float))
        sigma = source_sigma(random_state(np.random.default_rng(166), 2), np.eye(2)/2)
        weighted, ws = kernel(angles)
        gram = np.einsum('abij,ji->ab', weighted, sigma)
        np.testing.assert_allclose(gram, gram_from_sigma(sigma), atol=3e-16)
        _, c, d = schur_root(gram)
        analytic = angle_gradient(angles, ws, sigma, c, d)
        def value(theta):
            return schur_root(np.einsum('abij,ji->ab', kernel(theta)[0], sigma))[0]
        eps = 1e-5
        finite = [(value(angles+eps*v)-value(angles-eps*v))/(2*eps) for v in np.eye(4)]
        np.testing.assert_allclose(analytic, finite, atol=3e-10)

    def test_all_six_compiled_maps_and_gate_counts(self):
        counts = []
        for g in range(6):
            circuit = decoder_circuit(g)
            np.testing.assert_allclose(circuit_matrix(circuit, 4), decoder(g), atol=4e-15)
            counts.append(gate_counts(circuit)['all'])
        self.assertEqual(counts, [33, 21, 33, 23, 11, 23])
        self.assertEqual(368+12+max(counts), 413)

    def test_report_weights_and_angle_tolerance_are_unchanged(self):
        eps = .001
        for g in range(6):
            actual = circuit_matrix(decoder_circuit(g, angle_errors=(eps, -eps)), 4)
            self.assertLessEqual(np.linalg.norm(actual-decoder(g), 2), eps+5e-15)
        p = report_probabilities(BIT_ERROR)
        self.assertAlmostEqual(sum(p[g] for g in MODIFIED), (4+BIT_ERROR)/6)

    def test_conjugation_symmetry_and_parameter_guards(self):
        point = np.array([.2, -.3, .4, -.2, .1, .3]); other = point.copy(); other[[1, 4]] *= -1
        self.assertAlmostEqual(error_jet(point)['value'], error_jet(other)['value'], places=14)
        for bad in ((F(0),)*3, (F(0),)*3+(.1,), (F(0),)*3+(True,), (F(0),)*3+(F(-1),)):
            with self.assertRaises(ValueError): unitaries(bad)
        with self.assertRaises(ValueError): decoder_circuit(True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(FourAngleHistoryDecoderTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {'round':165, 'half_angle_tangents_exact':[str(t) for t in PARAMETERS],
              'angles_radians_diagnostic':(2*np.arctan(np.array(PARAMETERS, dtype=float))).tolist(),
              'report_order':['A_X', 'A_Z', 'B_X', 'B_Z'],
              'decoder_native_gate_counts':[33, 21, 33, 23, 11, 23],
              'new_native_gates_worst':413, 'noisy_reads':3, 'coherent_history_rebits':1,
              'expected_gate_count_same_as_round_164':True,
              'all_eight_raw_reports_retained':True, 'parameter_optimality_proved_here':False,
              'automated_checks':{'run':checks.testsRun,'failures':len(checks.failures),'errors':len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name('four_angle_history_decoder_results.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
