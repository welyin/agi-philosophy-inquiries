"""Round 168: a coherent relation correction outside the four-angle family.

One report-independent orthogonal correction acts on the retained negative
orientation sector, after the four-angle conditional inverse and before
the original intake inverse. All six reports, including Y, are affected.
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
from flagged_subject_joining import joining_kraus, random_state
from four_angle_history_decoder import (branch_errors as previous_errors,
    unitaries as previous_unitaries, recovery_kraus as previous_kraus)
from joint_history_channel import I2, X, Y, Z, BIT_ERROR, source_sigma
from minimal_coherent_joining_history import conditional_unitary
from noisy_coherent_history_recovery import confusion_matrix
from soft_history_decoder import rational_circle, branch_weights, schur_root


PARAMETER = F(3, 800)
J = np.array([[0, -1], [1, 0]])
G = np.kron(Z, np.kron(X, J))
H = np.kron(Z, np.kron(J, Z))


@lru_cache(None)
def correction(parameter=PARAMETER):
    c, s = rational_circle(parameter)
    return (c*np.eye(8, dtype=int)-s*G) @ (c*np.eye(8, dtype=int)+s*H)


@lru_cache(None)
def branch_errors(parameter=PARAMETER):
    e = correction(parameter)
    return tuple(e @ w for w in previous_errors())


def decoder(report, parameter=PARAMETER):
    if type(report) is not int or not 0 <= report < 6:
        raise ValueError('Use a report in 0..5.')
    negative = np.array(correction(parameter) @ previous_unitaries()[report].T, dtype=float)
    ks = joining_kraus((2, 2))
    return np.vstack((ks[(0, 0)], negative.T @ ks[(0, 1)])).T


def recovery_kraus(parameter=PARAMETER, error=BIT_ERROR):
    matrix = confusion_matrix(error)
    return tuple(np.sqrt(float(matrix[r][g])/6)*decoder(g, parameter) @ conditional_unitary(r+1)
                 for r in range(6) for g in range(6))


@lru_cache(None)
def gram_kernel(parameter=PARAMETER):
    ws = np.array(branch_errors(parameter)+(np.eye(8),), dtype=float)
    products = np.einsum('aki,bkj->abij', ws, ws).reshape(37, 37, 2, 4, 2, 4)
    partial = np.trace(products, axis1=2, axis2=4)
    return np.eye(4)[None, None]/2+(partial+partial.transpose(0, 1, 3, 2))/8


@lru_cache(None)
def gram_basis(parameter=PARAMETER):
    scale = np.r_[np.sqrt(np.array(branch_weights(), dtype=float)), 1.]
    kernel = gram_kernel(parameter)*scale[:, None, None, None]*scale[None, :, None, None]
    return np.array([[np.einsum('abij,ji->ab', kernel, np.kron(a, b.conj()).real/4)
                      for b in (I2, X, Y, Z)] for a in (I2, X, Y, Z)])


def gram_from_sigma(sigma, parameter=PARAMETER):
    scale = np.r_[np.sqrt(np.array(branch_weights(), dtype=float)), 1.]
    return np.einsum('abij,ji->ab', gram_kernel(parameter), sigma)*scale[:, None]*scale[None, :]


def error_jet(point, parameter=PARAMETER):
    point = np.asarray(point, dtype=float)
    if point.shape != (6,):
        raise ValueError('Use two three-coordinate Bloch vectors.')
    u, v = np.r_[1., point[:3]], np.r_[1., point[3:]]
    basis = gram_basis(parameter)
    gram = np.einsum('i,j,ijab->ab', u, v, basis)
    value, c, d = schur_root(gram); z = np.r_[-c, 1.]
    first = np.concatenate((np.einsum('j,ijab->iab', v, basis[1:]),
                            np.einsum('i,ijab->jab', u, basis[:, 1:])))
    return {'value':float(value), 'gradient':np.einsum('a,iab,b->i', z, first, z)/d,
            'coefficients':c*np.sqrt(np.array(branch_weights(), dtype=float))}


def correction_angle_derivative(sigma, parameter=PARAMETER):
    from soft_decoder_certificate import partial_reference
    value, c, d = schur_root(gram_from_sigma(sigma, parameter))
    a = c*np.sqrt(np.array(branch_weights(), dtype=float))
    ws = np.array(branch_errors(parameter), dtype=float)
    residual = np.eye(8)-np.einsum('j,jab->ab', a, ws)
    derivative = -np.einsum('j,jab->ab', a, np.einsum('ij,ajk->aik', H-G, ws))
    return float(np.trace(sigma @ partial_reference(residual.T @ derivative+derivative.T @ residual))/(4*d))


class RelationalHistoryDecoderTests(unittest.TestCase):
    def test_exact_generators_commute_and_correction_is_orthogonal(self):
        for generator in (G,H):
            np.testing.assert_array_equal(generator.T, -generator)
            np.testing.assert_array_equal(generator @ generator, -np.eye(8,dtype=int))
        np.testing.assert_array_equal(G @ H, H @ G)
        np.testing.assert_array_equal(G @ H, np.kron(I2,np.kron(Z,X)))
        for t in (F(0), PARAMETER, F(1,10)):
            e=correction(t)
            np.testing.assert_array_equal(e.T @ e, np.eye(8,dtype=int))

    def test_zero_correction_recovers_previous_channel(self):
        np.testing.assert_allclose(recovery_kraus(F(0)), previous_kraus(), atol=3e-16)
        ks=recovery_kraus()
        np.testing.assert_allclose(sum(k.T @ k for k in ks),np.eye(16),atol=3e-15)

    def test_gram_matches_complete_external_purifications(self):
        rng=np.random.default_rng(168)
        for _ in range(6):
            states=(random_state(rng,2),random_state(rng,2))
            actual=pure_extension_metrics(recovery_kraus(),purification(independent_encoding(states)))['trace_error']
            self.assertAlmostEqual(actual,schur_root(gram_from_sigma(source_sigma(*states)))[0],places=13)

    def test_source_gradient_matches_finite_differences(self):
        p=np.array([.2,.1,-.4,-.1,.3,.2]);eps=1e-5
        finite=[(error_jet(p+eps*v)['value']-error_jet(p-eps*v)['value'])/(2*eps) for v in np.eye(6)]
        np.testing.assert_allclose(error_jet(p)['gradient'],finite,atol=3e-10)

    def test_correction_angle_derivative_and_departure_direction(self):
        from four_angle_global_certificate import witness_data
        _,_,sigma=witness_data();sigma=np.array(sigma,dtype=float)
        self.assertLess(correction_angle_derivative(sigma,F(0)), -.06)
        step=F(1,10**6)
        t=PARAMETER
        delta=2*np.arctan(float(t+step))-2*np.arctan(float(t-step))
        finite=(schur_root(gram_from_sigma(sigma,t+step))[0]-schur_root(gram_from_sigma(sigma,t-step))[0])/delta
        self.assertAlmostEqual(correction_angle_derivative(sigma),finite,places=8)

    def test_all_reports_including_y_change_without_reading_coherent_history(self):
        for g in range(6):
            self.assertGreater(np.linalg.norm(decoder(g)-decoder(g,F(0))), .001)
        e=correction();full=np.zeros((16,16),dtype=object)
        full[:8,:8]=np.eye(8,dtype=int);full[8:,8:]=e
        np.testing.assert_array_equal(full.T @ full,np.eye(16,dtype=int))
        # Reference R remains block diagonal, so the sigma reduction survives.
        np.testing.assert_array_equal(e[:4,4:],np.zeros((4,4),dtype=int))

    def test_conjugation_symmetry_and_guards(self):
        p=np.array([.1,.3,-.4,.2,-.2,.4]);other=p.copy();other[[1,4]]*=-1
        self.assertAlmostEqual(error_jet(p)['value'],error_jet(other)['value'],places=14)
        with self.assertRaises(ValueError): decoder(True)
        with self.assertRaises(ValueError): correction(.01)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RelationalHistoryDecoderTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    from four_angle_global_certificate import witness_data
    _,_,sigma=witness_data()
    report={'round':168,'half_angle_tangent_exact':str(PARAMETER),
            'correction_angle_radians_diagnostic':float(2*np.arctan(float(PARAMETER))),
            'generators':['-Z_R X_A J_B','Z_R J_A Z_B'],
            'J_definition':[[0,-1],[1,0]], 'all_six_reports_modified':True,
            'old_witness_derivative_at_zero_diagnostic':correction_angle_derivative(np.array(sigma,float),F(0)),
            'finite_noise_readout_and_all_reports_unchanged':True,
            'coherent_history_rebits':1,'global_improvement_proved_here':False,
            'automated_checks':{'run':checks.testsRun,'failures':len(checks.failures),'errors':len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name('relational_history_decoder_results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__': main()
