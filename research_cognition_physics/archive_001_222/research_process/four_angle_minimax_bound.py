"""Round 167: a uniform lower bound for every four-angle recovery.

A single legal source and fixed external Rayleigh vector work for all angles.
Four two-dimensional quadratic upper certificates cover entire circles,
including angle choices outside the implementation's quarter-circle range.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import unittest

import numpy as np

from certified_intervals import Interval as I, SCALE
from four_angle_global_certificate import UPPER, residual_blocks, witness_data
from four_angle_history_decoder import MODIFIED, PARAMETERS, unitaries
from joint_history_certificate import positive_ldl
from joint_relation_records import interval_fields
from noisy_classical_history import majority_error
from noisy_coherent_history_recovery import confusion_matrix
from soft_decoder_certificate import block_expectation
from soft_history_decoder import soft_unitaries, schur_root


LOWER = F(4493947394, 10**10)
# Positive Lagrange multipliers are selected numerically, then fixed rationally.
# Any multipliers passing the interval positivity check give valid bounds.
MULTIPLIERS = (F(84323,10**7), F(112057,10**7), F(111976,10**7), F(83451,10**7))


@lru_cache(None)
def witness_components():
    data, _, sigma = witness_data()
    b, residual = residual_blocks(tuple(data['coefficient_numerators']))
    norm = block_expectation(sigma, b*b, residual.T @ residual)
    target = block_expectation(sigma, b, residual.T)
    if norm <= 0:
        raise ValueError('The fixed external test vector has zero norm.')
    probabilities = confusion_matrix(majority_error(1))
    true = soft_unitaries(F(0))
    zero, ninety = unitaries((F(0),)*4), unitaries((F(1),)*4)
    fixed = I.exact(0); reports = []
    for g in range(6):
        p = tuple(probabilities[r][g]/6 for r in range(6))
        if g not in MODIFIED:
            overlaps = [block_expectation(sigma, b, residual.T @ true[g].T @ true[r]) for r in range(6)]
            fixed += sum(weight*x*x for weight, x in zip(p, overlaps))
            continue
        # For v=(cos(theta),sin(theta)), each overlap is alpha + beta_r.v.
        alpha = b/2
        beta = tuple(tuple(block_expectation(sigma, F(0), residual.T @ w[g].T @ true[r])
                           for w in (zero, ninety)) for r in range(6))
        k = sum(p)*alpha*alpha
        linear = tuple(sum(p[r]*alpha*beta[r][i] for r in range(6)) for i in range(2))
        q = tuple(tuple(sum(p[r]*beta[r][i]*beta[r][j] for r in range(6)) for j in range(2)) for i in range(2))
        reports.append({'report':g, 'constant':k, 'linear':linear, 'quadratic':q,
                        'overlap_constant_exact':alpha, 'overlap_linear_exact':beta})
    return {'norm':norm, 'target':target, 'fixed':fixed, 'reports':reports}


def circle_upper(constant, linear, quadratic, multiplier):
    """f(v)<=k+lambda+l^T(lambda I-Q)^-1 l for every ||v||=1."""
    if type(multiplier) not in (int, F):
        raise ValueError('Use an exact rational multiplier.')
    a = [[I.exact(multiplier)*int(i == j)-quadratic[i][j] for j in range(2)] for i in range(2)]
    pivots = positive_ldl(a)
    determinant = a[0][0]*a[1][1]-a[0][1]*a[1][0]
    if determinant.lo <= 0:
        raise ValueError('The inverse denominator is not strictly positive.')
    l0, l1 = linear
    correction = (a[1][1]*l0*l0-(a[0][1]+a[1][0])*l0*l1+a[0][0]*l1*l1)/determinant
    return constant+multiplier+correction, pivots


@lru_cache(None)
def family_certificate():
    parts = witness_components(); upper_sum = parts['fixed']; details = []
    for item, multiplier in zip(parts['reports'], MULTIPLIERS):
        upper, pivots = circle_upper(item['constant'], item['linear'], item['quadratic'], multiplier)
        upper_sum += upper
        details.append({'report_zero_based':item['report'], 'multiplier_exact':str(multiplier),
                        'overlap_square_sum_upper':interval_fields(upper),
                        'positive_ldl_pivots':[interval_fields(p) for p in pivots],
                        'quadratic':[[interval_fields(x) for x in row] for row in item['quadratic']],
                        'linear':[interval_fields(x) for x in item['linear']],
                        'constant':interval_fields(item['constant'])})
    lower = (parts['target']**2-upper_sum)/parts['norm']
    if F(lower.lo, SCALE) <= LOWER:
        raise ValueError('The uniform family lower endpoint failed.')
    return {'round':167, 'strict_family_minimax_lower_exact':str(LOWER),
            'fixed_source_uniform_lower':interval_fields(lower),
            'available_candidate_upper_exact':str(UPPER),
            'certified_candidate_suboptimality_less_than_exact':str(UPPER-LOWER),
            'fixed_external_vector_norm_squared_exact':str(parts['norm']),
            'fixed_external_target_overlap_exact':str(parts['target']),
            'unmodified_Y_overlap_square_sum':interval_fields(parts['fixed']),
            'report_circle_certificates':details,
            'same_legal_source_and_external_test_vector_for_all_angles':True,
            'covers_all_four_angles_on_entire_circles':True,
            'covers_classically_randomized_angle_choices':True,
            'readout_confusion_and_two_Y_decoders_fixed':True,
            'all_CPTP_decoders_or_readout_designs_covered':False,
            'exact_optimizer_or_uniqueness_proved':False,
            'certificate_requires_scipy':False}


class FourAngleMinimaxBoundTests(unittest.TestCase):
    def test_uniform_lower_closes_the_family_minimax_gap(self):
        result = family_certificate()
        self.assertLess(F(result['certified_candidate_suboptimality_less_than_exact']), F(11,10**9))
        self.assertTrue(result['covers_all_four_angles_on_entire_circles'])
        self.assertFalse(result['all_CPTP_decoders_or_readout_designs_covered'])

    def test_full_circle_quadratics_obey_the_certificate(self):
        for item, multiplier in zip(witness_components()['reports'], MULTIPLIERS):
            upper, _ = circle_upper(item['constant'], item['linear'], item['quadratic'], multiplier)
            k = np.mean(item['constant'].floats())
            l = np.array([np.mean(x.floats()) for x in item['linear']])
            q = np.array([[np.mean(x.floats()) for x in row] for row in item['quadratic']])
            for theta in np.linspace(-np.pi, np.pi, 41):
                v = np.array([np.cos(theta),np.sin(theta)])
                self.assertLessEqual(k+2*l @ v+v @ q @ v, upper.floats()[1]+1e-16)

    def test_fixed_witness_lower_holds_for_unseen_recoveries(self):
        from search_four_angle_history import kernel
        _, _, exact_sigma = witness_data(); sigma = np.array(exact_sigma,dtype=float)
        rng = np.random.default_rng(167)
        for angles in [np.zeros(4), np.full(4,np.pi/2)]+[rng.uniform(-np.pi,np.pi,4) for _ in range(10)]:
            gram = np.einsum('abij,ji->ab',kernel(angles)[0],sigma)
            self.assertGreater(schur_root(gram)[0],float(LOWER))

    def test_quadratic_overlaps_match_direct_physical_operators(self):
        data, _, _ = witness_data()
        angles=np.array([.3,-.4,1.,2.])
        # The vector is defined using the fixed candidate's errors, whereas
        # the output overlaps use new errors. Compute these cross overlaps
        # directly in the physical 16-dimensional input coordinates.
        from flagged_subject_joining import joining_kraus
        from search_four_angle_history import unitaries as float_unitaries, TRUE
        from classical_joining_history import independent_encoding
        from joint_history_channel import bloch_state
        _, states, _ = witness_data(); omega=independent_encoding(tuple(bloch_state(np.array(s,float)) for s in states))
        ks=joining_kraus((2,2)); t=np.vstack((ks[(0,0)],ks[(0,1)]))
        b,residual=residual_blocks(tuple(data['coefficient_numerators']))
        block=np.zeros((16,16));block[:8,:8]=float(b)*np.eye(8);block[8:,8:]=np.array(residual,float)
        B=t.T@block@t; repair=float_unitaries(angles)
        parts=witness_components()
        for item in parts['reports']:
            g=item['report'];theta=angles[MODIFIED.index(g)];v=np.array([np.cos(theta),np.sin(theta)])
            for r in range(6):
                err=np.zeros((16,16));err[:8,:8]=np.eye(8);err[8:,8:]=repair[g].T@TRUE[r]
                expected=np.trace(omega @ B.T @ t.T @ err @ t)
                actual=float(item['overlap_constant_exact'])+np.array(item['overlap_linear_exact'][r],float)@v
                self.assertAlmostEqual(actual,expected.real,places=13)

    def test_bad_multiplier_is_rejected(self):
        item=witness_components()['reports'][0]
        with self.assertRaises(ValueError): circle_upper(item['constant'],item['linear'],item['quadratic'],F(0))
        with self.assertRaises(ValueError): circle_upper(item['constant'],item['linear'],item['quadratic'],.1)

    def test_positive_completion_identity_for_general_quadratic(self):
        q=np.array([[.2,.1],[.1,.4]]);l=np.array([.1,-.3]);k=.7;lam=.8;a=lam*np.eye(2)-q
        center=np.linalg.solve(a,l);upper=k+lam+l@center
        for theta in (.0,.7,2.,4.):
            v=np.array([np.cos(theta),np.sin(theta)])
            self.assertAlmostEqual(upper-(k+2*l@v+v@q@v),(v-center)@a@(v-center),places=14)

    def test_source_is_fixed_and_independent_of_chosen_angles(self):
        data, states, _=witness_data()
        self.assertTrue(all(sum(x*x for x in state)==1 for state in states))
        self.assertEqual(data['parameters_exact'],[str(t) for t in PARAMETERS])
        self.assertEqual(len(family_certificate()['report_circle_certificates']),4)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(FourAngleMinimaxBoundTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report=family_certificate()
    report['automated_checks']={'run':checks.testsRun,'failures':len(checks.failures),'errors':len(checks.errors)}
    if args.write_results:
        Path(__file__).with_name('four_angle_minimax_bound_results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__': main()
