"""Round 424: one contracting endpoint automorphism in place of all dilations.

Analytic group classification is cited in the note. These checks exercise the
strictly nonhomogeneous budget, exact finite certificates, and old coordinate
decoder under the new declared contract. No old test suite is run.
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
from direction_only_coordinate_audit import EffectOracle, DirectionChart, SIGMA, PREPS

OBS = {}
S = np.array([[1.2, .3, -.1], [0., .8, .2], [0., 0., 1.5]])
SI = np.linalg.inv(S)
Q = S.T @ S
ALPHA = SI @ np.diag([-.5, .25, .125]) @ S


def phi(r):
    return r + min(r, 1)


def budget(v):
    return phi(float(np.linalg.norm(S @ np.asarray(v))))


def unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def rotate(axis, angle):
    x, y, z = unit(axis)
    cross = np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])
    return np.eye(3) + math.sin(angle)*cross + (1-math.cos(angle))*cross@cross


def cubic_rotations():
    result = []
    for permutation in itertools.permutations(range(3)):
        p = np.eye(3)[:, permutation]
        for signs in itertools.product((-1., 1.), repeat=3):
            r = p @ np.diag(signs)
            if np.linalg.det(r) > 0:
                result.append(r)
    return result


def log_periodic(r):
    return 0. if r == 0 else r*(1+.5*math.sin(2*math.pi*math.log2(r)))


def root(left, right):
    fl = log_periodic(left)-1
    for _ in range(70):
        middle = (left+right)/2
        fm = log_periodic(middle)-1
        if fl*fm <= 0:
            right = middle
        else:
            left, fl = middle, fm
    return (left+right)/2


class Audit(unittest.TestCase):
    def test_01_nonlinear_budget_is_a_metric(self):
        rng = np.random.default_rng(42401)
        worst = -math.inf
        symmetry = 0.
        for _ in range(300):
            a, b = rng.normal(size=(2,3))*rng.uniform(.01,10.)
            worst = max(worst, budget(a+b)-budget(a)-budget(b))
            symmetry = max(symmetry, abs(budget(-a)-budget(a)))
            self.assertGreater(budget(a), 0.)
        self.assertLessEqual(worst, 1e-12)
        self.assertEqual(budget(np.zeros(3)), 0.)
        self.assertEqual(symmetry, 0.)
        OBS['metric_triangle_max_defect'] = worst
        OBS['metric_inverse_symmetry_error'] = symmetry

    def test_02_sharp_discrete_contraction_rational_certificate(self):
        rows = []
        for a in (Fraction(1,4), Fraction(1,2), Fraction(3,4)):
            q = 2*a/(1+a)
            radii = [Fraction(k,19) for k in range(1,160)] + [1/a]
            ratios = [phi(a*r)/phi(r) for r in radii]
            self.assertLessEqual(max(ratios), q)
            self.assertEqual(phi(a/a)/phi(1/a), q)
            # Independent piecewise upper certificates for all radii are in the note.
            rows.append(dict(linear_norm=str(a), sharp_ratio=str(q), attained_radius=str(1/a)))
        OBS['exact_contraction_certificates'] = rows

    def test_03_anisotropic_automorphism_and_iterates(self):
        rng = np.random.default_rng(42403)
        law_error = 0.
        ratio = 0.
        iterate_error = 0.
        for _ in range(70):
            a, b = rng.normal(size=(2,3))*4
            law_error = max(law_error, float(np.linalg.norm(ALPHA@(a+b)-ALPHA@a-ALPHA@b)))
            ratio = max(ratio, budget(ALPHA@a)/budget(a))
            current = a.copy()
            for k in range(1,16):
                current = ALPHA@current
                iterate_error = max(iterate_error, budget(current)-(2/3)**k*budget(a))
        self.assertLess(law_error, 1e-12)
        self.assertLessEqual(ratio, 2/3+1e-12)
        self.assertLess(iterate_error, 1e-12)
        singular = np.linalg.svd(S@ALPHA@SI, compute_uv=False)
        self.assertTrue(np.allclose(singular,[.5,.25,.125]))
        self.assertLess(np.linalg.det(ALPHA), 0.)
        OBS['endpoint_addition_error'] = law_error
        OBS['sampled_contraction_ratio'] = ratio
        OBS['alpha_metric_singular_values'] = singular.tolist()
        OBS['alpha_determinant'] = float(np.linalg.det(ALPHA))

    def test_04_no_nontrivial_exact_budget_homothety(self):
        # Small radii force |delta u|=s|u|; large radii then force 1=s.
        rows = []
        for s in (Fraction(1,2), Fraction(2), Fraction(3)):
            small = Fraction(1,10)
            large = Fraction(10)
            self.assertEqual(phi(s*small), s*phi(small))
            mismatch = phi(s*large)-s*phi(large)
            self.assertEqual(mismatch, 1-s)
            self.assertNotEqual(mismatch, 0)
            rows.append(dict(candidate_factor=str(s), exact_large_radius_mismatch=str(mismatch)))
        self.assertEqual(phi(Fraction(1,2)),1)
        self.assertEqual(phi(Fraction(1)),2)
        self.assertEqual(phi(Fraction(2)),3)
        OBS['homogeneity_obstruction_certificates'] = rows
        OBS['old_power_calibration_counterexample'] = dict(unit_budget=1,double_budget=2,
            quadruple_budget=3,power_prediction_from_first_two=4)

    def test_05_compact_average_and_full_shell(self):
        rotations = cubic_rotations()
        self.assertEqual(len(rotations),24)
        h = [SI@r@S for r in rotations]
        average = sum(a.T@a for a in h)/len(h)
        scale = np.trace(SI.T@SI)/3
        self.assertLess(np.linalg.norm(average-scale*Q),1e-12)
        worst_invariance = max(np.linalg.norm(a.T@average@a-average) for a in h)
        rng = np.random.default_rng(42405)
        shell_error = 0.
        for _ in range(50):
            v = SI@unit(rng.normal(size=3))/2
            shell_error = max(shell_error,abs(budget(v)-1))
            r = rotate(rng.normal(size=3),rng.uniform(-3,3))
            self.assertAlmostEqual(budget(SI@r@S@v),1)
        self.assertLess(worst_invariance,1e-12)
        self.assertLess(shell_error,1e-12)
        OBS['finite_quadratic_average_error'] = float(np.linalg.norm(average-scale*Q))
        OBS['averaged_form_invariance_error'] = float(worst_invariance)
        OBS['shell_budget_error'] = float(shell_error)
        # 24 rotations check the averaging formula; density is not inferred from a finite sample.

    def test_06_same_qubit_full_direction_covariance(self):
        rng = np.random.default_rng(42406)
        contrast = .8
        worst = 0.
        directions = np.vstack([np.eye(3),-np.eye(3)])
        effects = [np.eye(2)/2+contrast/2*np.einsum('i,ijk->jk',n,SIGMA) for n in directions]
        self.assertEqual(np.linalg.matrix_rank(directions[1:]-directions[0]),3)
        for _ in range(50):
            n, axis = unit(rng.normal(size=3)),unit(rng.normal(size=3))
            angle = rng.uniform(-3,3)
            r = rotate(axis,angle)
            generator = np.einsum('i,ijk->jk',axis,SIGMA)
            u = math.cos(angle/2)*np.eye(2)-1j*math.sin(angle/2)*generator
            e = np.eye(2)/2+contrast/2*np.einsum('i,ijk->jk',n,SIGMA)
            er = np.eye(2)/2+contrast/2*np.einsum('i,ijk->jk',r@n,SIGMA)
            worst = max(worst,float(np.linalg.norm(u@e@u.conj().T-er)))
        self.assertLess(worst,1e-12)
        self.assertTrue(all(np.linalg.eigvalsh(e).min() > 0 for e in effects))
        OBS['qubit_covariance_error'] = worst
        OBS['antipodal_operator_gap'] = float(np.linalg.norm(effects[0]-effects[3],2))

    def test_07_existing_coordinates_without_budget_homogeneity(self):
        # Reuse the frozen 413 decoder, only supplying token probabilities.
        oracle = EffectOracle(frame=S)
        inputs = [np.zeros(3),np.array([1.,0.,0.]),np.array([.2,.9,.3]),
                  np.array([.7,-.4,.8]),np.array([2.,0.,0.]),np.array([0.,0.,2.])]
        inputs += [ALPHA@v for v in inputs[3:]]
        tokens = [oracle.add(v) for v in inputs]
        chart = DirectionChart(oracle.probabilities,tokens[:3])
        baseline = np.linalg.norm(S@inputs[1])
        errors = [np.linalg.norm(chart.locate(t)-S@v/baseline) for t,v in zip(tokens,inputs)]
        maximum = float(max(errors))
        self.assertLess(maximum,1e-12)
        OBS['old_decoder_on_new_contract_error'] = maximum
        OBS['new_contract_coordinate_queries'] = oracle.query_count
        OBS['coordinate_budget_queries'] = 0

    def test_08_metric_is_not_its_path_length(self):
        direct = phi(Fraction(4))
        split = [n*phi(Fraction(4,n)) for n in (1,2,4,8,16)]
        self.assertEqual(direct,5)
        self.assertEqual(split[-1],8)
        self.assertEqual(split,[5,6,8,8,8])
        OBS['distance_of_four_unit_separation'] = float(direct)
        OBS['refined_straight_path_costs'] = [float(x) for x in split]

    def test_09_discrete_contraction_alone_does_not_make_one_shell(self):
        endpoints = [1.,math.sqrt(2),2**.75,2.]
        values = [log_periodic(r) for r in endpoints]
        self.assertAlmostEqual(values[0],1.)
        self.assertGreater(values[1],1.)
        self.assertLess(values[2],1.)
        self.assertGreater(values[3],1.)
        radii = [1.,root(endpoints[1],endpoints[2]),root(endpoints[2],endpoints[3])]
        self.assertTrue(radii[0] < radii[1] < radii[2])
        self.assertTrue(all(abs(log_periodic(r)-1)<1e-12 for r in radii))
        error = max(abs(log_periodic(r/2)-log_periodic(r)/2) for r in np.geomspace(.001,100,100))
        self.assertLess(error,1e-12)
        OBS['multiple_shell_radii_without_minimality'] = radii
        OBS['log_periodic_discrete_scaling_error'] = float(error)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise SystemExit(1)
    saved = dict(round=424,baseline_round=423,status='verified_conditional_premise_reduction',
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,observations=OBS,
        scope=dict(single_contracting_endpoint_automorphism_suffices=True,
            full_positive_scale_action_required=False,exact_budget_homogeneity_required=False,
            actual_endpoint_torsor_still_input=True,connected_position_sector_still_input=True,
            proper_continuous_budget_still_input=True,composition_preserving_minimal_reorientation_still_input=True,
            full_qubit_direction_interface_still_input=True,complete_displacement_effect_queries_still_input=True,
            old_budget_power_conclusion_retained=False,old_direction_decoder_reused=True,
            contraction_alone_makes_spherical_shell=False,budget_identified_as_path_length=False,
            all_contractions_claimed_physically_implemented=False,three_dimensional_space_unconditionally_derived=False,
            full_cognitive_countermodel_completed=False,full_GR_goal_completed=False,phase_closure_triggered=False))
    if not args.check:
        with Path(__file__).with_name('single_contraction_coordinate_audit_results.json').open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(saved,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    print(json.dumps(saved,ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()
