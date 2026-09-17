"""Round 169: certify improvement beyond the entire four-angle family.

Reuses round 161's ball/box inequalities, with new exact branch matrices.
The optional cover search never changes the round 161 cover or results.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
import time
import unittest

import numpy as np

from certified_intervals import Interval as I, SCALE
from relational_history_decoder import PARAMETER, branch_errors, error_jet, gram_from_sigma
from global_history_recovery import (ROOT, REAL_PRODUCTS, box_bound, outside_ball,
    split_box, maximum)
from joint_history_certificate import exact_product_sigma, stereographic_bloch
from joint_relation_records import interval_fields
from noisy_classical_history import majority_error
from soft_decoder_certificate import coefficients, DENOMINATOR, partial_reference, block_expectation
from four_angle_minimax_bound import LOWER as PREVIOUS_LOWER
from soft_history_decoder import branch_weights, schur_root


UPPER = F(44916, 10**5)
LOWER = F(4491536259, 10**10)
COVER_PATH = Path(__file__).with_name('relational_history_cover.json')
WITNESS_PATH = Path(__file__).with_name('relational_history_witness.json')


def residual_blocks(integers):
    a = coefficients(integers)
    residual = np.eye(8, dtype=object)
    for weight, w in zip(a, branch_errors()):
        residual = residual-weight*w
    return 1-sum(a), residual


@lru_cache(None)
def dual_data(integers):
    a = coefficients(integers)
    d = 1-sum(I.exact(x*x)/p for x, p in zip(a, branch_weights(majority_error(1))))
    if d.lo <= 0:
        raise ValueError('The positive completion denominator failed.')
    b, residual = residual_blocks(integers)
    numerator = b*b*np.eye(4, dtype=object)+partial_reference(residual.T @ residual)/2
    return d, numerator


@lru_cache(None)
def supporting_polynomial(integers):
    d, numerator = dual_data(integers)
    return tuple(tuple(sum(numerator[j, i]*int(REAL_PRODUCTS[r][s][i, j])
                           for i in range(4) for j in range(4))/(8*d)
                       for s in range(4)) for r in range(4))


def witness_data():
    data = json.loads(WITNESS_PATH.read_text(encoding='utf-8'))
    if data['parameter_exact'] != str(PARAMETER):
        raise ValueError('Wrong recovery parameters in the witness.')
    if data.get('coefficient_denominator') != DENOMINATOR:
        raise ValueError('Wrong coefficient denominator in the witness.')
    if (type(data.get('stereographic_denominator')) is not int or data['stereographic_denominator'] <= 0
            or len(data['stereographic_numerators']) != 4
            or any(type(n) is not int for n in data['stereographic_numerators'])):
        raise ValueError('Use four exact stereographic integers and a positive denominator.')
    p = [F(n, data['stereographic_denominator']) for n in data['stereographic_numerators']]
    states = (stereographic_bloch(*p[:2]), stereographic_bloch(*p[2:]))
    if any(sum(x*x for x in state) != 1 for state in states):
        raise ValueError('The witness must use two legal independent pure states.')
    return data, states, exact_product_sigma(*states)


def witness_lower():
    data, _, sigma = witness_data()
    b, residual = residual_blocks(tuple(data['coefficient_numerators']))
    norm = block_expectation(sigma, b*b, residual.T @ residual)
    target = block_expectation(sigma, b, residual.T)
    overlap = [block_expectation(sigma, b, residual.T @ w) for w in branch_errors()]
    if norm <= 0:
        raise ValueError('The Rayleigh vector has zero norm.')
    value = (target*target-sum(p*x*x for p, x in zip(branch_weights(majority_error(1)), overlap)))/norm
    if F(value.lo, SCALE) <= LOWER:
        raise ValueError('The legal lower witness failed.')
    return value


def verify_cover(data):
    if data.get('format') != 1 or data.get('root') != '[-1,1] x [0,1] x [-1,1]':
        raise ValueError('Unsupported cover format or domain.')
    if data.get('parameter_exact') != str(PARAMETER):
        raise ValueError('Wrong recovery parameters.')
    if data.get('coefficient_denominator') != DENOMINATOR or F(data['upper_exact']) != UPPER:
        raise ValueError('Wrong coefficient denominator or endpoint.')
    duals = [tuple(row) for row in data['dual_coefficient_numerators']]
    trie = {}
    for leaf in data['leaves']:
        path = leaf['path']
        if not isinstance(path, str) or any(x not in '01' for x in path):
            raise ValueError('Use binary partition paths.')
        node = trie
        for bit in path:
            if 'leaf' in node:
                raise ValueError('Overlapping leaves.')
            node = node.setdefault(bit, {})
        if node:
            raise ValueError('Overlapping or duplicate leaves.')
        node['leaf'] = leaf
    stack = [(trie, ROOT, 0)]
    values, outside, depth_max = [], 0, 0
    while stack:
        node, box, depth = stack.pop()
        depth_max = max(depth_max, depth)
        if 'leaf' in node:
            leaf = node['leaf']
            if leaf.get('outside') is True:
                if not outside_ball(box):
                    raise ValueError('An outside claim intersects the unit ball.')
                outside += 1
            else:
                index = leaf.get('dual')
                if type(index) is not int or not 0 <= index < len(duals):
                    raise ValueError('Invalid dual index.')
                bound = box_bound(supporting_polynomial(duals[index]), box)
                if F(bound.hi, SCALE) >= UPPER:
                    raise ValueError('A region does not certify the upper endpoint.')
                values.append(bound)
        else:
            if set(node) != {'0', '1'}:
                raise ValueError('An uncovered region remains.')
            left, right = split_box(box)
            stack.extend(((node['0'], left, depth+1), (node['1'], right, depth+1)))
    if not values:
        raise ValueError('No legal source region certified.')
    return {'round':169, 'parameter_exact':str(PARAMETER),
            'strict_lower_exact':str(LOWER), 'strict_upper_exact':str(UPPER),
            'legal_witness_lower':interval_fields(witness_lower()),
            'cover_upper':interval_fields(maximum(values)),
            'previous_four_angle_family_lower_exact':str(PREVIOUS_LOWER),
            'strict_improvement_over_entire_four_angle_family_at_least_exact':str(PREVIOUS_LOWER-UPPER),
            'certified_boxes':len(values), 'outside_boxes':outside, 'max_binary_depth':depth_max,
            'saved_supporting_matrices':len(duals), 'all_mixed_independent_sources_covered':True,
            'all_external_extensions_covered':True, 'three_original_reads_and_all_reports_retained':True,
            'coherent_history_rebits':1, 'certificate_requires_scipy':False,
            'new_native_gates_counted_separately_in_round_170':True,
            'general_recovery_optimality_proved':False}


def generate_cover():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'.research_runtime'))
    from scipy.optimize import minimize
    start = time.monotonic()
    data, states, _ = witness_data()
    duals, leaves = [tuple(data['coefficient_numerators'])], []
    stack = [('', ROOT, 0, np.array(states[1], dtype=float))]
    nodes = 0
    while stack:
        path, box, parent, initial = stack.pop(); nodes += 1
        if nodes > 30000:
            raise RuntimeError('Node budget exceeded; no incomplete certificate written.')
        if outside_ball(box):
            leaves.append({'path':path, 'outside':True}); continue
        if F(box_bound(supporting_polynomial(duals[parent]), box).hi, SCALE) < UPPER:
            leaves.append({'path':path, 'dual':parent}); continue
        center = np.array([(lo+hi)/2 for lo, hi in zip(*box)], dtype=float)
        center /= max(1., np.linalg.norm(center))
        def objective(v):
            norm = np.linalg.norm(v); legal = v/max(1., norm)
            jet = error_jet(np.r_[center, legal]); gradient = jet['gradient'][3:]
            if norm > 1:
                gradient = (gradient-legal*(legal @ gradient))/norm
            return -jet['value'], -gradient
        result = minimize(objective, initial, jac=True, method='SLSQP',
                          constraints={'type':'ineq','fun':lambda v:1-v @ v,'jac':lambda v:-2*v},
                          options={'ftol':1e-14,'maxiter':80})
        v = result.x/max(1., np.linalg.norm(result.x))
        jet = error_jet(np.r_[center, v])
        if jet['value'] > float(UPPER)+1e-13:
            raise RuntimeError('Numerical counterexample exceeds the proposed endpoint.')
        integers = tuple(int(x) for x in np.rint(jet['coefficients']*DENOMINATOR))
        duals.append(integers); index = len(duals)-1
        if F(box_bound(supporting_polynomial(integers), box).hi, SCALE) < UPPER:
            leaves.append({'path':path,'dual':index})
        else:
            left, right = split_box(box)
            stack.extend(((path+'1',right,index,v),(path+'0',left,index,v)))
        if nodes % 100 == 1:
            print(f'Cover: {nodes} nodes, {len(stack)} pending, {time.monotonic()-start:.1f}s', flush=True)
    used = sorted({leaf['dual'] for leaf in leaves if 'dual' in leaf})
    mapping = {old:new for new, old in enumerate(used)}
    for leaf in leaves:
        if 'dual' in leaf:
            leaf['dual'] = mapping[leaf['dual']]
    report = {'format':1,'root':'[-1,1] x [0,1] x [-1,1]',
              'parameter_exact':str(PARAMETER), 'upper_exact':str(UPPER),
              'coefficient_denominator':DENOMINATOR,
              'dual_coefficient_numerators':[duals[i] for i in used], 'leaves':leaves,
              'candidate_search_nodes':nodes}
    verify_cover(report)
    COVER_PATH.write_text(json.dumps(report,separators=(',',':'))+'\n',encoding='utf-8')


class RelationalHistoryCertificateTests(unittest.TestCase):
    def test_saved_cover_proves_strict_improvement(self):
        result = verify_cover(json.loads(COVER_PATH.read_text(encoding='utf-8')))
        self.assertGreater(F(result['strict_improvement_over_entire_four_angle_family_at_least_exact']), F(234,10**6))
        self.assertTrue(result['all_mixed_independent_sources_covered'])

    def test_legal_lower_and_full_external_gram(self):
        _, _, sigma = witness_data()
        self.assertAlmostEqual(witness_lower().floats()[0], schur_root(gram_from_sigma(np.array(sigma,dtype=float)))[0], places=12)

    def test_missing_region_is_rejected(self):
        data = json.loads(COVER_PATH.read_text(encoding='utf-8')); data['leaves'].pop()
        with self.assertRaisesRegex(ValueError,'uncovered'): verify_cover(data)

    def test_wrong_parameters_and_coefficients_are_rejected(self):
        data = json.loads(COVER_PATH.read_text(encoding='utf-8')); data['parameter_exact']='0'
        with self.assertRaisesRegex(ValueError,'parameters'): verify_cover(data)
        with self.assertRaises(ValueError): dual_data((2*DENOMINATOR,)+(0,)*35)

    def test_false_outside_and_overlaps_are_rejected(self):
        data = json.loads(COVER_PATH.read_text(encoding='utf-8')); data['leaves']=[{'path':'','outside':True}]
        with self.assertRaisesRegex(ValueError,'intersects'): verify_cover(data)
        data['leaves'].append({'path':'0','outside':True})
        with self.assertRaisesRegex(ValueError,'Overlapping'): verify_cover(data)

    def test_new_supports_bound_unseen_sources(self):
        data, _, _ = witness_data(); ints = tuple(data['coefficient_numerators'])
        poly = np.array([[np.mean(x.floats()) for x in row] for row in supporting_polynomial(ints)])
        rng = np.random.default_rng(269)
        for _ in range(20):
            u,v = rng.normal(size=(2,3)); u/=max(1.,np.linalg.norm(u)); v/=max(1.,np.linalg.norm(v))
            self.assertLessEqual(error_jet(np.r_[u,v])['value'], np.r_[1.,u] @ poly @ np.r_[1.,v]+2e-14)

    def test_uninformative_certificate_cannot_certify_small_error(self):
        data = json.loads(COVER_PATH.read_text(encoding='utf-8'))
        data['dual_coefficient_numerators']=[[0]*36 for _ in data['dual_coefficient_numerators']]
        with self.assertRaisesRegex(ValueError,'does not certify'): verify_cover(data)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate-cover',action='store_true')
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    if args.generate_cover: generate_cover()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RelationalHistoryCertificateTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report=verify_cover(json.loads(COVER_PATH.read_text(encoding='utf-8')))
    report['cover_sha256']=hashlib.sha256(COVER_PATH.read_bytes()).hexdigest()
    report['automated_checks']={'run':checks.testsRun,'failures':len(checks.failures),'errors':len(checks.errors)}
    if args.write_results:
        Path(__file__).with_name('relational_history_certificate_results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__ == '__main__': main()
