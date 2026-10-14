"""1207 exact certificates for first-hit amplification and native reference.

No large Hilbert-space array, prior author-code import, or controller build.
All primary arithmetic uses Fraction. Default mode is read-only.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import argparse
import hashlib
import json
import math

HERE = Path(__file__).resolve().parent


def decimal_bound(x, upper=False, digits=30):
    """Exact directed decimal rounding via integer division."""
    scale = 10**digits
    value = -((-x.numerator*scale)//x.denominator) if upper else x.numerator*scale//x.denominator
    absolute = abs(value)
    return ('-' if value < 0 else '') + str(absolute//scale) + '.' + str(absolute % scale).zfill(digits)


def bounds(x):
    return [decimal_bound(x), decimal_bound(x, True)]


def first_hit_examples():
    p = (F(1,10), F(1,5), F(7,10))
    energy = (F(1,1000), F(0), F(2))
    q = p[0]+p[1]
    checked = []
    for n in range(1, 7):
        strings = list(product(range(3), repeat=n))
        def perm(s):
            out = list(s)
            for j, value in enumerate(s):
                if value in (0, 1):
                    out[j] = 1-value
                    break
            return tuple(out)
        probability = {s: math.prod(p[i] for i in s) for s in strings}
        for s in strings:
            assert perm(perm(s)) == s
        norm = max(abs(sum(energy[i] for i in perm(s))-sum(energy[i] for i in s)) for s in strings)
        distance = sum(abs(probability[s]-probability[perm(s)]) for s in strings)/2
        claimed = (1-(1-q)**n)*abs(p[0]-p[1])/q
        assert norm == F(1,1000) and distance == claimed
        checked.append({'blocks': n, 'basis_strings': len(strings),
                        'energy_operator_norm_exact': str(norm),
                        'trace_distance_exact': str(distance)})
    return checked


def signed_interval(z, boxes):
    lo = sum(c*(a if c >= 0 else b) for c, (a,b) in zip(z, boxes))
    hi = sum(c*(b if c >= 0 else a) for c, (a,b) in zip(z, boxes))
    return lo, hi


def native_example():
    path = HERE.parent/'1202'/'results.json'
    old = json.loads(path.read_text(encoding='utf-8'))
    p = [F(x) for x in old['one_marginal_rational']]
    denominator = old['preparation']['root_denominator']
    boxes = [tuple(F(x, denominator) for x in ab) for ab in old['preparation']['root_boxes']]
    assert sum(p) == 1 and all(p[i] > p[i+1] > 0 for i in range(3))
    assert all(boxes[i][1] < boxes[i+1][0] for i in range(3))

    # A two-copy inversion certifies non-Gibbs, but is not a work apparatus.
    activation_z = (-1, 1, 1, -1)
    gap_lo, gap_hi = signed_interval(activation_z, boxes)
    probability_difference = p[1]*p[2]-p[0]*p[3]
    assert gap_lo > 0 and probability_difference > 0
    assert gap_lo*probability_difference > F('1.1434966974e-7')

    # The old 1206 occupation difference, stripped of common spectators.
    x = (0, 15, 0, 10)
    y = (14, 0, 11, 0)
    z = tuple(a-b for a,b in zip(x,y))
    assert z == (-14, 15, -11, 10) and sum(x) == sum(y) == 25
    px = math.prod(p[i]**x[i] for i in range(4))
    py = math.prod(p[i]**y[i] for i in range(4))
    q = px+py
    kappa = abs(px-py)/q
    assert px > py > 0 and q < 1
    assert kappa > F('0.4479')
    n_blocks = (q.denominator+q.numerator-1)//q.numerator
    assert n_blocks == 4375693525386771
    assert (n_blocks-1)*q < 1 <= n_blocks*q
    total_materials = 25*n_blocks
    assert total_materials == 109392338134669275
    delta_lo, delta_hi = signed_interval(z, boxes)
    assert F('.0004237212509') < delta_lo < delta_hi < F('.0004237212511')

    # (1-q)^N <= exp(-Nq) <= exp(-1) < 60/163,
    # because exp(1) > sum_{j=0}^5 1/j! = 163/60.
    exponential_partial_sum = sum(F(1, math.factorial(j)) for j in range(6))
    assert exponential_partial_sum == F(163,60)
    hit_lower = F(103,163)
    distance_lower = F('0.4479')*hit_lower
    assert distance_lower == F(461337,1630000) and distance_lower > F('.283')
    return {
        'reference': '1202 one-sample cold marginal, independent copies',
        'full_rank_and_single_copy_passive': True,
        'non_Gibbs_certificate': {
            'higher_energy_string': [1,2], 'lower_energy_string': [0,3],
            'energy_gap_bounds': [decimal_bound(gap_lo), decimal_bound(gap_hi, True)],
            'higher_minus_lower_probability_bounds': bounds(probability_difference),
            'cyclic_swap_mean_energy_drop_bounds': [
                decimal_bound(gap_lo*probability_difference),
                decimal_bound(gap_hi*probability_difference, True)],
            'battery_or_control_implementation_claimed': False
        },
        'first_hit_witness': {
            'block_length': 25, 'x_occupation': list(x), 'y_occupation': list(y),
            'occupation_difference': list(z),
            'P_x_bounds': bounds(px), 'P_y_bounds': bounds(py), 'q_bounds': bounds(q),
            'conditional_trace_distance_bounds': bounds(kappa),
            'N_blocks_exact_decimal_integer': str(n_blocks),
            'N_choice': 'ceil(1/q); (N-1)q<1<=Nq verified exactly',
            'total_materials_exact_decimal_integer': str(total_materials),
            'log2_Hilbert_dimension_exact_decimal_integer': str(2*total_materials),
            'energy_operator_norm_bounds': [decimal_bound(delta_lo), decimal_bound(delta_hi, True)],
            'hit_probability_strict_lower': str(hit_lower),
            'trace_distance_strict_lower_rational': str(distance_lower),
            'declared_trace_distance_strict_lower': '0.283',
            'whole_matrix_materialized': False,
            'current_preparation_readout_or_control_resources_certified': False
        }
    }


def main():
    names = ['1202/results.json', '1206/proof.md', '1206/results.json']
    hashes = {name: hashlib.sha256((HERE.parent/name).read_bytes()).hexdigest() for name in names}
    result = {
        'schema': 'composition_robust_reference_1207_v1',
        'status': 'AUTHOR_CHECKS_PASS_NOT_ACCEPTANCE',
        'source_sha256': hashes,
        'exact_first_hit_small_examples': first_hit_examples(),
        'native_finite_witness': native_example(),
        'analytic_claims_require_proof': [
            'Gibbs all real beta: R(delta)<=min(1,sqrt(abs(beta)*delta/2))',
            'full-rank stationary non-Gibbs: R(delta)=1 for every delta>0',
            'native h: R(0)=0 for every diagonal reference by 1206',
            'suprema need not be achieved at finite n',
            'full operation class is a candidate permission, not native QQ reachability'
        ],
        'external_reference_scope': 'Only the material marginal is tested; R is retained, not reset.',
        'new_cognitive_axioms_adopted': 0
    }
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    output = HERE/'results.json'
    if args.write:
        output.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    saved = json.loads(output.read_text(encoding='utf-8'))
    assert saved == result, 'Saved result differs; default mode makes no changes.'
    print(json.dumps({'status': 'AUTHOR_CHECKS_PASS', 'all_checks_passed': True,
                      'results_match': True, 'small_first_hit_cases': 6,
                      'native_trace_distance_strict_lower': '0.283'}))


if __name__ == '__main__':
    main()
