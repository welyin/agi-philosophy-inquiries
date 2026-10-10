"""Independent round-1061 statistical audit; no author calculation imports.

Exact rational half-plane vertices / closest segment points replace the author's
monotone-endpoint solver. Algebraic square roots have rational enclosures.
Gaussian tails use integrated Taylor bounds and Machin's formula, not erfc.
Default reads and compares the saved independent result; --write is exclusive.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
from itertools import combinations
from hashlib import sha256
import argparse
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PAD = F(1, 10**38)
SCALE = 10**80


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def decimal_text(x, upper=False, places=30):
    unit = 10**places
    n = x.numerator * unit
    d = x.denominator
    k = -((-n)//d) if upper else n//d
    sign = '-' if k < 0 else ''
    k = abs(k)
    return f'{sign}{k//unit}.{k%unit:0{places}d}'


def sqrt_bounds(x):
    assert x >= 0
    n = isqrt((x.numerator*SCALE*SCALE)//x.denominator)
    lo, hi = F(n, SCALE), F(n+1, SCALE)
    assert lo*lo <= x <= hi*hi
    return lo, hi


def atan_bounds(x, terms=100):
    value = sum(((-1)**n)*x**(2*n+1)/F(2*n+1)
                for n in range(terms))
    following = ((-1)**terms)*x**(2*terms+1)/F(2*terms+1)
    return min(value, value+following), max(value, value+following)


def tail_bounds(q):
    a0, a1 = atan_bounds(F(1, 5))
    b0, b1 = atan_bounds(F(1, 239))
    pi0, pi1 = 16*a0-4*b1, 16*a1-4*b0
    den0 = sqrt_bounds(2*pi0)[0]
    den1 = sqrt_bounds(2*pi1)[1]
    # Integral_0^q exp(-t^2/2) dt. Taylor's integral remainder has
    # alternating sign and magnitude bounded by the next integrated term.
    value = F(0)
    factorial = 1
    terms = 160
    for n in range(terms):
        if n:
            factorial *= n
        value += ((-1)**n)*q**(2*n+1)/F(2**n*factorial*(2*n+1))
    following = q**(2*terms+1)/F(2**terms*factorial*terms*(2*terms+1))
    int0, int1 = value, value+following  # terms is even: next term positive.
    lo, hi = F(1, 2)-int1/den0, F(1, 2)-int0/den1
    assert 0 < lo < hi < F(1, 2)
    assert hi-lo < F(1, 10**75)
    return lo, hi


def make_box(data, contract, rounded, mode, outer):
    answer = {}
    for k in ('x', 'y', 'r', 'z'):
        item = data[k]
        mean = F(item['mean'])
        half = F(item['rounding_half_unit']) if rounded else F(0)
        if 'sigma_components' in item:
            c = [F(s)+half for s in item['sigma_components']]
            sig0, sig1 = sqrt_bounds(sum(v*v for v in c))
            if k == 'z' and mode == 'linear_widening':
                sig0 = sig1 = sum(c)
        else:
            sig0 = sig1 = F(item['sigma'])+half
        if contract == 'supplement_three_two_sided' and k == 'r':
            answer[k] = (F(0), F(1))
            continue
        item_c = data['contracts'][contract]
        q = F(item_c['q']) if 'q' in item_c else F(
            item_c['ratio_two_sided_q'] if k == 'r' else item_c['upper_q'])
        # Both versions enclose the author's padded ideal irrational box:
        # inner uses sigma_lower; outer uses sigma_upper.
        sig = sig1 if outer else sig0
        answer[k] = (max(F(0), mean-half-q*sig-PAD),
                     min(F(1), mean+half+q*sig+PAD))
        if contract == 'supplement_directional' and k != 'r':
            answer[k] = (F(0), answer[k][1])
    return answer


def polygon_extrema(box):
    xl, xu = box['x']; yl, yu = box['y']; rl, ru = box['r']
    hs = [(F(-1), F(0), -xl), (F(1), F(0), xu),
          (F(0), F(-1), -yl), (F(0), F(1), yu),
          (-ru, F(1), F(0)), (rl, F(-1), F(0))]
    valid = lambda v: all(a*v[0]+b*v[1] <= c for a,b,c in hs)
    vertices = set()
    for (a,b,c),(d,e,f) in combinations(hs, 2):
        det = a*e-b*d
        if det:
            v = ((c*e-b*f)/det, (a*f-c*d)/det)
            if valid(v):
                vertices.add(v)
    if not vertices:
        return None
    candidates = set(vertices)
    if valid((F(0),F(0))):
        candidates.add((F(0),F(0)))
    # In 2D the nearest boundary point lies on an edge. All vertex-pair
    # segments include every edge and stay inside this convex polygon.
    for u,v in combinations(vertices, 2):
        d = (v[0]-u[0],v[1]-u[1])
        dd = d[0]*d[0]+d[1]*d[1]
        if dd:
            t = max(F(0),min(F(1),-(u[0]*d[0]+u[1]*d[1])/dd))
            candidates.add((u[0]+t*d[0],u[1]+t*d[1]))
    norm = lambda v: v[0]*v[0]+v[1]*v[1]
    vlo = min(candidates,key=norm); vhi = max(vertices,key=norm)
    lo = norm(vlo)+box['z'][0]**2
    hi = norm(vhi)+box['z'][1]**2
    return lo, hi, len(vertices), vlo, vhi


def run():
    inputs = json.loads((HERE/'public_inputs.json').read_text(encoding='utf-8'))
    author = json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    receipt = json.loads((HERE/'research_round_1061_checks.json').read_text(encoding='utf-8'))
    for rel,h in receipt['owned_sha256'].items():
        assert sha(HERE/rel) == h, rel
    for rel,h in receipt['historical_sha256'].items():
        assert sha(ROOT/rel) == h, rel

    tails = {q:tail_bounds(F(q)) for q in ('2.50','2.40','2.25')}
    allocations = {
        'main_four_two_sided': [('2.50',8)],
        'supplement_three_two_sided': [('2.40',6)],
        'supplement_directional': [('2.25',3),('2.50',2)]}
    coverage = {}
    for key, terms in allocations.items():
        lo = 1-sum(n*tails[q][1] for q,n in terms)
        hi = 1-sum(n*tails[q][0] for q,n in terms)
        assert lo > F(95,100)
        saved = F(author['coverage_lower_bounds_under_adopted_marginals'][key])
        assert lo-F(1,10**15) < saved < hi+F(1,10**15)
        coverage[key] = [decimal_text(lo),decimal_text(hi,True)]

    reports = {}
    witness_errors = []
    for mode in inputs['z_sigma_modes']:
        for rounded in (False,True):
            tag = f'zsigma_{mode}_publication_rounding_{str(rounded).lower()}'
            reports[tag] = {}
            for contract in inputs['contracts']:
                inner = make_box(inputs,contract,rounded,mode,False)
                outer = make_box(inputs,contract,rounded,mode,True)
                a = polygon_extrema(inner); b = polygon_extrema(outer)
                assert a is not None and b is not None
                # The two rational sandwiches give the same strict verdict.
                feasible_inner = a[0] <= 1 <= a[1]
                feasible_outer = b[0] <= 1 <= b[1]
                assert feasible_inner == feasible_outer
                saved = author['results'][tag][contract]
                assert feasible_outer == saved['feasible']
                assert F(saved['norm_min']) <= b[0]
                assert F(saved['norm_max']) >= b[1]
                assert abs(F(saved['norm_min'])-b[0]) < F(2,10**24)
                assert abs(F(saved['norm_max'])-b[1]) < F(2,10**24)
                for k in outer:
                    l,u = map(F,saved['box'][k])
                    assert l <= outer[k][0] and u >= outer[k][1]
                    assert abs(l-outer[k][0]) < F(2,10**24)
                    assert abs(u-outer[k][1]) < F(2,10**24)
                if 'witness' in saved:
                    w = {k:F(v) for k,v in saved['witness'].items()}
                    tol = F(3,10**23)
                    assert all(l-tol <= w[k] <= u+tol for k,(l,u) in inner.items())
                    e1 = abs(w['y']-w['x']*w['r'])
                    e2 = abs(w['x']**2+w['y']**2+w['z']**2-1)
                    assert max(e1,e2) < tol
                    witness_errors.append(max(e1,e2))
                reports[tag][contract] = {
                    'feasible':feasible_outer,'vertices':b[2],
                    'norm_min_enclosure':[decimal_text(b[0]),decimal_text(a[0],True)],
                    'norm_max_enclosure':[decimal_text(a[1]),decimal_text(b[1],True)],
                    'directional_ratio_nonbinding_at_max': (
                        b[4] == (outer['x'][1],outer['y'][1])
                        if contract == 'supplement_directional' else None)}
    return {
        'round':1061,'independent_checks_passed':True,
        'new_scientific_calibration_groups':0,
        'independent_method':'exact rational half-plane vertices and segment projections; rational sqrt and Gaussian-tail enclosures',
        'author_calculation_imported':False,'contracts_checked':12,
        'author_assets_checked':len(receipt['owned_sha256']),
        'historical_assets_checked':len(receipt['historical_sha256']),
        'author_receipt_sha256':sha(HERE/'research_round_1061_checks.json'),
        'author_assets_sha256':receipt['owned_sha256'],
        'coverage_enclosures_conditional_on_adopted_marginals':coverage,
        'all_published_endpoints_enclose_the_exact_target_to_display_precision':True,
        'max_printed_witness_residual_upper':decimal_text(max(witness_errors),True,30),
        'results':reports,
        'directional_selection_adjusted':False,
        'physical_joint_likelihood_or_cognitive_generation_validated':False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args = parser.parse_args()
    result = run()
    target = HERE/'independent_checks.json'
    if args.write:
        with target.open('x',encoding='utf-8') as out:
            json.dump(result,out,ensure_ascii=False,indent=2);out.write('\n')
    else:
        assert result == json.loads(target.read_text(encoding='utf-8'))
    print(json.dumps({'round':1061,'independent_checks_passed':True,
                      'contracts':12,'assets':result['author_assets_checked'],
                      'new_scientific_calibration_groups':0}))
