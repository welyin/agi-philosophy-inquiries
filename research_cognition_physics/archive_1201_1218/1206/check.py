"""1206: exact permissions for common-reference selection on the native h.

Primary certificates use integers/Fraction. NumPy is a diagnostic only.
No import or execution of earlier author code. Earlier cold population data
are read with their exact rational representation and source hash recorded.
Default mode compares the saved result; --write creates/updates this round only.
"""
from fractions import Fraction as F
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent


def polynomial(x):
    return 32000*x**3 - 31200*x*x - 2095*x - 17


def out_decimal(x, upper=False):
    with localcontext() as ctx:
        ctx.prec = 100
        y = Decimal(x.numerator)/Decimal(x.denominator)
        return str(y.quantize(Decimal('1e-30'), rounding=ROUND_CEILING if upper else ROUND_FLOOR))


def root_boxes():
    out = []
    for a, b in [('-0.054072', '-0.054071'), ('-0.009461', '-0.009460'),
                 ('1.038532', '1.038533')]:
        a, b = F(a), F(b)
        assert polynomial(a)*polynomial(b) < 0
        for _ in range(180):
            m = (a+b)/2
            if polynomial(a)*polynomial(m) < 0:
                b = m
            else:
                a = m
        assert b-a < F('1e-59')
        out.append((a, b))
    return out[:2] + [(F(1), F(1))] + out[2:]


def rank(rows):
    a = [[F(x) for x in r] for r in rows]
    k = 0
    for j in range(len(a[0])):
        p = next((i for i in range(k, len(a)) if a[i][j]), None)
        if p is None:
            continue
        a[k], a[p] = a[p], a[k]
        s = a[k][j]
        a[k] = [x/s for x in a[k]]
        for i in range(len(a)):
            if i != k:
                s = a[i][j]
                a[i] = [x-s*y for x, y in zip(a[i], a[k])]
        k += 1
        if k == len(a):
            break
    return k


def native_graph():
    c = 2
    edges = []
    for a in range(16):
        i, j = divmod(a, 4)
        for b in range(a+1, 16):
            k, l = divmod(b, 4)
            if ((i == c) != (k == c)) and ((j == c) != (l == c)):
                edges.append((a, b))
    remaining = set(range(16))
    components = []
    while remaining:
        component = {min(remaining)}
        changed = True
        while changed:
            old = len(component)
            for a, b in edges:
                if a in component or b in component:
                    component.update((a, b))
            changed = len(component) != old
        remaining -= component
        components.append(sorted(component))
    rows = []
    for a, b in edges:
        row = [0]*16
        row[a], row[b] = 1, -1
        rows.append(row)
    assert sorted(map(len, components)) == [6, 10]
    assert rank(rows) == 14
    return {'edges': len(edges), 'components': components, 'diagonal_stationarity_rank': 14}


def minimum_gap(boxes, n=30):
    entries = []
    for a in range(n+1):
        for b in range(n+1-a):
            for c in range(n+1-a-b):
                v = (a, b, c, n-a-b-c)
                lo = sum(x*t[0] for x, t in zip(v, boxes))
                hi = sum(x*t[1] for x, t in zip(v, boxes))
                entries.append(((lo+hi)/2, v, lo, hi))
    entries.sort()
    gaps = []
    for x, y in zip(entries, entries[1:]):
        lo, hi = y[2]-x[3], y[3]-x[2]
        assert 0 < lo <= hi
        gaps.append((lo, hi, x[1], y[1]))
    lower = min(x[0] for x in gaps)
    witness = min(gaps, key=lambda x: x[1])
    upper = witness[1]
    assert upper-lower < F('1e-55')
    assert F('.0004237212509') < lower < upper < F('.0004237212510')
    z = [b-a for a, b in zip(witness[2], witness[3])]
    assert z == [-14, 15, -11, 10]
    count = math.comb(n+3, 3)
    assert len(entries) == count == 5456
    return {'n': n, 'occupation_types': count,
            'minimum_nonzero_gap_lower': out_decimal(lower),
            'minimum_nonzero_gap_upper': out_decimal(upper, True),
            'witness_lower_energy_counts': list(witness[2]),
            'witness_upper_energy_counts': list(witness[3]),
            'occupation_difference': z,
            'pigeonhole_bound_using_width_1_094': str(F(n)*F('1.094')/(count-1)),
            'all_adjacent_interval_orderings_certified': True,
            'physical_swap_implementation_claimed': False}


def cold_witness(boxes):
    path = HERE.parent/'1202'/'results.json'
    raw = path.read_bytes()
    old = json.loads(raw)
    r = [F(x) for x in old['one_marginal_rational']]
    assert len(r) == 4 and sum(r) == 1 and min(r) > 0
    a, b = boxes[0]
    # q_0^2 is the squared native |1> component of the e_0 eigenvector.
    aa = (1/(25*(1-a)**2), 1/(25*(1-b)**2))
    bb = (3/(6400*(a+F(1,40))**2), 3/(6400*(b+F(1,40))**2))
    assert 0 < aa[0] < aa[1] and 0 < bb[0] < bb[1]
    q2 = (aa[0]/(1+aa[1]+bb[1]), aa[1]/(1+aa[0]+bb[0]))
    difference = r[0]**2-r[2]**2
    bounds = [F(1,10)*q*difference for q in q2]
    assert bounds[0] > F('.0003202089318')
    assert bounds[1] < F('.0003202089319')
    return {'preparation': 'two independent cold pools; one sample from each',
            'source': '../1202/results.json', 'source_sha256': hashlib.sha256(raw).hexdigest(),
            'cold_population_exact': list(map(str, r)),
            'entry': '<00|[H0+0.1 QQ,r_C tensor r_C]|cc>; c=2',
            'entry_absolute_lower': out_decimal(bounds[0]),
            'entry_absolute_upper': out_decimal(bounds[1], True),
            'declared_strict_lower': '0.0003202089318',
            'finite_time_readout_difference_claimed': False}


def calibration_example():
    rows = []
    for a, b in [(0, 1), (1, 2)]:
        for i in [1, 2]:
            row = [0]*9
            row[3*a+i] += 1
            row[3*a] -= 1
            row[3*b+i] -= 1
            row[3*b] += 1
            rows.append(row)
    conversion = [1, -2, 1]+[0]*6
    assert rank(rows) == 4
    assert rank(rows+[conversion]) == 5
    assert rank(rows[:2]+[conversion]) == 3
    r = [F(1, 2), F(1, 3), F(1, 6)]
    diff = r[1]**2-r[0]*r[2]
    assert diff == F(1, 36)
    sv = np.linalg.svd(np.array(rows+[conversion], dtype=float), compute_uv=False)
    return {'spectrum': [0, 1, 2], 'subjects': 3, 'edges': [[0, 1], [1, 2]],
            'conversion': conversion, 'exchange_only_rank': 4,
            'with_one_conversion_rank': 5, 'delete_last_edge_rank': 3,
            'counterexample_probability': list(map(str, r)),
            'conversion_population_imbalance': str(diff),
            'smallest_nonzero_singular_value_diagnostic': round(float(min(sv)), 12),
            'native_material_permission': False, 'beta_sign_selected': False}


def main():
    # Direct determinant formula on the native 3x3 block, checked as polynomials.
    for x in map(F, [-2, -1, 0, 1, 2]):
        determinant = x*(x-1)*(x+F(1,40)) - F(1,25)*(x+F(1,40)) - F(3,6400)*(x-1)
        assert 32000*determinant == polynomial(x)
    residues = [polynomial(x) % 11 for x in range(11)]
    assert residues == [5, 8, 9, 3, 7, 5, 3, 7, 1, 2, 5]
    assert 32000 % 11 and 0 not in residues
    a, b, c, d = 32000, -31200, -2095, -17
    delta = b*b*c*c-4*a*c**3-4*b**3*d-27*a*a*d*d+18*a*b*c*d
    k = math.isqrt(delta)
    assert delta == 2736128576000000 and k == 52308016
    assert k*k < delta < (k+1)**2
    assert F(3)-F(39,40) == F(81,40)
    assert polynomial(F(1)) != 0 and polynomial(F(-1,40)) != 0
    boxes = root_boxes()
    source_names = ['1201/proof.md', '1201/results.json', '1202/results.json',
                    '1204/results.json', '1205/results.json']
    source_hashes = {name: hashlib.sha256((HERE.parent/name).read_bytes()).hexdigest()
                     for name in source_names}
    result = {'schema': 'native_reference_permission_1206_v1',
              'status': 'AUTHOR_CHECKS_PASS_NOT_ACCEPTANCE',
              'source_sha256': source_hashes,
              'algebraic_certificate': {'polynomial_descending': [a,b,c,d],
                  'irreducibility_prime': 11, 'mod_11_values': residues,
                  'discriminant': delta, 'sqrt_floor': k,
                  'galois_group': 'S3', 'root_sum': '39/40',
                  'number_energy_2x2_determinant': '81/40',
                  'all_finite_n_claim_requires_proof_not_enumeration': True},
              'native_QQ_graph': native_graph(),
              'original_cold_reference_witness': cold_witness(boxes),
              'finite_resolution_scope': minimum_gap(boxes),
              'alternative_calibration_example': calibration_example(),
              'limits': ['marginal invariance need not preserve correlations with R',
                         'bare energy preserving is not all cyclic work extraction',
                         'all-time stationary is not one-time zero directed flow',
                         'full-rank diagonal product is not correlated interacting Gibbs',
                         'no new physical apparatus, axiom, or geometry identity']}
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    path = HERE/'results.json'
    if args.write:
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    saved = json.loads(path.read_text(encoding='utf-8'))
    assert result == saved, 'Saved certificate differs; no default writes were performed.'
    print(json.dumps({'status': 'AUTHOR_CHECKS_PASS', 'all_checks_passed': True, 'results_match': True,
                      'exact_algebraic_certificate': True,
                      'native_cold_commutator_lower': result['original_cold_reference_witness']['declared_strict_lower'],
                      'occupation_types_checked': 5456}))


if __name__ == '__main__':
    main()
