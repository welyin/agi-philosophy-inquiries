"""Finite algebra and symbol checks only; default is read-only replay."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json

VARS = 'alpha eps w wp f fp q r p dp X dX y K g div nu A TA'.split()
D = len(VARS)
ZERO = (0,) * D

def c(value):
    return {ZERO: F(value)} if value else {}

def v(name):
    e = [0] * D
    e[VARS.index(name)] = 1
    return {tuple(e): F(1)}

def add(*ps):
    out = {}
    for p in ps:
        for key, value in p.items():
            out[key] = out.get(key, F(0)) + value
    return {key: value for key, value in out.items() if value}

def mul(*ps):
    out = c(1)
    for p in ps:
        new = {}
        for a, av in out.items():
            for b, bv in p.items():
                key = tuple(x + y for x, y in zip(a, b))
                new[key] = new.get(key, F(0)) + av * bv
        out = {key: value for key, value in new.items() if value}
    return out

def term(*names, sign=1):
    return mul(c(sign), *(v(name) for name in names))

def algebra():
    jq = add(term('alpha', 'eps', 'wp', 'q', 'p', sign=-1),
             term('g', 'fp', 'X', 'y', sign=-1))
    dq = v('r')
    dr = add(term('q', sign=-1), term('alpha', 'eps', 'w', 'p', sign=-1))
    dy = v('K')
    dK = term('g', 'f', 'X', sign=-1)
    derivative = add(jq, mul(v('q'), dq), mul(v('r'), dr), mul(v('K'), dK),
        term('alpha', 'eps', 'wp', 'q', 'p'),
        mul(term('alpha', 'eps', 'w', 'p'), dq),
        term('alpha', 'eps', 'w', 'q', 'dp'),
        term('g', 'fp', 'X', 'y'), term('g', 'f', 'dX', 'y'),
        mul(term('g', 'f', 'X'), dy))
    target = add(term('alpha', 'eps', 'w', 'q', 'dp'), term('g', 'f', 'y', 'dX'))
    assert derivative == target
    dnuA = add(term('nu', 'div', 'A', sign=-1), term('nu', 'TA'))
    split_rhs = mul(v('nu'), add(v('TA'), term('A', 'div', sign=-1)))
    assert dnuA == split_rhs
    return {'mass_derivative_residual_terms': 0, 'source_transport_residual_terms': 0,
            'independent_formal_symbols': VARS}

def deriv(poly):
    return [F(i) * poly[i] for i in range(1, len(poly))]

def ev(poly, x):
    out = F(0)
    for a in reversed(poly):
        out = out * x + a
    return out

def pulses():
    p = list(map(F, [0, 0, 0, 64, -192, 192, -64]))
    ends = []
    for order in range(4):
        pair = [ev(p, F(0)), ev(p, F(1))]
        ends.append([str(a) for a in pair])
        if order < 3:
            assert pair == [F(0), F(0)]
        else:
            assert pair == [F(384), F(-384)]
        p = deriv(p)
    return {'w_endpoint_derivatives_orders_0_to_3': ends,
            'f_is_translation_of_w': True, 'zero_extension_is_C2_not_C3': True}

def ledger():
    ell = {'Y0': 3, 'Yc': 3, 'Yr': 3, 'Uc': 2, 'Ur': 2,
           'zc': 2, 'zr': 2, 'nc': 1, 'nr': 1}
    rows = {
        'Y0': (1, {'Y0': 2, 'Yc': 1, 'Yr': 1}),
        'Yc': (0, {'Yc': 3, 'Y0': 1, 'Yr': 1, 'Uc': 1, 'zc': 0, 'nc': 0}),
        'Yr': (0, {'Yr': 3, 'Y0': 1, 'Yc': 1, 'Ur': 1, 'zr': 0, 'nr': 0}),
        'Uc': (1, {'Uc': 1, 'Y0': 1, 'Yc': 1, 'Yr': 1, 'zc': 0}),
        'Ur': (1, {'Ur': 1, 'Y0': 1, 'Yc': 1, 'Yr': 1, 'zr': 0}),
        'zc': (1, {'zc': 1, 'Y0': 1, 'Yc': 1, 'Yr': 1, 'Uc': 0}),
        'zr': (1, {'zr': 1, 'Y0': 1, 'Yc': 1, 'Yr': 1, 'Ur': 0}),
        'nc': (0, {'nc': 1, 'Uc': 1, 'Y0': 1, 'Yc': 1, 'Yr': 1}),
        'nr': (0, {'nr': 1, 'Ur': 1, 'Y0': 1, 'Yc': 1, 'Yr': 1}),
    }
    tested = 0
    for name, (index, entries) in rows.items():
        for col, order in entries.items():
            allowed = ell[col] - index
            assert order == allowed if col == name else order < allowed
            tested += 1
    return {'unknown_indices': ell, 'rows': rows,
            'derivative_entries_checked': tested, 'all_off_diagonal_strict': True,
            'absent_entries_are_not_zero_order_dependencies': True}

def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), F(0))

def mv(a, x):
    return [dot(row, x) for row in a]

def symbols():
    xi = list(map(F, [1, 2, 2]))
    k2 = dot(xi, xi)
    assert k2 == 9
    cases = [('equal', [F(1,10), F(1,5), F(0)], [F(1,10), F(1,5), F(0)]),
             ('opposed', [F(1,3), F(0), F(0)], [F(-1,3), F(0), F(0)]),
             ('distinct', [F(0), F(1,4), F(1,5)], [F(1,5), F(-1,4), F(0)])]
    out = []
    for name, vc, vr in cases:
        blocks = []
        for vel in [vc, vr]:
            assert dot(vel, vel) < 1
            a = dot(vel, xi)
            roots = [F(-3), a, F(3)]
            comp = [[F(0), F(1), F(0)], [F(0), F(0), F(1)], [-a*k2, k2, a]]
            for lam in roots:
                eig = [F(1), lam, lam*lam]
                assert mv(comp, eig) == [lam*x for x in eig]
            vandermonde = (roots[1]-roots[0])*(roots[2]-roots[0])*(roots[2]-roots[1])
            assert vandermonde != 0
            gap = min(roots[1]-roots[0], roots[2]-roots[1])
            assert gap > 0
            blocks.append({'v': [str(x) for x in vel], 'roots': [str(x) for x in roots],
                           'vandermonde': str(vandermonde), 'light_transport_gap': str(gap)})
        out.append({'case': name, 'independent_blocks': blocks})
    return {'frequency': [str(x) for x in xi], 'cases': out,
            'no_cross_velocity_gap_assumed': True}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = {'status': 'PASS_FINITE_ALGEBRA_AND_SYMBOL_ONLY',
              'does_not_prove': ['PDE existence', 'continuous dependence', 'geometric uniqueness',
                                'common writing window', 'quantum causality', 'actual menu admission'],
              'algebra': algebra(), 'pulse': pulses(), 'leray_ledger': ledger(), 'symbol_samples': symbols()}
    out = Path(__file__).with_name('results.json')
    canonical = json.loads(json.dumps(result, ensure_ascii=False))
    if args.write:
        with out.open('x', encoding='utf-8', newline='\n') as f:
            json.dump(canonical, f, indent=2, ensure_ascii=False)
            f.write('\n')
    else:
        saved = json.loads(out.read_text(encoding='utf-8'))
        assert saved == canonical, 'saved evidence differs from recomputation'
    print(json.dumps({'status': canonical['status'], 'replay_matches': not args.write,
                      'results_sha256': hashlib.sha256(out.read_bytes()).hexdigest()}, ensure_ascii=False))

if __name__ == '__main__':
    main()
