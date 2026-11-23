"""Response of full Yukawa matrices and exact neutral finite-order competitor."""
import argparse
import json
import math
from fractions import Fraction as Q
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent
OLD = runpy.run_path(str(HERE.parent / '1678/independent_check.py'))


def at_root(j):
    # 1, z, z² = -1-z, z³=1; work in Q[z]/(z²+z+1).
    return ((1, 0), (0, 1), (-1, -1))[j % 3]


def solve_rational(matrix, target):
    a = [[Q(x) for x in row + [b]] for row, b in zip(matrix, target)]
    n = len(a)
    for i in range(n):
        k = next(j for j in range(i, n) if a[j][i])
        a[k], a[i] = a[i], a[k]
        a[i] = [v / a[i][i] for v in a[i]]
        for j in range(n):
            if j != i:
                m = a[j][i]
                a[j] = [x - m*y for x, y in zip(a[j], a[i])]
    return [row[-1] for row in a]


def neutral_jets(order):
    n = 2 * (order + 1)
    def coefficient(j, r, c):
        return (4*j)**r * at_root(j)[c]
    matrix = [[coefficient(j,r,c) for j in range(n)] for r in range(order+1) for c in (0,1)]
    w = solve_rational(matrix, [2,1] * (order+1))
    assert all(sum(a*b for a,b in zip(row,w)) == q for row,q in zip(matrix,[2,1]*(order+1)))
    r = order + 1
    mismatch = [sum(coefficient(j,r,c)*x for j,x in enumerate(w))-q for c,q in enumerate((2,1))]
    assert w[-1] and mismatch != [0,0]
    return {'degree':n-1, 'coefficients_before_division_by_sqrt3':list(map(str,w)),
            'next_euler_derivative_difference_in_root_basis':list(map(str,mismatch))}


def diagonal_response(y, dy):
    u, v, vh = np.linalg.svd(y)
    p = v*v / (v@v)
    slope = np.diag(u.conj().T @ dy @ vh.conj().T) / v
    return (p-1/3) @ slope


def run():
    def read(path):
        return json.loads((HERE/path).read_bytes())
    witness = read('../1678/independent_results.json')['joint_witness']
    d = read('../1397/inputs.json')['high_scale']
    center, widths = np.array(d['down'] + d['angles_theta12_theta13_theta23_delta'][:3]).T
    l, b, _, _ = OLD['get_matrices'](center+widths*np.array(witness['quark_pulls']),63.2)
    ld = l.copy()
    ld[1,1] /= 3
    yd = ld + b
    t = witness['common_abs_t'] * np.exp(-1j*np.pi/6)
    ye = (l + t*b).T
    wd = np.linalg.norm(b)**2/np.linalg.norm(yd)**2
    we = np.linalg.norm(t*b)**2/np.linalg.norm(ye)**2
    # Both slopes have independent complex common terms; the two complex edges vary freely.
    db = np.zeros_like(b)
    db[1,2], db[2,1] = .39-.13j,-.27+.15j
    output = {}
    for name, eta in [('one_power',1j),('two_paths',4/math.sqrt(3))]:
        dd = (.32+.48j)*yd + db
        de = (-.19+.25j)*ye + (t*db + eta*t*b).T
        c = diagonal_response(ye,de) - (we/wd)*diagonal_response(yd,dd)
        np.testing.assert_allclose(c, (we-1/3)*eta, rtol=2e-11, atol=1e-12)
        output[name] = [round(c.real,10), round(c.imag,10)]
    return {'third_support_weights':[round(we,12),round(wd,12)], 'signed_diagonal_observable':output,
            'neutral_interpolation':{str(n):neutral_jets(n) for n in range(4)}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    result = run()
    path = HERE/'independent_results.json'
    if parser.parse_args().check:
        assert result == json.loads(path.read_bytes())
        print('PASS: six mass-basis vertices and exact first four neutral interpolation budgets')
    else:
        with path.open('x',encoding='utf-8') as f:
            f.write(json.dumps(result,indent=2)+'\n')
