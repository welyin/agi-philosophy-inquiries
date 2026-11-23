import argparse
import importlib.util
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial

HERE = Path(__file__).resolve().parent


def read(rel):
    return json.loads((HERE / rel).read_bytes())


def make_texture(data):
    # Apply the same mass-basis reconstruction as 1489 at the chosen scale.
    spec = importlib.util.spec_from_file_location('texture1489', HERE.parent/'1489/check.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    angles = data['angles_theta12_theta13_theta23_delta']
    sample = {'yukawa': {'down': data['down']},
              'quark_sines_and_phase': [[math.sin(x[0])] for x in angles[:3]] + [angles[-1]]}
    mat = module.construction(sample)
    return [*mat['normalized_real_edges_A_B_C_E_F'], *mat['Z_real_imag']]


def matrix(edges, factors, t):
    a, b, c, e, f, zr, zi = edges
    ca, cb, cc = factors
    return np.array([[0, cb*a, 0], [ca*b, cc*c, t*complex(zr, zi)], [0, t*e, t*f]])


def shape_polynomials(edges, factors):
    a, b, c, e, f, zr, zi = edges
    ca, cb, cc = factors
    X, Y, W = cb*a, ca*b, cc*c
    zero = Fraction(0)
    T = Polynomial([X**2+Y**2+W**2, zero, zr**2+zi**2+e**2+f**2])
    U = Polynomial([(X*Y)**2, zero,
                    X**2*(zr**2+zi**2+f**2)+Y**2*(e**2+f**2)+(W*f)**2,
                    -2*W*f*e*zr, e**2*(zr**2+zi**2)])
    V = Polynomial([zero, zero, (X*Y*f)**2])
    return T, U, V


def elimination(T, U, V, k):
    # q is the sum of a squared mass ratio and its reciprocal.
    q = k + 1/k
    return (V*V*q**3 - (T*U-3*V)*V*q**2 + (T*T*T*V+U*U*U-5*T*U*V+3*V*V)*q
            - (T*T-2*U)*(U*U-2*T*V) + V*V)


def exact_real_roots(coeff):
    # Descartes after t = sign*u/(1-u). Dyadic Bernstein subdivision
    # uses common-denominator integer arithmetic on both real half axes.
    assert coeff[0] and coeff[-1] and all(isinstance(x, Fraction) for x in coeff)
    degree = len(coeff)-1
    intervals = []
    for sign in (-1, 1):
        bern = [x*sign**j/math.comb(degree, j) for j, x in enumerate(coeff)]
        den = math.lcm(*(x.denominator for x in bern))
        initial = [int(x*den) for x in bern]
        queue = [(0, 1, initial)]
        while queue:
            i, n, values = queue.pop()
            nonzero = [v for v in values if v]
            variations = sum((a>0) != (b>0) for a, b in zip(nonzero, nonzero[1:]))
            if not variations:
                continue
            if variations == 1 and i+1<n:
                endpoints = sorted([Fraction(sign*i, n-i), Fraction(sign*(i+1), n-i-1)])
                if endpoints[1]-endpoints[0] < Fraction(1, 10**14):
                    intervals.append(endpoints)
                    continue
            assert n < 2**90
            row, left, right = values, [], []
            for j in range(degree+1):
                left.append(row[0] << (degree-j))
                right.append(row[-1] << (degree-j))
                row = [a+b for a, b in zip(row, row[1:])]
            assert left[-1] != 0, 'exact dyadic root needs separate endpoint handling'
            queue.extend([(2*i, 2*n, left), (2*i+1, 2*n, right[::-1])])
    return sorted(intervals)


def classify(edges, factors, intervals, leptons):
    electron, mu, tau = leptons
    ratio = tau/mu
    solutions = []
    tagged = {'tau_mu': 0, 'mu_e': 0, 'tau_e': 0}
    for low, high in intervals:
        t = float((low+high)/2)
        raw = matrix(edges, factors, t)
        s = np.linalg.svd(raw, compute_uv=False)[::-1]
        candidates = np.array([s[2]/s[1], s[1]/s[0], s[2]/s[0]])
        index = int(np.argmin(abs(candidates/ratio-1)))
        assert abs(candidates[index]/ratio-1)<1e-6
        assert sorted(abs(candidates/ratio-1))[1]>.01
        tag = list(tagged)[index]
        tagged[tag] += 1
        if index == 0:
            prediction = tau*s[0]/s[-1]
            solutions.append({'t': t, 'k_times_yb': tau/s[-1],
                              'predicted_ye': prediction, 'predicted_over_measured':prediction/electron})
    assert len(intervals) == sum(tagged.values())
    return {'isolated_ratio_roots_by_pair':tagged, 'solutions_fitting_mu_tau':solutions}


def compute():
    cfg = read('inputs.json')
    data = read(cfg['quarks'])['high_scale']
    inputs = read(cfg['masses'])['scales'][cfg['scale']]
    leps = [x[0] for x in inputs[3:]]
    q = Fraction(str(leps[2]))/Fraction(str(leps[1]))
    raw = make_texture(data)
    vals = [Fraction(f'{a:.{cfg["round_digits"]}f}') for a in raw]
    edges = list(map(float, vals))
    assert max(abs(a-b) for a,b in zip(raw,edges)) <= 5.1e-13
    models = {}
    for label, row in read(cfg['clebsches'])['clebsch'].items():
        factors = [Fraction(str(x)) for x in row]
        T, U, V = shape_polynomials(vals, factors)
        elim = elimination(T, U, V, q*q)
        for probe in (.2, -1, 4):
            singular = np.linalg.svd(matrix(edges,row,probe), compute_uv=False)
            sq = singular**2
            np.testing.assert_allclose([float(T(probe)), float(U(probe)), float(V(probe))],
                                       [sq.sum(), sum(sq[i]*sq[j] for i,j in [(0,1),(1,2),(2,0)]), sq.prod()],
                                       rtol=2e-12, atol=1e-23)
        intervals = exact_real_roots(elim.coef)
        models[label] = classify(edges, row, intervals, leps)
    return {'edges_A_B_C_E_F_ReZ_ImZ':edges, 'lepton_reference':leps, 'models':models}


def tidy(data):
    if isinstance(data, dict):
        return {k:tidy(v) for k,v in data.items()}
    if isinstance(data, list):
        return [tidy(v) for v in data]
    if isinstance(data, (float, np.floating)):
        return float(f'{data:.10g}')
    return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args = parser.parse_args()
    result = tidy(compute())
    if args.check:
        assert read('results.json') == result
        print('PASS: full signed third-family response, complete roots, mu/tau calibration, electron holdout')
    else:
        with (HERE / 'results.json').open('x',encoding='utf-8') as out:
            out.write(json.dumps(result,indent=2)+'\n')
