"""Test the leading-order three-flavour topology inversion at the weak mass-ratio points."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial


HERE = Path(__file__).resolve().parent


def observables(masses):
    weights = 1 / np.array(masses, dtype=float)
    weights /= sum(weights)
    p3, p5, p7 = (sum(weights ** n) for n in (3, 5, 7))
    cumulants = np.array([-p3, 10*p3**2 - 9*p5, -280*p3**3 + 504*p3*p5 - 225*p7])
    # Independently compose the truncated inverse vacuum equation as a polynomial.
    inverse = Polynomial([0, 1, 0, cumulants[0]/6, 0, cumulants[1]/120,
                          0, cumulants[2]/5040])
    theta = inverse + p3*inverse**3/6 + 3*p5*inverse**5/40 + 5*p7*inverse**7/112
    np.testing.assert_allclose(theta.coef[:8], [0,1,0,0,0,0,0,0], atol=2e-14)
    return weights, cumulants


def recover(c4, c6):
    p3 = -c4
    p5 = (10*p3**2 - c6)/9
    d = (1-p3)/3
    if d <= 0:
        raise ValueError("the inverse requires three positive masses")
    e2 = (p5-1+5*d)/(5*d)
    # A companion-matrix solve can turn a triple real root into a complex pair.
    spread = 1/9 - e2/3
    if spread < -2e-15:
        raise ValueError("outside the positive three-mass domain")
    if spread <= 2e-15:
        roots = np.full(3, 1/3)
    else:
        offset = ((e2-d) - e2/3 + 2/27)/2
        scale = spread**1.5
        if abs(offset) > scale + 2e-15:
            raise ValueError("outside the real three-mass domain")
        angle = math.acos(np.clip(offset/scale, -1, 1))/3
        roots = np.array([1/3 + 2*math.sqrt(spread)*math.cos(angle+2*k*math.pi/3)
                          for k in range(3)])
    if min(roots) <= 0 or not np.allclose([sum(roots**3), sum(roots**5)], [p3, p5],
                                          rtol=0, atol=2e-13):
        raise ValueError("outside the positive three-mass domain")
    return np.sort(roots)


def spectra(y, q):
    roots = np.roots([4*q*q-1, -2, 4*y*y-1-4*q*q])
    return [(float(x.real), 1, y) for x in roots if abs(x.imag)<1e-10 and 0<x.real<1]


def solve():
    for masses in ((1,1,1), (.99999999,1,1), (1,1,3)):
        weights, cumulants = observables(masses)
        np.testing.assert_allclose(recover(*cumulants[:2]), np.sort(weights), atol=1e-7)
    family = json.loads((HERE/'results.json').read_bytes())['charged_family']
    result = {}
    for name in ('updated', '2024'):
        y = family[name]['sd']
        cases = {}
        for q in (21.4, 22.1, 22.8):
            masses, = spectra(y, q)
            weights, cumulants = observables(masses)
            np.testing.assert_allclose(recover(*cumulants[:2]), np.sort(weights), atol=2e-13)
            for perm in itertools.permutations(masses):
                np.testing.assert_allclose(observables(perm)[1], cumulants, atol=8e-14)
            cases[str(q)] = dict(mu_over_md=masses[0], C4=cumulants[0], C6=cumulants[1],
                                 C8=cumulants[2], b2=cumulants[0]/12)
        result[name] = dict(ms_over_md=y, Q_cases=cases)

    # Distinct nearby mass shape on the same empirical Q ellipse, with two allowed local terms.
    y = family['updated']['sd']
    _, target = observables(spectra(y, 22.1)[0])
    _, other = observables(spectra(1.001*y, 22.1)[0])
    c4, c6 = target[:2]
    matrix = [[1+c4,16+4*c4], [1-c6,64-4*c6]]
    a,b = np.linalg.solve(matrix,[other[0]-c4, c6-other[1]])
    chi = 1+a+4*b
    new = (other + np.array([-a-16*b, a+64*b, -a-256*b]))/chi
    np.testing.assert_allclose(new[:2], target[:2], atol=2e-14)
    assert chi>0 and abs(new[2]-target[2])>1e-6
    result['higher_order_witness'] = dict(shifted_ms_over_md=1.001*y, a=a, b=b,
                                          relative_chi=chi, held_out_C8_shift=new[2]-target[2])
    return json.loads(json.dumps(result), parse_float=lambda v:float(f'{float(v):.10g}'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    path = HERE/'topology_results.json'
    result = solve()
    if args.check:
        assert result == json.loads(path.read_bytes())
        print("1766 topology: cubic inverse, eighth derivative, higher-order counterexample")
    else:
        with path.open('x', encoding='utf8') as output:
            json.dump(result, output, indent=2, allow_nan=False)
            output.write('\n')
