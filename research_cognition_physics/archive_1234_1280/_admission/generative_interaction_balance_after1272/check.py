"""Read-only replay of bounded perturbative diagnostics; --write creates once.

Independent proper-time and momentum radial integrals verify the same tadpole.
No all-orders remainder, finite-coupling solution or physical apparatus simulated.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent


def integrals(t, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    x, w = (nodes + 1) / 2, weights / 2
    proper = float(w @ np.exp(-t / x))
    z = x / (1 - x)
    radial = float(w @ (np.exp(-t-z) * z / (z+t) / (1-x)**2))
    return proper, radial


def pairings(items):
    if not items:
        return 1
    return sum(pairings(items[1:j] + items[j+1:]) for j in range(1, len(items)))


def report():
    wick = F(pairings(tuple(range(4))), math.factorial(4))
    assert wick == F(1, 8)
    # Keep I1 and pi symbolic in these rational coefficient checks.
    v_factor = wick / 16**2
    k_factor = F(1, 32)
    slope_factor = 2 * v_factor / k_factor
    assert v_factor == F(1, 2048) and slope_factor == F(1, 32)
    # Independent constant metric scaling: W(z)=V*z^4-K*R*z^2.
    # At lambda order: 4v1-2K0 R1=0, regardless of K'(0).
    for i1 in (F(1, 7), F(2, 3), F(4, 5)):
        v1 = v_factor * i1**2
        k0 = k_factor * i1
        r1 = slope_factor * i1
        assert 4*v1 - 2*k0*r1 == 0
    rows = []
    for t in (0.01, 0.1, 1.0, 4.0, 25.0):
        p0, r0 = integrals(t, 256)
        p, r = integrals(t, 512)
        relerr = abs(p-r) / p
        drift = max(abs(p-p0), abs(r-r0)) / p
        assert 0 < p < 1
        assert relerr < 1e-9 and drift < 1e-8, (t, relerr, drift)
        slope = p/(32*math.pi**2)
        assert 0 < slope < 1/(32*math.pi**2)
        rows.append({'mass_squared_times_a': t,
                     'a_I1_proper_time': p, 'a_I1_momentum_radial': r,
                     'integral_relative_difference': relerr,
                     'quadrature_relative_change': drift,
                     'dimensionless_vacuum_derivative': p*p/(2048*math.pi**4),
                     'truncated_dimensionless_curvature_slope': slope})
    return {'status': 'PASS_INTERACTION_BALANCE_DIAGNOSTICS',
            'scope': 'Specified finite per-line regulator; first perturbative coefficient and local 0+2 derivative formal curvature tangent only.',
            'exact_factors': {'wick': str(wick), 'v1_I1_squared_pi_minus4': str(v_factor),
                             'Rprime_I1_pi_minus2': str(slope_factor)},
            'mass_rows': rows,
            'curvature_slope_upper_limit': 1/(32*math.pi**2),
            'full_model_or_remainder_certified': False,
            'new_scientific_rounds': 0}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = report()
    path = HERE / 'results.json'
    if args.write:
        assert not path.exists(), 'Use default replay; preserve saved results.'
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    else:
        assert json.loads(path.read_text(encoding='utf-8')) == result
    print(json.dumps({'status': result['status'], 'mass_cases': len(result['mass_rows']),
                      'wick': result['exact_factors']['wick'], 'write': args.write}))


if __name__ == '__main__':
    main()
