"""Finite diagnostics for a conditional induced-vacuum admission test.

Analytic inequalities are proved in proof.md; quadrature is not an all-mass proof.
Default: recompute and compare without writes. --write creates results once.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent


def moments(t, order):
    # x = a/s. Factor exp(-t) from both integrals for stable heavy-mass ratios.
    nodes, weights = np.polynomial.legendre.leggauss(order)
    x = (nodes + 1.0) / 2.0
    w = weights / 2.0 * np.exp(-t * (1.0 / x - 1.0))
    return float(w @ x), float(w.sum())


def report():
    rows = []
    for t in (0.0, 0.01, 0.1, 1.0, 4.0, 25.0, 100.0):
        ac, bc = moments(t, 256)
        aa, bb = moments(t, 512)
        ratio = aa / bb
        error = max(abs(aa-ac), abs(bb-bc))
        # Independently obtained integration-by-parts identity after exp(t) scaling.
        recurrence = abs(2.0 * aa + t * bb - 1.0)
        assert error < 1e-9, (t, error)
        assert recurrence < 1e-9, (t, recurrence)
        assert 0.5 - 1e-13 <= ratio < 1.0
        rows.append({'mass_squared_times_a': t, 'scaled_I0_a2': aa,
                     'scaled_I1_a': bb, 'a_I0_over_I1': ratio,
                     'a_R_truncated': -12.0 * ratio,
                     'quadrature_change': error, 'parts_residual': recurrence})
    spectra = [[(0, 1)], [(0, 1000000)], [(0, 2), (2, 7), (4, 3)],
               [(1, 3), (3, 4), (5, 2), (6, 1)]]
    mixtures = []
    for spectrum in spectra:
        a_sum = sum(n * math.exp(-rows[i]['mass_squared_times_a']) * rows[i]['scaled_I0_a2'] for i, n in spectrum)
        b_sum = sum(n * math.exp(-rows[i]['mass_squared_times_a']) * rows[i]['scaled_I1_a'] for i, n in spectrum)
        ratio = a_sum / b_sum
        assert 0.5 - 1e-13 <= ratio < 1.0
        mixtures.append({'spectrum_row_multiplicities': [[i, n] for i, n in spectrum],
                         'a_A_over_B': ratio, 'a_R_truncated': -12.0*ratio})
    eta = F(1, 10)
    required = 1 - eta / 6
    assert required == F(59, 60)
    assert 6 * (1 - F(1, 10)) > eta
    # Constant-curvature Weyl scaling: sqrt(g)->z^4 sqrt(g), R->R/z^2.
    # Derivative of A*z^4 + B*R*z^2/6 at z=1 is 4A + B*R/3.
    for a, n in ((F(1), 1), (F(1, 100), 3), (F(7, 5), 1000000)):
        aa, bb = F(n, 2)/a**2, F(n)/a
        rr = -12*aa/bb
        assert 4*aa + bb*rr/3 == 0
        assert a*rr == -6
    return {'status': 'PASS_FINITE_DIAGNOSTICS',
            'scope': 'Positive minimal scalar, zero independent bare terms, common finite proper-time regulator. No full determinant remainder or Lorentzian observable certificate.',
            'mass_rows': rows, 'mixtures': mixtures,
            'exact_budget': {'eta': str(eta), 'necessary_delta': str(required),
                            'massless_a_R': '-6'},
            'source_of_all_mass_bound': 'Analytic covariance proof, not quadrature.',
            'new_physical_rounds': 0}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = report()
    out = HERE / 'results.json'
    if args.write:
        assert not out.exists(), 'Preserve saved results; use default read-only replay.'
        out.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    else:
        saved = json.loads(out.read_text(encoding='utf-8'))
        assert saved == result
    print(json.dumps({'status': result['status'], 'mass_cases': len(result['mass_rows']),
                      'mixtures': len(result['mixtures']),
                      'required_delta': result['exact_budget']['necessary_delta'],
                      'write': args.write}))


if __name__ == '__main__':
    main()
