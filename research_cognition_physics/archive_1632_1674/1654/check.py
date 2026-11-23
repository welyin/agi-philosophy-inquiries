"""Extract scalar cubic interactions from the defining potential as polynomials."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial

HERE = Path(__file__).resolve().parent


def evaluate():
    cfg = json.loads((HERE / 'inputs.json').read_bytes())
    ref = json.loads((HERE / cfg['previous']).read_bytes())
    c, v, m, M = [cfg[k] for k in ('c', 'v_GeV', 'mh_GeV', 'M_GeV')]
    s = math.sqrt(1 - c*c)
    O = np.array([[c, s], [-s, c]])
    mass = O.T @ np.diag([m*m, M*M]) @ O
    scale = M*M + 2*m*m
    xz = (cfg['mZ_GeV']/M)**2
    phase_Z = math.sqrt(1 - 4*xz) * (1 - 4*xz + 12*xz*xz)
    phase_h = math.sqrt(1 - (2*m/M)**2)
    actual, point = [], cfg['xi_trial']

    def couplings(f, deformations=(0., 0.)):
        quartic = mass / np.outer([v, f], [v, f]) / 2
        assert np.min(np.linalg.eigvalsh(quartic)) > 0

        def cubic_along(a, b):
            x, u = O.T @ [a, b]
            X, L = Polynomial([0, v*x, x*x/2]), Polynomial([0, f*u, u*u/2])
            potential = quartic[0, 0]*X**2 + 2*quartic[0, 1]*X*L + quartic[1, 1]*L**2
            potential += deformations[0]*(X/v)**3 + deformations[1]*(X/v)**2*(L/f)
            np.testing.assert_allclose(2*potential.coef[2], m*m*a*a + M*M*b*b)
            return 6*potential.coef[3]

        hhh = cubic_along(1, 0)
        Hhh = (cubic_along(1, 1) - cubic_along(1, -1) - 2*cubic_along(0, 1))/6
        return np.array([hhh, Hhh])

    def report(xi):
        f = v*s/xi
        hhh, Hhh = couplings(f)
        kappa = v*hhh/(3*m*m)
        A = -v*Hhh/(s*scale)
        ratio = Hhh**2*v*v*phase_h/(s*s*M**4*phase_Z)
        unsigned = math.sqrt(ratio*phase_Z/phase_h)*M*M/scale
        np.testing.assert_allclose([A, kappa, unsigned], [c*(c-xi), c**3+s*s*xi, abs(A)])
        np.testing.assert_allclose((kappa-c**3)/(s*s), xi)
        np.testing.assert_allclose(kappa, c - s*s/c*A)
        return {'xi': xi, 'f_GeV': f, 'A': A, 'kappa_lambda': kappa,
                'Gamma_Hhh_over_HZZ': ratio, 'flavour_z1_and_z2': cfg['N']*xi,
                'inclusive_FC_over_SMbb_lower': ref['FC_lower_coefficient_on_equal_contrasts']*(cfg['N']*xi)**2}

    for xi in [point, 2*c-point]:
        actual.append(report(xi))
    np.testing.assert_allclose(actual[0]['Gamma_Hhh_over_HZZ'], actual[1]['Gamma_Hhh_over_HZZ'])
    assert abs(actual[0]['kappa_lambda'] - actual[1]['kappa_lambda']) > 1e-3
    f1, f2 = [item['f_GeV'] for item in actual]
    base, target = couplings(f2), couplings(f1)
    response = np.column_stack([couplings(f2, d) - base for d in [(1, 0), (0, 1)]])
    np.testing.assert_allclose(np.linalg.det(response), 12*c**4)
    correction = np.linalg.solve(response, target-base)
    np.testing.assert_allclose(couplings(f2, correction), target, atol=1e-10)
    return {'two_width_solutions': actual, 'cubic_response_determinant': np.linalg.det(response),
            'EFT_counterterms_GeV_matching_both_signed_cubics': correction.tolist(),
            'absolute_error_gain_xi_from_kappa_lambda_at_fixed_c': 1/(1-c*c)}


def rounded(item):
    if isinstance(item, dict):
        return {key: rounded(val) for key, val in item.items()}
    if isinstance(item, list):
        return [rounded(val) for val in item]
    if isinstance(item, (float, np.floating)):
        return float(f'{item:.10g}')
    return item


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = rounded(evaluate())
    output = HERE / 'results.json'
    if args.check:
        assert json.loads(output.read_bytes()) == result
        print('PASS: both scalar widths, signed cubics, two isolated solutions, EFT competitor')
    else:
        with output.open('x', encoding='utf-8') as stream:
            stream.write(json.dumps(result, indent=2) + '\n')
