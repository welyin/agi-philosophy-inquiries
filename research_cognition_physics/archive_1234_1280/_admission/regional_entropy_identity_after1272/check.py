"""Same normal pulse: regional relative entropy versus physical energy.

Default execution only compares with results.json. --write creates it once.
Three coordinate integrals test the analytic identities in README.md.
Quadrature refinement is not a certified apparatus or field-equation error.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import math
from fractions import Fraction
import numpy as np

HERE = Path(__file__).resolve().parent
OLD = HERE.parent/'finite_horizon_transfer_after1272/check.py'
spec = importlib.util.spec_from_file_location('old_pulse', OLD)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
H, A, B, L = old.H, old.A, old.B, old.L
PI = math.pi
TIMES = (0.0, 0.5, 1.0, 1.5, 2.0)


def weight(x):
    return 2*PI/H*(math.cosh(H*L)-np.cosh(H*x))/math.sinh(H*L)


def row(t, n):
    lo, hi = A+t, min(B+t, L)
    out = {'time': t}
    keys = ('relative_entropy_x', 'relative_entropy_minkowski',
            'physical_weighted_charge_r', 'positive_volume_correction',
            'physical_ball_energy_r', 'ball_energy_reduced_boundary',
            'outward_power', 'entropy_time_derivative',
            'minkowski_identity_residual', 'weighted_source_identity_residual',
            'energy_identity_residual')
    if hi <= lo:
        out.update({k: 0.0 for k in keys})
        return out
    x, wx = old.nodes(lo, hi, n)
    f, fp = old.profile(x-t)
    sd = 4*PI*float(wx @ (weight(x)*fp**2))
    deriv = -8*PI**2/math.sinh(H*L)*float(wx @ (np.sinh(H*x)*fp**2))
    corr = 4*PI**2*H/math.sinh(H*L)*float(wx @ (f**2/np.cosh(H*x)))

    # Independent conformal Minkowski ball integral, including its Born term.
    ell = math.tanh(H*L/2)/H
    rho, wrho = old.nodes(math.tanh(H*lo/2)/H, math.tanh(H*hi/2)/H, n)
    xm = 2*np.arctanh(H*rho)/H
    fm, fpm = old.profile(xm-t)
    C = np.cosh(H*xm)
    phi_s = -(1+C)*fpm/rho
    phi_rho = (1+C)*fpm/rho-fm/rho**2
    parabolic = (ell**2-rho**2)/ell
    s_flat = float(wrho @ (2*PI**2*parabolic*rho**2*(phi_s**2+phi_rho**2)
                          +4*PI**2/ell*fm**2))

    # Independent original minimal Hilbert source in areal radius.
    r, wr = old.nodes(math.tanh(H*lo)/H, math.tanh(H*hi)/H, n)
    Fr = 1-H**2*r**2
    xr = np.arctanh(H*r)/H
    fr, fpr = old.profile(xr-t)
    phi, phi_t = fr/r, -fpr/r
    phi_r = fpr/(r*Fr)-fr/r**2
    density = 2*PI*(r**2/Fr*phi_t**2+r**2*Fr*phi_r**2+2*H**2*r**2*phi**2)
    q_weight = float(wr @ (weight(xr)*density))
    e_ball = float(wr @ density)
    R = math.tanh(H*L)/H
    FR = 1-H**2*R**2
    fb, fpb = (float(v) for v in old.profile(np.array(L-t)))
    e_boundary = 4*PI*float(wx @ fp**2)-2*PI*FR/R*fb**2
    power = 4*PI*(fpb**2-FR/R*fb*fpb)
    out.update(dict(zip(keys, (sd, s_flat, q_weight, corr, e_ball, e_boundary,
                    power, deriv, abs(sd-s_flat), abs(sd-q_weight-corr),
                    abs(e_ball-e_boundary)))))
    return out


def calculate():
    runs = [{'nodes': n, 'rows': [row(t, n) for t in TIMES]} for n in (96, 192, 384)]
    end = runs[-1]['rows']
    x, wx = old.nodes(A, B, 384)
    _, fp = old.profile(x)
    Q = 4*PI*float(wx @ fp**2)
    whole = 2*PI/H*Q
    residual_keys = ('minkowski_identity_residual', 'weighted_source_identity_residual',
                     'energy_identity_residual')
    max_residual = max(r[k] for r in end for k in residual_keys)
    refinement = max(abs(r[k]-s[k]) for r, s in zip(end, runs[-2]['rows'])
                     for k in r if k != 'time')
    gap = end[0]['relative_entropy_x']-end[1]['relative_entropy_x']
    analytic_gap = 20*PI**2*math.exp(-2/3)/math.sinh(3)
    rational_gap = Fraction(40, 9)
    assert max_residual < 1e-8
    assert refinement < 1e-8
    assert gap > analytic_gap > float(rational_gap)
    for r in end[:3]:
        assert abs(r['physical_ball_energy_r']-Q) < 1e-8
        assert r['outward_power'] == 0
        assert r['entropy_time_derivative'] < 0
    assert all(end[i]['relative_entropy_x'] > end[i+1]['relative_entropy_x']
               for i in range(len(end)-1))
    assert all(0 <= r['relative_entropy_x'] < whole for r in end)
    assert end[-1]['relative_entropy_x'] == 0
    assert end[3]['positive_volume_correction'] > 0
    return {'status': 'PASS_REGIONAL_IDENTITY_SCOPE',
        'parameters': {'H': H, 'minimal_mass_squared': 2*H**2, 'xi': 0,
                       'support': [A, B], 'L': L},
        'whole_patch_energy': Q, 'whole_patch_relative_entropy_constant': whole,
        'integrals': runs, 'max_independent_identity_residual': max_residual,
        'max_last_refinement_difference': refinement,
        'early_entropy_drop': gap, 'analytic_drop_lower_bound': analytic_gap,
        'strict_rational_drop_lower_bound': str(rational_gap),
        'scope': {'actual_readout_certified': False,
                  'coupled_quantum_geometry_certified': False,
                  'geometric_entropy_identified': False,
                  'governance_conjecture_refuted': False,
                  'new_numbered_rounds': 0, 'new_axioms': 0}}


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for k in a: compare(a[k], b[k])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b): compare(x, y)
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=2e-11, abs_tol=2e-10), (a, b)
    else:
        assert a == b, (a, b)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = calculate()
    dest = HERE/'results.json'
    if args.write:
        with dest.open('x', encoding='utf-8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
        mode = 'created_once'
    else:
        compare(result, json.loads(dest.read_text(encoding='utf-8')))
        mode = 'read_only_match'
    print(json.dumps({'status': result['status'], 'mode': mode,
        'early_entropy_drop': result['early_entropy_drop'],
        'max_identity_residual': result['max_independent_identity_residual'],
        'rows': result['integrals'][-1]['rows']}, ensure_ascii=False))


if __name__ == '__main__':
    main()
