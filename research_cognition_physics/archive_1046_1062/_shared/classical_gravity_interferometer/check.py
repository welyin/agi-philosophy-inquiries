"""Algebra and rounded-parameter diagnostics, not an experimental reanalysis.

Default is read-only. --write creates results.json exclusively on first save.
No data fit, visibility confidence bound, pulse solver, or new scientific group.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import math


def compute():
    # AM-GM identity for the adopted symmetric pair generator.
    h, k, d = F(7, 5), F(-9, 4), F(2, 3)
    q = 1/(4*d) + k*k*d/(4*h*h)
    excess = (h-abs(k)*d)**2/(4*d*h*h)
    assert q-abs(k)/(2*h) == excess >= 0
    dopt = h/abs(k)
    assert 1/(4*dopt)+k*k*dopt/(4*h*h) == abs(k)/(2*h)

    # Existing 1033 asymmetric boundary: not a universal single-end minimum.
    g, a = F(1), F(1, 100)
    gamma_a, gamma_b = g*a, g/a
    assert gamma_a*gamma_b == g*g
    assert gamma_a < g < gamma_b

    # Rounded model inputs, not exact measured parameters.
    G, hbar = F('6.67430e-11'), F('1.054571817e-34')
    mass_earth, radius, mass = F('6e24'), F('6e6'), F('1.4e-25')
    recoil, T, N = F('0.0058'), F('1.04'), 45
    v = 2*N*recoil
    delta_max = v*T
    area = F(2, 3)*v*v*T**3
    assert area == F(2, 3)*delta_max**2*T
    assert delta_max == F('0.54288')
    # Independent Simpson integral is exact for each quadratic half-path.
    simpson = 2*(T/6)*(0 + 4*(v*T/2)**2 + (v*T)**2)
    assert simpson == area

    rows = []
    for C in [F('1'), F('0.47'), F('0.1')]:
        D = C*G*mass_earth*mass/(hbar*radius**3)
        exponent = D*area
        # The printed Eq.(7) form, eliminating k via recoil=hbar*k/m.
        wave_number = recoil*mass/hbar
        alternate = F(2, 3)*C*G*hbar*mass_earth/(mass*radius**3)*(2*N*wave_number)**2*T**3
        assert exponent == alternate > 0
        rows.append({
            'C': str(C),
            'D_per_m2_s': format(float(D), '.15g'),
            'minus_log_model_visibility': format(float(exponent), '.15g'),
            'model_visibility_diagnostic': format(math.exp(-float(exponent)), '.15g'),
        })
    return {
        'scope': 'symmetric_pairwise_KTM_only_algebra_and_rounded_inputs_no_data_fit',
        'amgm_identity_residual': '0',
        'old_1033_asymmetric_boundary': {'gamma_A': str(gamma_a), 'gamma_B': str(gamma_b), 'product': str(g*g)},
        'half_time_s': str(T), 'full_sequence_time_s': str(2*T),
        'momentum_separation_hbar_k': 2*N,
        'maximum_path_separation_m': str(delta_max),
        'integral_path_separation_squared_m2_s': str(area),
        'two_exponent_formula_residual': '0',
        'rows': rows,
        'no_experimental_uncertainty_or_sigma_recomputed': True,
        'all_assertions_passed': True,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = compute()
    path = Path(__file__).with_name('results.json')
    if args.write:
        with path.open('x', encoding='utf-8') as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
        print('First-save algebra/rounded-parameter diagnostics passed.')
    else:
        assert json.loads(path.read_text(encoding='utf-8')) == result
        print('Read-only checks passed; no experimental significance or general LOCC exclusion claimed.')
