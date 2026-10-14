"""Read-only, exact-arithmetic screening; --write creates results.json once.

This substitutes a known Einstein-aether family into published low-energy
formulae. It does not derive those field equations or construct a quantum model.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json


def bracket_root(value, degree, bits=100):
    lo, hi = F(0), F(2)
    while hi ** degree < value:
        hi *= 2
    for _ in range(bits):
        mid = (lo + hi) / 2
        if mid ** degree <= value:
            lo = mid
        else:
            hi = mid
    assert lo ** degree <= value <= hi ** degree
    return lo, hi


def published_formulae(t):
    assert 0 < t < F(1, 2)
    c1, c2, c3, c4 = t, t / (1 - 2*t), -t, F(0)
    c13, c14, c123 = c1 + c3, c1 + c4, c1 + c2 + c3
    numerator1 = 2*c1 - c1*c1 + c3*c3
    speeds = [1/(1-c13), numerator1/(2*c14*(1-c13)),
              c123*(2-c14)/(c14*(1-c13)*(2+c13+3*c2))]
    alpha1 = -8*(c3*c3+c1*c4)/numerator1
    alpha2 = alpha1/2 - ((c1+2*c3-c4)*(2*c1+3*c2+c3+c4)
                         / (c123*(2-c14)))
    gn_over_g = 1/(1-c14/2)
    gc_over_g = 1/(1+(c13+3*c2)/2)
    ratio = gc_over_g/gn_over_g
    signs = [numerator1/(1-c13), c14*(2-c14)]
    assert speeds == [1, 1, 1]
    assert alpha1 == -4*t and alpha2 == 0
    assert ratio == 1-2*t
    assert signs == [2*t, t*(2-t)] and min(signs) > 0
    return dict(t=str(t), ci=list(map(str, (c1, c2, c3, c4))),
                speed_squared=list(map(str, speeds)),
                vector_scalar_energy_sign_factors=list(map(str, signs)),
                alpha1=str(alpha1), alpha2=str(alpha2),
                GN_over_G=str(gn_over_g), Gcosmo_over_G=str(gc_over_g),
                Gcosmo_over_GN=str(ratio))


def run():
    # These exact witnesses check implementation. The all-t identities are
    # proved algebraically in README.md, not inferred from sampling.
    rows = [published_formulae(t) for t in (F(1,1000), F(1,10), F(2,5))]
    # Common GN, rho0, a(0)=1; x=sqrt(8*pi*GN*rho0/3)*proper_time.
    # Compare at x=2/3, so a_GR^3=4, a_t^3=(1+sqrt(4/5))^2.
    ratio = F(4,5)
    slo, shi = bracket_root(ratio, 2)
    ag_lo, ag_hi = bracket_root(F(4), 3)
    at_lo, _ = bracket_root((1+slo)**2, 3)
    _, at_hi = bracket_root((1+shi)**2, 3)
    delta_lo, delta_hi = ag_lo-at_hi, ag_hi-at_lo
    assert delta_lo > F(5,100)
    assert delta_hi-delta_lo < F(1,10**25)
    # First integral and continuity imply a''=-ratio/(2*a^2)<0.
    # Algebraic identity for f=a^(3/2): f'=3*sqrt(ratio)/2;
    # then (a')^2=ratio/a and rho/rho0=a^-3.
    assert F(9,4)*ratio == F(9,5)
    return {
        'schema': 'reference_selection_screen_after1198_v1',
        'status': 'PASS_SCOPED_REUSE',
        'published_formula_substitution': rows,
        'finite_classical_FLRW_comparison': {
            't': '1/10', 'x': '2/3', 'ratio': '4/5',
            'a_GR_interval': list(map(str, (ag_lo,ag_hi))),
            'a_aether_interval': list(map(str, (at_lo,at_hi))),
            'scale_difference_interval': list(map(str, (delta_lo,delta_hi))),
            'difference_decimal_diagnostic': float((delta_lo+delta_hi)/2),
            'certified_difference_greater_than': '1/20',
            'both_solutions_decelerate': True,
            'common_initial_data': ['a0=1', 'positive dust density rho0', 'measured GN'],
            'not_common_initial_data': ['H0, which is fixed by each model constraint'],
        },
        'scope': ('Known classical Einstein-aether formulae, flat-constant-vector '
                  'linear modes and separate flat-FLRW dust solutions. No all-background '
                  'stability, full quantum process, actual instruments, cognitive-axiom '
                  'countermodel or dark-sector identification.'),
        'new_numbered_rounds': 0,
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    rendered = json.dumps(result, ensure_ascii=False, indent=2) + '\n'
    if args.write:
        with Path(__file__).with_name('results.json').open('x', encoding='utf-8') as out:
            out.write(rendered)
    print(rendered, end='')
