"""Rational bounds for a proposed Hartree response screen, not an acceptance.

This checks constants conditional on the analytic response bounds in README.
It does not certify the full finite-K spectrum, a detector, or a GR no-go.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json


def run():
    eps, nu, g, om, m, ell = map(F, ('3/20', '7/20', '1/5', '1/2', '2', '16/65'))
    zlo, freq = F(49, 50), F(1, 4)
    sqrt3_up, xs_up = F(17321, 10000), F(199, 1000)
    assert sqrt3_up**2 > 3 and xs_up**2 > 1-zlo*zlo
    nlo, vlo = 1+eps*zlo, 1+nu*zlo
    bs = 2*eps**2*m/nlo+2*nu**2*sqrt3_up/(vlo*ell)
    bt = 2*g**2*sqrt3_up/(vlo*ell)
    cs, ct = xs_up**3/om, 1/(om*nlo)
    static = cs*bs+ct*bt
    rs2lo, rt2lo = om**2/(1-zlo*zlo), (om*nlo)**2
    gaplo = 2*m*nlo
    finite_frequency = (cs*bs/(1-freq**2/(4*rs2lo))
                        +ct*bt/(1-freq**2/(4*rt2lo)))/(1-freq/gaplo)
    assert static < F(753, 1000)
    assert finite_frequency < F(835, 1000)
    def out(x):
        return {'exact': str(x), 'decimal_diagnostic': float(x)}
    return {
        'status': 'PASS_CONDITIONAL_ARITHMETIC_ONLY',
        'accepted_round': None,
        'reference_condition': 'stationary isotropic Hartree filled band; z_s >= 49/50',
        'fixed_lattice_spacing': str(ell),
        'Bs_upper': out(bs), 'Bt_upper': out(bt),
        'static_feedback_norm_upper': out(static),
        'frequency_interval': ['-1/4', '1/4'],
        'source_pair_gap_lower': out(gaplo),
        'low_frequency_feedback_norm_upper': out(finite_frequency),
        'conditional_inverse_margin_lower': out(1-finite_frequency),
        'analytic_spectral_chain_independently_reviewed': False,
        'full_finite_K_gap_certified': False,
        'dark_sector_identity_or_GR_exclusion': False,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    result['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    output = json.dumps(result, ensure_ascii=False, indent=2)+'\n'
    if args.write:
        with Path(__file__).with_name('results.json').open('x', encoding='utf-8') as stream:
            stream.write(output)
    print(output, end='')
