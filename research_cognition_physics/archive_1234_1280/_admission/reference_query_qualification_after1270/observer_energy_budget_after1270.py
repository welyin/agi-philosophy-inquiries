"""Exact finite-prefix budget bounds for the stated quadratic product reference.

No field integration, truncation approximation or physical source certification.
Default: read-only replay. --write: create the first result exclusively.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

SCRIPT = Path(__file__).resolve()
HERE = SCRIPT.parent
ROOT = next(p for p in (HERE, *HERE.parents)
            if (p / 'research_cognition_physics').is_dir()
            and (p / 'scripts').is_dir())
MEMO = HERE / 'observer_physical_window_after1270.md'
RESULT = SCRIPT.with_suffix('.json')
OLD = ROOT / 'research_cognition_physics/archive_1234_/_admission/bounded_reference_query_after1270/observer_finite_query_gate_after1270.md'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational(value):
    return {'numerator': str(value.numerator), 'denominator': str(value.denominator)}


def output():
    records = []
    bits = 64
    for energy in (1, 10, 100):
        ell = 1
        while (ell * (ell + 2)) // 2 < 4 * energy + bits:
            ell += 1
        modes = ell * (ell + 2)
        assert modes == sum(2 * j + 1 for j in range(1, ell + 1))
        assert ((ell - 1) * (ell + 1)) // 2 < 4 * energy + bits
        # t=2: P(H_BE,N <= E) <= exp(2E) * 2**(-N/2).
        # E>0 and e<3 give a strict rational majorant; flooring N/2 is conservative.
        upper = Fraction(3 ** (2 * energy), 2 ** (modes // 2))
        target = Fraction(1, 2 ** bits)
        assert upper < target
        records.append({'energy_budget_dimensionless': energy,
                        'complete_multiplets_through_ell': ell,
                        'number_of_real_charge_modes': modes,
                        'strict_probability_upper': rational(upper),
                        'requested_upper': rational(target),
                        'mean_binding_energy_dimensionless': rational(Fraction(modes, 4)),
                        'strictly_below_requested_upper': True})
    paths = {'proof': MEMO, 'code': SCRIPT, 'previous_reference_query': OLD}
    return {'status': 'PASS_FINITE_PREFIX_BOUND_ONLY',
            'model_scope': {
                'joint_reference': 'independent Gaussian charge spectra from neutral Lebesgue spectral weights',
                'mode_energy': 'H_j=c_j*q_j^2, c_j>0; reference proportional to exp(-H_j)',
                'binding_energy': 'H_BE,N=(1/2)*sum_j H_j in fixed dimensionless units',
                'multiplets': 'all m for each ell=1,...,L; rotationally invariant coefficients',
                'physical_cutoff_or_minimum_length_assumed': False,
                'finite_effective_models_ruled_out': False,
                'paper_or_entire_programme_refuted': False,
                'complete_source_certified': False,
                'new_scientific_rounds': 0, 'new_science_groups': 0},
            'analytic_identity_used': 'E[exp(-t*H_BE,N)]=(1+t/2)^(-N/2), t>0',
            'analytic_consequence': 'For each finite E, probability of the infinite positive sum <= E is zero in this normal product-reference sector.',
            'cases': records,
            'source_bindings': {k: {'path': str(v.relative_to(ROOT)).replace('\\', '/'),
                                    'sha256': digest(v)} for k, v in paths.items()}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = output()
    if args.write:
        with RESULT.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, sort_keys=True)
            stream.write('\n')
        print('FIRST_WRITE_FINITE_PREFIX_BOUND')
    else:
        assert json.loads(RESULT.read_text('utf-8')) == result, 'Saved result differs.'
        print('READ_ONLY_MATCH_FINITE_PREFIX_BOUND')
    for case in result['cases']:
        print('E=', case['energy_budget_dimensionless'],
              'L=', case['complete_multiplets_through_ell'],
              'N=', case['number_of_real_charge_modes'], ': probability < 2^-64')
    print('No full-source or actual readout certificate; no field integration.')


if __name__ == '__main__':
    main()
