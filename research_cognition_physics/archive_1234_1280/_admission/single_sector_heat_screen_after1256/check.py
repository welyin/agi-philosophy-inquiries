"""Read-only rational consequence of the frozen 1250 source and null tube.

No new evolution, quantum reference, instrument, or numbered science group.
Use --write-results only for the first creation of this screen's result.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
STAGE = HERE.parents[1]
BINDINGS = {
    '1247/proof.md': 'f3a79f8f43aa235d2fa13cc87a232e2b64643fc1bade375bfee96c204ab4bf3d',
    '1250/proof.md': '5e9f5e1495156f6bcd81c48a5046ac538ea72b2823a6cb7a9e3b7a59ff499d9f',
    '1250/results.json': '84fa683b1e343f94ba2f5a4f2e26653dfa49c7d6c04d3833a436daa8cb880284',
}


def compute():
    for name, digest in BINDINGS.items():
        assert hashlib.sha256((STAGE / name).read_bytes()).hexdigest() == digest, name
    original = json.loads((STAGE / '1250/results.json').read_text('utf-8-sig'))
    assert original['status'] == 'PASS'
    assert original['fixed']['time_width'] == '1/64'
    eta = Q(original['fixed']['error_per_menu'])
    branches = []
    for branch in original['group_same_history_null']['branches']:
        tube = dict(zip(branch['null_names'], branch['null_tube']['all_time_tube']))
        u, p, V = [[Q(x) for x in tube[n]] for n in ['u', 'p', 'V']]
        mc = Q(branch['frozen_t2']['M_c'][0])
        fmin = 1 + p[0] / 5
        dust_min = mc + Q(1, 50) + p[0] / 100
        assert u[0] > 0 and fmin > 0 and dust_min > 0
        x_max = max(V[0] ** 2, V[1] ** 2)
        other_min = fmin * u[0] ** 2 / 10 + dust_min * u[0] ** 3
        # Dropping v^2 only weakens the lower bound on the non-X source.
        ratio = x_max / (x_max + other_min)
        assert 0 < ratio < Q(1, 25)
        boost_lower = Q(branch['screen']['boost_full'][0])
        gap_lower = (1 - ratio) * boost_lower
        memory_lower = Q(branch['screen']['positive_memory_gap'][0])
        assert gap_lower > Q(19, 10**6) > 7 * eta
        assert gap_lower > memory_lower > 7 * eta
        branches.append({
            'alpha': branch['alpha'],
            'X_null_source_upper': str(x_max),
            'non_X_null_source_lower_without_v_squared': str(other_min),
            'X_fraction_upper': str(ratio),
            'absolute_X_only_heat_omission_lower': str(gap_lower),
            'old_memory_omission_lower': str(memory_lower),
            'conditional_margin_over_old_seven_eta': str(gap_lower - 7 * eta),
            'display_X_fraction_upper': float(ratio),
            'display_heat_omission_lower': float(gap_lower),
        })
    return {
        'status': 'PASS_ABSOLUTE_SOURCE_SCOPE_ONLY',
        'source_sha256': BINDINGS,
        'branches': branches,
        'uniform_X_fraction_strict_upper': '1/25',
        'uniform_absolute_heat_omission_strict_lower': '19/1000000',
        'new_scientific_rounds': 0,
        'new_science_groups': 0,
        'quantum_reference_certified': False,
        'reference_subtracted_difference_covered': False,
        'new_actual_X_readout_certified': False,
        'scope': 'Original 1250 absolute total source, geometry and positive null weight. Conditional use of its old read-error contract only. No change of reference, geometry or renormalization.',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    result = compute()
    path = HERE / 'results.json'
    if args.write_results:
        with path.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
        action = 'FIRST_WRITE'
    else:
        assert json.loads(path.read_text('utf-8')) == result
        action = 'READ_ONLY_MATCH'
    print(json.dumps({'status': result['status'], 'action': action,
                      'X_fraction_strict_upper': '1/25',
                      'heat_omission_strict_lower': '19/1000000',
                      'new_rounds': 0}))


if __name__ == '__main__':
    main()
