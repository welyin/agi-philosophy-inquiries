"""Exact budget replay. No sampling, image, or external dependencies."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent

def calculate():
    b0, upper = 4 * F(499, 500)**2, F(72)
    d = F(591628623, 124500500000)
    e = F(61, 20000)
    sigma, cutoff, step = F(1, 100000), F(1, 1000), F(1, 10000)
    failure = 3 * sigma**2 / cutoff**2
    precision_margin = (F(1, 100) - d) * b0 - 2 * e
    geometry_gap = b0 / 4 - F(15, 4) * e
    thermal_gap = (F(1, 4) - 2*d) * b0 - F(15, 4) * e
    # One overflow code, plus both ends and each step on the interval.
    number_codes = int(200 / step) + 2
    bits_per_channel = (number_codes - 1).bit_length()
    initial = dict(clock_support=F(2), control=F(2), store=F(4),
                   pointer=F(1), record=F(1), outgoing=F(0))
    transfers = [('store', 'control', F(3, 10)),
                 ('control', 'pointer', F(1, 10)),
                 ('control', 'record', F(1, 5)),
                 ('control', 'outgoing', F(1, 10))]
    final = initial.copy()
    for sender, receiver, value in transfers:
        final[sender] -= value
        final[receiver] += value
    assert all(x >= 0 for x in final.values())
    assert sum(initial.values()) == sum(final.values()) == 10
    assert F(0) < failure < F(1, 1000)
    assert (1+d)*upper + e < 100
    assert bits_per_channel == 21
    assert precision_margin > 0 and geometry_gap > 0 and thermal_gap > F(94, 100)
    assert (F(1, 100)-d)*b0-2 < 0
    q = lambda x: str(x)
    return {
        'status': 'PASS_EFFECTIVE_BUDGET_ONLY',
        'bounds': {'B_min': q(b0), 'B_max_strict': q(upper),
                   'relative_reference_bound': q(d), 'channel_radius': q(e),
                   'joint_failure_probability_upper': q(failure),
                   'one_percent_margin': q(precision_margin),
                   'geometry_omission_gap': q(geometry_gap),
                   'thermal_omission_gap': q(thermal_gap)},
        'finite_record': {'codes_per_channel': number_codes,
                          'bits_per_channel': bits_per_channel,
                          'total_bits': 3*bits_per_channel},
        'energy_ledger': {'initial': {k:q(v) for k,v in initial.items()},
                          'final': {k:q(v) for k,v in final.items()},
                          'sum_each': '10'},
        'removal_check': {'channel_radius_1_fails_one_percent_certificate': True},
        'scope': {'new_science_round': False, 'new_axiom': False,
                  'microscopic_instrument_proved': False,
                  'source_and_backreaction_bounds_are_effective_inputs': True,
                  'whole_goal_completed': False},
        'fixed_spec_sha256': hashlib.sha256((HERE/'fixed_spec.md').read_bytes()).hexdigest(),
    }

if __name__ == '__main__':
    data = calculate()
    dest = HERE/'results.json'
    if '--write' in sys.argv:
        with dest.open('x', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
            f.write('\n')
    else:
        assert json.loads(dest.read_text(encoding='utf-8')) == data
    print(json.dumps(data, ensure_ascii=False, indent=2))
