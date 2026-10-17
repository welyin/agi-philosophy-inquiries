"""Read-only direction checks; --write-once creates the initial result only."""
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent


def rat(x):
    return str(x.numerator) + '/' + str(x.denominator)


def calculate():
    mu, radius, arm = F(10**24), F(10**30), F(10**11)
    # Strict bounds 26 < m0 < 27, independent of any actual optical apparatus.
    near_ratio_lower = F(26, 1) / mu * (radius / arm)**3
    near_ratio_upper = F(27, 1) / mu * (radius / arm)**3
    assert near_ratio_lower == 26 * 10**33
    assert near_ratio_upper == 27 * 10**33
    s = F(10**40)
    assert near_ratio_upper / s < F(1, 100000)
    # Continuous finite bound uses the polynomial identity, not these samples.
    # (1-2x)^3-(1-3x)^2 = 3x^2-8x^3, whose factors are >=0 on [0,1/8].
    def mul(a, b):
        out = [F(0)] * (len(a) + len(b) - 1)
        for i, ai in enumerate(a):
            for j, bj in enumerate(b):
                out[i+j] += ai*bj
        return out
    cube = mul(mul([F(1), F(-2)], [F(1), F(-2)]), [F(1), F(-2)])
    square = mul([F(1), F(-3)], [F(1), F(-3)]) + [F(0)]
    lhs_coefficients = [a-b for a, b in zip(cube, square)]
    assert lhs_coefficients == [F(0), F(0), F(3), F(-8)]
    assert 3 - 8 * F(1, 8) == 2
    sampled = []
    with localcontext() as ctx:
        ctx.prec = 80
        for x in [F(0), F(1, 10**6), F(1, 100), F(1, 8)]:
            xd = Decimal(x.numerator) / Decimal(x.denominator)
            value = 1 - (1 - 3 * xd) / (1 - 2 * xd).sqrt()
            assert 2 * xd <= value <= 3 * xd
            sampled.append(dict(x=rat(x), full_sector_per_radian=str(value)))
    # Interface terms cancel exactly, without bounding them separately.
    for slope_inner, slope_outer in [(F(3, 4), F(999,1000)), (F(11,10), F(1,2))]:
        assert (1-slope_inner)+(slope_inner-slope_outer) == 1-slope_outer
    # Illustrative supplied constants only: no common-model qualification.
    signal, fixed, a, b = F(1,100), F(1,1000), F(2), F(3)
    delta = signal/2-fixed
    selected_s = F(10**8)
    bound = a/F(10**4)+b/selected_s
    threshold = max(F(1), (2*a/delta)**2, 2*b/delta)
    assert selected_s > threshold and fixed+bound < signal/2
    return dict(status='PASS_DIRECTION_ARITHMETIC_ONLY', new_scientific_rounds=0,
        near_field_diagnostic=dict(lower_strict=rat(near_ratio_lower), upper_strict=rat(near_ratio_upper),
            conditions='Weak monopoles with sampling distance comparable to arm; not a vertex error lower bound.',
            nonhomothetic_ratio_upper_at_s_1e40=rat(near_ratio_upper/s)),
        regular_sector=dict(analytic_lower='2x', analytic_upper='3x', x_domain='[0,1/8]',
            proof='Positive square-root comparison and x^2(3-8x)>=0; matched full sector only.', sampled=sampled,
            interface_cancellation='Exact algebra; no multiple-body geometry inferred.'),
        supplied_budget_example=dict(signal=rat(signal), fixed_error=rat(fixed), A=rat(a), B=rat(b),
            sufficient_threshold=rat(threshold), selected_s=rat(selected_s), actual_bound=rat(fixed+bound)),
        scope=dict(actual_null_triangle=False, source_or_instrument_qualified=False,
            full_quantum_process_certified=False, shared_finite_parameter_certificate=False,
            hbar_rescaled=False, CT_adopted=False, GR_completed=False))


def main():
    result = calculate()
    result['check_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['note_sha256'] = hashlib.sha256((HERE/'README.md').read_bytes()).hexdigest()
    target = HERE/'results.json'
    if sys.argv[1:] == ['--write-once']:
        with target.open('x', encoding='utf-8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write('\n')
        print('WROTE_DIRECTION_RESULT_NO_NEW_ROUND')
    else:
        assert not sys.argv[1:]
        assert json.loads(target.read_text('utf-8')) == result
        print('READ_ONLY_MATCH_DIRECTION_ARITHMETIC_NO_NEW_ROUND')


if __name__ == '__main__':
    main()
