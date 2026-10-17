"""Exact initial-constraint diagnostic, not a new numbered scientific group.

No ADM evolution, actual preparation, finite-record discrimination, or EFT
no-go is claimed. Default is read-only; --write exclusively creates this
nonnumbered diagnostic's result on its first run.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUTPUT = HERE / 'finite_ct_continuation_results_after1245.json'


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p)


def add(*polys):
    n = max(map(len, polys))
    out = [F(0)] * n
    for p in polys:
        for i, c in enumerate(p):
            out[i] += c
    return trim(out)


def scale(c, p):
    return trim([c * x for x in p])


def mul(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] += a*b
    return trim(out)


def derivative(p):
    return trim([i*p[i] for i in range(1, len(p))] or [F(0)])


def value(p, x):
    out = F(0)
    for c in reversed(p):
        out = out*x + c
    return out


def average_sine_polynomial(p, epsilon):
    """Exact normalized integral for p(epsilon*sin(x)); no quadrature."""
    out = F(0)
    for power, c in enumerate(p):
        if power % 2 == 0:
            k = power // 2
            out += c * epsilon**power * F(comb(2*k, k), 4**k)
    return out


def pack_poly(p):
    return [str(c) for c in p]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    one, p = (F(1),), (F(0), F(1))
    eps = F(1, 10)
    a, k2 = F(6), F(1, 10)
    p2 = mul(p, p)
    grad_p_squared = add((eps**2,), scale(-1, p2))
    f = add(one, scale(F(1, 5), p))
    mc = F(27, 10)
    mr = add((F(1, 5),), scale(F(1, 10), p))
    nr = F(1, 10)

    pi_diag = [add(one, scale(F(1, 2), p)), one, one]
    trace = add(*pi_diag)
    norm2 = add(*(mul(item, item) for item in pi_diag))
    hg = scale(a, add(norm2, scale(-F(1, 2), mul(trace, trace))))
    hf = add((F(1, 2),), scale(F(1, 2), grad_p_squared),
             scale(F(1, 2), p2), scale(3*k2/2, f))
    hr = scale(nr, mr)
    # Fixed positive clock preparation, independent of Q=ab.
    nc = scale(1/mc, (F(333, 40), F(74, 25), -F(3, 4)))
    hc = scale(mc, nc)
    total_h = add(hg, hf, hr, hc)
    dg_factor = scale(-2, derivative(pi_diag[0]))
    dm_factor = one  # P_p=1, all rod and mechanical momenta vanish.
    total_d_factor = add(dg_factor, dm_factor)

    assert hg == (-F(9), -F(3), F(3, 4))
    assert hf == (F(131, 200), F(3, 100))
    assert hr == (F(1, 50), F(1, 100))
    assert hc == (F(333, 40), F(74, 25), -F(3, 4))
    assert total_h == (0,) and total_d_factor == (0,)
    assert derivative(trace) == (F(1, 2),)

    # These are exact all-p interval bounds, not a sample grid.
    nc_derivative = derivative(nc)
    nc_derivative_min = min(value(nc_derivative, -eps),
                            value(nc_derivative, eps))
    assert nc_derivative_min > 0  # derivative is affine on [-eps, eps].
    nc_bounds = [value(nc, -eps), value(nc, eps)]
    f_bounds = [value(f, -eps), value(f, eps)]
    mr_bounds = [value(mr, -eps), value(mr, eps)]
    action_rate_bounds = [(1-eps)/10, (1+eps)/10]
    assert nc_bounds == [F(16043, 5400), F(17227, 5400)]
    assert min(nc_bounds + f_bounds + mr_bounds + action_rate_bounds) > 0
    assert f_bounds[0] == F(49, 50)
    assert mr_bounds[0] == F(19, 100)

    lapse_specs = [
        ('unit', one, F(1)),
        ('fixed_wave', add(one, scale(F(1, 4)/eps, p)), F(3, 4)),
    ]
    integrals = []
    for name, lapse, n0 in lapse_specs:
        square_norm = mul(mul(lapse, lapse), grad_p_squared)
        integral = average_sine_polynomial(square_norm, eps)
        guaranteed = n0*n0*eps*eps/2
        assert integral >= guaranteed > 0
        # Exact divergence integration by parts: average[p*d_x(N^2*p')].
        # d_x(A(p)*p') = A'(p)*(eps^2-p^2) - A(p)*p.
        lapse2 = mul(lapse, lapse)
        divergence = add(mul(derivative(lapse2), grad_p_squared),
                         scale(-1, mul(lapse2, p)))
        ibp = average_sine_polynomial(mul(p, divergence), eps)
        assert ibp == -integral
        integrals.append(dict(name=name, lapse_poly=pack_poly(lapse),
                              minimum_lapse=str(n0),
                              weighted_gradient_norm=str(integral),
                              analytic_positive_lower_bound=str(guaranteed),
                              integration_by_parts=str(ibp)))
    assert integrals[0]['weighted_gradient_norm'] == '1/200'
    assert integrals[1]['weighted_gradient_norm'] == '13/2560'

    cases = []
    for q in (F(1), F(1001, 1000)):
        b = q/a
        assert b > 0 and a*b == q
        w_factor = 1-q  # lambda=1/2; geometric trace term is exactly zero.
        # At unit lapse dot(H density)=d_x W=(q-1)*p.
        unit_derivative = scale(-w_factor, p)
        smear_derivative = average_sine_polynomial(mul(p, unit_derivative), eps)
        assert smear_derivative == (q-1)/200
        fixed_wave_derivative = (q-1)*F(13, 2560)
        cases.append(dict(q=str(q), a=str(a), b=str(b),
                          same_initial_H='0', same_initial_D='0',
                          W_gradient_factor=str(w_factor),
                          unit_lapse_Hdot_poly=pack_poly(unit_derivative),
                          unit_lapse_H_p_derivative=str(smear_derivative),
                          fixed_wave_H_Np_derivative=str(fixed_wave_derivative),
                          passes_this_initial_preservation_condition=(q == 1)))
    assert cases[0]['unit_lapse_H_p_derivative'] == '0'
    assert cases[1]['unit_lapse_H_p_derivative'] == '1/200000'
    assert cases[1]['fixed_wave_H_Np_derivative'] == '13/2560000'

    sources = [
        'research_cognition_physics/archive_1234_/1245/proof.md',
        'research_cognition_physics/archive_1234_/1235/proof.md',
        'research_cognition_physics/archive_1161_1200/1183/proof.md',
    ]
    return dict(
        status='PASS_EXACT_INITIAL_CONSTRAINT_SCOPE',
        scope='Same-material mathematical initial data only; no physical '
              'preparation, future record prediction, generic PDE existence, '
              'or EFT no-go certification.',
        lambda_H='1/2', fixed_spatial_wave_number=1, epsilon=str(eps),
        initial_data_shared_across_q=True,
        original_1242_initial_preparation=False,
        fixed_label_measure_reprepared=True,
        hamiltonian_polynomials=dict(gravity=pack_poly(hg), fields=pack_poly(hf),
                                    recorder=pack_poly(hr), clock=pack_poly(hc),
                                    total=pack_poly(total_h)),
        momentum_gradient_factors=dict(gravity=pack_poly(dg_factor),
                                      matter=pack_poly(dm_factor),
                                      total=pack_poly(total_d_factor)),
        positivity=dict(clock_mass=str(mc), recorder_number_density=str(nr),
                        clock_density=[str(v) for v in nc_bounds],
                        clock_density_derivative_lower=str(nc_derivative_min),
                        F=[str(v) for v in f_bounds],
                        recorder_mass=[str(v) for v in mr_bounds],
                        recorder_action_rate=[str(v) for v in action_rate_bounds]),
        fixed_low_wave_number_integrals=integrals,
        coefficient_cases=cases,
        zero_deviation_scope='Necessary preservation obstruction vanishes; '
                             'this script does not prove local existence.',
        nonzero_deviation_scope='The analytic positive-lapse integral excludes '
                                'this preparation for the exact parent '
                                'constraints; it does not exclude all data '
                                'or any consistent EFT completion.',
        new_numbered_rounds=0,
        new_formal_scientific_groups=0,
        source_sha256={rel: sha(ROOT/rel) for rel in sources},
        code_sha256=sha(Path(__file__)),
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write:
        with OUTPUT.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    else:
        stored = json.loads(OUTPUT.read_text('utf-8'))
        assert result == stored, 'Read-only result or source hash mismatch'
    print(json.dumps(dict(status=result['status'],
                          mode='INITIAL_SAVE' if args.write else 'READ_ONLY_MATCH',
                          exact_initial_constraints=True,
                          positive_lapse_integrals=[r['weighted_gradient_norm']
                                                   for r in result['fixed_low_wave_number_integrals']],
                          unit_lapse_contrast=result['coefficient_cases'][1]['unit_lapse_H_p_derivative'],
                          new_numbered_rounds=0, new_formal_scientific_groups=0)))
