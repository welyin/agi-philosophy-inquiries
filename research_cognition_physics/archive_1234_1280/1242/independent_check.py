"""Independent 1242 audit: exact polynomial identities + a-variable Decimal RK4.

The Decimal integrations are diagnostics, never the rigorous time-window proof.
The small-window certificate uses exact rational inequalities. Author interval
Taylor code is not imported. Default execution is read-only; --write creates once.
"""
from __future__ import annotations

import argparse
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
NVAR = 5  # u, H, p, v, constant dust coefficient D0
ZERO = (0,) * NVAR


def const(c):
    c = F(c)
    return {ZERO: c} if c else {}


def variable(j):
    e = list(ZERO)
    e[j] = 1
    return {tuple(e): F(1)}


def add(*ps):
    result = {}
    for p in ps:
        for e, c in p.items():
            result[e] = result.get(e, F(0)) + c
    return {e: c for e, c in result.items() if c}


def scale(p, a):
    return {e: c * a for e, c in p.items() if c * a}


def mul(*ps):
    result = const(1)
    for p in ps:
        nxt = {}
        for e, c in result.items():
            for f, b in p.items():
                ef = tuple(x + y for x, y in zip(e, f))
                nxt[ef] = nxt.get(ef, F(0)) + c * b
        result = {e: c for e, c in nxt.items() if c}
    return result


def derivative(p, j):
    result = {}
    for e, c in p.items():
        if e[j]:
            f = list(e)
            f[j] -= 1
            result[tuple(f)] = c * e[j]
    return result


def lie(p, fields):
    return add(*(mul(derivative(p, j), fields[j]) for j in range(NVAR)))


def evaluate(p, point):
    return sum((c * product(x ** n for x, n in zip(point, e))
                for e, c in p.items()), F(0))


def product(xs):
    r = F(1)
    for x in xs:
        r *= x
    return r


def exact_checks():
    u, H, p, v, dust0 = [variable(j) for j in range(NVAR)]
    u2, u3 = mul(u, u), mul(u, u, u)
    rod = scale(mul(add(const(1), scale(p, F(1, 5))), u2), F(3, 20))
    dust = mul(add(dust0, scale(p, F(1, 100))), u3)
    kinetic, potential = scale(mul(v, v), F(1, 2)), scale(mul(p, p), F(1, 2))
    rho = add(rod, dust, kinetic, potential)
    pressure = add(scale(rod, -F(1, 3)), kinetic, scale(potential, -1))
    rhoP = add(rho, pressure)
    fields = [scale(mul(H, u), -1), scale(rhoP, -F(3, 2)), v,
              add(scale(mul(H, v), -3), scale(p, -1),
                  scale(u2, -F(3, 100)), scale(u3, -F(1, 100))), {}]
    continuity = add(lie(rho, fields), scale(mul(H, rhoP), 3))
    constraint = add(mul(H, H), scale(rho, -1))
    constraint_dot = lie(constraint, fields)
    assert not continuity and not constraint_dot

    # Derive both deletions, rather than inserting their anticipated signs.
    no_force = list(fields)
    no_force[3] = add(no_force[3], scale(u3, F(1, 100)))
    force_defect = add(lie(rho, no_force), scale(mul(H, rhoP), 3))
    target = scale(mul(v, u3), F(1, 100))
    assert force_defect == target
    no_stress_rho = add(rho, scale(mul(p, u3), -F(1, 100)))
    stress_defect = add(lie(no_stress_rho, fields),
                        scale(mul(H, add(no_stress_rho, pressure)), 3))
    assert stress_defect == scale(target, -1)

    # Independent exact Lie derivative for the initial fourth record jet.
    z4 = scale(lie(lie(v, fields), fields), F(1, 10))
    for d0 in (F(51, 50), F(41, 50)):
        for h0 in (F(1), F(11, 10)):
            assert evaluate(z4, [F(1), h0, F(0), F(0), d0]) == F(21, 1000) * h0

    # Hand-derived common tube at T=1/100, with exact rational arithmetic.
    # u in [49/50,1], H in [9/10,11/10], |p|,|v| <= 1/1000.
    T, umax, hmax, pmax, vmax = F(1, 100), F(1), F(11, 10), F(1, 1000), F(1, 1000)
    b, Q = F(3, 100), F(1, 100)
    Fmax, dustmax = 1 + F(1, 5) * pmax, F(51, 50) + Q * pmax
    h1 = F(3, 2) * (vmax ** 2 + F(1, 10) * Fmax + dustmax)
    p2 = 3 * hmax * vmax + pmax + b + Q
    assert h1 < F(17, 10) and p2 == F(443, 10000)
    assert F(98, 100) ** 2 < F(97, 100)
    assert F(11, 10) ** 2 > F(117, 100)
    assert F(98, 100) - T * F(17, 10) > F(9, 10)
    assert 1 - hmax * T > F(49, 50)
    assert T * vmax < pmax and T * p2 < vmax
    h2 = F(3, 2) * (2 * vmax * p2 + Q * vmax
                    + 3 * hmax * dustmax + F(1, 50) * vmax
                    + F(1, 5) * hmax * Fmax)
    p3 = 3 * F(17, 10) * vmax + 3 * hmax * p2 + vmax + 2 * b * hmax + 3 * Q * hmax
    assert h2 < F(27, 5) and p3 < F(63, 250)
    p4 = (3 * F(27, 5) * vmax + 6 * F(17, 10) * p2
          + 3 * hmax * F(63, 250) + p2
          + 2 * b * (F(17, 10) + 2 * hmax ** 2)
          + 3 * Q * (F(17, 10) + 3 * hmax ** 2))
    assert p4 == F(87553, 50000) and p4 < 2
    # z = theta - omega*t, so |z^(5)| < 1/5 for each history.
    # H+ - H- > (1/5)/(21/10)=2/21, strictly.
    low = T ** 4 / 12000 - T ** 5 / 300
    assert low == F(1, 2 * 10 ** 12)
    assert 1 - F(1, 5) * pmax > 0
    assert F(1, 5) - F(1, 10) * pmax > 0
    # Independently bound the author common-domain Jacobian in the infinity norm.
    # On u<=1.05,H<=1.2,|p|<=.02,|v|<=.05, both D values share this bound.
    U, HH, PP, VV = F(21, 20), F(6, 5), F(1, 50), F(1, 20)
    duH = F(3, 2) * (3 * (F(51, 50) + Q * PP) * U ** 2
                       + F(1, 5) * (1 + F(1, 5) * PP) * U)
    dpH = F(3, 2) * (Q * U ** 3 + F(1, 50) * U ** 2)
    lipschitz = max(HH + U, duH + dpH + 3 * VV, F(1),
                    2 * b * U + 3 * Q * U ** 2 + 3 * VV + 1 + 3 * HH,
                    F(1, 10))
    assert lipschitz < 6 and lipschitz / 100 < F(3, 50)
    assert F(1, 10) - F(1, 10) * PP > 0  # positive record frequency, not only mass
    # The fixed-geometry subsystem p,v,z contains no dust0, hence same IVP.
    assert all(e[4] == 0 for e in fields[3])
    return {
        'status': 'PASS',
        'identity_method': 'exact multivariate rational polynomial Lie derivative; all coefficients vanish',
        'continuity_polynomial_terms': len(continuity),
        'constraint_derivative_polynomial_terms': len(constraint_dot),
        'deleting_pointer_force': '+Q*v*u^3',
        'deleting_interaction_stress': '-Q*v*u^3',
        'z_fourth_initial_derivative': '(21/1000)*H0',
        'finite_short_window': {
            'T': str(T), 'abs_Hprime_bound': str(h1),
            'abs_p_second_bound': str(p2), 'abs_H_second_bound': str(h2),
            'abs_p_third_bound': str(p3), 'abs_p_fourth_bound': str(p4),
            'strict_record_gap_lower_bound': str(low),
            'method': 'first-exit common tube plus Taylor theorem and rational strict inequalities',
        },
        'fixed_geometry_record_gap': '0 by same scalar IVP; counterfactual only',
        'author_common_domain_Picard_audit': {
            'infinity_norm_Jacobian_bound': str(lipschitz),
            'step_contraction_bound': str(lipschitz / 100),
            'strictly_below_one': True,
            'role': 'combined with author all-step tube self-inclusion, certifies unique whole-step solution',
        },
    }


def energy(a, p, v, dust0):
    inv = 1 / a
    rod = D('.15') * (1 + D('.2') * p) * inv ** 2
    dust = (dust0 + D('.01') * p) * inv ** 3
    return dust + rod + (v * v + p * p) / 2


def evolution(y, dust0):
    a, p, v, z = y
    H = energy(a, p, v, dust0).sqrt()
    return [a * H, v, -3 * H * v - p - D('.03') / a ** 2 - D('.01') / a ** 3,
            D('.1') * p]


def rk4(dust0, steps, final=D('.2')):
    y = [D(1), D(0), D(0), D(0)]
    dt = final / steps
    for _ in range(steps):
        k1 = evolution(y, dust0)
        k2 = evolution([x + dt * k / 2 for x, k in zip(y, k1)], dust0)
        k3 = evolution([x + dt * k / 2 for x, k in zip(y, k2)], dust0)
        k4 = evolution([x + dt * k for x, k in zip(y, k3)], dust0)
        y = [x + dt * (q1 + 2 * q2 + 2 * q3 + q4) / 6
             for x, q1, q2, q3, q4 in zip(y, k1, k2, k3, k4)]
    a, p, v, z = y
    H = energy(a, p, v, dust0).sqrt()
    return [1 / a, H, p, v, z]


def diagnostics(author):
    rows = []
    finals = []
    with localcontext() as context:
        context.prec = 70
        for label, dust0 in [('high', D('1.02')), ('low', D('.82'))]:
            runs = {n: rk4(dust0, n) for n in [1024, 2048, 4096]}
            diff1 = max(abs(x - y) for x, y in zip(runs[1024], runs[2048]))
            diff2 = max(abs(x - y) for x, y in zip(runs[2048], runs[4096]))
            assert diff2 < diff1 / 14
            assert diff2 < D('1e-17')
            intervals = author['group_finite_records'][label]['state']
            margins = []
            for x, pair in zip(runs[4096], intervals):
                lower, upper = map(D, pair)
                assert lower < x < upper, (label, x, pair)
                margins.append(min(x - lower, upper - x))
            finals.append(runs[4096])
            rows.append({
                'branch': label, 'dust_D': str(dust0),
                'state_u_H_p_v_z': [str(x) for x in runs[4096]],
                'step_max_difference_1024_2048': str(diff1),
                'step_max_difference_2048_4096': str(diff2),
                'observed_refinement_ratio': str(diff1 / diff2),
                'all_components_strictly_inside_author_interval': True,
                'minimum_author_interval_margin': str(min(margins)),
            })
        gap = finals[0][4] - finals[1][4]
        gl, gu = map(D, author['group_finite_records']['gap'])
        assert gl < gap < gu
        assert gl - D('2e-8') > D('6.47e-8')
        Bgap = D('.1') * (finals[1][0] ** 2 - finals[0][0] ** 2)
        bl, bu = map(D, author['group_finite_records']['B_separation'])
        assert bl < Bgap < bu
        # Full theta stays in a single chart during the declared window.
        theta = [D('.02') + a[4] for a in finals]
        assert all(D(0) < x < D('.021') for x in theta)
        return {
            'status': 'PASS', 'method': '70-digit Decimal RK4 with a variable and H=sqrt(full rho)',
            'author_algorithm_imported': False,
            'precision_claim': 'diagnostic independent convergence, not a rigorous ODE error certificate',
            'branches': rows, 'record_gap': str(gap), 'B_separation': str(Bgap),
            'contract_error_per_history': '1e-8',
            'author_rigorous_gap_minus_two_contract_errors_lower': str(gl - D('2e-8')),
            'actual_quantum_readout_built': False,
            'readout_contract_inherited_from_1241': False,
        }


def run():
    author = json.loads((HERE / 'results.json').read_text('utf-8'))
    exact = exact_checks()
    numerical = diagnostics(author)
    bindings = {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
                for name in ['admission.md', 'proof.md', 'check.py', 'results.json', 'NEXT.md']}
    return {
        'round': 1242, 'status': 'PASS', 'science_groups': 2,
        'group_exact_source_and_short_window': exact,
        'group_independent_finite_history_diagnostics': numerical,
        'author_sha256': bindings,
        'scope': 'conditional finite homogeneous classical histories and terminal phase-readout contract',
        'CT_adopted': False, 'whole_GR_proved': False,
        'frozen_geometry_is_another_Einstein_solution': False,
        'preparation_or_conversion_protocol_certified': False,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    output = HERE / 'independent_results.json'
    if args.write:
        with output.open('x', encoding='utf-8') as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
        mode = 'CREATED_ONCE'
    else:
        assert json.loads(output.read_text('utf-8')) == result, 'independent result mismatch'
        mode = 'READ_ONLY_MATCH'
    print(json.dumps({'status': 'PASS', 'mode': mode, 'groups': 2,
                      'gap': result['group_independent_finite_history_diagnostics']['record_gap']},
                     ensure_ascii=False))
