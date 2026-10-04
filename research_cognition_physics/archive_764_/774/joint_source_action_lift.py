"""774: exact diagnostics for a source-preserving formal action lift.

These are finite algebra checks, not computations of the original loop anomaly.
The continuum statement has the explicitly listed local-QAP premises in the note.
"""
import argparse
from fractions import Fraction as Q
from math import comb, prod
from pathlib import Path
import hashlib
import itertools
import json
import numpy as np
import joint_material_local_brst as old

HERE = Path(__file__).resolve().parent
TARGET = HERE/'joint_source_action_lift_results.json'
v, add, scale, mul = old.var, old.add, old.scale, old.mul
one = old.ONE


def degree(k):
    return sum(k[0]) + k[1].bit_count()


def cut(p, n=7):
    return {k: x for k, x in p.items() if degree(k) <= n}


def part(p, n):
    return {k: x for k, x in p.items() if degree(k) == n}


def power(p, n):
    out = one
    for _ in range(n):
        out = mul(out, p)
    return out


def inverse_t(n, order=9):
    return add(*(scale(power(v('q0'), j), Q((-1)**j*comb(n+j-1, j)))
                 for j in range(order+1)))


# Exact formal gauge model: t=1+q, y=r/t, S=y^2/2; sq=t*c, sr=r*c.
# p,a,z are the gauge, physical and ghost antifields, with the 773 convention.
def s(p):
    t = add(one, v('q0'))
    rules = {'q0': mul(t, v('c0')), 'r': mul(v('r'), v('c0')),
             'p': add(scale(mul(power(v('r'), 2), inverse_t(3)), -1),
                      mul(v('p'), v('c0'))),
             'a': add(mul(v('r'), inverse_t(2)), mul(v('a'), v('c0'))),
             'z': add(mul(t, v('p')), mul(v('r'), v('a')))}
    return cut(old.derive(p, rules, 1))


# Adding the invariant insertion zeta*y changes only the antifield equations.
def insertion_d(p):
    return cut(old.derive(p, {'p': scale(mul(v('r'), inverse_t(2)), -1),
                             'a': inverse_t(1)}, 1))


def solve_closed(rhs, prescribed=None, order=6):
    result = prescribed or {}
    for n in range(1, order+1):
        residual = part(add(rhs, scale(s(result), -1)), n)
        assert not old.s(residual), ('not free-closed', n)
        correction = old.homotopy(residual)
        assert old.s(correction) == residual
        result = add(result, correction)
    assert not cut(add(rhs, scale(s(result), -1)), order)
    return result


def first_jet():
    # Minimal BV projection includes a one-loop transformation term off shell.
    E = (Q(0), Q(0))
    R = ((Q(2), Q(1)), (Q(-1), Q(3)))
    j = (Q(5, 7), Q(-4, 9))
    r1 = ((Q(3, 5), Q(7, 8)), (Q(-2, 7), Q(1, 3)))
    ward = tuple(sum(j[i]*R[i][a]+E[i]*r1[i][a] for i in range(2))
                 for a in range(2))
    linear = tuple(sum(j[i]*R[i][a] for i in range(2)) for a in range(2))
    assert ward == linear
    Eoff = (Q(1, 4), Q(-3, 7))
    omitted = tuple(sum(Eoff[i]*r1[i][a] for i in range(2)) for a in range(2))
    assert all(x != 0 for x in omitted)
    # Reuse the external-leg identity as an audit, not a new graph theorem.
    allowed = []
    for count in range(1, 5):
        for valences in itertools.combinations_with_replacement(range(3, 8), count):
            if sum(n-2 for n in valences) == 1:
                allowed.append(list(valences))
    assert allowed == [[3]]
    return dict(on_shell_projection_exact=True,
                projected_ward=[str(x) for x in ward],
                omitted_off_shell_term=[str(x) for x in omitted],
                one_loop_one_leg_vertex_lists=allowed,
                scope='Finite BV first-jet projection and reused valence counting; no original continuum loop coefficients evaluated.')


def simultaneous_lift():
    y = cut(mul(v('r'), inverse_t(1)))
    for name in ('q0', 'r', 'c0', 'p', 'a', 'z'):
        z = v(name)
        assert not cut(s(s(z)), 6)
        assert not cut(add(s(insertion_d(z)), insertion_d(s(z))), 6)
        assert not cut(insertion_d(insertion_d(z)), 6)
    assert not cut(s(y), 6)
    seeds = [add(scale(v('q0'), Q(2, 3)), scale(v('r'), Q(4, 5)),
                 old.product(['q0', 'q0', 'r']), old.product(['p', 'c0', 'r'])),
             add(old.product(['q0', 'r']), old.product(['a', 'c0', 'q0'])),
             old.product(['q0', 'q0', 'r'])]
    anomalies = []
    for k in range(3):
        anomalies.append(add(s(seeds[k]), insertion_d(seeds[k-1]) if k else {}))
    desired = add(scale(v('q0'), Q(2, 3)), scale(v('r'), Q(7, 11)))
    # Retain an independently prescribed physical first jet, unlike the seed.
    result = []
    for k in range(3):
        rhs = add(anomalies[k], scale(insertion_d(result[k-1]), -1) if k else {})
        # insertion_d lowers field degree: reserve two extra orders in C_0,
        # one extra in C_1, and verify the common final order four.
        result.append(solve_closed(rhs, desired if k == 0 else None, 6-k))
    assert part(result[0], 1) == desired
    assert part(result[0], 1) != part(seeds[0], 1)
    for k in range(3):
        residual = add(s(result[k]), insertion_d(result[k-1]) if k else {},
                       scale(anomalies[k], -1))
        assert not cut(residual, 4)
    # Counterterms with the negative of these coefficients cancel the defect.
    return dict(source_orders_checked=3, common_field_degree_verified=4,
                internal_computation_degrees=[6, 5, 4],
                action_counterterm_terms=[len(x) for x in result],
                prescribed_first_jet_exact=True, nonlinear_nilpotency_exact=True,
                invariant_observable_exact_to_checked_order=True,
                joint_action_single_and_double_insertion_residual_zero=True,
                original_loop_anomaly_computed=False,
                scope='Exact rational nonlinear minimal-BV diagnostic with antifields; anomalies are manufactured coboundaries and do not measure the original theory.')


def finite_parts():
    # Centered Gaussian Wick powers: hbar*C_d first contributes to a mean
    # force at hbar^(1+(d-1)/2), only for odd d. No state subtraction is used.
    orders = {}
    for d in range(1, 8):
        n = d-1
        if n % 2:
            orders[str(d)] = None
        else:
            moment = prod(range(1, n, 2)) if n else 1
            orders[str(d)] = {'hbar_order': 1+n//2, 'coefficient': d*moment}
    assert orders['1']['hbar_order'] == 1
    assert all(x is None or x['hbar_order'] >= 2 for k, x in orders.items() if k != '1')
    H = np.array([[2., .3], [.3, 1.5]])
    B = np.array([[.7, -.2], [-.2, .4]])
    C = np.linalg.inv(H)
    correction = -C@B@C
    assert np.linalg.norm(correction) > .1
    errors = []
    for h in (1e-3, 5e-4, 2.5e-4):
        actual = h*np.linalg.inv(H+h*B)
        errors.append(float(np.max(np.abs(actual-h*C-h*h*correction))/h**2))
    assert errors[2] < .27*errors[0]
    # A source-dependent action provides symmetric companion derivatives.
    p = add(old.product(['q0', 'q0', 'r']), old.product(['q0', 'r', 'r']))
    dq = lambda f: old.derive(f, {'q0': one}, 0)
    dr = lambda f: old.derive(f, {'r': one}, 0)
    assert dq(dr(p)) == dr(dq(p))
    return dict(centered_first_mean_orders=orders,
                quadratic_finite_term_changes_higher_correlations=True,
                covariance_hbar2_coefficient=correction.tolist(),
                scaled_expansion_errors=errors,
                companion_derivatives_commute=True,
                scope='Finite Gaussian order diagnostic, not original interacting positivity or an instrument; higher finite action terms can have physical effects.')


def run():
    dependencies = ['research_note_770.md', 'research_note_771.md',
                    'research_note_772.md', 'research_note_773.md',
                    'joint_material_local_brst.py', 'round774_drafts/STATUS.md']
    return dict(round=774, tests_run=3, failures=0, errors=0,
                first_jet=first_jet(), joint_lift=simultaneous_lift(),
                finite_parts=finite_parts(),
                dependency_hashes={p: hashlib.sha256((HERE/p).read_bytes()).hexdigest()
                                   for p in dependencies},
                scope='Conditional local one-loop action/source/insertion lift in the regular reference patch and enlarged formal coefficient class; no global quantum theory, interacting positive state, actual instrument, continuum equivalence or cognition-only GR derivation.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write:
        with TARGET.open('x', encoding='utf8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
    else:
        assert result == json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result, ensure_ascii=False, indent=2))
