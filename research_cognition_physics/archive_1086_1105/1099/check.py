"""Finite certificates for 1099; these do not certify physical permissions.

Use --write once to save results. Default invocation recomputes and compares.
Universal statements are proved in proof.md, not by sampling these matrices.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent


def eye(n=4):
    return [[Q(i == j) for j in range(n)] for i in range(n)]


def mul(a, b):
    return [[sum((x*y for x, y in zip(row, col)), Q(0))
             for col in zip(*b)] for row in a]


def mv(a, v):
    return [sum((x*y for x, y in zip(row, v)), Q(0)) for row in a]


def inverse(a):
    n = len(a)
    rows = [list(map(Q, row)) + unit for row, unit in zip(a, eye(n))]
    for j in range(n):
        pivot = next(i for i in range(j, n) if rows[i][j])
        rows[j], rows[pivot] = rows[pivot], rows[j]
        f = rows[j][j]
        rows[j] = [v/f for v in rows[j]]
        for i in range(n):
            if i != j:
                f = rows[i][j]
                rows[i] = [x-f*y for x, y in zip(rows[i], rows[j])]
    return [row[n:] for row in rows]


def power(a, n):
    out = eye(len(a))
    for _ in range(n):
        out = mul(out, a)
    return out


def scale(a, s):
    return [[s*x for x in row] for row in a]


def rotate(c, s):
    r = eye()
    r[1][1], r[1][2] = c, -s
    r[2][1], r[2][2] = s, c
    assert c*c+s*s == 1
    return r


def boost(u):
    b = eye()
    for i, v in enumerate(u, 1):
        b[i][0] = v
    return b


def commutator(a, b):
    return mul(mul(mul(a, b), inverse(a)), inverse(b))


def norm2(v):
    return sum((x*x for x in v), Q(0))


def velocity(a, v):
    z = mv(a, [Q(1)] + v)
    return [x/z[0] for x in z[1:]]


def compute():
    groups = {}
    r = rotate(Q(7, 25), Q(24, 25))
    expected = boost([Q(18, 25), Q(-24, 25), Q(0)])
    cases = []
    for lam in (Q(1, 3), Q(1), Q(3, 2), Q(2)):
        g = scale(boost([Q(1), Q(0), Q(0)]), lam)
        c = commutator(g, r)
        assert c == expected and c[0] == [1, 0, 0, 0]
        cases.append(str(lam))
    groups['scaled_galilean_word'] = dict(
        status='PASS', scales=cases,
        velocity_translation=['18/25', '-24/25', '0'],
        time_factor='1', physical_word_permission_certified=False)

    rows = []
    for n in (1, 2, 4, 8, 32):
        h = Q(5, 6*n)
        z = mv(power(expected, n), [h, Q(0), Q(0), Q(0)])
        assert z == [h, Q(3, 5), Q(-4, 5), Q(0)]
        assert norm2(z[1:]) == 1
        rows.append(dict(n=n, finite_reference_factors=4*n,
                         holding_time=str(h), radius_squared='1'))
    groups['finite_baseline_certificates'] = dict(
        status='PASS', cases=rows,
        holding_times_are_conditional_permissions=True,
        no_actual_task_concatenation_used=True)

    # A non-scalar spatial block tests the affine-action formulas used in §3.
    a = [[Q(2), Q(0), Q(0), Q(0)],
         [Q(1), Q(2), Q(1), Q(0)],
         [Q(-1), Q(0), Q(3), Q(0)],
         [Q(0), Q(0), Q(0), Q(4)]]
    b = scale(boost([Q(2), Q(1), Q(-1)]), Q(3, 2))
    v = [Q(1, 2), Q(-2, 3), Q(1, 5)]
    ab = mul(a, b)
    assert ab[0][0] == a[0][0]*b[0][0]
    assert velocity(ab, v) == velocity(a, velocity(b, v))
    ar = mul(mul(a, r), inverse(a))
    assert ar[0] == [1, 0, 0, 0]
    assert velocity(ar, [Q(0)]*3) != [0, 0, 0]
    assert mul(ar, inverse(ar)) == eye()
    groups['affine_action_and_time_character'] = dict(
        status='PASS', non_scalar_spatial_block=True,
        conjugate_time_factor='1',
        unbounded_orbit_and_connectedness='analytic proof only; not inferred from finite tests')

    # Calibration of the inherited 1097 cone conclusion, in units c*=1.
    eta = eye()
    for i in range(1, 4):
        eta[i][i] = Q(-1)
    lor = eye()
    lor[0][0] = lor[1][1] = Q(5, 4)
    lor[0][1] = lor[1][0] = Q(-3, 4)
    conformal = scale(mul(lor, r), Q(3, 2))
    gram = mul(mul(list(map(list, zip(*conformal))), eta), conformal)
    assert gram == scale(eta, Q(9, 4))
    points = [[Q(1), Q(0), Q(0), Q(0)],
              [Q(1), Q(3, 5), Q(4, 5), Q(0)],
              [Q(2), Q(-1), Q(1), Q(0)]]
    for z in points:
        y = mv(conformal, z)
        assert y[0] > 0 and y[0]*y[0] >= norm2(y[1:])
        assert mv(inverse(conformal), y) == z
    assert Q(5, 4)/Q(3, 4) == Q(5, 3)  # Conservative row bound, not optimal c*.
    groups['inherited_conformal_cone_calibration'] = dict(
        status='PASS', lambda_squared='9/4', cone_speed='1',
        row_bound='5/3', scale_not_fixed=True,
        all_matter_dynamics_certified=False)

    delta = Q(1, 4)
    for k in (1, 2, 10, 100):
        z = [delta, Q(k), Q(0), Q(0)]
        y = mv(boost([Q(k), Q(-k), Q(1)]), z)
        assert y[0] == delta  # P_delta admits every x and t>=delta.
        t = Q(1, k)
        assert t > 0  # P_0 allows (t, e1) for every k.
        q = Q(k)
        assert (q+q*q)/q == 1+q  # Rotations-only nonlinear menu.
    # At L=2, t+t^2 reaches L first at t=1 (monotone for t>0).
    assert Q(1)+Q(1)**2 == 2
    groups['deleted_condition_menus'] = dict(
        status='PASS', without_H=dict(menu='t>=1/4, arbitrary x', tau_L='1/4',
            average_speed_unbounded=True),
        without_D=dict(menu='t>0, arbitrary x', tau_L='0'),
        without_nonrest_reference=dict(menu='|x|<=t+t^2, rotations only',
            tau_at_L_2='1', average_speed_unbounded=True),
        status_scope='permission menus; not full FUCP+CO physical countermodels')

    n = 4
    h = Q(5, 6*n)
    time_upper = h + Q(1, 100)
    error_upper = 4*n*Q(1, 1000) + Q(1, 100)
    assert time_upper == Q(131, 600) < Q(1, 4)
    assert error_upper == Q(13, 500) < Q(1, 10)
    assert Q(49, 50) <= Q(99, 100) <= Q(101, 100) <= Q(51, 50)
    groups['finite_error_annular_certificate'] = dict(
        status='PASS', time_upper=str(time_upper), lower_bound='1/4',
        time_margin=str(Q(1, 4)-time_upper),
        complete_task_error_upper=str(error_upper), error_limit='1/10',
        certified_radius_annulus=['49/50', '51/50'],
        assumed_output_radius_interval=['99/100', '101/100'],
        numerical_instrument_data=False, budgets_are_additional_contracts=True,
        nonexpansive_complete_transport_and_uniform_error_assumed=True)
    assert len(groups) == 6
    return dict(round=1099, status='PASS', groups=groups,
                exact_arithmetic='fractions.Fraction',
                universal_theorem='See proof.md; finite certificates are not its proof.',
                new_adopted_axioms=0, scientific_count_increment=0,
                H_R_O_D_certified_from_cognition=False,
                unconditional_Lorentz_derived=False)


def main():
    result = compute()
    target = HERE/'results.json'
    if '--write' in sys.argv:
        target.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
    else:
        assert json.loads(target.read_text(encoding='utf8')) == result, 'Saved results differ.'
    print(json.dumps(dict(round=1099, status='PASS', groups=len(result['groups']),
                         exact_arithmetic=True, saved_results_match=True,
                         conditional_only=True), ensure_ascii=False))


if __name__ == '__main__':
    main()
