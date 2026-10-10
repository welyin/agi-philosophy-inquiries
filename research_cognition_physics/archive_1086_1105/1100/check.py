"""1100 finite certificates; reuses frozen 1099 exact matrix arithmetic.

The universal FD theorem is analytic. Menu examples are not physical models.
Default invocation recomputes and compares results; --write saves explicitly.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
import runpy
import sys

HERE = Path(__file__).resolve().parent
M = runpy.run_path(str(HERE.parent/'1099/check.py'), run_name='frozen_round1099')
eye, mul, mv = (M[k] for k in ('eye', 'mul', 'mv'))
inverse, boost, scale = (M[k] for k in ('inverse', 'boost', 'scale'))
commutator, norm2 = (M[k] for k in ('commutator', 'norm2'))


def menu_lorentz(z, delta=Q(1)):
    return z[0] > 0 and z[0]*z[0]-norm2(z[1:]) >= delta*delta


def relation_probability(x1):
    return (1+x1/(1+abs(x1)))/2


def compute():
    groups = {}
    rot = eye()
    rot[1][1] = rot[2][2] = Q(-1)
    rows = []
    # Each chosen translation is a finite actual-reference word only if R holds.
    for t, start, target in (
        (Q(1), [Q(0), Q(0), Q(0)], [Q(1), Q(0), Q(0)]),
        (Q(2), [Q(1), Q(-1), Q(0)], [Q(-2), Q(3), Q(0)]),
        (Q(3, 2), [Q(-1), Q(2), Q(1)], [Q(5), Q(-4), Q(1)]),
    ):
        u = [(y-x)/(2*t) for x, y in zip(start, target)]
        g = scale(boost(u), Q(2))
        word = commutator(g, rot)
        assert word[0] == [1, 0, 0, 0]
        result = mv(word, [t]+start)
        assert result == [t]+target
        rows.append(dict(time=str(t), start=list(map(str, start)),
                         target=list(map(str, target)), reference_factors=4))
    groups['fixed_time_transport_certificate'] = dict(status='PASS', cases=rows,
        stationary_task_not_required=True,
        general_transitivity='1099 orbit lemma plus proof.md; not a sampling inference')

    lor = eye()
    lor[0][0] = lor[1][1] = Q(5, 4)
    lor[0][1] = lor[1][0] = Q(-3, 4)
    rest = [Q(1), Q(0), Q(0), Q(0)]
    output = mv(lor, rest)
    assert output == [Q(5, 4), Q(-3, 4), Q(0), Q(0)]
    assert menu_lorentz(rest) and menu_lorentz(output)
    assert mv(inverse(lor), output) == rest
    for h in (Q(1, 2), Q(1, 4), Q(1, 100)):
        assert not menu_lorentz([h, Q(0), Q(0), Q(0)])
    deadline = Q(5, 4)
    assert deadline*deadline-1 == Q(9, 16)
    assert menu_lorentz([deadline, Q(3, 4), Q(0), Q(0)])
    assert not menu_lorentz([deadline, Q(1), Q(0), Q(0)])
    # For t<=deadline the target x=e1 has still smaller quadratic form.
    assert deadline*deadline-1 < 1
    for s in (Q(3, 5), Q(4, 5), Q(99, 100)):
        # t=1/(1-s^2) is finite and gives q=t^2(1-s^2)>=1.
        t = 1/(1-s*s)
        assert menu_lorentz([t, s*t, Q(0), Q(0)])
    groups['positive_proper_latency_menu'] = dict(status='PASS', c='1', delta='1',
        deadline='5/4', reachable_radius_squared='9/16',
        excluded_target=['1', '0', '0'], H_holds=False, FD_holds=True,
        scope='Lorentz-invariant permission menu; no physical implementation or minimum-time axiom')

    def resource_menu(z, n):
        return z[0] >= 1 and norm2(z[1:]) <= n*n*z[0]*z[0]

    resource_rows = []
    for n in (1, 2, 4, 16):
        z = [Q(1), Q(n+1), Q(0), Q(0)]
        assert not resource_menu(z, n) and resource_menu(z, n+1)
        resource_rows.append(dict(budget=n, outside_target=n+1,
                                  sufficient_larger_budget=n+1))
    # T<delta has no completed task, so the empty-window modulus can be zero.
    assert all(not resource_menu([Q(1, 2), Q(0), Q(0), Q(0)], n)
               for n in (1, 2, 4, 16))
    groups['resource_union_and_empty_window'] = dict(status='PASS', cases=resource_rows,
        analytic_union='t>=1, arbitrary x',
        per_budget_FD=True, union_FD=False, short_empty_window_modulus='0',
        scope='Quantifier diagnostic; not an internal-resource physical countermodel')

    # Frozen 1093 permission example: FD alone does not select a finite speed.
    nonlinear_rows = []
    for t in (Q(1), Q(2), Q(10), Q(100)):
        radius = t+t*t
        assert radius/t == 1+t
        nonlinear_rows.append(dict(deadline=str(t), max_radius=str(radius),
                                   boundary_average_speed=str(radius/t)))
    groups['missing_actual_motion'] = dict(status='PASS', cases=nonlinear_rows,
        FD=True, finite_global_speed=False, actual_references='spatial rotations only',
        source='1093 permission menu, reused; O fails')

    g = scale(boost([Q(1, 2), Q(0), Q(0)]), Q(2))
    word = commutator(g, rot)
    assert mv(word, rest) == [Q(1), Q(1), Q(0), Q(0)]
    time_upper = Q(1)+Q(1, 10)
    assert time_upper < deadline
    assert Q(1, 10) < Q(1, 5)
    total_error = Q(1, 100)+4*Q(1, 100)
    assert total_error == Q(1, 20) < Q(1, 10)
    anchor_p = relation_probability(Q(0))
    target_gap = relation_probability(Q(1))-anchor_p
    nearby_gap = relation_probability(Q(9, 10))-anchor_p
    assert anchor_p == Q(1, 2) and target_gap == Q(1, 4)
    assert nearby_gap == Q(9, 38) > Q(1, 10)
    groups['finite_error_and_original_relation_witness'] = dict(status='PASS',
        task_time='1', deadline=str(deadline), time_upper=str(time_upper),
        time_margin=str(deadline-time_upper), location_error='1/10', excluded_radius='1/5',
        total_complete_error=str(total_error), error_limit='1/10',
        target_relation_gap=str(target_gap), near_target_relation_gap=str(nearby_gap),
        assumed_all_outputs_relation_gap_bound='1/10',
        relation_margin=str(nearby_gap-Q(1, 10)),
        uniform_nonexpansive_complete_transport_assumed=True,
        exclusion_is_additional_contract=True, real_instrument_data=False)
    assert len(groups) == 5
    return dict(round=1100, status='PASS', groups=groups,
        theorem='FD iff positive finite task speed, conditional on actual R/O and affine positive time',
        new_adopted_axioms=0, scientific_count_increment=0,
        FD_derived_from_cognition=False, unconditional_Lorentz_derived=False,
        previous_orbit_and_cone_proofs_reused=True)


def main():
    data = compute()
    target = HERE/'results.json'
    if '--write' in sys.argv:
        target.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
    else:
        assert json.loads(target.read_text(encoding='utf8')) == data, 'Saved results differ.'
    print(json.dumps(dict(round=1100, status='PASS', groups=len(data['groups']),
        exact_arithmetic=True, saved_results_match=True, conditional_only=True), ensure_ascii=False))


if __name__ == '__main__':
    main()
