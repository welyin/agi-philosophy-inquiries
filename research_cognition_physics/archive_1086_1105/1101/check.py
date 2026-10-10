"""1101 exact finite certificates; reuses frozen 1099 Fraction arithmetic.

These checks do not establish actual reference or task permissions.
Universal claims and hypotheses belong to proof.md. Default invocation is
read-only and compares saved results; --write explicitly saves this file's result.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
import runpy
import sys

HERE = Path(__file__).resolve().parent
M = runpy.run_path(str(HERE.parent / "1099/check.py"),
                   run_name="frozen_round1099")
eye, mul, mv = (M[k] for k in ("eye", "mul", "mv"))
inverse, power, scale = (M[k] for k in ("inverse", "power", "scale"))
rotate, commutator, norm2 = (M[k] for k in ("rotate", "commutator", "norm2"))


def transpose(a):
    return list(map(list, zip(*a)))


def quadratic(z):
    return z[0] * z[0] - norm2(z[1:])


def string_matrix(a):
    return [list(map(str, row)) for row in a]


def det3(a):
    return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))


def canonical_boost(z):
    gamma, u = z[0], z[1:]
    assert gamma >= 1 and quadratic(z) == 1
    b = eye()
    b[0][0] = gamma
    for i in range(3):
        b[0][i + 1] = b[i + 1][0] = u[i]
        for j in range(3):
            b[i + 1][j + 1] += u[i] * u[j] / (gamma + 1)
    return b


def compute():
    groups = {}
    eta = eye()
    for i in range(1, 4):
        eta[i][i] = Q(-1)
    e0 = [Q(1), Q(0), Q(0), Q(0)]
    lor = eye()
    lor[0][0] = lor[1][1] = Q(5, 4)
    lor[0][1] = lor[1][0] = Q(3, 4)
    assert mul(mul(transpose(lor), eta), lor) == eta
    flip = rotate(Q(-1), Q(0))
    expected = eye()
    expected[0][0] = expected[1][1] = Q(17, 8)
    expected[0][1] = expected[1][0] = Q(15, 8)
    cases = []
    for lam in (Q(1, 3), Q(1), Q(2)):
        g = scale(lor, lam)
        word = commutator(g, flip)
        assert word == expected == power(lor, 2)
        assert mul(mul(transpose(word), eta), word) == eta
        cases.append(dict(initial_scale=str(lam),
                          result_cosh=str(word[0][0]),
                          result_sinh=str(word[0][1]),
                          result_scale_squared="1",
                          reference_factors=4))
    groups["conjugated_flip_cancels_scale"] = dict(
        status="PASS", cases=cases, exact_word="g R g^-1 R^-1",
        scope="Finite matrix identity conditional on actual inverses and rotations.",
        all_reference_scales_eliminated_by_this_identity=False)

    quarter = rotate(Q(0), Q(1))
    target_q = [Q(25, 16), Q(15, 16), Q(-12, 16), Q(0)]
    b = canonical_boost(target_q)
    second_cases = []
    for lam in (Q(1, 3), Q(1), Q(2)):
        g = scale(lor, lam)
        h = mul(mul(g, quarter), inverse(g))
        assert mv(h, e0) == target_q
        k = mul(inverse(b), h)
        assert mv(k, e0) == e0 and k[0] == e0
        spatial_k = [row[1:] for row in k[1:]]
        assert mul(transpose(spatial_k), spatial_k) == eye(3)
        assert det3(spatial_k) == 1
        assert mul(h, inverse(k)) == b
        assert mul(mul(transpose(b), eta), b) == eta
        assert mv(b, e0) == target_q
        second_cases.append(dict(initial_scale=str(lam),
                                 rotation_determinant="1",
                                 reference_word_factors=4))
    groups["canonical_boost_from_actual_rotation_word"] = dict(
        status="PASS", cases=second_cases,
        h_e0=list(map(str, target_q)), canonical_boost=string_matrix(b),
        spatial_rotation=string_matrix(spatial_k),
        identity="h k^-1 = B(h e0)",
        scope=("k is certified algebraically to be SO(3); its actual availability "
               "uses O. B is realized conditionally as the finite word, "
               "not added as a free boost permission."),
        universal_boost_family_certified_by_samples=False)

    # Choose the same elementary boost lor as B in A=(1/2)B.
    a = scale(lor, Q(1, 2))
    start = [Q(2), Q(0), Q(0), Q(0)]
    end = mv(power(a, 2), start)
    assert end == [Q(17, 16), Q(15, 16), Q(0), Q(0)]
    assert quadratic(start) == 4 and quadratic(end) == Q(1, 4)
    d = Q(1)
    coordinate_error = Q(1, 100)
    # For every |e_mu|<=eps, |q(z+e)-q(z)| <= 2 eps ||z||_1+4 eps^2.
    q_error = (2 * coordinate_error * sum(map(abs, end), Q(0))
               + 4 * coordinate_error ** 2)
    q_lower, q_upper = quadratic(end) - q_error, quadratic(end) + q_error
    assert q_error == Q(101, 2500)
    assert q_lower == Q(131, 625) > 0
    assert q_upper == Q(363, 1250) < d * d
    assert end[0] - coordinate_error == Q(421, 400) > 0
    groups["finite_scale_contradiction_with_coordinate_margin"] = dict(
        status="PASS", lower_bound_d="1", initial_tau="2",
        scale="1/2", iterations=2, result_tau="1/2",
        result_coordinates=list(map(str, end)),
        assumed_coordinate_error_per_component=str(coordinate_error),
        quadratic_error_bound=str(q_error), quadratic_lower=str(q_lower),
        quadratic_upper=str(q_upper),
        strict_quadratic_margin=str(d * d - q_upper),
        future_time_lower=str(end[0] - coordinate_error),
        scope=("Conditional finite certificate for an already qualified complete "
               "task. Coordinate errors are assumed certificate bounds, not "
               "measured instrument data or proof of full-process fidelity."),
        whole_menu_positive_gap_is_additional=True,
        actual_reference_word_and_complete_task_qualification_assumed=True,
        real_instrument_data=False)

    def budget_menu(z, n):
        assert isinstance(n, int) and n >= 1
        return z[0] > 0 and quadratic(z) >= Q(1, n * n)

    union_cases = []
    for n in (1, 2, 4, 16, 100):
        z = mv(lor, [Q(1, n), Q(0), Q(0), Q(0)])
        doubled = [2 * a for a in z]
        halved = [a / 2 for a in z]
        assert budget_menu(z, n)
        assert budget_menu(doubled, n)
        assert not budget_menu(halved, n) and budget_menu(halved, 2 * n)
        assert quadratic(z) == Q(1, n * n)
        assert quadratic(doubled) == 4 * quadratic(z)
        assert quadratic(halved) == quadratic(z) / 4
        union_cases.append(dict(budget=n, threshold_squared=str(Q(1, n * n)),
                                inverse_scale_needed_budget=2 * n,
                                finite_coordinates=list(map(str, z))))
    deadline = Q(1)
    successful = [Q(1, 2), Q(0), Q(0), Q(0)]
    excluded = [deadline, Q(2), Q(0), Q(0)]
    assert budget_menu(successful, 2)
    assert quadratic(excluded) == -3
    assert all(not budget_menu(excluded, n) for n in (1, 2, 4, 16, 100))
    groups["positive_budget_gaps_disappear_under_union"] = dict(
        status="PASS", cases=union_cases,
        analytic_menu="P_n = {t>0, q(t,x)>=1/n^2}",
        analytic_union="P = {t>0, q(t,x)>0}",
        union_infimum_tau="0", union_admits_scale_and_inverse=["2", "1/2"],
        FD=dict(deadline="1", completed_example=list(map(str, successful)),
                excluded_target=["2", "0", "0"],
                all_deadline_outputs_have_norm_less_than="1",
                analytic_exclusion="For t<=1 at x=(2,0,0), q<=-3."),
        scope=("Finite samples check analytic witnesses only. The union and "
               "all-budget statements require the proof; these are permission "
               "menus, not a full FUCP+CO+six-protocol physical countermodel."))
    assert len(groups) == 4
    return dict(
        round=1101, status="PASS", groups=groups,
        exact_arithmetic="fractions.Fraction",
        reused_arithmetic="../1099/check.py (runpy, main not executed)",
        universal_theorem="See proof.md; finite matrix tests do not prove permissions.",
        new_adopted_axioms=0, scientific_count_increment=0,
        positive_task_gap_derived_from_cognition=False,
        full_actual_clock_proper_time_identification=False,
        all_matter_Lorentz_certified=False,
        unconditional_Lorentz_derived=False)


def main():
    result = compute()
    target = HERE / "results.json"
    if "--write" in sys.argv:
        target.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n",
                          encoding="utf8")
    else:
        saved = json.loads(target.read_text(encoding="utf8"))
        assert saved == result, "Saved results differ."
    print(json.dumps(dict(round=1101, status="PASS", groups=len(result["groups"]),
                          exact_arithmetic=True, saved_results_match=True,
                          conditional_only=True), ensure_ascii=False))


if __name__ == "__main__":
    main()
