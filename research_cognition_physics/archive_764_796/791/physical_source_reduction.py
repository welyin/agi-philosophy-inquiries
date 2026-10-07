"""791: fixed-Wick physical tadpole reduction, with all gauge partners.

Exact rational diagnostics; finite differences here test jet and grading
identities, not the existence of the original continuum Hadamard state.
"""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import json
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('bv791', HERE.parent/'777/joint_auxiliary_bv_reduction.py')
bv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bv)
N = 3
bv.EVEN = tuple(f'{name}{i}' for name in ('w', 'q', 'b') for i in range(N))
bv.ODD = tuple(f'{name}{i}' for name in ('c', 'h', 'psi', 'barpsi') for i in range(N))
bv.NAMES = bv.EVEN+bv.ODD
bv.ZERO = (0,)*len(bv.NAMES)
bv.ONE = {bv.ZERO: Q(1)}
V = {name: bv.var(name) for name in bv.NAMES}
DIFF = [[Q(-1), Q(1), Q(0)], [Q(0), Q(-1), Q(1)], [Q(1), Q(0), Q(-1)]]


def gamma(poly):
    return bv.add(*(bv.add(bv.mul(V[f'c{i}'], bv.diff(poly, f'q{i}')),
                           bv.mul(V[f'b{i}'], bv.diff(poly, f'h{i}'))) for i in range(N)))


def wick(poly, kernel):
    return bv.add(*(bv.scale(bv.diff(bv.diff(poly, a), b), value/2)
                    for (a, b), value in kernel.items()))


def put(kernel, a, b, value):
    value = Q(value)
    kernel[(a, b)] = value
    kernel[(b, a)] = -value if a in bv.ODD and b in bv.ODD else value


def kernel(case, physical_scale=Q(1)):
    out = {}
    for i in range(N):
        for j in range(N):
            put(out, f'w{i}', f'w{j}', physical_scale*(Q(2) if i == j else Q(1, 3)))
            put(out, f'q{i}', f'q{j}', Q(case+2, i+j+2))
            put(out, f'q{i}', f'w{j}', Q((i+1)*(case+1), j+4))
            # Ordered, generally nonsymmetric mixed blocks. gamma(h_i q_j)
            # fixes the matching ghost covariance including its slot order.
            cross = Q((i+1)*(case+2)+j, i+j+3)
            put(out, f'q{j}', f'b{i}', cross)
            put(out, f'h{i}', f'c{j}', cross)
            put(out, f'barpsi{i}', f'psi{j}', Q(3, 7) if i == j else Q(1, 11))
    return out


def fields_and_action():
    dw = [bv.add(*(bv.scale(V[f'w{j}'], DIFF[i][j]) for j in range(N))) for i in range(N)]
    dq = [bv.add(*(bv.scale(V[f'q{j}'], DIFF[i][j]) for j in range(N))) for i in range(N)]
    psi_terms, red_terms = [], []
    for i in range(N):
        w, q, b, h = [V[f'{name}{i}'] for name in ('w', 'q', 'b', 'h')]
        gauge_function = bv.add(bv.scale(bv.mul(dq[i], dw[i]), Q(2, 3)),
                                bv.scale(bv.mul(w, w), Q(5, 14)),
                                bv.scale(bv.mul(q, b), Q(3, 5)),
                                bv.scale(bv.mul(w, b), Q(-2, 7)),
                                bv.scale(bv.mul(q, q), Q(4, 9)))
        psi_terms.append(bv.mul(h, gauge_function))
        red_terms += [bv.scale(bv.prod(w, w, w), Q(5, 42)),
                      bv.scale(bv.prod(w, dw[i], dw[i]), Q(2, 9)),
                      bv.scale(bv.prod(w, V[f'barpsi{i}'], V[f'psi{i}']), Q(7, 10))]
    return bv.add(*psi_terms), bv.add(*red_terms)


def gauge_fiber_tadpole():
    psi, red = fields_and_action()
    gf = gamma(psi)
    assert not gamma(gf) and not gamma(red)
    rows = []
    source_rows = []
    for case in (0, 1, 3):
        c = kernel(case)
        # Every quadratic monomial: sufficient to test covariance of an even
        # constant second-order contraction and a linear field differential.
        checks = 0
        for i, a in enumerate(bv.NAMES):
            for b in bv.NAMES[i:]:
                mono = bv.mul(V[a], V[b])
                if not mono:
                    continue
                assert wick(gamma(mono), c) == gamma(wick(mono, c))
                checks += 1
        assert wick(gf, c) == gamma(wick(psi, c))
        gauge_linear = wick(gf, c)
        assert gauge_linear  # The whole gauge tadpole is NOT identically zero.
        assert all(not bv.diff(gauge_linear, f'w{i}') for i in range(N))
        assert all(not bv.diff(gauge_linear, f'q{i}') for i in range(N))
        physical = [wick(bv.diff(red, f'w{i}'), c) for i in range(N)]
        full = [wick(bv.diff(bv.add(red, gf), f'w{i}'), c) for i in range(N)]
        assert full == physical and any(full)
        missing_ghost = {key: value for key, value in c.items() if not any(x.startswith(('h', 'c')) for x in key)}
        wrong = [wick(bv.diff(gf, f'w{i}'), missing_ghost) for i in range(N)]
        assert any(wrong)
        no_fermion = {key: value for key, value in c.items() if not any(x.startswith(('psi', 'barpsi')) for x in key)}
        bad_physical = [wick(bv.diff(red, f'w{i}'), no_fermion) for i in range(N)]
        assert physical != bad_physical
        rows.append(dict(case=case, covariance_checks=checks,
                         physical_raw_source=[bv.display(x) for x in physical],
                         surviving_auxiliary_source=bv.display(gauge_linear),
                         omitted_ghost_error=[bv.display(x) for x in wrong],
                         omitted_fermion_error=[bv.display(bv.add(a, bv.scale(b, -1))) for a, b in zip(bad_physical, physical)]))
        source_rows.append(physical)
    assert source_rows[0] == source_rows[1] == source_rows[2]
    other_state = kernel(0, Q(7, 6))
    new_physical = [wick(bv.diff(red, f'w{i}'), other_state) for i in range(N)]
    assert new_physical != source_rows[0]
    ell = bv.add(*(bv.scale(V[f'w{i}'], Q(i+1, 13)) for i in range(N)),
                 *(bv.scale(V[f'b{i}'], Q(i+2, 17)) for i in range(N)))
    assert not gamma(ell)
    wrong_ell = bv.add(ell, bv.scale(V['q1'], Q(2, 5)))
    assert gamma(wrong_ell) == bv.scale(V['c1'], Q(2, 5))
    return dict(cases=rows, same_physical_contractions_same_source=True,
                physical_state_dependence_retained=True,
                field_independent_source_repair_preserved=True,
                bad_q_counterterm_defect=bv.display(gamma(wrong_ell)),
                same_field_BRST_not_full_BV=True, continuum_state_simulated=False)


def solve(matrix, rhs):
    a = [[Q(x) for x in row]+[Q(y)] for row, y in zip(matrix, rhs)]
    size = len(a)
    for i in range(size):
        pivot = next(j for j in range(i, size) if a[j][i])
        a[i], a[pivot] = a[pivot], a[i]
        d = a[i][i]
        a[i] = [x/d for x in a[i]]
        for j in range(size):
            if j != i:
                d = a[j][i]
                a[j] = [x-d*y for x, y in zip(a[j], a[i])]
    return [row[-1] for row in a]


def source_response():
    # This repeats no claim of a new Green theorem: 785/787 already supply it.
    # The check connects the NEW no-q source condition to the mixed response.
    p = [[Q(3), Q(1, 2)], [Q(1, 2), Q(2)]]
    j_w = [Q(2, 3), Q(-5, 7)]
    expected = solve(p, [-x for x in j_w])
    rows = []
    for case in (1, 2, 4):
        a = [[Q(case+1), Q(1, 3)], [Q(-2, 5), Q(case+2)]]
        b = [[Q(2, 7), Q(-3, 5)], [Q(4, 9), Q(case, 6)]]
        n = [[Q(case, 3), Q(1, 4)], [Q(1, 4), Q(-2, 7)]]
        matrix = [[Q(0) for _ in range(6)] for _ in range(6)]
        for i in range(2):
            for j in range(2):
                matrix[i][j] = p[i][j]
                matrix[i][j+4] = b[j][i]
                matrix[i+4][j] = b[i][j]
                matrix[i+2][j+4] = a[j][i]
                matrix[i+4][j+2] = a[i][j]
                matrix[i+4][j+4] = n[i][j]
        j_b = [Q(case, 5), Q(-case, 11)]
        full = solve(matrix, [-x for x in j_w+[Q(0), Q(0)]+j_b])
        assert full[:2] == expected and full[4:] == [Q(0), Q(0)]
        wrong = solve(matrix, [-x for x in j_w+[Q(1, 3), Q(0)]+j_b])
        assert wrong[:2] != expected
        rows.append(dict(case=case, physical_mean=[str(x) for x in full[:2]],
                         auxiliary_mean=[str(x) for x in full[4:]],
                         illegal_q_source_error=[str(x-y) for x, y in zip(wrong[:2], expected)]))
    return dict(cases=rows, inherited_free_inverse_used=True,
                no_gauge_endpoint_limit=True, arithmetic='Fraction',
                scope='Finite mixed-block response to a Ward-closed linear source; causal continuum kernel and compatible initial data are inherited, not numerically simulated.')


def run():
    return dict(round=791, fresh_test_groups=2,
                fixed_wick_gauge_fiber=gauge_fiber_tadpole(), source_response=source_response(),
                all_checks_passed=True,
                original_one_loop_physical_source_reduction_proven=True,
                original_physical_finite_counterterm_retained=True,
                full_interacting_gauge_fiber_decoupling_proven=False,
                original_interacting_positive_state_proven=False,
                all_order_source_dictionary_computed=False)


if __name__ == '__main__':
    results = run()
    HERE.joinpath('physical_source_reduction_results.json').write_text(
        json.dumps(results, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(results, ensure_ascii=False, indent=2))
