"""783: exact checks of the local, formal all-order repair interface.

The finite BV Laplacian below is a consistency diagnostic, NOT a replacement
for the renormalized continuum anomaly or a calculation of its coefficients.
Reuses 777's exact polynomial algebra. Only Python standard library is needed.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from pathlib import Path
import argparse
import importlib.util
import json
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('bv777', HERE.parent/'777/joint_auxiliary_bv_reduction.py')
bv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bv)
add, mul, scale, diff, bracket = bv.add, bv.mul, bv.scale, bv.diff, bv.bracket
V = {name: bv.var(name) for name in bv.NAMES}


def lap(poly):
    return add(*(scale(diff(diff(poly, anti), field), -1 if field in bv.ODD else 1)
                 for field, anti in bv.PAIRS))


def homotopy(poly):
    # q,c and c*,q* are the contractible pairs. z denotes c* in 777.
    out = {}
    for key, coefficient in poly.items():
        degree = sum(key[bv.NAMES.index(x)] for x in ('q', 'c', 'p', 'z'))
        if degree:
            term = {key: coefficient}
            out = add(out, scale(add(mul(V['q'], diff(term, 'c')),
                                     mul(V['z'], diff(term, 'p'))), Q(1, degree)))
    return out


def defect_at(actions, order):
    terms = [scale(bracket(actions[k], actions[order-k]), Q(1, 2))
             for k in range(order+1) if k < len(actions) and order-k < len(actions)]
    if order and order-1 < len(actions):
        terms.append(scale(lap(actions[order-1]), -1))
    return add(*terms)


def all_order_source_diagnostic():
    q, r, p, t, c = [V[x] for x in ('q', 'r', 'p', 't', 'c')]
    # The unused even symbol b is ONLY an external source zeta in this check.
    # Its partner rho never occurs; b is not an auxiliary quantum field here.
    zeta = V['b']
    classical = add(scale(mul(r, r), Q(1, 2)), mul(p, c), mul(zeta, r))
    generator = add(bv.prod(p, q, q, r), bv.prod(t, q, r, r),
                    bv.prod(zeta, p, q, r), scale(bv.prod(zeta, t, r, r), Q(1, 2)))
    first = add(bracket(classical, generator), scale(r, Q(7, 11)))
    assert not bracket(classical, classical) and not lap(classical)
    assert not bracket(classical, first)
    assert lap(first), 'Choose a diagnostic with a real quantum contact term.'
    actions = [classical, first]
    two_loop_defect = defect_at(actions, 2)
    assert two_loop_defect and bracket(first, first)
    records = []
    for order in range(2, 7):
        residual = defect_at(actions, order)
        assert residual
        assert not bracket(classical, residual), ('consistency', order)
        correction = scale(homotopy(residual), -1)
        assert not add(bracket(classical, correction), residual), ('primitive', order)
        actions.append(correction)
        assert not defect_at(actions, order), ('QME', order)
        contacts = diff(diff(correction, 'b'), 'b')
        records.append(dict(order=order, defect_terms=len(residual),
                            correction_terms=len(correction),
                            second_source_derivative_terms=len(contacts),
                            coefficient_equation_exact=True))
    for order in range(7):
        assert not defect_at(actions, order)
    # All contact derivatives are obtained from one action, with its partners.
    derivative_checks = 0
    for correction in actions[1:]:
        for name in ('q', 'r', 'p', 't', 'c'):
            assert diff(diff(correction, 'b'), name) == diff(diff(correction, name), 'b')
            derivative_checks += 1
    assert any(row['second_source_derivative_terms'] for row in records)
    # Extract BEFORE setting fields to zero.
    physical_first = bv.substitute(diff(first, 'r'), {x: {} for x in bv.NAMES})
    assert physical_first == scale(bv.ONE, Q(7, 11))
    return dict(qme_coefficients_checked_through_order=6,
                equation='1/2 (S_epsilon,S_epsilon) - epsilon Delta S_epsilon = 0',
                epsilon_is_a_formal_algebraic_loop_symbol_not_physical_hbar=True,
                first_physical_linear_coefficient='7/11', first_correction_never_replaced=True,
                unrepaired_second_order_defect=bv.display(two_loop_defect),
                quantum_laplacian_first_correction=bv.display(lap(first)),
                second_order_correction=bv.display(actions[2]),
                order_checks=records, mixed_source_derivatives_checked=derivative_checks,
                continuum_anomaly_coefficients_computed=False,
                source_symbol='b in inherited engine denotes external zeta, not a quantum field')


def filtration_checks():
    q, r, p, c, t = [V[x] for x in ('q', 'r', 'p', 'c', 't')]
    free = add(scale(mul(r, r), Q(1, 2)), mul(p, c))
    s0 = lambda f: bracket(free, f)
    checks = 0
    names = ('q', 'r', 'p', 't', 'c', 'z')
    for degree in range(5):
        for word in combinations_with_replacement(names, degree):
            f = bv.prod(*(V[x] for x in word))
            if not f:
                continue
            projected = bv.substitute(f, {x: {} for x in ('q', 'c', 'p', 'z')})
            assert not add(s0(homotopy(f)), homotopy(s0(f)), scale(f, -1), projected)
            assert not lap(lap(f))
            checks += 1
    # Enumerate valences and loop labels, not individual Feynman amplitudes.
    types = [(0, d) for d in range(3, 7)] + [(ell, d) for ell in (1, 2) for d in range(1, 5)]
    graph_counts = 0
    maximum_vertices = 0
    for count in range(1, 6):
        for indices in combinations_with_replacement(range(len(types)), count):
            vertices = [types[i] for i in indices]
            ends = sum(d for ell, d in vertices)
            inserted_loops = sum(ell for ell, d in vertices)
            for internal in range(count-1, ends//2+1):
                external = ends-2*internal
                total_loop = internal-count+1+inserted_loops
                if external < 1 or total_loop > 4:
                    continue
                excess = sum(d+2*ell-2 for ell, d in vertices)
                assert excess == external+2*total_loop-2
                assert count <= excess
                graph_counts += 1
                maximum_vertices = max(maximum_vertices, count)
    # Differential contractions lower field order: fixed-degree truncation is unsafe.
    omitted = bv.prod(p, q, q, q)
    assert lap(omitted) == scale(mul(q, q), 3)
    truncated = {key: value for key, value in omitted.items() if sum(key) <= 2}
    assert not lap(truncated) and lap(omitted)
    # Finite Taylor truncations of genuinely positive density matrices can be negative.
    eps = Q(1, 10)
    exact_det = (1/(1+eps**2))*(eps**2/(1+eps**2))-(eps/(1+eps**2))**2
    first_order_det = -eps**2
    assert exact_det == 0 and first_order_det < 0
    return dict(homotopy_and_laplacian_monomials=checks,
                vertex_valence_loop_count_cases=graph_counts,
                maximum_vertices_enumerated=maximum_vertices,
                identity='sum(d_v+2*r_v-2) = external_legs+2*total_loops-2',
                types_have_strictly_positive_excess=True,
                omitted_degree_four_contact_at_degree_two=bv.display(lap(omitted)),
                exact_positive_rank_one_density_determinant=str(exact_det),
                first_order_density_determinant_at_epsilon_one_tenth=str(first_order_det),
                enumeration_is_not_a_continuum_Feynman_integral=True)


def run():
    return dict(round=783, test_groups=2, exact_arithmetic=True,
                all_order_source_repair=all_order_source_diagnostic(),
                completion_and_scope=filtration_checks(),
                original_physical_anomaly_computed=False,
                interacting_charge_or_positive_state_constructed=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    target = HERE/'all_order_local_repair_results.json'
    if args.write:
        target.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    if args.check:
        assert result == json.loads(target.read_text(encoding='utf-8'))
    print(json.dumps(result, ensure_ascii=False, indent=2))
