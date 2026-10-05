"""779: stationary auxiliary contact prescription, exact finite diagnostics.

Tree enumeration and formal stationary series are computed independently.
Neither is a continuum measure or a numerical computation of the original QAP.
"""
from fractions import Fraction as Q
from itertools import product
from math import factorial
from pathlib import Path
import argparse
import json

DIM = 4  # e0, e1, propagating spectator x, independent insertion source j
ZERO = (0,)*DIM
ONE = {ZERO: Q(1)}


def add(*values):
    out = {}
    for value in values:
        for key, coeff in value.items():
            out[key] = out.get(key, Q(0))+coeff
    return {key: coeff for key, coeff in out.items() if coeff}


def scale(value, coeff):
    return {key: coeff*c for key, c in value.items() if coeff*c}


def mul(a, b):
    out = {}
    for key, coeff in a.items():
        for other, value in b.items():
            exponent = tuple(i+j for i, j in zip(key, other))
            out[exponent] = out.get(exponent, Q(0))+coeff*value
    return {key: coeff for key, coeff in out.items() if coeff}


def power(value, exponent):
    out = ONE
    for _ in range(exponent):
        out = mul(out, value)
    return out


def variable(axis):
    key = list(ZERO)
    key[axis] = 1
    return {tuple(key): Q(1)}


def derivative(value, axis):
    out = {}
    for key, coeff in value.items():
        if key[axis]:
            exponent = list(key)
            exponent[axis] -= 1
            out[tuple(exponent)] = key[axis]*coeff
    return out


def multiply_series(a, b, order):
    out = [{} for _ in range(order+1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i+j <= order:
                out[i+j] = add(out[i+j], mul(x, y))
    return out


def substitute_series(value, shifts, order):
    replacements = {axis: [add(variable(axis), series[0])]+series[1:]
                    for axis, series in shifts.items()}
    out = [{} for _ in range(order+1)]
    for exponent, coeff in value.items():
        term = [scale(ONE, coeff)]+[{} for _ in range(order)]
        for axis, degree in enumerate(exponent):
            replacement = replacements.get(axis, [variable(axis)])
            for _ in range(degree):
                term = multiply_series(term, replacement, order)
        out = [add(a, b) for a, b in zip(out, term)]
    return out


def stationary(vertex, masses, order):
    """For V=t*vertex, solve N u=V'(e+u), return u and R through order."""
    u = {axis: [{} for _ in range(order+1)] for axis in masses}
    gradients = {axis: derivative(vertex, axis) for axis in masses}
    for degree in range(1, order+1):
        # All new coefficients use only lower orders, simultaneously.
        updates = {axis: scale(substitute_series(gradients[axis], u, degree-1)[degree-1], 1/n)
                   for axis, n in masses.items()}
        for axis in masses:
            u[axis][degree] = updates[axis]
    evaluated = substitute_series(vertex, u, order)
    result = [{}]+evaluated[:-1]
    for axis, n in masses.items():
        square = multiply_series(u[axis], u[axis], order)
        result = [add(a, scale(b, -n/2)) for a, b in zip(result, square)]
    return u, result


def show(value):
    return [{'powers': list(k), 'coefficient': str(v)} for k, v in sorted(value.items())]


def tree_comparison():
    e = variable(0)
    vertex = add(scale(e, Q(2, 5)), scale(power(e, 2), Q(-1, 6)),
                 scale(power(e, 3), Q(1, 9)), scale(power(e, 4), Q(1, 24)))
    n_aux, max_order = Q(3, 2), 6
    u, series = stationary(vertex, {0: n_aux}, max_order)
    derivatives = [vertex]
    for _ in range(max_order):
        derivatives.append(derivative(derivatives[-1], 0))
    cases = []
    for count in range(1, max_order+1):
        tree_sum, total = {}, 0
        if count == 1:
            tree_sum, total = vertex, 1
        else:
            for prufer in product(range(count), repeat=count-2):
                degrees = [1+prufer.count(i) for i in range(count)]
                term = ONE
                for degree in degrees:
                    term = mul(term, derivatives[degree])
                tree_sum = add(tree_sum, term)
                total += 1
            tree_sum = scale(tree_sum, Q(1, factorial(count))/n_aux**(count-1))
        assert tree_sum == series[count]
        cases.append(dict(vertices=count, labelled_trees=total,
                          local_polynomial_terms=len(tree_sum)))
    return dict(compared_orders=max_order, independently_enumerated_trees=sum(x['labelled_trees'] for x in cases),
                all_tree_coefficients_exact=True, cases=cases,
                second_order_coefficient=show(series[2]))


def envelope_checks():
    e0, e1, x, j = [variable(i) for i in range(DIM)]
    vertex = add(scale(mul(x, power(add(e1, scale(e0, -1)), 2)), Q(1, 2)),
                 scale(power(e0, 3), Q(1, 6)), mul(j, e0))
    masses, order = {0: Q(2), 1: Q(3)}, 5
    u, series = stationary(vertex, masses, order)
    checks = 0
    for axis in range(DIM):
        expected = [{}]+substitute_series(derivative(vertex, axis), u, order)[:-1]
        for degree in range(order+1):
            assert derivative(series[degree], axis) == expected[degree]
            checks += 1
    for axis, n in masses.items():
        for degree in range(1, order+1):
            assert derivative(series[degree], axis) == scale(u[axis][degree], n)
            checks += 1
    # Derivative with respect to the independent linear e0 source j.
    expected_j = [{}]+[e0]+u[0][1:order]
    for degree in range(order+1):
        assert derivative(series[degree], 3) == expected_j[degree]
        checks += 1
    # Omitting the stationary cost violates the envelope identity at order two.
    wrong = [{}]+substitute_series(vertex, u, order)[:-1]
    expected = [{}]+substitute_series(derivative(vertex, 2), u, order)[:-1]
    defect = add(derivative(wrong[2], 2), scale(expected[2], -1))
    assert defect
    return dict(order=order, sites=2, background_weights=['2','3'],
                derivative_source_checks=checks, includes_difference_stencil=True,
                stationary_cost_omission_defect=show(defect),
                all_envelope_and_basic_insertion_identities_exact=True)


def mixed_menu_and_loop_choice():
    e, x, j = variable(0), variable(2), variable(3)
    n = Q(3, 2)
    _, mixed = stationary(mul(mul(j, e), x), {0: n}, 4)
    assert mixed[1] == mul(mul(j, e), x)
    assert mixed[2] == scale(mul(power(j, 2), power(x, 2)), 1/(2*n))
    assert not any(mixed[3:])
    _, pure = stationary(add(power(x, 3), mul(j, x)), {0: n}, 4)
    assert not any(pure[2:])
    _, quadratic_aux = stationary(scale(power(e, 2), Q(1, 2)), {0: n}, 4)
    assert all(not any(key[0] == 0 for key in p) for p in quadratic_aux)
    # A finite ordinary Gaussian DOES have a two-edge loop at two e^2 vertices.
    # Our prescription is a choice of continuum contact normalization, not that measure.
    finite_gaussian_double_edge = Q(2)/n**2
    assert finite_gaussian_double_edge == Q(8, 9)
    return dict(mixed_e_times_field_second_contact=show(mixed[2]),
                arbitrary_physical_vertex_unchanged=True,
                pure_auxiliary_closed_cycle_coefficient_chosen_zero=True,
                finite_gaussian_double_edge_coefficient=str(finite_gaussian_double_edge),
                equals_unsubtracted_finite_gaussian=False,
                no_global_or_interacting_positive_measure_claim=True)


def run():
    return dict(round=779, groups={
        'independent_tree_enumeration': tree_comparison(),
        'stationary_envelope_and_linear_sources': envelope_checks(),
        'mixed_menu_and_loop_normalization': mixed_menu_and_loop_choice()},
        all_checks_passed=True, finite_tests_prove_continuum_prescription=False,
        dynamic_frame_quantum_map_proven=False, original_full_N1_N2_proven=False,
        original_loop_anomaly_computed=False, actual_instrument_proven=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    path = Path(__file__).with_name('auxiliary_tree_normalization_results.json')
    if args.write:
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    if args.check:
        assert result == json.loads(path.read_text(encoding='utf-8'))
    print(json.dumps(result, ensure_ascii=False, indent=2))
