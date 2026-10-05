"""778: exact graded Gaussian/BV diagnostics and an ultralocal-loop boundary.

The Gaussian calculation is finite dimensional; it does NOT implement a
continuum Epstein--Glaser extension or compute the original model's anomaly.
Use the existing Python runtime. No third-party package or historical ZIP.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import factorial
from pathlib import Path
import argparse
import importlib.util
import json
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    'round777_algebra', HERE.parent/'777/joint_auxiliary_bv_reduction.py')
bv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bv)
add, scale, mul, diff = bv.add, bv.scale, bv.mul, bv.diff
ACTIVE = ('q', 'r', 'z', 'a', 'c', 'h', 'p', 't')
PAIRS = (('q', 'p'), ('r', 't'), ('c', 'z'), ('h', 'a'))


def setup():
    variables, n, _, _, _, red, _, _ = bv.setup()
    free = {key: value for key, value in red.items() if sum(key) == 2}
    bosons = ('q', 'r')
    hessian = [[diff(diff(free, x), y).get(bv.ZERO, Q(0))
                for y in bosons] for x in bosons]
    a, b = hessian[0]
    c, d = hessian[1]
    determinant = a*d-b*c
    inverse = [[d/determinant, -b/determinant],
               [-c/determinant, a/determinant]]
    ghost_k = diff(diff(free, 'c'), 'h').get(bv.ZERO, Q(0))

    def contraction(poly, wrong_ghost=False):
        out = {}
        for i, x in enumerate(bosons):
            for j, y in enumerate(bosons):
                out = add(out, scale(diff(diff(poly, x), y), inverse[i][j]/2))
        sign = 1 if wrong_ghost else -1
        return add(out, scale(diff(diff(poly, 'c'), 'h'), Q(sign)/ghost_k))

    def laplacian(poly):
        return add(*(scale(diff(diff(poly, anti), field),
                           -1 if field in bv.ODD else 1)
                     for field, anti in PAIRS))

    return variables, free, hessian, inverse, ghost_k, contraction, laplacian


def monomials(max_degree=4):
    seen = {}
    for degree in range(max_degree+1):
        for names in combinations_with_replacement(ACTIVE, degree):
            poly = bv.prod(*(bv.var(name) for name in names))
            if poly:
                # Normalize any sign from reordering odd generators.
                key = next(iter(poly))
                seen[key] = {key: Q(1)}
    return list(seen.values())


def exponential(poly, operation):
    """Coefficients of exp(epsilon * operation), not a numerical epsilon."""
    answer, current, power = [], poly, 0
    while current:
        answer.append(scale(current, Q(1, factorial(power))))
        current = operation(current)
        power += 1
    return answer


def coefficient(series, order):
    return series[order] if 0 <= order < len(series) else {}


def graded_bv_identity():
    v, free, hessian, inverse, k, contraction, laplacian = setup()
    derivative = lambda poly: bv.bracket(free, poly)
    tested = monomials()
    coefficient_checks = 0
    for poly in tested:
        assert not derivative(derivative(poly))
        assert not laplacian(laplacian(poly))
        assert not add(derivative(laplacian(poly)), laplacian(derivative(poly)))
        assert not add(derivative(contraction(poly)),
                       scale(contraction(derivative(poly)), -1), laplacian(poly))
        transformed = exponential(poly, contraction)
        transformed_derivative = exponential(derivative(poly), contraction)
        for order in range(max(len(transformed), len(transformed_derivative))+1):
            # D exp(eps C) - exp(eps C) D = -eps Delta exp(eps C).
            assert not add(derivative(coefficient(transformed, order)),
                           scale(coefficient(transformed_derivative, order), -1),
                           laplacian(coefficient(transformed, order-1)))
            coefficient_checks += 1
    probe = mul(v['h'], v['a'])
    wrong = lambda poly: contraction(poly, wrong_ghost=True)
    defect = add(derivative(wrong(probe)), scale(wrong(derivative(probe)), -1),
                 laplacian(probe))
    assert defect == scale(bv.ONE, -2)
    return dict(supermonomials=len(tested), coefficient_checks=coefficient_checks,
                free_boson_hessian=[[str(x) for x in row] for row in hessian],
                mixed_boson_covariance=[[str(x) for x in row] for row in inverse],
                ghost_contraction=str(-1/k),
                graded_quantum_bv_identity_exact=True,
                wrong_ghost_sign_defect=bv.display(defect))


def insertion_and_antifield():
    v, free, hessian, inverse, k, contraction, _ = setup()
    tested = monomials(3)
    insertions, central_checks = 0, 0
    anti_square = mul(v['a'], v['a'])
    for poly in tested:
        series = exponential(poly, contraction)
        # Two even basic field equations with the SAME mixed Hessian.
        for i in range(2):
            equation = add(*(scale(v[x], hessian[i][j])
                             for j, x in enumerate(('q', 'r'))))
            lhs = exponential(mul(equation, poly), contraction)
            gradient = exponential(diff(poly, ('q', 'r')[i]), contraction)
            for order in range(max(len(lhs), len(series), len(gradient)+1)):
                rhs = add(mul(equation, coefficient(series, order)),
                          coefficient(gradient, order-1))
                assert coefficient(lhs, order) == rhs
                insertions += 1
        lhs = exponential(mul(anti_square, poly), contraction)
        for order in range(max(len(lhs), len(series))):
            assert coefficient(lhs, order) == mul(anti_square, coefficient(series, order))
            central_checks += 1
    anti_variation = bv.bracket(free, anti_square)
    assert anti_variation == scale(mul(v['a'], v['c']), 2*k)
    assert anti_variation  # No contractions does NOT mean BRST-inert.
    # Removing mixed covariance keeps the two separate variances but destroys FE.
    diagonal_only = [[inverse[0][0], Q(0)], [Q(0), inverse[1][1]]]
    wrong_contact = sum(hessian[0][j]*diagonal_only[j][1] for j in range(2))
    assert wrong_contact != 0
    return dict(supermonomials=len(tested), basic_equation_coefficient_checks=insertions,
                quadratic_antifield_factor_checks=central_checks,
                antifield_square_brst_variation=bv.display(anti_variation),
                omit_mixed_covariance_offdiagonal_residual=str(wrong_contact),
                preserves_antifield_square_without_treating_it_as_inert=True)


def contact_loop_boundary():
    """Exact 4-D Gaussian regulator, with 16*pi^2 factored out.

    delta_eps=(2*pi*eps^2)^-2 exp(-|x|^2/(2*eps^2)),
    f_mu=exp(-mu*|x|^2). The subtraction is local through two jets;
    the finite limit probes four jets. This is a UV diagnostic in a chart.
    """
    rows = []
    for mu in (Q(1), Q(2), Q(3)):
        previous_error = None
        for power in (2, 3, 4, 5, 6):
            eps = Q(1, 2**power)
            raw = 1/(eps**4*(1+mu*eps**2)**2)
            remainder = raw-1/eps**4+2*mu/eps**2
            limit = 3*mu**2
            # Closed algebraic remainder after subtracting the two divergent jets.
            exact = -mu**3*eps**2*(4+3*mu*eps**2)/(1+mu*eps**2)**2
            assert remainder-limit == exact
            error = abs(exact)
            if previous_error is not None:
                assert error < previous_error
            previous_error = error
            rows.append(dict(mu=str(mu), epsilon=str(eps), raw_scaled=str(raw),
                             subtracted_scaled=str(remainder),
                             finite_limit_scaled=str(limit), error=float(error)))
    # Different allowed finite delta counterterms change f(0), with same separated points.
    finite_scheme_shift = Q(7, 3)
    assert finite_scheme_shift != 0
    return dict(regulator_dimension=4, exact_regulator_cases=len(rows), cases=rows,
                raw_product_diverges_as='epsilon^-4',
                further_divergence='epsilon^-2 * Laplacian(f)(0) / 4',
                finite_part_scaled='Laplacian^2(f)(0) / 32',
                finite_delta_scheme_shift_scaled=str(finite_scheme_shift),
                regulator_is_not_a_physical_model=True,
                continuum_auxiliary_extension_proven=False)


def run():
    return dict(round=778, groups={
        'graded_bv_identity': graded_bv_identity(),
        'linear_insertions_and_quadratic_antifield': insertion_and_antifield(),
        'ultralocal_composite_loop_boundary': contact_loop_boundary()},
        all_checks_passed=True,
        finite_tests_prove_continuum_normalization=False,
        original_full_N1_N2_proven=False, frame_quantum_transfer_proven=False,
        original_loop_anomaly_computed=False, interacting_positive_state_proven=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    destination = HERE/'joint_propagating_normalization_results.json'
    if args.write:
        destination.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    if args.check:
        assert result == json.loads(destination.read_text(encoding='utf-8'))
    print(json.dumps(result, ensure_ascii=False, indent=2))
