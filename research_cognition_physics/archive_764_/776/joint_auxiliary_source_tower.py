"""776: exact diagnostics for auxiliary-source fusion and local contact jets.

Only the Python standard library is needed. These finite rational identities
check the analytic source map; they are not a continuum loop calculation.
Run with --check (default) to compare the saved result, or --write to save it.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import comb, factorial
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
RESULT = HERE / 'joint_auxiliary_source_tower_results.json'


def transpose(a):
    return tuple(zip(*a))


def mm(a, b):
    return tuple(tuple(sum(x*y for x, y in zip(row, col))
                       for col in transpose(b)) for row in a)


def addm(a, b):
    return tuple(tuple(x+y for x, y in zip(r, s)) for r, s in zip(a, b))


def inv(a):
    n = len(a)
    b = [[Q(x) for x in row] + [Q(i == j) for j in range(n)]
         for i, row in enumerate(a)]
    for i in range(n):
        pivot = next(j for j in range(i, n) if b[j][i])
        b[i], b[pivot] = b[pivot], b[i]
        scale = b[i][i]
        b[i] = [x/scale for x in b[i]]
        for j in range(n):
            if j != i:
                scale = b[j][i]
                b[j] = [x-scale*y for x, y in zip(b[j], b[i])]
    return tuple(tuple(row[n:]) for row in b)


def exponents(d, degree):
    if d == 1:
        yield (degree,)
    else:
        for i in range(degree+1):
            for tail in exponents(d-1, degree-i):
                yield (i,)+tail


def pmul(a, b):
    out = {}
    for e, c in a.items():
        for f, d in b.items():
            k = tuple(x+y for x, y in zip(e, f))
            out[k] = out.get(k, Q(0)) + c*d
    return {e: c for e, c in out.items() if c}


def padd(a, b):
    out = dict(a)
    for e, c in b.items():
        out[e] = out.get(e, Q(0)) + c
    return {e: c for e, c in out.items() if c}


def moment_function(cov):
    @lru_cache(None)
    def moment(e):
        if not any(e):
            return Q(1)
        if sum(e) % 2:
            return Q(0)
        i = next(i for i, n in enumerate(e) if n)
        a = list(e)
        a[i] -= 1
        value = Q(0)
        for j, n in enumerate(a):
            if n and cov[i][j]:
                b = a.copy()
                b[j] -= 1
                value += n*cov[i][j]*moment(tuple(b))
        return value
    return moment


def source_tower():
    # A rational finite symbol with two gauge directions and one physical one.
    # PK=0 is exact. The covariances below have their common i*hbar stripped.
    p = ((0, 0, 0), (0, 0, 0), (0, 0, 3))
    k = ((1, 1), (0, 1), (0, 0))
    n = ((2, 1), (1, 3))
    ni = inv(n)
    u = mm(ni, transpose(k))
    sigma = inv(addm(p, mm(k, u)))
    ell = tuple(tuple(p[i])+tuple(k[i]) for i in range(3)) + tuple(
        tuple(transpose(k)[i])+tuple(-Q(x) for x in n[i]) for i in range(2))
    full = inv(ell)
    assert all(full[i][j] == 0 for i in (3, 4) for j in (3, 4))
    direct = moment_function(full)
    diagonal = tuple(tuple(sigma[i])+(Q(0), Q(0)) for i in range(3)) + tuple(
        (Q(0),)*3+tuple(-x for x in ni[i]) for i in range(2))
    reduced = moment_function(diagonal)
    beta_forms = []
    for j in range(2):
        form = {}
        for i, value in enumerate(u[j]):
            if value:
                e = tuple(int(a == i) for a in range(5))
                form[e] = value
        form[tuple(int(a == 3+j) for a in range(5))] = Q(1)
        beta_forms.append(form)
    vertex = {(2, 1, 0): Q(2, 7), (0, 0, 4): Q(1, 11),
              (1, 1, 1): Q(3, 5)}
    powers = [{(0, 0, 0): Q(1)}]
    for _ in range(2):
        powers.append(pmul(powers[-1], vertex))
    checks = 0
    nonzero = 0
    for total in range(7):
        for e in exponents(5, total):
            expanded = {e[:3]+(0, 0): Q(1)}
            for j in range(2):
                for _ in range(e[3+j]):
                    expanded = pmul(expanded, beta_forms[j])
            for vertices in powers:
                # Compare every power of the common i*hbar separately: the
                # cubic and quartic vertices must not mask an order mismatch.
                lhs, rhs = {}, {}
                for ve, c in vertices.items():
                    contractions = (sum(e)+sum(ve))//2
                    left = c*direct(tuple(e[i]+ve[i] for i in range(3))+e[3:])
                    right = sum(c*d*reduced(tuple(f[i]+ve[i] for i in range(3))+f[3:])
                                for f, d in expanded.items())
                    lhs[contractions] = lhs.get(contractions, Q(0))+left
                    rhs[contractions] = rhs.get(contractions, Q(0))+right
                assert lhs == rhs, (e, vertices, lhs, rhs)
                checks += 1
                nonzero += int(any(lhs.values()))
    # Omitting the contact incorrectly leaves a nonzero pure b covariance.
    wrong_bb = mm(mm(u, sigma), transpose(u))
    assert wrong_bb == ni
    controls = {str(r): str(factorial(r)//(2**(r//2)*factorial(r//2))
                            * wrong_bb[0][0]**(r//2)) for r in (2, 4, 6)}
    mixed = direct((2, 0, 0, 2, 0))
    assert mixed != 0
    return dict(exact_mixed_source_and_vertex_checks=checks,
                nonzero_cases=nonzero, largest_source_degree=6,
                interaction_orders=[0, 1, 2], rational_identity_failures=0,
                common_i_hbar_powers_compared_separately=True,
                pure_b_covariance='0',
                omitted_contact_pure_b_moments=controls,
                mixed_b0_squared_v0_squared=str(mixed),
                scope='Finite 3+2 rational calibration of the all-order analytic source map; not the original 111-field continuum loop amplitudes.')


def deriv(p, axis):
    out = {}
    for e, c in p.items():
        if e[axis]:
            f = list(e)
            f[axis] -= 1
            out[tuple(f)] = c*e[axis]
    return out


def wick(p, cov):
    # Variables x, beta, external antifield a, hbar. No a contractions.
    term, out, order = dict(p), dict(p), 0
    while term:
        order += 1
        nxt = {}
        for i, j in product(range(3), repeat=2):
            if not cov[i][j]:
                continue
            for e, c in deriv(deriv(term, i), j).items():
                f = e[:3]+(e[3]+1,)
                nxt[f] = nxt.get(f, Q(0))-cov[i][j]*c/(2*order)
        term = {e: c for e, c in nxt.items() if c}
        out = padd(out, term)
    return out


def wick_transport():
    hb = ((Q(4, 9), 0, 0), (0, 0, 0), (0, 0, 0))
    smooth = ((0, Q(2, 7), 0), (Q(2, 7), Q(3, 5), 0), (0, 0, 0))
    h = addm(hb, smooth)
    monomials = derivatives = 0
    for total in range(7):
        for e in exponents(3, total):
            f = {e+(0,): Q(1)}
            assert wick(f, h) == wick(wick(f, smooth), hb)
            monomials += 1
            for axis in range(3):
                assert deriv(wick(f, smooth), axis) == wick(deriv(f, axis), smooth)
                derivatives += 1
    beta2 = {(0, 2, 0, 0): Q(1)}
    betaxx = {(2, 1, 0, 0): Q(1)}
    assert wick(beta2, smooth) == {**beta2, (0, 0, 0, 1): Q(-3, 5)}
    assert wick(betaxx, smooth) == {**betaxx, (1, 0, 0, 1): Q(-4, 7)}
    cubic_and_external = {(3, 0, 0, 0): Q(2, 5), (0, 1, 1, 0): Q(1)}
    assert wick(cubic_and_external, smooth) == cubic_and_external
    return dict(exact_monomial_transport_checks=monomials,
                fixed_background_field_derivative_checks=derivatives,
                source_cubic_and_external_antifield_unchanged=True,
                lost_contact_if_transport_omitted={
                    'beta^2': '-3*hbar/5', 'beta*x^2': '-4*hbar*x/7'},
                scope='Finite smooth Wick transport; beta-composite insertions require this transport, not a zero-contact convention.')


def regular_coefficient_contacts():
    # Inverse Gram coefficient of N(x)=[[2+x,1+2x],[1+2x,3-x]].
    # det N=5-3x-5x^2; positive on [-1/10,1/10].
    denominator, numerator = (Q(5), Q(-3), Q(-5)), (Q(3), Q(-1))
    coeff = []
    for j in range(7):
        rhs = numerator[j] if j < len(numerator) else Q(0)
        rhs -= sum(denominator[k]*coeff[j-k] for k in range(1, min(2, j)+1))
        coeff.append(rhs/denominator[0])
    exact_checks = 0
    for r in range(7):
        for q in range(9):
            # Test on x**q: LHS applies delta^(r) to a(x)*x**q.
            lhs = (-1)**r*factorial(r)*coeff[r-q] if q <= r else Q(0)
            rhs = Q(0)
            for j in range(r+1):
                delta_order = r-j
                delta_on_test = (-1)**q*factorial(q) if delta_order == q else 0
                rhs += (-1)**j*comb(r, j)*factorial(j)*coeff[j]*delta_on_test
            assert lhs == rhs
            exact_checks += 1
    missing = 2*coeff[2]
    assert missing != 0
    return dict(inverse_gram_taylor_coefficients=[str(x) for x in coeff],
                delta_derivative_product_checks=exact_checks,
                largest_delta_derivative_order=6,
                frozen_coefficient_wrong_answer_on_constant_test='0',
                correct_second_derivative_contact_on_constant_test=str(missing),
                scope='Exact contact-jet identity on a regular coefficient patch; finite scaling-degree extension is an analytic cited theorem, not tested numerically.')


def run():
    return {'round': 776, 'test_groups': 3, 'all_passed': True,
            'auxiliary_source_tower': source_tower(),
            'smooth_wick_transport': wick_transport(),
            'regular_coefficient_contacts': regular_coefficient_contacts(),
            'full_N1_N2_realization_proven': False,
            'interacting_positive_state_or_instrument_proven': False,
            'continuum_loop_anomaly_coefficients_computed': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write:
        RESULT.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    else:
        assert json.loads(RESULT.read_text(encoding='utf-8')) == result, 'Saved result mismatch'
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
