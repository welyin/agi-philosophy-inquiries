"""785 diagnostics: constrained projection and formal relational slicing.

Group 1 is a finite matrix calibration, not a spacetime Green kernel or a
Hadamard-state existence test. Group 2 uses exact rational finite jets.
The continuous statements and their hypotheses are in research_note_785.md.
"""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json
import numpy as np


def matrix_projection():
    n = 7
    eye = np.eye(n)
    derivative = (np.roll(eye, 1, axis=1)-np.roll(eye, -1, axis=1))/2
    derivative = np.diag(np.linspace(.6, 1.4, n)) @ derivative
    derivative += np.diag(np.linspace(-.2, .3, n))
    k = np.vstack((derivative, eye))
    f = np.hstack((eye, -derivative))
    l0 = np.hstack((np.zeros((n, n)), eye))
    t = .17*eye + .09*derivative
    l = l0+t@f
    pi = np.eye(2*n)-k@l
    h = 2*eye+derivative.T@derivative
    p = f.T@h@f
    d1 = p-k@k.T
    d0 = -k.T@k
    g1 = np.linalg.inv(d1)
    g0 = np.linalg.inv(d0)
    e = pi@g1@pi.T
    residuals = {
        'left_inverse': l@k-eye,
        'projector': pi@pi-pi,
        'kills_gauge': pi@k,
        'slice_constraint': l@pi,
        'noether': p@k,
        'ward': g1@k-k@g0,
        'dual_ward': k.T@g1-g0@k.T,
        'source_inverse': p@e-pi.T,
        'field_inverse_on_slice': e@p@pi-pi,
        'source_redundancy': pi.T@l.T,
    }
    errors = {name: float(np.max(np.abs(value))) for name, value in residuals.items()}
    assert max(errors.values()) < 3e-12, errors
    physical_lift = f.T@np.linalg.inv(f@f.T)
    sigma = eye+.13*derivative.T@derivative
    c = physical_lift@sigma@physical_lift.T - .4*k@k.T
    a = np.arange(n*2*n, dtype=float).reshape(n, 2*n)/93
    change = k@a + a.T@k.T
    reduced = pi@c@pi.T
    projected_change = pi@change@pi.T
    assert np.max(np.abs(projected_change)) < 3e-12
    eig = np.linalg.eigvalsh(reduced)
    assert eig.min() > -3e-12
    assert np.linalg.eigvalsh(c).min() < -.1
    rng = np.random.default_rng(785)
    tests = rng.normal(size=(2*n, 31))
    physical_tests = pi.T@tests
    lhs = tests.T@reduced@tests
    rhs = physical_tests.T@c@physical_tests
    assert np.max(np.abs(lhs-rhs)) < 3e-12
    # Pi is oblique: replacing its dual by itself fails the source equation.
    wrong = p@(pi@g1@pi)-pi.T
    wrong_residual = float(np.max(np.abs(wrong)))
    assert wrong_residual > .1
    return {
        'sites': n, 'field_dimension': 2*n, 'gauge_dimension': n,
        'projector_rank': int(np.linalg.matrix_rank(pi)),
        'max_identity_residual': float(max(errors.values())),
        'identity_residuals': errors,
        'unprojected_covariance_min_eigenvalue': float(np.linalg.eigvalsh(c).min()),
        'projected_covariance_min_eigenvalue': float(eig.min()),
        'projected_covariance_positive_eigenvalues': int(np.count_nonzero(eig > 1e-10)),
        'pure_gauge_extension_change_residual': float(np.max(np.abs(projected_change))),
        'wrong_dual_projection_residual': wrong_residual,
        'scope': 'Finite Ward/projection/covariance calibration; not a causal or on-shell covariance model.'
    }


ORDER = 4
JETS = 10
NVAR = 3*JETS
ZERO = (0,)*NVAR


def clean(p):
    return {key: Q(value) for key, value in p.items() if value and sum(key) <= ORDER}


def add(*terms):
    result = {}
    for term in terms:
        for key, value in term.items():
            result[key] = result.get(key, Q(0))+value
    return clean(result)


def scale(p, a):
    return clean({k: a*v for k, v in p.items()})


def mul(p, q):
    result = {}
    for a, av in p.items():
        for b, bv in q.items():
            key = tuple(x+y for x, y in zip(a, b))
            if sum(key) <= ORDER:
                result[key] = result.get(key, Q(0))+av*bv
    return clean(result)


def power(p, degree):
    result = {ZERO: Q(1)}
    for _ in range(degree):
        result = mul(result, p)
    return result


def jet(family, degree):
    key = [0]*NVAR
    key[family*JETS+degree] = 1
    return {tuple(key): Q(1)}


def deriv(p, times=1):
    for _ in range(times):
        result = {}
        for key, value in p.items():
            for index, exponent in enumerate(key):
                if exponent:
                    assert index % JETS < JETS-1, 'Insufficient jet order'
                    new = list(key)
                    new[index] -= 1
                    new[index+1] += 1
                    new = tuple(new)
                    result[new] = result.get(new, Q(0))+value*exponent
        p = clean(result)
    return p


def homogeneous(p, degree):
    return {k: v for k, v in p.items() if sum(k) == degree}


def shift(p, displacement):
    return add(*(scale(mul(power(displacement, j), deriv(p, j)), Q(1, factorial(j)))
                 for j in range(ORDER)))


def inverse_reference(q):
    # Lagrange-Buermann inversion of chi+q(x+chi)=0 at fixed x.
    return add(*(scale(deriv(power(q, n), n-1), Q((-1)**n, factorial(n)))
                 for n in range(1, ORDER+1)))


def substitute(p, replacements):
    result = {}
    cache = {}
    for key, value in p.items():
        term = {ZERO: value}
        for index, exponent in enumerate(key):
            if exponent:
                cache_key = (index, exponent)
                if cache_key not in cache:
                    cache[cache_key] = power(replacements.get(index, jet(index//JETS, index%JETS)), exponent)
                term = mul(term, cache[cache_key])
        result = add(result, term)
    return result


def exact_relational_slice():
    q, r, alpha = jet(0, 0), jet(1, 0), jet(2, 0)
    slope = Q(2, 3)
    chi = inverse_reference(q)
    residual = add(chi, shift(q, chi))
    assert not residual
    w = add(scale(chi, slope), shift(r, chi))
    q_transformed = add(alpha, shift(q, alpha))
    r_transformed = add(scale(alpha, slope), shift(r, alpha))
    replacements = {}
    for j in range(ORDER):
        replacements[j] = deriv(q_transformed, j)
        replacements[JETS+j] = deriv(r_transformed, j)
    w_transformed = substitute(w, replacements)
    assert w_transformed == w
    chi_transformed = substitute(chi, replacements)
    assert not add(chi_transformed, shift(q_transformed, chi_transformed))
    # Independent recursive solution, without the inversion formula.
    recursive = {}
    for degree in range(1, ORDER+1):
        defect = homogeneous(add(recursive, shift(q, recursive)), degree)
        recursive = add(recursive, scale(defect, -1))
    assert recursive == chi
    w1 = add(r, scale(q, -slope))
    w2 = add(scale(mul(q, deriv(q)), slope), scale(mul(q, deriv(r)), -1))
    w3 = add(scale(mul(q, power(deriv(q), 2)), -slope),
             scale(mul(power(q, 2), deriv(q, 2)), -slope/2),
             mul(mul(q, deriv(q)), deriv(r)),
             scale(mul(power(q, 2), deriv(r, 2)), Q(1, 2)))
    assert homogeneous(w, 1) == w1
    assert homogeneous(w, 2) == w2
    assert homogeneous(w, 3) == w3
    linear_failure = add(substitute(w1, replacements), scale(w1, -1))
    expected = mul(alpha, add(deriv(r), scale(deriv(q), -slope)))
    assert homogeneous(linear_failure, 2) == expected
    assert expected
    return {
        'coefficient_arithmetic': 'fractions.Fraction', 'field_degree_through': ORDER,
        'background_slope': str(slope),
        'inverse_terms': len(chi), 'relational_field_terms': len(w),
        'terms_by_degree': {str(n): len(homogeneous(w, n)) for n in range(1, ORDER+1)},
        'slice_equation_residual_terms': len(residual),
        'gauge_transformation_residual_terms': len(add(w_transformed, scale(w, -1))),
        'recursive_and_Lagrange_inverses_agree': recursive == chi,
        'linear_projection_noninvariance_degree': 2,
        'linear_projection_noninvariance_leading_terms': len(expected),
        'scope': 'Formal one-dimensional reference/relational-field jet diagnostic, not a spacetime dimension or a quantum matching proof.'
    }


def run():
    return {'round': 785, 'fresh_test_groups': 2,
            'matrix_projection': matrix_projection(),
            'exact_relational_slice': exact_relational_slice(),
            'all_checks_passed': True,
            'original_interacting_quantum_dictionary_proven': False,
            'original_interacting_positive_state_proven': False}


if __name__ == '__main__':
    result = run()
    path = Path(__file__).with_name('local_physical_slice_results.json')
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
