"""777: exact BV canonical splitting of the original quadratic auxiliary sector.

The finite nonlinear gauge model is a sign/closure diagnostic, not the original
continuum interaction or a quantum anomaly calculation. Python standard library
only. Use --write once, then --check to reproduce the stored result.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from pathlib import Path
import argparse
import json

# Same symbols are reused for old coordinates and new coordinates; substitutions
# below always state which chart is meant. b means the new e in the split chart.
EVEN = ('q', 'r', 'b', 'z', 'a')  # z=c*, a=bar c*
ODD = ('c', 'h', 'p', 't', 'rho')  # h=bar c (or eta), p=q*, t=r*, rho=b*
NAMES = EVEN+ODD
PAIRS = (('q', 'p'), ('r', 't'), ('b', 'rho'), ('c', 'z'), ('h', 'a'))
ZERO = (0,)*len(NAMES)
ONE = {ZERO: Q(1)}


def add(*polys):
    out = {}
    for poly in polys:
        for key, value in poly.items():
            out[key] = out.get(key, Q(0))+value
    return {key: value for key, value in out.items() if value}


def scale(poly, value):
    return {key: value*c for key, c in poly.items() if value*c}


def var(name):
    e = list(ZERO)
    e[NAMES.index(name)] = 1
    return {tuple(e): Q(1)}


def mul(a, b):
    out = {}
    start = len(EVEN)
    for e, c in a.items():
        for f, d in b.items():
            if any(e[i] and f[i] for i in range(start, len(NAMES))):
                continue
            sign = (-1)**sum(e[i]*f[j] for i in range(start, len(NAMES))
                            for j in range(start, i))
            key = tuple(x+y for x, y in zip(e, f))
            out[key] = out.get(key, Q(0))+sign*c*d
    return {key: value for key, value in out.items() if value}


def prod(*polys):
    out = ONE
    for poly in polys:
        out = mul(out, poly)
    return out


def diff(poly, name, right=False):
    axis = NAMES.index(name)
    out = {}
    odd = name in ODD
    for e, c in poly.items():
        if not e[axis]:
            continue
        sign = (-1)**sum(e[len(EVEN):axis]) if odd else 1
        if right and odd:
            sign *= (-1)**(sum(e[len(EVEN):])+1)
        f = list(e)
        f[axis] -= 1
        out[tuple(f)] = sign*e[axis]*c
    return out


def bracket(f, g):
    out = {}
    for field, anti in PAIRS:
        out = add(out, mul(diff(f, field, True), diff(g, anti)),
                  scale(mul(diff(f, anti, True), diff(g, field)), -1))
    return out


def substitute(poly, mapping):
    out = {}
    for e, c in poly.items():
        term = scale(ONE, c)
        for name, power in zip(NAMES, e):
            for _ in range(power):
                term = mul(term, mapping.get(name, var(name)))
        out = add(out, term)
    return out


def display(poly):
    if not poly:
        return '0'
    terms = []
    for e, c in sorted(poly.items()):
        monomial = '*'.join(name+(f'^{power}' if power > 1 else '')
                            for name, power in zip(NAMES, e) if power)
        terms.append(str(c)+(('*'+monomial) if monomial else ''))
    return ' + '.join(terms)


def setup():
    v = {name: var(name) for name in NAMES}
    q, r, b, a, c, h, p, t, rho = [v[x] for x in ('q','r','b','a','c','h','p','t','rho')]
    n, gq, gr = Q(3, 2), Q(2, 3), Q(-1, 5)
    f = add(scale(q, gq), scale(r, gr))
    yq, yr = c, mul(q, c)
    gy = add(scale(yq, gq), scale(yr, gr))
    invariant = add(r, scale(mul(q, q), Q(-1, 2)))
    minimal = add(scale(mul(invariant, invariant), Q(1, 2)), mul(p, yq), mul(t, yr))
    full = add(minimal, mul(b, f), scale(mul(b, b), -n/2),
               scale(mul(h, gy), -1), scale(mul(a, b), -1))
    red = add(minimal, scale(mul(f, f), 1/(2*n)), scale(mul(h, gy), -1),
              scale(mul(a, f), -1/n), scale(mul(a, a), 1/(2*n)))
    old_in_new = {'b': add(b, scale(f, 1/n), scale(a, -1/n)),
                  'h': add(h, scale(rho, -1/n)),
                  'p': add(p, scale(rho, -gq/n)),
                  't': add(t, scale(rho, -gr/n))}
    new_in_old = {'b': add(b, scale(f, -1/n), scale(a, 1/n)),
                  'h': add(h, scale(rho, 1/n)),
                  'p': add(p, scale(rho, gq/n)),
                  't': add(t, scale(rho, gr/n))}
    return v, n, f, gy, full, red, old_in_new, new_in_old


def canonical_and_master():
    v, n, f, gy, full, red, old_in_new, new_in_old = setup()
    assert not bracket(full, full)
    assert not bracket(red, red)
    split = add(red, scale(mul(v['b'], v['b']), -n/2))
    assert substitute(full, old_in_new) == split
    pairs = 0
    for x in NAMES:
        assert substitute(substitute(v[x], old_in_new), new_in_old) == v[x]
        for y in NAMES:
            assert bracket(substitute(v[x], old_in_new), substitute(v[y], old_in_new)) == bracket(v[x], v[y])
            pairs += 1
    expected_h = add(scale(f, 1/n), scale(v['a'], -1/n))
    assert bracket(red, v['h']) == expected_h
    assert bracket(red, v['a']) == gy
    omitted = add(red, scale(mul(v['a'], v['a']), -1/(2*n)))
    lost_nilpotency = bracket(omitted, bracket(omitted, v['h']))
    assert lost_nilpotency == scale(gy, 1/n) and lost_nilpotency
    # The 776 source variable has a_src=-a in these canonical conventions.
    beta = add(v['b'], scale(f, -1/n))
    assert bracket(full, beta) == scale(gy, -1/n)
    raw_beta_ideal_defect = substitute(bracket(full, beta), {'b': scale(f, 1/n)})
    assert raw_beta_ideal_defect
    wrong = dict(old_in_new)
    wrong['p'], wrong['t'] = v['p'], v['t']
    missing_antifield_map = add(substitute(full, wrong), scale(split, -1))
    assert missing_antifield_map
    return dict(canonical_generator_brackets=pairs,
                nonlinear_full_and_reduced_classical_master_equations_exact=True,
                transformed_action_exactly_decoupled=True,
                source_antifield_sign_dictionary='a_src(round776)=-a_canonical',
                reduced_s_antighost=display(expected_h),
                omit_antifield_square_nilpotency_defect=display(lost_nilpotency),
                set_beta_zero_off_shell_ideal_defect=display(raw_beta_ideal_defect),
                omit_antifield_shift_action_defect=display(missing_antifield_map))


def homotopy_and_menu():
    v, n, f, gy, full, red, old_in_new, _ = setup()
    split = add(red, scale(mul(v['b'], v['b']), -n/2))
    s = lambda p: bracket(split, p)
    projection = lambda p: substitute(p, {'b': {}, 'rho': {}})
    def h(poly):
        out = {}
        for e, c in poly.items():
            degree = e[NAMES.index('b')]+e[NAMES.index('rho')]
            if degree:
                term = mul(v['rho'], diff({e: c}, 'b'))
                out = add(out, scale(term, -1/(n*degree)))
        return out
    count = 0
    for degree in range(4):
        for names in combinations_with_replacement(NAMES, degree):
            poly = prod(*(v[name] for name in names))
            if not poly:
                continue
            assert not s(s(poly))
            assert add(s(h(poly)), h(s(poly))) == add(poly, scale(projection(poly), -1))
            assert projection(s(poly)) == bracket(red, projection(poly))
            count += 1
    # Whole local insertion algebra, not just linear external sources.
    auxiliary_square = bracket(full, mul(v['h'], v['b']))
    assert auxiliary_square == mul(v['b'], v['b'])
    # The projection is a chain map but is NOT a Poisson map on arbitrary inputs.
    poisson_defect = add(projection(bracket(v['b'], v['rho'])),
                        scale(bracket(projection(v['b']), projection(v['rho'])), -1))
    assert poisson_defect == ONE
    return dict(exact_supermonomials=count,
                contraction_nilpotency_and_chain_map_exact=True,
                actual_menu_s_antighost_times_b=display(auxiliary_square),
                projection_is_not_full_poisson_map_defect=display(poisson_defect),
                contraction='h=-(rho/n)*partial_e/N_aux',
                scope='Finite nonlinear BV algebra check; continuum finite-jet cohomology statement is proved in the report.')


def finite_derivative_adjoint():
    # Two sites with distinct pointwise N. A finite stencil tests the transpose
    # required by the continuum formal adjoint; it is not a new lattice model.
    from sys import path
    path.insert(0, str(Path(__file__).resolve().parents[1]/'776'))
    from joint_auxiliary_source_tower import mm, transpose
    ni = ((Q(1, 2), 0), (0, Q(1, 3)))
    g = ((Q(-1), Q(1)), (Q(-2), Q(2)))
    u = mm(ni, g)
    size = 6
    identity = lambda: [[Q(i == j) for j in range(size)] for i in range(size)]
    ae, bo = identity(), identity()
    # Even: (v,e,a); odd: (v*,rho,eta).
    for i in range(2):
        for j in range(2):
            ae[2+i][j] = -u[i][j]
            ae[2+i][4+j] = ni[i][j]
            bo[i][2+j] = u[j][i]
            bo[4+i][2+j] = ni[i][j]
    pairing = [[Q(0) for _ in range(size)] for _ in range(size)]
    for i in range(size):
        pairing[i][i] = Q(1 if i < 4 else -1)
    expected = tuple(tuple(row) for row in pairing)
    assert mm(mm(ae, pairing), transpose(bo)) == expected
    wrong = [row[:] for row in bo]
    for i in range(2):
        for j in range(2):
            wrong[i][2+j] = u[i][j]
    bad = mm(mm(ae, pairing), transpose(wrong))
    defect = max(abs(bad[i][j]-expected[i][j]) for i in range(size) for j in range(size))
    assert defect
    return dict(two_site_canonical_pairing_exact=True,
                wrong_formal_adjoint_max_defect=str(defect),
                pointwise_N=['2','3'],
                scope='Finite stencil checks the adjoint and variable coefficient order, not continuum propagation.')


def run():
    return dict(round=777, test_groups=3, all_passed=True,
                canonical=canonical_and_master(),
                reduction=homotopy_and_menu(), adjoint=finite_derivative_adjoint(),
                continuum_quantum_pushforward_proven=False,
                full_N1_N2_realization_proven=False,
                interacting_positive_state_or_instrument_proven=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--write', action='store_true')
    modes.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    target = Path(__file__).with_name('joint_auxiliary_bv_reduction_results.json')
    if args.write:
        target.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    else:
        assert json.loads(target.read_text(encoding='utf-8')) == result
    print(json.dumps(result, ensure_ascii=False, indent=2))
