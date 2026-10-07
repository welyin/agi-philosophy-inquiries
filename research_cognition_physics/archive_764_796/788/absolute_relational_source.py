"""788: absolute one-loop coordinate-source and relational-mean diagnostics.

Exact rational calibrations of identities proved in research_note_788.md.
Finite matrices/kernels are NOT the original spacetime state or a proof of
endpoint equivalence, a continuous measure, or interacting positivity.
"""
from fractions import Fraction as Q
from itertools import product
from math import factorial
from pathlib import Path
import importlib.util
import json
import random
import sys

sys.dont_write_bytecode = True
N = 4
ZERO = (0,)*N


def add(*polys):
    out = {}
    for p in polys:
        for key, value in p.items():
            out[key] = out.get(key, Q(0))+value
    return {k: v for k, v in out.items() if v}


def scale(p, a):
    return {k: a*v for k, v in p.items() if a*v}


def mul(p, q):
    out = {}
    for a, av in p.items():
        for b, bv in q.items():
            key = tuple(x+y for x, y in zip(a, b))
            if sum(key) <= 3:
                out[key] = out.get(key, Q(0))+av*bv
    return {k: v for k, v in out.items() if v}


def var(i):
    key = [0]*N
    key[i] = 1
    return {tuple(key): Q(1)}


def jet(p, indices):
    key = tuple(indices.count(i) for i in range(N))
    return p.get(key, Q(0))*product_factorials(key)


def product_factorials(key):
    value = 1
    for k in key:
        value *= factorial(k)
    return value


def transpose(a):
    return [list(row) for row in zip(*a)]


def mm(a, b):
    return [[sum((x*y for x, y in zip(row, col)), Q(0))
             for col in zip(*b)] for row in a]


def mv(a, v):
    return [sum((x*y for x, y in zip(row, v)), Q(0)) for row in a]


def plus(*vectors):
    return [sum(col, Q(0)) for col in zip(*vectors)]


def identity():
    return [[Q(i == j) for j in range(N)] for i in range(N)]


def contract(tensor, c):
    return [sum((tensor[k][i][j]*c[i][j]/2 for i, j in product(range(N), repeat=2)), Q(0))
            for k in range(N)]


def pullback_action(f, d, v):
    # Expand S_Y(f(x)) as a polynomial; extract derivatives afterwards.
    terms = []
    for a, b in product(range(N), repeat=2):
        terms.append(scale(mul(f[a], f[b]), d[a][b]/2))
    for a, b, c in product(range(N), repeat=3):
        terms.append(scale(mul(mul(f[a], f[b]), f[c]), v[a][b][c]/6))
    return add(*terms)


def absolute_source_identity():
    rows = []
    for seed in range(788, 791):
        rng = random.Random(seed)
        t, ti = identity(), identity()
        # An invertible triangular stencil, with a non-self-adjoint dual.
        for i, j in ((0, 2), (0, 3), (1, 2), (1, 3)):
            t[i][j] = Q(rng.randint(-3, 3), 5)
            ti[i][j] = -t[i][j]
        assert mm(t, ti) == identity()
        dy = [[Q(0) for _ in range(N)] for _ in range(N)]
        dy[0][0], dy[1][1] = Q(2), Q(3)
        dx = mm(mm(transpose(t), dy), t)
        u = [[[Q(0) for _ in range(N)] for _ in range(N)] for _ in range(N)]
        for a in range(N):
            for i in range(N):
                for j in range(i, N):
                    u[a][i][j] = u[a][j][i] = Q(rng.randint(-4, 4), 7)
        v_by_key = {tuple(sorted(idx)): Q(rng.randint(-3, 3), 5)
                    for idx in product(range(N), repeat=3)}
        v = [[[v_by_key[tuple(sorted((a, b, c))) ] for c in range(N)]
              for b in range(N)] for a in range(N)]
        f = [add(*(scale(var(i), t[a][i]) for i in range(N)),
                 *(scale(mul(var(i), var(j)), u[a][i][j]/2)
                   for i, j in product(range(N), repeat=2))) for a in range(N)]
        sx = pullback_action(f, dy, v)
        assert [[jet(sx, [i, j]) for j in range(N)] for i in range(N)] == dx
        gx = [[[jet(sx, [i, j, k]) for k in range(N)] for j in range(N)] for i in range(N)]
        h = [[Q((i+1)*(j+1), 11)+Q(i == j, 3) for j in range(N)] for i in range(N)]
        contacts, means = [], []
        for amplitude in (Q(1), Q(5, 2)):
            wy = [[Q(0) for _ in range(N)] for _ in range(N)]
            wy[2][2], wy[3][3] = amplitude, amplitude+1
            w = mm(mm(ti, wy), transpose(ti))
            assert not any(x for row in mm(dx, w) for x in row)
            c = [[w[i][j]-h[i][j] for j in range(N)] for i in range(N)]
            cy = mm(mm(t, c), transpose(t))
            k = contract(u, c)
            jx, jy = contract(gx, c), contract(v, cy)
            # Cross term computed from the original operators, independently
            # of the expanded action polynomial.
            rho = [sum((c[i][j]*u[a][ell][i]*dy[a][b]*t[b][j]
                        for i, j, a, b in product(range(N), repeat=4)), Q(0))
                   for ell in range(N)]
            predicted = plus(mv(transpose(t), plus(jy, mv(dy, k))), rho)
            assert jx == predicted
            # Replacing C by -H yields precisely the same local remainder.
            local = [-sum((h[i][j]*u[a][ell][i]*dy[a][b]*t[b][j]
                           for i, j, a, b in product(range(N), repeat=4)), Q(0))
                     for ell in range(N)]
            assert rho == local
            lx = [Q(i+1, 13) for i in range(N)]
            ly = mv(transpose(ti), plus(lx, rho))
            mux = [Q(i-1, 17) for i in range(N)]
            muy = plus(mv(t, mux), k)
            residual_x = plus(mv(dx, mux), jx, lx)
            residual_y = mv(transpose(t), plus(mv(dy, muy), jy, ly))
            assert residual_x == residual_y
            # A nonlinear measured quantity is expanded separately.
            oa = [Q(i+1, 3) for i in range(N)]
            ob = [[Q(i+j+1, 9) for j in range(N)] for i in range(N)]
            ox = add(*(scale(f[a], oa[a]) for a in range(N)),
                     *(scale(mul(f[a], f[b]), ob[a][b]/2)
                       for a, b in product(range(N), repeat=2)))
            ox_mean = sum((jet(ox, [i])*mux[i] for i in range(N)), Q(0))
            ox_mean += sum((jet(ox, [i, j])*c[i][j]/2 for i, j in product(range(N), repeat=2)), Q(0))
            oy_mean = sum((oa[a]*muy[a] for a in range(N)), Q(0))
            oy_mean += sum((ob[a][b]*cy[a][b]/2 for a, b in product(range(N), repeat=2)), Q(0))
            assert ox_mean == oy_mean
            assert any(rho) and any(mv(transpose(t), mv(dy, k)))
            contacts.append(rho)
            means.append(k)
        assert contacts[0] == contacts[1] and means[0] != means[1]
        rows.append(dict(seed=seed, states=2, exact_source_and_mean_residual='0',
                         remainder_unchanged_by_state=True,
                         omitted_contact_max_error=str(max(map(abs, contacts[0]))),
                         transported_mean_changes_with_state=True))
    return dict(cases=rows, arithmetic='Fraction', action_derivatives='independent polynomial expansion',
                scope='Finite algebra with a nonzero on-shell kernel; not a discretized Hadamard state.')


def derivative_relational_mean():
    path = Path(__file__).resolve().parents[1]/'785/local_physical_slice.py'
    spec = importlib.util.spec_from_file_location('round785_reused', path)
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    q, r = old.jet(0, 0), old.jet(1, 0)
    slope = Q(2, 3)
    chi = old.inverse_reference(q)
    second = old.homogeneous(old.add(old.scale(chi, slope), old.shift(r, chi)), 2)
    rows = []
    for theta in (Q(0), Q(1, 3), Q(-4, 5)):
        # At x=0: q=z1, q'=-theta*z2, r=z2, r'=theta*z1.
        # q(x),r(x)=R(theta*x)z is an exact smooth positive rank-two kernel.
        vectors = {0: (Q(1), Q(0)), 1: (Q(0), -theta),
                   old.JETS: (Q(0), Q(1)), old.JETS+1: (theta, Q(0))}
        wick = Q(0)
        for key, value in second.items():
            inds = [i for i, power in enumerate(key) for _ in range(power)]
            assert len(inds) == 2
            wick += value*sum((a*b for a, b in zip(vectors[inds[0]], vectors[inds[1]])), Q(0))
        assert wick == -theta
        rows.append(dict(theta=str(theta), same_diagonal_covariance='identity_2',
                         cross_right_derivative=str(theta), relational_wick_mean=str(wick)))
    return dict(cases=rows, nonzero_derivative_contacts_despite_identical_diagonal=True,
                scope='Smooth positive finite-rank kernel calibration; no field equations or quantum CCR claimed.')


def interacting_finite_gauge_source():
    mass, cubic, curvature = Q(3, 2), Q(5, 7), Q(4, 3)
    rows = []
    for a, b in ((Q(1), Q(0)), (Q(2), Q(1, 3)), (Q(3, 2), Q(-2, 5))):
        cqq = (mass+b*b)/(a*a*mass)
        cqr, crr = -b/(a*mass), 1/mass
        jq_b = -mass*curvature*cqr
        jq_gh = -b*curvature/a
        jr = -mass*curvature*cqq/2+cubic*crr/2
        # Independent off-shell normal derivative of the complete Hessian det.
        det0 = a*a*mass
        det_r = -mass*curvature*(mass+b*b)+a*a*cubic
        assert jr == det_r/(2*det0) and jq_b+jq_gh == 0
        mean_r = -crr*jr
        mean_w = mean_r-curvature*cqq/2
        assert mean_w == -cubic/(2*mass*mass)
        omitted_ghost_mean = mean_w-cqr*jq_b
        rows.append(dict(a=str(a), b=str(b), bosonic_source_q=str(jq_b), ghost_source_q=str(jq_gh),
                         source_r=str(jr), relational_mean=str(mean_w),
                         wrong_mean_without_ghost=str(omitted_ghost_mean)))
    assert len({r['source_r'] for r in rows}) == 3
    assert all(r['relational_mean'] == str(-cubic/(2*mass*mass)) for r in rows)
    assert rows[1]['wrong_mean_without_ghost'] != rows[1]['relational_mean']
    return dict(mass=str(mass), physical_cubic=str(cubic), slice_curvature=str(curvature), cases=rows,
                nonzero_physical_mean_preserved=True,
                scope='Finite local Euclidean perturbative gauge diagnostic, not the original prepared Lorentzian state.')


def run():
    return dict(round=788, fresh_test_groups=3,
                absolute_source_identity=absolute_source_identity(),
                derivative_relational_mean=derivative_relational_mean(),
                interacting_finite_gauge_source=interacting_finite_gauge_source(),
                all_checks_passed=True, absolute_fixed_gauge_one_loop_insertion_dictionary_proven=True,
                auxiliary_endpoint_quantum_matching_proven=False,
                common_action_integrability_of_transported_contacts_proven=False,
                original_interacting_positive_state_proven=False)


if __name__ == '__main__':
    result = run()
    Path(__file__).with_name('absolute_relational_source_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
