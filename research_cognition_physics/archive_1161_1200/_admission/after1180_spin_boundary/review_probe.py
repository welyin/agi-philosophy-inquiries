"""Independent fixed-spin EPRL screen; exact magnetic CG contraction.

This does not import or read the author's probe.  No SymPy is required.
All algebraic certificates use rational square-root arithmetic. NumPy only
supplies diagnostic norms for polar maps. Default is read-only comparison;
--write creates the result once. This is admission support, not a new round.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np


@lru_cache(None)
def squarefree(n):
    outside, inside, p = 1, 1, 2
    while p * p <= n:
        power = 0
        while n % p == 0:
            n //= p
            power += 1
        outside *= p ** (power // 2)
        if power % 2:
            inside *= p
        p += 1
    return outside, inside * n


class R:
    """Real multiquadratic numbers: sum q_d sqrt(d), d squarefree."""
    def __init__(self, value=0):
        if isinstance(value, R):
            self.a = value.a.copy()
        elif isinstance(value, dict):
            self.a = {d: Q(c) for d, c in value.items() if c}
        else:
            self.a = {1: Q(value)} if value else {}

    def __add__(self, other):
        other, out = R(other), self.a.copy()
        for d, c in other.a.items():
            out[d] = out.get(d, Q(0)) + c
        return R(out)

    __radd__ = __add__

    def __neg__(self):
        return R({d: -c for d, c in self.a.items()})

    def __sub__(self, other):
        return self + (-R(other))

    def __rsub__(self, other):
        return R(other) - self

    def __mul__(self, other):
        out = {}
        for d, c in self.a.items():
            for e, f in R(other).a.items():
                outside, inside = squarefree(d * e)
                out[inside] = out.get(inside, Q(0)) + c * f * outside
        return R(out)

    __rmul__ = __mul__

    def __truediv__(self, rational):
        return self * (1 / Q(rational))

    def __bool__(self):
        return bool(self.a)

    def __eq__(self, other):
        return self.a == R(other).a

    def __float__(self):
        return sum((float(c) * math.sqrt(d) for d, c in self.a.items()), 0.0)

    def text(self):
        if not self.a:
            return "0"
        return " + ".join(str(c) if d == 1 else f"({c})*sqrt({d})"
                          for d, c in sorted(self.a.items()))


def sqrtq(q):
    q = Q(q)
    if not q:
        return R()
    assert q > 0
    out, inside = squarefree(q.numerator * q.denominator)
    return R({inside: Q(out, q.denominator)})


@lru_cache(None)
def cg(j1, m1, j2, m2, j, m):
    """Condon--Shortley CG; arguments are twice spin/magnetic number."""
    if m1 + m2 != m or abs(m) > j or j < abs(j1-j2) or j > j1+j2:
        return R()
    vals = (j+j1-j2, j-j1+j2, j1+j2-j, j1+j2+j+2,
            j+m, j-m, j1-m1, j1+m1, j2-m2, j2+m2)
    if any(v < 0 or v % 2 for v in vals):
        return R()
    a, b, c, den, *rest = [v//2 for v in vals]
    pref = Q((j+1)*math.factorial(a)*math.factorial(b)*math.factorial(c),
             math.factorial(den))
    for n in rest:
        pref *= math.factorial(n)
    total = Q(0)
    for k in range((j1+j2+j)//2 + 2):
        args2 = (2*k, j1+j2-j-2*k, j1-m1-2*k, j2+m2-2*k,
                 j-j2+m1+2*k, j-j1-m2+2*k)
        if any(v < 0 or v % 2 for v in args2):
            continue
        fac = math.prod(math.factorial(v//2) for v in args2)
        total += Q((-1)**k, fac)
    return sqrtq(pref) * total


def invariant(j, ell):
    out = {}
    for ms in itertools.product(range(-j, j+1, 2), repeat=4):
        if sum(ms):
            continue
        m = ms[0] + ms[1]
        coeff = (cg(j, ms[0], j, ms[1], 2*ell, m)
                 * cg(j, ms[2], j, ms[3], 2*ell, -m)
                 * cg(2*ell, m, 2*ell, -m, 0, 0))
        if coeff:
            out[ms] = coeff
    return out


def identity(n):
    return [[R(int(i == j)) for j in range(n)] for i in range(n)]


def transpose(a):
    return [list(r) for r in zip(*a)]


def multiply(a, b):
    return [[sum((x*y for x, y in zip(row, col)), R())
             for col in zip(*b)] for row in a]


def scale(a, s):
    return [[x*s for x in row] for row in a]


def plus(a, b):
    return [[x+y for x, y in zip(r, s)] for r, s in zip(a, b)]


def minus(a, b):
    return plus(a, scale(b, -1))


def zero(a):
    return all(not x for row in a for x in row)


def trace(a):
    return sum((a[i][i] for i in range(len(a))), R())


def numeric(a):
    return np.array([[float(x) for x in row] for row in a])


def texts(a):
    return [[x.text() for x in row] for row in a]


def inner(a, b):
    return sum((x*b.get(k, R()) for k, x in a.items()), R())


def source_maps():
    ein = [invariant(3, ell) for ell in range(4)]
    ep = [invariant(2, ell) for ell in range(3)]
    em = [invariant(1, ell) for ell in range(2)]
    for basis in (ein, ep, em):
        assert [[inner(x, y) for y in basis] for x in basis] == identity(len(basis))
    f = []
    for plus_state in ep:
        for minus_state in em:
            row = []
            for original in ein:
                val = R()
                for mp, cp in plus_state.items():
                    for mm, cm in minus_state.items():
                        ms = tuple(a+b for a, b in zip(mp, mm))
                        if ms not in original:
                            continue
                        y = R(1)
                        for a, b, c in zip(mp, mm, ms):
                            y *= cg(2, a, 1, b, 3, c)
                        val += cp*cm*y*original[ms]
                row.append(val)
            f.append(row)
    # Swap legs 2 and 3: independent exact recoupling on the same invariant space.
    swapped = [{(m[0], m[2], m[1], m[3]): c for m, c in e.items()} for e in ein]
    recoupling = [[inner(x, y) for y in swapped] for x in ein]
    assert multiply(transpose(recoupling), recoupling) == identity(4)
    return f, recoupling


def spectral_projectors(g, values):
    out = {}
    for val in values:
        p = identity(4)
        for other in values:
            if other != val:
                p = scale(multiply(p, minus(g, scale(identity(4), other))), 1/(val-other))
        assert multiply(p, p) == p
        assert multiply(g, p) == scale(p, val)
        out[val] = p
    assert sum_matrix(out.values(), 4) == identity(4)
    return out


def sum_matrix(seq, n):
    result = scale(identity(n), 0)
    for a in seq:
        result = plus(result, a)
    return result


def run():
    f, u = source_maps()
    ft, g = transpose(f), multiply(transpose(f), f)
    values = (Q(1, 3), Q(35, 81), Q(13, 27))
    projectors = spectral_projectors(g, values)
    multiplicities = {str(v): projectors[v] and trace(projectors[v]).text() for v in values}
    assert multiplicities == {'1/3': '1', '35/81': '1', '13/27': '2'}
    ginv = sum_matrix([scale(projectors[v], 1/v) for v in values], 4)
    invsqrt = sum_matrix([scale(projectors[v], sqrtq(1/v)) for v in values], 4)
    sqrtg = sum_matrix([scale(projectors[v], sqrtq(v)) for v in values], 4)
    vmap = multiply(f, invsqrt)
    assert multiply(transpose(vmap), vmap) == identity(4)
    edge = multiply(multiply(f, ginv), ft)
    assert multiply(edge, edge) == edge and edge == transpose(edge)
    assert multiply(edge, f) == f
    assert edge == multiply(vmap, transpose(vmap))
    a12 = [[R(i*(i+1) if i == j else 0) for j in range(4)] for i in range(4)]
    a13 = multiply(multiply(u, a12), transpose(u))
    comm12 = minus(multiply(g, a12), multiply(a12, g))
    comm13 = minus(multiply(g, a13), multiply(a13, g))
    # A12 has simple spectrum. A commuting matrix must be diagonal; the
    # nonzero entries of A13 connect its four diagonal coordinates.
    links = [(i, j) for i in range(4) for j in range(i+1, 4) if a13[i][j]]
    reached = {0}
    for _ in range(4):
        for i, j in links:
            if i in reached or j in reached:
                reached.update((i, j))
    assert len(reached) == 4
    assert not zero(comm12) and not zero(comm13)
    # In the polar-transported fixed A12 dictionary, raw successful filtering
    # can take an original A12=0 state into the L=2 (A12=6) outcome.
    mixing = sqrtg[2][0]*sqrtg[2][0] / Q(37, 81)
    # A retained coherent success/failure flag is an isometry. This is a
    # mathematical dilation certificate, not an implemented boundary time.
    fail = sum_matrix([scale(projectors[v], sqrtq(1-v)) for v in values], 4)
    assert plus(g, multiply(fail, fail)) == identity(4)
    # Extreme-eigenvalue Bell pair witness on an arbitrary reference.
    lo, hi = min(values), max(values)
    weights = (lo/(lo+hi), hi/(lo+hi))
    scalar_opt = 2/(lo+hi)
    scalar_error = (hi-lo)/(hi+lo)
    assert scalar_error == Q(2, 11)
    assert weights == (Q(9, 22), Q(13, 22))
    diagnostics = {
        'gram_eigenvalues': np.linalg.eigvalsh(numeric(g)).tolist(),
        'commutator_A12_operator_norm': float(np.linalg.norm(numeric(comm12), 2)),
        'commutator_A13_operator_norm': float(np.linalg.norm(numeric(comm13), 2)),
        'raw_A12_0_to_6_conditional_probability': float(mixing),
        'polar_isometry_floating_error': float(np.linalg.norm(numeric(vmap).T@numeric(vmap)-np.eye(4), 2)),
    }
    return {
        'status': 'PASS_FIXED_SPIN_ADMISSION_SCREEN',
        'not_a_numbered_research_round': True,
        'algorithm': 'Exact factorial CG sums and magnetic-index contraction; rational square-root arithmetic; no author code imported/read.',
        'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'parameters': {'gamma': '1/3', 'four_k': '3/2', 'j_plus': '1', 'j_minus': '1/2'},
        'basis': 'rows (Lplus,Lminus)=(0,0),(0,1),(1,0),(1,1),(2,0),(2,1); columns L=0,1,2,3; exact Condon--Shortley pairing convention',
        'exact': {
            'F_6_by_4': texts(f), 'G': texts(g), 'G_inverse': texts(ginv),
            'eigenvalue_multiplicities': multiplicities,
            'recoupling_swap_23': texts(u), 'A12': texts(a12), 'A13': texts(a13),
            'G_A12_commutator': texts(comm12), 'G_A13_commutator': texts(comm13),
            'joint_commutant_complex_dimension': 1,
            'A13_connectivity_edges': links,
            'raw_A12_0_to_6_conditional_probability': mixing.text(),
            'optimal_scalar_edge_weight': str(scalar_opt),
            'optimal_scalar_image_identity_operator_error': str(scalar_error),
            'extreme_Bell_success_probability': str((lo+hi)/2),
            'extreme_Bell_conditional_reference_weights': list(map(str, weights)),
            'extreme_Bell_reference_total_variation': str(abs(weights[0]-Q(1,2))),
            'checks': {
                'all_three_intertwiner_bases_orthonormal': True,
                'spectral_projectors_exact': True,
                'recoupling_orthogonal_exact': True,
                'polar_isometry_exact': True,
                'inverse_Gram_edge_projector_and_image_identity_exact': True,
                'coherent_failure_completion_isometry_exact': True,
            },
        },
        'floating_diagnostics_not_the_exact_certificate': diagnostics,
        'scope': [
            'Known EPRL nonisometry and inverse-Gram correction are reused, not discoveries.',
            'P_W F=F with injective F forces W G=I, hence W=G^-1; keeping all image states rules out reduced projectors.',
            'Joint A12/A13 selfadjointness in a changed positive input metric forces that metric to be scalar; G is not scalar.',
            'The scalar 2/11 error is an operator identity defect, not yet a probability or diamond distance.',
            'The Bell and A12 witnesses concern the declared successful filter with retained failure outcome; postselection is not deterministic transfer.',
            'Polar V transports every original finite input/effect and arbitrary reference mathematically; target geometric effects must be stated as V O Vdagger.',
            'No clock, causal boundary implementation, vertex dynamics, local source, GR selection or natural instrument is constructed.',
        ],
        'primary_sources': {
            'EPRL_map_and_inverse_Gram': 'https://arxiv.org/html/0912.0540v2#S2.SS5',
            'injectivity_admissibility': 'https://arxiv.org/abs/1109.5023',
        },
    }


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for k in a:
            compare(a[k], b[k])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            compare(x, y)
    elif isinstance(a, float):
        assert abs(a-b) <= 2e-13
    else:
        assert a == b, (a, b)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    path = Path(__file__).with_name('review_probe_results.json')
    result = run()
    # JSON round trip converts tuple-only presentation fields to lists.
    result = json.loads(json.dumps(result))
    if args.write:
        with path.open('x', encoding='utf-8', newline='\n') as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write('\n')
        print('WRITE_ONCE_PASS')
    else:
        compare(result, json.loads(path.read_text(encoding='utf-8')))
        print('READ_ONLY_MATCH')
    print(json.dumps({'eigenvalues': result['exact']['eigenvalue_multiplicities'],
                      'joint_commutant': 1, 'scalar_identity_error': '2/11'}, ensure_ascii=False))


if __name__ == '__main__':
    main()
