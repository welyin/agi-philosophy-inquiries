"""Exact independent rational-CG library reused for round 1182.

Copied from the frozen 1181 independent algebra, SHA256
2faa02174ac81202eac342c335e854c3925c8acec2d1f74d0c34a52f2688e275.
The 1182 Peter--Weyl construction is in reviewer_check.py.
No author code is imported; this library performs no writes.
"""
from __future__ import annotations

from fractions import Fraction as Q
from functools import lru_cache
import itertools
import math
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

