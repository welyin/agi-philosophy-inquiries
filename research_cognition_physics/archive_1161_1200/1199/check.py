"""Exact finite source-identity audit of the Gu--Wen L-type lattice.

Not a spectral/continuum/Einstein reconstruction. Uses finite-field stabilizer
traces, never constructs a 5**384 matrix. Default prints; --write saves only a
new result, refusing to replace a different existing result.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np


def lattice(side=4):
    assert side >= 4 and side % 2 == 0
    sites = list(itertools.product(range(side), repeat=3))
    cells = {r: i for i, r in enumerate(sites)}
    pairs = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)]
    types = {ab: k for k, ab in enumerate(pairs)}
    n = 6 * len(sites)
    S, V, R, T, rweights, tweights = [], [], [], [], [], []

    def shift(r, *axes):
        out = list(r)
        for axis, sign in axes:
            out[axis] += sign
        return tuple(x % side for x in out)

    def add(row, r, a, b, c):
        row[6 * cells[tuple(x % side for x in r)] + types[tuple(sorted((a, b)))]] += c

    for i in sites:
        # Scalar stencil: diagonal transverse second sums minus face sums.
        s = np.zeros(n, dtype=np.int64)
        for a in range(3):
            for b in range(3):
                if a != b:
                    add(s, shift(i, (a, 1)), b, b, 1)
                    add(s, shift(i, (a, -1)), b, b, 1)
                    add(s, i, b, b, 2)
        for a, b in pairs[3:]:
            for da, db in itertools.product((0, -1), repeat=2):
                add(s, shift(i, (a, da), (b, db)), a, b, -1)
        S.append(s)
        for a in range(3):
            v = np.zeros(n, dtype=np.int64)
            add(v, i, a, a, 1)
            add(v, shift(i, (a, 1)), a, a, 1)
            for b in range(3):
                if a != b:
                    add(v, i, a, b, 1)
                    add(v, shift(i, (b, -1)), a, b, 1)
            V.append(v)
        for a, b, c in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
            # R_aa, cyclic version of the original R_zz formula.
            rd = np.zeros(n, dtype=np.int64)
            for u, v in ((b, c), (c, b)):
                for d, k in ((1, 1), (-1, 1), (0, 2)):
                    add(rd, shift(i, (v, d)), u, u, k)
            for db, dc in itertools.product((0, -1), repeat=2):
                add(rd, shift(i, (b, db), (c, dc)), b, c, -1)
            R.append(rd)
            rweights.append(F(-1, 2))
            # R_bc, cyclic version of original equation (47).
            ro = np.zeros(n, dtype=np.int64)
            for db, dc in itertools.product((0, 1), repeat=2):
                add(ro, shift(i, (b, db), (c, dc)), a, a, 2)
            for d in (0, -1):
                for e in (0, 1):
                    add(ro, shift(i, (a, d), (c, e)), a, b, -1)
                    add(ro, shift(i, (a, d), (b, e)), a, c, -1)
            for d, k in ((1, 1), (-1, 1), (0, 2)):
                add(ro, shift(i, (a, d)), b, c, k)
            R.append(ro)
            rweights.append(F(-1, 4))
            # sigma^a_a, sigma^a_b, sigma^a_c: equations (49)--(51).
            td = np.zeros(n, dtype=np.int64)
            add(td, i, a, c, -1)
            add(td, shift(i, (b, 1)), a, c, -1)
            add(td, i, a, b, 1)
            add(td, shift(i, (c, 1)), a, b, 1)
            tb = np.zeros(n, dtype=np.int64)
            add(tb, i, b, c, 2)
            add(tb, shift(i, (b, -1)), b, c, 2)
            for d in (0, 1):
                for k, weight in ((b, -1), (a, 1), (c, 1)):
                    add(tb, shift(i, (c, d)), k, k, weight)
            tc = np.zeros(n, dtype=np.int64)
            for d in (0, 1):
                for k, weight in ((c, 1), (a, -1), (b, -1)):
                    add(tc, shift(i, (b, d)), k, k, weight)
            add(tc, i, b, c, -2)
            add(tc, shift(i, (c, -1)), b, c, -2)
            T.extend((td, tb, tc))
            tweights.extend((F(-1, 2), F(-1, 8), F(-1, 8)))
    return tuple(np.array(x, dtype=np.int64) for x in (S, V, R, T)), rweights, tweights


def rref(mat, p):
    a = mat.copy() % p
    rank, pivots = 0, []
    for col in range(a.shape[1]):
        nz = np.flatnonzero(a[rank:, col])
        if not len(nz):
            continue
        row = rank + int(nz[0])
        a[[rank, row]] = a[[row, rank]]
        a[rank] = a[rank] * pow(int(a[rank, col]), -1, p) % p
        coeff = a[:, col].copy()
        coeff[rank] = 0
        a = (a - coeff[:, None] * a[rank]) % p
        pivots.append(col)
        rank += 1
        if rank == len(a):
            break
    return a[:rank], pivots


def quotient(rows, basis, pivots, p):
    out = rows.copy() % p
    for row, col in zip(basis, pivots):
        out = (out - out[:, col:col+1] * row) % p
    return out


def moments(rows, weights, basis, pivots, p):
    reduced = quotient(rows, basis, pivots, p)
    grouped = defaultdict(F)
    for v, w in zip(reduced, weights):
        # Each row denotes one unitary, with its Hermitian conjugate included.
        for sign in (-1, 1):
            grouped[tuple(int(x) for x in sign * v % p)] += p * w
    zero = (0,) * rows.shape[1]
    mean = grouped.pop(zero, F(0))
    assert all(grouped[tuple((-np.array(k)) % p)] == val for k, val in list(grouped.items()))
    variance = sum((c*c for c in grouped.values()), F(0))
    return mean, variance, grouped, reduced


def check(side=4, p=5):
    assert p > 2 and all(p % d for d in range(2, int(p**0.5)+1))
    (S, V, R, T), rw, tw = lattice(side)
    n = S.shape[1]
    # Exact INTEGER zero is stronger than zero modulo the qudit prime.
    commutators = {
        'S_V': int(np.max(np.abs(S @ V.T))),
        'R_V': int(np.max(np.abs(R @ V.T))),
        'T_S': int(np.max(np.abs(T @ S.T))),
    }
    assert set(commutators.values()) == {0}
    assert np.array_equal(R[0::6] + R[2::6] + R[4::6], S)
    sx, sp = rref(S, p)
    vz, vp = rref(V, p)
    mg, vg, cg, rg = moments(R, rw, sx, sp, p)
    mj, vj, cj, tj = moments(T, tw, vz, vp, p)
    assert vg > 0 and vj > 0
    assert mg == mj == 0
    # Pure X and pure Z nonidentity logical Weyls are trace-orthogonal.
    # Thus full-H variance is vg*g**2 + vj*J**2, with no cross term.
    variance = vg + vj  # certificate instance J=g=1; arbitrary U1,U2>0
    assert F(77) ** 2 < variance
    dynamic_norm_bound = 2*p*(sum(abs(w) for w in rw) + sum(abs(w) for w in tw))
    assert dynamic_norm_bound == 2880
    # For p=5, min nonzero HU = 2000*5*(1-cos(2*pi/5)).
    # cos(2*pi/5)=(sqrt(5)-1)/4 and sqrt(5)<9/4.
    assert p == 5 and F(9, 4)**2 > 5
    penalty_gap_lower = 2000*p*(1-(F(9, 4)-1)/4)
    isolated_code_gap_lower = penalty_gap_lower - 2*dynamic_norm_bound
    assert isolated_code_gap_lower == 1115
    # Independent commutator witness: two actual Hamiltonian terms do not commute.
    nz = np.argwhere((R @ T.T) % p != 0)
    assert len(nz)
    a, b = map(int, nz[0])
    assert np.any(rg[a]) and np.any(tj[b])
    return {
        'schema': 'finite_l_type_complete_energy_source_v1',
        'status': 'PASS_EXACT_FINITE_ALGEBRA',
        'box': [side]*3, 'qudit_prime': p, 'qudits': n,
        'full_hilbert_dimension': f'{p}^{n}',
        'constraint_generators': {'scalar': len(S), 'vector': len(V)},
        'constraint_ranks': {'scalar': len(sx), 'vector': len(vz)},
        'code_dimension': f'{p}^{n-len(sx)-len(vz)}',
        'integer_commutators_max_abs': commutators,
        'hamiltonian_rows_without_hc': {'R': len(R), 'T': len(T)},
        'nonidentity_cosets': {'R': len(cg), 'T': len(cj)},
        'normalized_code_trace_mean_Hg_per_g': str(mg),
        'normalized_code_trace_mean_HJ_per_J': str(mj),
        'variance_g_squared_coefficient': str(vg),
        'variance_J_squared_coefficient': str(vj),
        'mixed_covariance': '0',
        'J_g_one_variance': str(variance),
        'J_g_one_error_strict_lower': '77',
        'positive_shifted_finite_band': {
            'J': '1', 'g': '1', 'U1': '2000', 'U2': '2000',
            'dynamic_operator_norm_upper': str(dynamic_norm_bound),
            'constant_positive_shift': str(dynamic_norm_bound),
            'whole_code_energy_interval': ['0', str(2*dynamic_norm_bound)],
            'outside_code_energy_strict_lower': str(penalty_gap_lower),
            'band_separation_strict_lower': str(isolated_code_gap_lower),
            'qualifier': 'The whole protected band is isolated by a chosen finite penalty; '
                         'no variance bound is claimed for every narrower energy window.'},
        'logical_noncommutation_witness': {'R_row': a, 'T_row': b,
             'pairing_mod_p': int(R[a] @ T[b] % p)},
        'stencils_sha256': {key: hashlib.sha256(arr.astype('<i8').tobytes()).hexdigest()
                           for key, arr in zip(('S','V','R','T'), (S,V,R,T))},
        'check_groups': 2,
        'scope': 'Full finite L-type Hamiltonian on its invariant defect-free sector, '
                 'with one finite-penalty isolated-band example; '
                 'not a bound on every narrower low-energy window, graviton dispersion certificate, '
                 'physical instrument, gravitational source dictionary or dark-sector identification.',
        'physical_inputs': {'lattice_and_rule': 'Gu-Wen L type; cubic periodic box',
                            'couplings': 'J,g positive; U1,U2 arbitrary positive',
                            'energy_unit': 'J=g=1 only for the reported numerical threshold'},
        'claim': 'Any self-adjoint source that is a function only of all S,V labels '
                 'is scalar on the code. For any nonzero real a, its best additive-constant '
                 'operator-norm approximation to a*H has error at least abs(a)*sqrt(variance).',
    }


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--write', action='store_true')
    args = ap.parse_args()
    result = check()
    raw = json.dumps(result, ensure_ascii=False, indent=2) + '\n'
    if args.write:
        path = Path(__file__).with_name('results.json')
        if path.exists() and json.loads(path.read_text(encoding='utf-8')) != result:
            raise RuntimeError('Refusing to replace a different recorded result')
        if not path.exists():
            path.write_text(raw, encoding='utf-8')
    print(raw, end='')
