"""Independent 1063 review checks; default prints results without writing files.

This preserves the integer-epsilon, hook/Weyl, and independent Clebsch-basis
algorithms actually run during mathematical review. It imports no author code
and adds no scientific calibration group. --save-exclusive creates only its own
result file and refuses to overwrite an existing file.
"""
import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np


def signed_perms(d):
    for p in itertools.permutations(range(d)):
        yield p, (-1) ** sum(
            p[i] > p[j] for i in range(d) for j in range(i + 1, d)
        )


def run():
    exact = []
    for d in range(2, 6):
        e = {}
        for r in range(d + 1):
            for g in range(d):
                v = {}
                for p, sgn in signed_perms(d):
                    x = list(p)
                    x.insert(r, g)
                    v[tuple(x)] = (-1) ** r * sgn
                e[g, r] = v
        pairs = 0
        for g, h in itertools.product(range(d), repeat=2):
            for r, s in itertools.product(range(d + 1), repeat=2):
                dot = sum(
                    value * e[h, s].get(x, 0) for x, value in e[g, r].items()
                )
                want = (
                    math.factorial(d) if r == s else -math.factorial(d - 1)
                ) if g == h else 0
                assert dot == want, (d, g, h, r, s, dot, want)
                pairs += 1
        swaps = 0
        for a, b in itertools.combinations(range(d + 1), 2):
            for (g, r), v in e.items():
                moved = {}
                for x, coef in v.items():
                    y = list(x)
                    y[a], y[b] = y[b], y[a]
                    moved[tuple(y)] = coef
                rp = b if r == a else a if r == b else r
                assert moved == {x: -coef for x, coef in e[g, rp].items()}
                swaps += 1
        exact.append(dict(
            d=d,
            integer_gram_comparisons=pairs,
            signed_swap_comparisons=swaps,
        ))

    for d in range(2, 51):
        lam = [2] + [1] * (d - 1)
        hook = 1
        for i in range(d):
            for j in range(lam[i]):
                hook *= lam[i] - j + sum(
                    lam[k] > j for k in range(i + 1, d)
                )
        specht = math.factorial(d + 1) // hook
        weyl = Fraction(1)
        for i in range(d):
            for j in range(i + 1, d):
                weyl *= Fraction(lam[i] - lam[j] + j - i, j - i)
        assert specht == d and weyl == d

    # Independent Clebsch basis: singlet-pair and triplet-pair coupling.
    w = np.zeros((8, 4))
    w[2, 0] = 1 / math.sqrt(2)
    w[4, 0] = -1 / math.sqrt(2)
    w[1, 1] = 2 / math.sqrt(6)
    w[2, 1] = w[4, 1] = -1 / math.sqrt(6)
    w[3, 2] = 1 / math.sqrt(2)
    w[5, 2] = -1 / math.sqrt(2)
    w[3, 3] = w[5, 3] = 1 / math.sqrt(6)
    w[6, 3] = -2 / math.sqrt(6)
    child = w.reshape(8, 2, 2)
    parent = w.reshape(2, 2, 2, 2, 2)
    w2 = np.einsum(
        'xia,yjb,zkc,ijkgm->xyzgmabc', child, child, child, parent
    ).reshape(512, 32)
    iso = float(np.max(np.abs(w2.T @ w2 - np.eye(32))))
    tensor = w2.reshape((2,) * 9 + (32,))
    rawswap = np.swapaxes(tensor, 0, 3).reshape(512, 32)
    leak = rawswap - w2 @ (w2.T @ rawswap)
    leak_norm = float(np.linalg.norm(leak))
    uniform = sum(
        np.swapaxes(tensor, i, j) for i in range(3) for j in range(3, 6)
    ).reshape(512, 32)
    uniform_leak = float(np.max(np.abs(uniform - w2 @ (w2.T @ uniform))))
    inner = np.swapaxes(tensor, 0, 1).reshape(512, 32)
    inner_leak = float(np.max(np.abs(inner - w2 @ (w2.T @ inner))))
    wp = np.kron(w2, np.eye(2))
    tp = wp.reshape((2,) * 10 + (2,) * 6)
    lhs = sum(np.swapaxes(tp, i, 9) for i in range(9))
    rhs = 4 * tp + np.swapaxes(tp, 10, 15)
    external_residual = float(np.max(np.abs(lhs - rhs)))
    assert max(iso, uniform_leak, inner_leak, external_residual) < 2e-12
    assert leak_norm > 1
    return dict(
        exact_epsilon_checks=exact,
        hook_and_weyl_dimensions_verified_through=50,
        recursive_qubit_depth2=dict(
            shape=list(w2.shape),
            isometry=iso,
            allowed_uniform_child_contact_leakage=uniform_leak,
            allowed_child_internal_swap_leakage=inner_leak,
            root_to_primitive_intertwiner=external_residual,
            excluded_single_cross_child_swap_leakage_frobenius=leak_norm,
        ),
        passed=True,
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--save-exclusive', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.save_exclusive:
        destination = Path(__file__).with_name('independent_results.json')
        with destination.open('x', encoding='utf8') as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
    print(json.dumps(result, ensure_ascii=False))
