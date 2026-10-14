"""1180: exact finite algebra and outward rational bounds. Default is read-only.

No author/reviewer module imports, oscillator truncation, or floating eigensolver.
The proof supplies the all-jet/all-band/all-input quantifiers; finite checks audit
the exact identities and explicit constants used by that proof.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement, product
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def exact(q):
    q = F(q)
    return f"{q.numerator}/{q.denominator}"


def float_down(q):
    q = F(q)
    x = float(q)
    while F.from_float(x) > q:
        x = math.nextafter(x, -math.inf)
    return x


def float_up(q):
    q = F(q)
    x = float(q)
    while F.from_float(x) < q:
        x = math.nextafter(x, math.inf)
    return x


def interval(lo, hi=None):
    hi = lo if hi is None else hi
    assert F(lo) <= F(hi)
    return {"lo_exact": exact(lo), "hi_exact": exact(hi),
            "lo_display": float_down(lo), "hi_display": float_up(hi)}


def sqrt_bounds(q, digits=45):
    q = F(q)
    assert q >= 0
    scale = 10**digits
    a = math.isqrt(q.numerator * scale**2 // q.denominator)
    lo, hi = F(a, scale), F(a + 1, scale)
    if lo * lo == q:
        hi = lo
    assert lo * lo <= q <= hi * hi
    return lo, hi


def rank(rows):
    a = [[F(x) for x in row] for row in rows]
    r = 0
    for c in range(len(a[0])):
        piv = next((j for j in range(r, len(a)) if a[j][c]), None)
        if piv is None:
            continue
        a[r], a[piv] = a[piv], a[r]
        d = a[r][c]
        a[r] = [x / d for x in a[r]]
        for j in range(len(a)):
            if j != r and a[j][c]:
                d = a[j][c]
                a[j] = [x - d * y for x, y in zip(a[j], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def matmul(a, b):
    return [[sum(x * y for x, y in zip(row, col))
             for col in zip(*b)] for row in a]


def jet_group():
    pair = list(combinations_with_replacement(range(3), 2))
    zeroth_map = []
    for i, j in pair:
        row = [0] * 9
        row[3*i+j] += 1
        row[3*j+i] += 1
        zeroth_map.append(row)
    assert rank(zeroth_map) == 6
    first = [(k, i, j) for k in range(3) for i, j in pair]
    b_basis = [(i, j, k) for i in range(3) for j, k in pair]
    b_index = {b: n for n, b in enumerate(b_basis)}
    first_map = []
    for k, i, j in first:
        row = [0] * 18
        row[b_index[(j, *sorted((k, i)))]] += 1
        row[b_index[(i, *sorted((k, j)))]] += 1
        first_map.append(row)
    fi = {key: n for n, key in enumerate(first)}
    inverse = []
    for i, j, k in b_basis:
        row = [F(0)] * 18
        for key, coeff in [((j, *sorted((i, k))), F(1, 2)),
                           ((k, *sorted((i, j))), F(1, 2)),
                           ((i, *sorted((j, k))), F(-1, 2))]:
            row[fi[key]] += coeff
        inverse.append(row)
    assert matmul(inverse, first_map) == [[int(i == j) for j in range(18)] for i in range(18)]

    second = [(a, b, i, j) for a, b in pair for i, j in pair]
    third = [(i, *jkl) for i in range(3)
             for jkl in combinations_with_replacement(range(3), 3)]
    ti = {key: n for n, key in enumerate(third)}
    second_map = []
    for a, b, i, j in second:
        row = [0] * 30
        row[ti[(j, *sorted((a, b, i)))]] += 1
        row[ti[(i, *sorted((a, b, j)))]] += 1
        second_map.append(row)
    si = {key: n for n, key in enumerate(second)}
    ricci2 = []
    for i, j in pair:
        row = [0] * 36
        for l in range(3):
            for a, b, c, d, sign in [(l, i, l, j, 1), (l, j, l, i, 1),
                                      (l, l, i, j, -1), (i, j, l, l, -1)]:
                row[si[(*sorted((a, b)), *sorted((c, d)))]] += sign
        ricci2.append(row)
    assert matmul(ricci2, second_map) == [[0] * 30 for _ in range(6)]
    ranks = [rank(first_map), rank(second_map), rank(ricci2)]
    assert ranks == [18, 30, 6]
    t = F(2, 3)
    best_error = max(abs(t - 1), abs(2 * t - 1))
    assert best_error == F(1, 3)
    return {"first_jet_gauge_rank": 18, "first_jet_inverse_exact": True,
            "zeroth_jet_gauge_rank": rank(zeroth_map),
            "second_jet_dimension": 36, "third_gauge_jet_rank": 30,
            "twice_Ricci_rank": 6, "curvature_times_gauge_zero": True,
            "TT_frequency_ratio_at_2k": 4,
            "optimal_bk_over_c": exact(t), "best_two_point_relative_error": exact(best_error)}


def escape_group():
    # c=1, k0=1. C=2 Ric^(1), hence C_TT(k)=k^2 h_TT(k).
    def omega_squared(s):
        return F(9, 7) * s * s * ((1 - s / 6)**2 + F(1, 12))
    w1, w2 = omega_squared(F(1)), omega_squared(F(4))
    assert (w1, w2) == (1, 4)
    return {"first_square_multiplier": "1 + Delta/6",
            "positive_base": exact(F(1, 12)), "overall_potential_factor": exact(F(9, 7)),
            "omega_squared_k1": exact(w1), "omega_squared_k2": exact(w2),
            "relative_error_at_both_points": exact(0),
            "strict_TT_positive_for_every_nonzero_k": True,
            "cross_derivative_orders_within_one_invariant_square": True}


def band_group():
    degree, delta, nmax, time = 32, F(1, 10**12), 1, F(1)
    coeff = [F(math.comb(2*j, j), 4**j) for j in range(degree + 1)]
    assert all(0 < c <= 1 for c in coeff)
    for j in range(degree):
        assert coeff[j+1] == coeff[j] * F(2*j+1, 2*j+2)
    eps_poly = F(16, 5) * F(3, 5)**(degree + 1)
    eps = eps_poly + 3 * delta
    assert eps < F(1, 10**6)
    # The added base shifts sqrt(r^2+s delta) by <=2 delta/(1-eps_poly)<3 delta.
    assert 2 / (1 - eps_poly) < 3
    shells = Counter(sum(x*x for x in k) for k in product(range(-2, 3), repeat=3)
                     if 0 < sum(x*x for x in k) <= 4)
    assert shells == Counter({1: 6, 2: 12, 3: 8, 4: 6})
    real_modes = 2 * sum(shells.values())
    assert real_modes == 64
    sum_lo = sum_hi = F(0)
    shell_audits = []
    for s, number_of_vectors in sorted(shells.items()):
        x = 1 - F(2*s, 5)
        poly = sum(c * x**j for j, c in enumerate(coeff))
        ratio_sq = s * (F(2, 5) * poly * poly + delta)
        assert (1 - eps)**2 <= ratio_sq <= (1 + eps)**2
        k_lo, k_hi = sqrt_bounds(s)
        sum_lo += 2 * number_of_vectors * k_lo
        sum_hi += 2 * number_of_vectors * k_hi
        rlo, rhi = sqrt_bounds(ratio_sq)
        shell_audits.append({"k_squared": s, "integer_vectors": number_of_vectors,
                             "real_TT_oscillators": 2*number_of_vectors,
                             "omega_over_k": interval(rlo, rhi)})
    # Public certificate uses the uniform 64 modes * omega<=2 bound, not an optimized sum.
    assert sum_hi < 128
    eta_hi = time * (eps + eps*eps/2) * (2*nmax + 3) * 128
    assert eta_hi < F(1, 10000)
    # Q^2 creates two extra quanta: this certificate keeps them, it never compresses H.
    root6lo, root6hi = sqrt_bounds(6)
    q2_triangle_hi = (root6hi + 3) / 2
    assert q2_triangle_hi < 2*nmax + 3
    return {"band_k": [1, 2], "degree_in_k_squared": degree,
            "maximum_h_derivative_order_in_square": 2 + 2*degree,
            "delta": exact(delta), "coefficients_S32": [exact(c) for c in coeff],
            "uniform_relative_frequency_error": interval(0, eps),
            "polynomial_error_bound": interval(0, eps_poly),
            "all_band_bound_is_analytic_not_sampling": True,
            "box_side": "2*pi", "cutoff": "0 < |k| <= 2", "shells": shell_audits,
            "real_TT_modes": real_modes, "reference_nmax_each_mode": nmax,
            "input_code_dimension": str(2**real_modes), "time_limit": exact(time),
            "sum_reference_frequencies": interval(sum_lo, sum_hi),
            "Q2_n1_triangle_bound_times_omega": interval((root6lo+3)/2, q2_triangle_hi),
            "Q2_declared_bound_times_omega": 2*nmax+3,
            "public_frequency_sum_upper": 128,
            "prepared_full_Fock_distance_bound": interval(0, eta_hi),
            "public_full_Fock_distance_upper": exact(F(1, 10000)),
            "same_canonical_Q_P_and_reference_input_code": True,
            "pair_creation_and_all_Fock_output_retained": True,
            "arbitrary_reference_and_common_output_instrument": True,
            "zero_modes_if_present_common_free_spectator": True,
            "strict_lattice_locality_or_infinite_UV_claim": False}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def calculate():
    text_assets = [ROOT / 'research_note_1180.md', HERE / 'admission.md',
                   HERE / 'proof.md', HERE / 'NEXT.md', HERE / 'check.py']
    history = [ROOT.parent / 'archive_342_369' / 'research_note_343.md',
               ROOT.parent / 'archive_342_369' / 'research_note_358.md',
               ROOT.parent / 'archive_702_741' / 'research_note_739.md']
    return {"round": 1180, "status": "PASS_WITH_SCOPE", "groups": 3,
            "method": "exact rational jet algebra and analytic band/prepared-Duhamel bounds",
            "author_asset_sha256_excluding_results":
                {str(p.relative_to(ROOT)).replace('\\', '/'): digest(p) for p in text_assets},
            "history_attribution_sha256":
                {str(p.relative_to(ROOT.parent)).replace('\\', '/'): digest(p) for p in history},
            "primary_source": {"url": "https://arxiv.org/html/0907.1203v3",
                               "sections": "IV.2-IV.4, VI-VII",
                               "role": "mature context; not a proof of the present finite Fock certificate"},
            "g1_low_order_class": jet_group(), "g2_strict_positive_two_point_escape": escape_group(),
            "g3_uniform_band_same_canonical_quantum_process": band_group()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true', help='write results once; never overwrite')
    args = parser.parse_args()
    result = calculate()
    target = HERE / 'results.json'
    if args.write:
        with target.open('x', encoding='utf-8', newline='\n') as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write('\n')
        label = 'WRITE_ONCE_PASS'
    else:
        old = json.loads(target.read_text(encoding='utf-8'))
        assert old == result, 'read-only replay differs from frozen result'
        label = 'READ_ONLY_MATCH'
    eta = result['g3_uniform_band_same_canonical_quantum_process']['prepared_full_Fock_distance_bound']['hi_display']
    print(f'{label}: 3 groups; full_Fock_distance_upper={eta:.17g}')


if __name__ == '__main__':
    main()
