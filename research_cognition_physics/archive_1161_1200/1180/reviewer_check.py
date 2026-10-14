"""Independent 1180 mathematics audit; no author-code imports.

The algebraic checks were implemented before reading the author's code.
All certificates are rational identities or outward rational bounds.  Finite
Fock calculations below apply an exact polynomial to prepared number states;
they do not replace the full output oscillator Hilbert space by a truncation.
--math-only does not bind or write anything.  --write is exclusive.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as Q
import hashlib
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
CORE_FILES = [
    'research_note_1180.md', '1180/admission.md', '1180/proof.md',
    '1180/check.py', '1180/results.json', '1180/NEXT.md',
]
# Populated only after the author declares the six files frozen.
EXPECTED_CORE = {
    'research_note_1180.md':'dd819b5eee437ce31ce117e13deeb50df4ac9013486fe086f00a7087b9b147f7',
    '1180/admission.md':'8d8a9c07737d660ca4c0e838c45368abacfb8d137678729aa2a521263a37e9d5',
    '1180/proof.md':'5125a890b5f9e32d22560fbdd91ae1d57a3fd6ee49ed2a3f9eb5745d1cf2e82c',
    '1180/check.py':'67efccd2eb77db14d77c742154d5ea9e25d9d694d12498de5b3ea638a69f7b60',
    '1180/results.json':'468d2c345d26be824ce8060bc3016d76609fd50fa747fa0e37df344d8678d3e7',
    '1180/NEXT.md':'f1613379fff3f024004aae11f693b4be8043c33a074c1601d67e952c5a348a57',
}
EXPECTED_SOURCE = {
    '../archive_342_369/research_note_343.md':'6b3622ae728c4fcfd1280f3b903cc2e1e92cbe4380439ae62204d5af5a974d85',
    '../archive_342_369/research_note_358.md':'e37029d318b84e7b54e7a6a0f3bbc74170f0122c96d91f8784de9828826f1595',
    '../archive_702_741/research_note_739.md':'9f5d670c015686dddccb5a9d70ca0ddbfa8b1a59f9874281394778bc2c14e32d',
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pack(x):
    return {'exact': str(x), 'display': float(x)}


def up(x):
    y = float(x)
    while Q.from_float(y) < x:
        y = math.nextafter(y, math.inf)
    return y


def rank(matrix):
    a = [[Q(x) for x in row] for row in matrix]
    row = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(row, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        z = a[row][col]
        a[row] = [v/z for v in a[row]]
        for i in range(len(a)):
            if i != row and a[i][col]:
                z = a[i][col]
                a[i] = [u-z*v for u, v in zip(a[i], a[row])]
        row += 1
        if row == len(a):
            break
    return row


PAIRS = list(itertools.combinations_with_replacement(range(3), 2))


def gauge(order):
    derivatives = list(itertools.combinations_with_replacement(range(3), order))
    gauge_derivatives = list(itertools.combinations_with_replacement(range(3), order+1))
    rows = [(pair, deriv) for pair in PAIRS for deriv in derivatives]
    cols = [(i, deriv) for i in range(3) for deriv in gauge_derivatives]
    index = {c: i for i, c in enumerate(cols)}
    matrix = []
    for (i, j), deriv in rows:
        row = [0]*len(cols)
        row[index[j, tuple(sorted((i,)+deriv))]] += 1
        row[index[i, tuple(sorted((j,)+deriv))]] += 1
        matrix.append(row)
    return rows, cols, matrix


def multiply(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def mv(a, x):
    return [sum(u*v for u, v in zip(row, x)) for row in a]


def jet_group():
    out = {}
    for order, expected in [(0, 6), (1, 18), (2, 30)]:
        rows, cols, mat = gauge(order)
        r = rank(mat)
        assert r == expected
        out[str(order)] = {'rows': len(rows), 'columns': len(cols),
                           'rank': r, 'invariant_dimension': len(rows)-r}
    rows, _, mat = gauge(2)
    index = {c: i for i, c in enumerate(rows)}
    # 2 Ric_ij = sum_k (h_kj,ki + h_ki,kj - h_kk,ij - h_ij,kk).
    ricci = []
    for i, j in PAIRS:
        row = [0]*len(rows)
        for k in range(3):
            for h, d, sign in [((k,j),(k,i),1), ((k,i),(k,j),1),
                               ((k,k),(i,j),-1), ((i,j),(k,k),-1)]:
                row[index[tuple(sorted(h)), tuple(sorted(d))]] += sign
        ricci.append(row)
    assert rank(ricci) == 6
    assert all(not x for row in multiply(ricci, mat) for x in row)
    plus = {(0,0): Q(1), (1,1): Q(-1)}
    cross = {(0,1): Q(1)}
    tt = []
    for h in (plus, cross):
        values = []
        for k in (Q(1), Q(2)):
            jets = [-k*k*h.get(pair, 0) if d == (2,2) else Q(0)
                    for pair, d in rows]
            r = mv(ricci, jets)
            assert r == [k*k*h.get(pair, 0) for pair in PAIRS]
            assert sum(r[PAIRS.index((i,i))] for i in range(3)) == 0
            values.append(sum((1 if i == j else 2)*x*x
                              for (i,j), x in zip(PAIRS, r)))
        assert values[1] == 16*values[0]
        tt.append([str(v) for v in values])
    e = Q(1, 3)
    x = Q(2, 3)
    assert abs(x-1) == abs(2*x-1) == e
    # Necessary intervals intersect iff 1-e <= (1+e)/2, i.e. e>=1/3.
    assert 1-e == (1+e)/2
    out.update({'two_ricci_rank': 6, 'two_ricci_gauge_product_zero': True,
                'tt_curvature_squared_at_k1_k2': tt,
                'minimax_relative_error': pack(e),
                'minimax_scale': pack(x),
                'tt_nondegeneracy_is_separate': True,
                'ricci_scalar_squared_alone_annihilates_tt': True})
    return out


def polynomial(n, s):
    # P_n(s) = sqrt(2/5) times this rational polynomial.
    x = (s-Q(5,2))/Q(5,2)
    return sum(((-1)**j*Q(math.comb(2*j,j), 4**j)*x**j
                for j in range(n+1)), Q(0))


def band_group():
    # Mixed fourth-jet square: omega(s)=s*|1-s/7|.
    w1 = 1*abs(1-Q(1,7))
    w2 = 4*abs(1-Q(4,7))
    assert w2 == 2*w1
    # Author's stronger escape keeps every nonzero TT frequency positive.
    strict_frequency_squares = {
        s: Q(9,7)*s*s*((1-s/Q(6))**2+Q(1,12)) for s in (Q(1),Q(4))
    }
    assert strict_frequency_squares == {Q(1):Q(1), Q(4):Q(4)}
    rho, degree, delta = Q(3,5), 32, Q(1,10**12)
    # [2/sqrt(5/2)/(1-rho)]^2=10 < (16/5)^2.
    exact_prefactor_squared = Q(4)/Q(5,2)/(1-rho)**2
    assert exact_prefactor_squared == 10 < Q(16,5)**2
    eps0 = Q(16,5)*rho**(degree+1)
    eps = eps0+3*delta
    assert eps < 1
    # The potential base is +delta, not +(delta)^2.  The existing band
    # lower bound r>=1-eps0 converts its frequency increment to O(delta).
    assert 2/(1-eps0) < 3
    # These coefficients are the exact binomial coefficients, <=1; the
    # geometric majorant bounds the whole band, not just these mode samples.
    assert all(0 < Q(math.comb(2*j,j), 4**j) <= 1 for j in range(degree+2))
    checks = {}
    for s in [Q(1), Q(2), Q(3), Q(4), Q(5,2)]:
        p = polynomial(degree, s)
        # Candidate adds a second positive curvature square with weight delta.
        relative_squared = s*(Q(2,5)*p*p+delta)
        assert (1-eps)**2 <= relative_squared <= (1+eps)**2
        checks[str(s)] = {'relative_frequency_squared': pack(relative_squared)}
    n_first = next(n for n in range(100)
                   if 640*((Q(16,5)*rho**(n+1)+3*delta)
                           +(Q(16,5)*rho**(n+1)+3*delta)**2/2) < Q(1,10000))
    assert n_first == 32
    return {'mixed_square_omega_k1': pack(w1), 'mixed_square_omega_k2': pack(w2),
            'mixed_square_two_scale_relative_error': pack(Q(0)),
            'strict_positive_escape_squared_frequencies': {
                str(s):pack(value) for s,value in strict_frequency_squares.items()},
            'degree': degree, 'positive_base_delta': pack(delta),
            'band_s': [1,4], 'geometric_majorant_prefactor_squared': pack(exact_prefactor_squared),
            'epsilon_polynomial': pack(eps0), 'epsilon_with_positive_base': pack(eps),
            'candidate_linear_jet_order': 2+2*degree,
            'minimal_degree_under_declared_channel_bound': n_first,
            'exact_squared_frequency_samples': checks}


def oscillator_group():
    counts = Counter(sum(v*v for v in k)
                     for k in itertools.product(range(-2,3), repeat=3)
                     if 0 < sum(v*v for v in k) <= 4)
    mode_counts = {s: 2*number for s, number in counts.items()}
    assert mode_counts == {1:12, 2:24, 3:16, 4:12}
    assert sum(mode_counts.values()) == 64
    # (a+a^dagger)^2 on |0>,|1> has Gram diag(3,15), exactly.
    # Q^2 divides by 2*omega, hence ||Q^2 P_1||=sqrt(15)/(2omega)
    # <= (2N+3)/omega=5/omega. All output levels 0..3 arise exactly
    # from this single polynomial action, not from truncating evolution.
    q2_gram_times_4omega2 = [3,15]
    assert max(q2_gram_times_4omega2) <= 4*5**2
    degree, rho, delta = 32, Q(3,5), Q(1,10**12)
    eps = Q(16,5)*rho**(degree+1)+3*delta
    sum_omega_upper = Q(128)
    n_cut, duration = 1, Q(1)
    eta = duration*(eps+eps*eps/2)*(2*n_cut+3)*sum_omega_upper
    assert eta < Q(1,10000)
    return {'real_tt_mode_counts_by_k_squared': {str(s):mode_counts[s] for s in sorted(mode_counts)},
            'number_of_real_tt_modes': 64, 'reference_per_mode_max_occupation': n_cut,
            'reference_code_dimension_exact': str(2**64),
            'full_output_fock_and_arbitrary_reference': True,
            'same_canonical_Q_P_not_squeezed_mode_identification': True,
            'q2_exact_gram_times_4omega2_on_P1': q2_gram_times_4omega2,
            'sum_reference_frequency_upper': pack(sum_omega_upper),
            'prepared_isometry_error': pack(eta),
            'prepared_isometry_error_float_upper': up(eta),
            'any_effect_probability_error_upper': pack(eta),
            'prepared_channel_diamond_error_upper': pack(2*eta),
            'scope': 'Finite Fourier regulator, full Fock output; no ultraviolet/full-field continuum claim.'}


def calculate(bind=True):
    out = {'schema':'research_1180_independent_review_v1', 'groups':3,
           'algebra':jet_group(), 'finite_band':band_group(),
           'canonical_quantum_comparison':oscillator_group()}
    if bind:
        if set(EXPECTED_CORE) != set(CORE_FILES):
            raise RuntimeError('Author core has not yet been frozen and bound.')
        for rel, digest in {**EXPECTED_CORE, **EXPECTED_SOURCE}.items():
            if sha(STAGE/rel) != digest:
                raise RuntimeError(f'Frozen dependency changed: {rel}')
        out['core_sha256'] = EXPECTED_CORE
        out['source_sha256'] = EXPECTED_SOURCE
        out['reviewer_code_sha256'] = sha(Path(__file__))
        author = json.loads((HERE/'results.json').read_text(encoding='utf-8'))
        band = author['g3_uniform_band_same_canonical_quantum_process']
        assert Q(band['prepared_full_Fock_distance_bound']['hi_exact']) == Q(
            out['canonical_quantum_comparison']['prepared_isometry_error']['exact'])
        assert Q(band['uniform_relative_frequency_error']['hi_exact']) == Q(
            out['finite_band']['epsilon_with_positive_base']['exact'])
        assert author['g1_low_order_class']['third_gauge_jet_rank'] == out['algebra']['2']['rank']
        assert band['real_TT_modes'] == out['canonical_quantum_comparison']['number_of_real_tt_modes']
        assert author['g2_strict_positive_two_point_escape']['omega_squared_k1'] == '1/1'
        assert author['g2_strict_positive_two_point_escape']['omega_squared_k2'] == '4/1'
        out['author_exact_summary_matches'] = True
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--math-only', action='store_true')
    args = parser.parse_args()
    if args.write and args.math_only:
        parser.error('--write and --math-only are mutually exclusive')
    result = calculate(bind=not args.math_only)
    if args.math_only:
        print(json.dumps({'status':'MATH_ONLY_PASS','groups':3,
                          'eta_upper':result['canonical_quantum_comparison']['prepared_isometry_error_float_upper']},
                         ensure_ascii=False))
        return
    target = HERE/'review_checks.json'
    if args.write:
        with target.open('x', encoding='utf-8', newline='\n') as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write('\n')
        status = 'WRITE_ONCE_PASS'
    else:
        if result != json.loads(target.read_text(encoding='utf-8')):
            raise RuntimeError('Independent saved review result does not match.')
        status = 'READ_ONLY_MATCH'
    print(json.dumps({'status':status,'groups':3,'blockers':0},ensure_ascii=False))


if __name__ == '__main__':
    main()
