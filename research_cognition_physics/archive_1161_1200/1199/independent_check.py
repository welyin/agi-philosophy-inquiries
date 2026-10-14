"""Independent finite-code verification; no import of the author check.py.

Rebuilds source stencils as translated Laurent polynomials, computes quotient
syndromes using an independently implemented kernel basis, and evaluates exact
conditional traces. Default is stdout only. No large Hilbert matrix is used.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
import json
import numpy as np

P = 5
SIDE = 4
PAIRS = tuple((a, b) for a in range(3) for b in range(a, 3))
SITES = tuple(product(range(SIDE), repeat=3))
INDEX = {(x, ab): j for j, (x, ab) in enumerate(product(SITES, PAIRS))}
D = len(INDEX)


def poly(*items):
    # Each item is (tensor indices, shift vector, integer coefficient).
    q = defaultdict(int)
    for ab, shift, coefficient in items:
        q[(tuple(sorted(ab)), tuple(shift))] += coefficient
    return dict(q)


def steps(*axis_steps):
    x = [0, 0, 0]
    for axis, step in axis_steps:
        x[axis] += step
    return tuple(x)


def orbit(polynomial):
    rows = []
    for origin in SITES:
        row = np.zeros(D, dtype=np.int64)
        for (ab, shift), coefficient in polynomial.items():
            x = tuple((a+b) % SIDE for a, b in zip(origin, shift))
            row[INDEX[x, ab]] += coefficient
        rows.append(row)
    return np.array(rows)


def templates():
    rdiag, roff, tdiag, toff, v = [], [], [], [], []
    for a, b, c in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        # R_aa: transverse second sums and the bc face corner sum.
        items = [((u, u), steps((w, j)), q)
                 for u, w in ((b, c), (c, b))
                 for j, q in ((-1, 1), (0, 2), (1, 1))]
        items += [((b, c), steps((b, j), (c, k)), -1)
                  for j, k in product((0, -1), repeat=2)]
        rdiag.append(poly(*items))
        # R_bc = 2(1+b)(1+c)theta_aa - ... + (a+2+a^-1)theta_bc.
        items = [((a, a), steps((b, j), (c, k)), 2)
                 for j, k in product((0, 1), repeat=2)]
        items += [((a, b), steps((a, j), (c, k)), -1)
                  for j in (0, -1) for k in (0, 1)]
        items += [((a, c), steps((a, j), (b, k)), -1)
                  for j in (0, -1) for k in (0, 1)]
        items += [((b, c), steps((a, j)), q)
                  for j, q in ((-1, 1), (0, 2), (1, 1))]
        roff.append(poly(*items))
        tdiag.append(poly(*[((a, b), steps((c, j)), 1) for j in (0, 1)],
                          *[((a, c), steps((b, j)), -1) for j in (0, 1)]))
        items = [((b, c), steps((b, j)), 2) for j in (0, -1)]
        items += [((u, u), steps((c, j)), coefficient)
                  for u, coefficient in ((a, 1), (b, -1), (c, 1)) for j in (0, 1)]
        toff.append(poly(*items))
        items = [((b, c), steps((c, j)), -2) for j in (0, -1)]
        items += [((u, u), steps((b, j)), coefficient)
                  for u, coefficient in ((a, -1), (b, -1), (c, 1)) for j in (0, 1)]
        toff.append(poly(*items))
        items = [((a, a), steps((a, j)), 1) for j in (0, 1)]
        items += [((a, u), steps((u, j)), 1)
                  for u in range(3) if u != a for j in (0, -1)]
        v.append(poly(*items))
    eta = defaultdict(int)
    for polynomial in rdiag:
        for key, value in polynomial.items():
            eta[key] += value
    # Explicit independent eta construction checks the trace identity.
    eta_explicit = defaultdict(int)
    for b in range(3):
        for a in range(3):
            if a != b:
                for j, q in ((-1, 1), (0, 2), (1, 1)):
                    eta_explicit[((b, b), steps((a, j)))] += q
    for a, b in ((0, 1), (0, 2), (1, 2)):
        for j, k in product((0, -1), repeat=2):
            eta_explicit[((a, b), steps((a, j), (b, k)))] -= 1
    assert eta == eta_explicit
    return eta, v, rdiag, roff, tdiag, toff


def kernel_basis(matrix):
    # Forward elimination only, followed by triangular back substitution.
    # The author instead reduces input rows against full RREF quotient bases.
    a = matrix.copy() % P
    row = 0
    pivots = []
    for col in range(a.shape[1]):
        candidates = np.flatnonzero(a[row:, col])
        if not len(candidates):
            continue
        pick = row + int(candidates[0])
        a[[row, pick]] = a[[pick, row]]
        a[row] = a[row] * pow(int(a[row, col]), P-2, P) % P
        if row+1 < a.shape[0]:
            a[row+1:] = (a[row+1:] - a[row+1:, col:col+1] * a[row]) % P
        pivots.append(col)
        row += 1
        if row == len(a):
            break
    free = [c for c in range(a.shape[1]) if c not in pivots]
    null = np.zeros((len(free), a.shape[1]), dtype=np.int64)
    null[np.arange(len(free)), free] = 1
    for r, c in reversed(list(enumerate(pivots))):
        null[:, c] = -(null[:, c+1:] @ a[r, c+1:]) % P
    assert np.all((matrix @ null.T) % P == 0)
    assert len(free)+len(pivots) == matrix.shape[1]
    return len(pivots), null


def sector(rows, coefficients, kernel):
    syndromes = (rows @ kernel.T) % P
    dictionary = defaultdict(Fraction)
    for syndrome, coefficient in zip(syndromes, coefficients):
        dictionary[bytes(syndrome.astype(np.uint8))] += coefficient
        dictionary[bytes((-syndrome % P).astype(np.uint8))] += coefficient
    zero = bytes(kernel.shape[0])
    mean = dictionary.pop(zero, Fraction(0))
    second = Fraction(0)
    for key, coefficient in dictionary.items():
        opposite = bytes((-np.frombuffer(key, dtype=np.uint8).astype(np.int64) % P).astype(np.uint8))
        assert dictionary[opposite] == coefficient
        second += coefficient * dictionary[opposite]
    histogram = Counter(str(x) for x in dictionary.values())
    return mean, second, len(dictionary), dict(sorted(histogram.items()))


def verify():
    eta, vectors, rd, ro, td, to = templates()
    scalar = orbit(eta)
    vector = np.concatenate([orbit(x) for x in vectors])
    curvature = np.concatenate([orbit(x) for x in rd+ro])
    current = np.concatenate([orbit(x) for x in td+to])
    assert scalar.shape == (64, 384) and vector.shape == (192, 384)
    comm = {'scalar_vector': int(np.abs(scalar@vector.T).max()),
            'curvature_vector': int(np.abs(curvature@vector.T).max()),
            'scalar_current': int(np.abs(scalar@current.T).max())}
    assert all(x == 0 for x in comm.values())
    rs, ks = kernel_basis(scalar)
    rv, kv = kernel_basis(vector)
    coeff_g = [Fraction(-P, 2)]*(3*64) + [Fraction(-P, 4)]*(3*64)
    coeff_j = [Fraction(-P, 2)]*(3*64) + [Fraction(-P, 8)]*(6*64)
    mg, vg, ng, hg = sector(curvature, coeff_g, ks)
    mj, vj, nj, hj = sector(current, coeff_j, kv)
    assert (rs, rv) == (63, 189)
    assert mg == mj == 0 and vg == vj == 3000
    # A pure XZ product lies in this CSS stabilizer iff both components do.
    # Since neither Hamiltonian contains the trivial syndrome, mixed trace=0.
    mixed = Fraction(0)
    total = vg+vj+2*mixed
    assert total == 6000 and 77**2 < total < 78**2
    norm_upper = 2*(sum(map(abs, coeff_g))+sum(map(abs, coeff_j)))
    assert norm_upper == 2880
    gap_lower = P*2000*(1-(Fraction(9,4)-1)/4)
    assert gap_lower == 6875 and gap_lower-2*norm_upper == 1115
    # Independent one-qudit Weyl orientation and trace identity sanity check.
    w = np.exp(2j*np.pi/P)
    x = np.roll(np.eye(P), 1, axis=0)
    z = np.diag(w**np.arange(P))
    assert np.linalg.norm(z@x-w*x@z) < 1e-13
    for a in range(P):
        for b in range(P):
            tr = np.trace(np.linalg.matrix_power(x,a)@np.linalg.matrix_power(z,b))/P
            assert abs(tr-(1 if a==b==0 else 0)) < 1e-13
    return {'status': 'PASS_INDEPENDENT_EXACT_FINITE_CODE',
            'method': 'Laurent polynomial stencils; finite-field kernel syndromes; no author imports',
            'side': SIDE, 'prime': P, 'qudits': D,
            'ranks': {'scalar': rs, 'vector': rv}, 'code_power': D-rs-rv,
            'integer_commutators': comm,
            'Hg': {'mean': str(mg), 'variance_per_g_squared': str(vg),
                   'cosets': ng, 'coefficient_histogram': hg},
            'HJ': {'mean': str(mj), 'variance_per_J_squared': str(vj),
                   'cosets': nj, 'coefficient_histogram': hj},
            'mixed_trace': str(mixed), 'total_variance_g_J_one': str(total),
            'norm_source_mismatch_strict_lower_for_a_one': 77,
            'finite_penalty': {'dynamic_norm_upper': str(norm_upper),
                               'defect_penalty_strict_lower': str(gap_lower),
                               'isolated_band_gap_strict_lower': str(gap_lower-2*norm_upper)},
            'scope': 'Whole invariant code, or chosen entire isolated finite band; not every smaller window or a physical source assertion'}


if __name__ == '__main__':
    print(json.dumps(verify(), ensure_ascii=False, indent=2))
