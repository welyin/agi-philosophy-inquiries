"""Finite-query matrix diagnostics; no physical menu permissions are inferred.

Default: read-only recomputation. --write-once creates results.json exclusively.
"""
from pathlib import Path
import argparse
import json
import math
from fractions import Fraction
import numpy as np

HERE = Path(__file__).resolve().parent
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
PAULI = [X, Y, Z]


def tensor(*matrices):
    out = np.ones((1, 1), dtype=complex)
    for matrix in matrices:
        out = np.kron(out, matrix)
    return out


def reduced(rho, dims, keep):
    """Trace complement, preserving the explicit order in keep."""
    rest = [k for k in range(len(dims)) if k not in keep]
    order = list(keep) + rest
    q = len(dims)
    dk = math.prod(dims[k] for k in keep)
    dr = math.prod(dims[k] for k in rest)
    arr = rho.reshape(tuple(dims) * 2).transpose(order + [k + q for k in order])
    return np.einsum('arbr->ab', arr.reshape(dk, dr, dk, dr))


def remainders(u, g, m):
    vec = u.reshape(-1) / math.sqrt(g * m)
    rho = np.outer(vec, vec.conj())
    dims = [g, m, g, m]  # Gout, Mout, Gin, Min
    xx = reduced(rho, dims, [0, 2])
    yy = reduced(rho, dims, [1, 3])
    rb = reduced(rho, dims, [0, 2, 3]) - tensor(xx, np.eye(m) / m)
    rf = reduced(rho, dims, [1, 3, 2]) - tensor(yy, np.eye(g) / g)
    return rb, rf, 1 - float(np.trace(xx @ xx).real)


def probability(u, gp, mp, effect, output):
    rho = u @ tensor(gp, mp) @ u.conj().T
    ef = tensor(effect, I) if output == 'G' else tensor(I, effect)
    return float(np.trace(ef @ rho).real)


def calculate():
    rng = np.random.default_rng(1271001)
    identities = []
    for g, m in [(2, 2), (2, 3), (3, 2)]:
        for sample in range(2):
            d = g * m
            a = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
            u, _ = np.linalg.qr(a)
            rb, rf, p = remainders(u, g, m)
            values = [m * np.linalg.norm(rb)**2, g * np.linalg.norm(rf)**2, p]
            error = max(values) - min(values)
            assert error < 2e-12
            identities.append({'dimensions': [g, m], 'sample': sample,
                               'three_values': values, 'max_difference': error})

    u = np.array([[1, 0, 0, 0], [0, 1, 0, 0],
                  [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex)
    effects = [(I + p) / 2 for p in PAULI]
    preparations = [I / 2] + effects
    pairs = [((I + p) / 2, (I - p) / 2) for p in PAULI]
    traceless = [p / math.sqrt(2) for p in PAULI]
    full = [I / math.sqrt(2)] + traceless
    basis = [tensor(a, b, c) for a in traceless for b in full for c in traceless]
    rb, rf, impurity = remainders(u, 2, 2)
    tests, gaps = [], []
    for effect in effects:
        for prep in preparations:
            for plus, minus in pairs:
                t = tensor(effect, prep.T, (plus - minus).T)
                tests.append([float(np.trace(t @ b).real) for b in basis])
                gap = probability(u, prep, plus, effect, 'G') - probability(u, prep, minus, effect, 'G')
                assert abs(gap - 4 * np.trace(rb @ t).real) < 1e-12
                gaps.append(gap)
    frame = np.asarray(tests)
    singular = np.linalg.svd(frame, compute_uv=False)
    s = float(singular[-1])
    expected_s2 = (5 - math.sqrt(21)) / 4
    assert abs(s*s - expected_s2) < 1e-12
    assert frame.shape == (36, 36) and np.linalg.matrix_rank(frame) == 36
    zero, one = (I + Z) / 2, (I - Z) / 2
    delta = abs(probability(u, zero, zero, zero, 'M') - probability(u, one, zero, zero, 'M'))
    ft = tensor(zero, zero.T, (zero - one).T)  # Mout, Min, Gin
    assert abs(delta - abs(4 * np.trace(rf @ ft).real)) < 1e-12
    af = float(np.linalg.norm(zero) * np.linalg.norm(zero - one) * np.linalg.norm(zero))
    constant = s / (math.sqrt(36) * af)
    back = max(abs(x) for x in gaps)
    assert abs(delta-1) < 1e-12 and abs(back-1) < 1e-12
    assert back >= constant*delta and constant > 1/30
    assert abs(2*np.linalg.norm(rb)**2 - 0.5) < 1e-12
    assert abs(impurity-0.5) < 1e-12

    # Exact sufficient rational certificate, not a floating-point proof.
    assert Fraction(23, 5)**2 > 21
    s2_lower = (5-Fraction(23, 5))/4
    c2_lower = s2_lower/72
    assert c2_lower > Fraction(1, 30)**2
    error = Fraction(1, 1000)
    lower = (1-2*error)/30-2*error
    assert lower == Fraction(469, 15000) > 0

    # Necessity examples: full mathematical unitary but deficient actual menu.
    z_read_gaps = []
    diagonal_prep_gaps = []
    for prep in preparations:
        for plus, minus in pairs:
            z_read_gaps.append(probability(u, prep, plus, zero, 'G')-probability(u, prep, minus, zero, 'G'))
    for prep in [zero, one]:
        for effect in effects:
            for plus, minus in pairs:
                diagonal_prep_gaps.append(probability(u, prep, plus, effect, 'G')-probability(u, prep, minus, effect, 'G'))
    assert max(map(abs,z_read_gaps)) < 1e-12
    assert max(map(abs,diagonal_prep_gaps)) < 1e-12
    assert np.linalg.matrix_rank(frame[24:36]) == 12
    # Reduced one-way channel: measure G in Z and prepare the copied Z value.
    # G output = diag(G input), independent of M input, while M learns G.
    kraus = []
    for a in range(2):
        for b in range(2):
            k = np.zeros((4,4),dtype=complex)
            k[2*a+a,2*a+b] = 1
            kraus.append(k)
    assert np.linalg.norm(sum(k.conj().T@k for k in kraus)-np.eye(4)) < 1e-12
    def channel(gp, mp):
        rho = tensor(gp,mp)
        return sum(k@rho@k.conj().T for k in kraus)
    forward = abs(np.trace(tensor(I,zero)@(channel(zero,zero)-channel(one,zero))).real)
    backward = max(abs(np.trace(tensor(e,I)@(channel(p,a)-channel(p,b))).real)
                   for e in effects for p in preparations for a,b in pairs)
    assert forward == 1 and backward < 1e-12
    return {
        'status':'PASS_CONDITIONAL_MATRIX_LEMMA_ONLY',
        'seed':1271001, 'closed_unitary_identity_diagnostics':identities,
        'frame':{'rows':36,'columns':36,'rank':36,'s_min':s,'s_min_squared_exact':'(5-sqrt(21))/4'},
        'CNOT':{'forward_gap':delta,'forward_test_norm':af,'backward_gap':back,
                'bound_constant':constant,'choi_impurity':impurity},
        'rational_certificate':{'sqrt21_upper':'23/5','s_min_squared_lower':str(s2_lower),
            'constant_squared_lower':str(c2_lower),'conservative_constant':'1/30',
            'per_branch_error':'1/1000','observed_forward_gap':'1',
            'robust_bound':str(lower)},
        'deletion_checks':{'Z_only_backward_gap':max(map(abs,z_read_gaps)),
            'diagonal_preparations_backward_gap':max(map(abs,diagonal_prep_gaps)),
            'Z_only_frame_rank':12,'one_way_reduced_channel_forward':float(forward),
            'one_way_reduced_channel_backward':float(backward)},
        'physical_geometry_frame_certified':False,'new_adopted_axioms':0,
        'scope':'Finite algebra and error certificate. Random tests are diagnostics, not a proof. No apparatus, geometric identity, curvature, Einstein selection, or exclusion of 1271 is claimed.'}


def same(a,b):
    if isinstance(a,dict):
        return isinstance(b,dict) and a.keys()==b.keys() and all(same(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return isinstance(b,list) and len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    if isinstance(a,float):
        return isinstance(b,(int,float)) and math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-12)
    return a==b


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-once',action='store_true')
    args=parser.parse_args()
    result=calculate()
    path=HERE/'results.json'
    if args.write_once:
        with path.open('x',encoding='utf-8') as output:
            json.dump(result,output,ensure_ascii=False,indent=2)
            output.write('\n')
        print('WRITE_ONCE_PASS_CONDITIONAL_MATRIX_LEMMA_ONLY')
    else:
        saved=json.loads(path.read_text(encoding='utf-8'))
        assert same(result,saved), 'saved result differs'
        print('READ_ONLY_MATCH_CONDITIONAL_MATRIX_LEMMA_ONLY')
