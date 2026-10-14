"""1197 exact source/response certificate. Default read-only.

Reuses the frozen 1196 rational interval implementation. The Pauli calculation
uses exact Gaussian rational coefficients; it tests a representative full
finite-copy operator identity, not a simulation of the spatial parent.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement
import argparse, hashlib, importlib.util, json

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('cert1196', HERE.parent/'1196/certificate.py')
iv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(iv)
pt, add, sub, mul, sc, div, sq, sqrt = (getattr(iv, k) for k in
    ('pt', 'add', 'sub', 'mul', 'sc', 'div', 'sq', 'sqrt'))


def gaussmul(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def opadd(*ops):
    out = {}
    for op in ops:
        for word, c in op.items():
            old = out.get(word, (F(0), F(0)))
            new = old[0]+c[0], old[1]+c[1]
            if new != (0, 0):
                out[word] = new
            else:
                out.pop(word, None)
    return out


def scale(op, a, b=0):
    z = (F(a), F(b))
    return {w: gaussmul(c, z) for w, c in op.items() if gaussmul(c, z) != (0, 0)}


def product(a, b):
    out = {}
    cyclic = {(1, 2): 3, (2, 3): 1, (3, 1): 2}
    for wa, ca in a.items():
        for wb, cb in b.items():
            word, phase = [], (F(1), F(0))
            for x, y in zip(wa, wb):
                if x == 0:
                    word.append(y)
                elif y == 0:
                    word.append(x)
                elif x == y:
                    word.append(0)
                elif (x, y) in cyclic:
                    word.append(cyclic[x, y]); phase = gaussmul(phase, (0, 1))
                else:
                    word.append(cyclic[y, x]); phase = gaussmul(phase, (0, -1))
            out = opadd(out, {tuple(word): gaussmul(gaussmul(ca, cb), phase)})
    return out


def comm(a, b):
    return opadd(product(a, b), scale(product(b, a), -1))


def pauli_identity(K):
    # Per copy: source S, relation s, relation t; plus one probe factor.
    size = 3*K+1
    def one(index=None, axis=0):
        w = [0]*size
        if index is not None:
            w[index] = axis
        return {tuple(w): (F(1), F(0))}
    I = one(); om, eps, nu = F(1, 2), F(3, 20), F(7, 20)
    H = {}
    fs, ft, bt = [], [], []
    for j in range(K):
        S, s, t = 3*j, 3*j+1, 3*j+2
        H = opadd(H, scale(one(S, 3), 2), one(S, 1),
                  scale(one(s, 1), -om), scale(one(t, 1), -om))
        fs.append(opadd(scale(one(S, 3), 2*eps), scale(one(S, 1), nu)))
        ft.append(scale(one(S, 2), F(1, 5)))
        bt.append(scale(one(t, 1), -eps*om))
    M = opadd(scale(I, F(2, 5)), scale(one(size-1, 3), -F(1, 10)))
    fp = scale(M, eps)
    H = opadd(H, M)
    for j in range(K):
        H = opadd(H, scale(product(one(3*j+1, 3), fp), F(1, K)))
        for k in range(K):
            if j != k:
                H = opadd(H, scale(product(one(3*j+1, 3), opadd(fs[k], bt[k])), F(1, K)),
                          scale(product(one(3*j+2, 3), ft[k]), F(1, K)))
    Z, X = one(1, 3), one(1, 1)
    A = scale(fp, F(1, K))
    for k in range(1, K):
        A = opadd(A, scale(opadd(fs[k], bt[k]), F(1, K)))
    lhs = scale(comm(H, comm(H, Z)), -1)
    rhs = opadd(scale(Z, -4*om**2), scale(product(X, A), -4*om))
    difference = opadd(lhs, scale(rhs, -1))
    assert not difference
    return {'copies': K, 'Hamiltonian_terms': len(H),
            'exact_double_commutator_difference_terms': len(difference)}


def calculate():
    old = iv.run()  # Recheck the frozen enclosure, do not trust a pasted decimal.
    z = tuple(F(old['fixed_point']['root_enclosure'][s]) for s in ('lo', 'hi'))
    ell, eps, nu, om, m = F(16, 65), F(3, 20), F(7, 20), F(1, 2), F(2)
    N, v = add(pt(1), sc(z, eps)), add(pt(1), sc(z, nu))
    k2 = [sc(sq(iv.trig(sc(iv.PI, F(2*j, 65)))), 1/ell**2) for j in range(33)]
    invsum, hopsum, weights = pt(0), pt(0), 0
    for i, j, k in combinations_with_replacement(range(33), 3):
        perm = 1 if i == k else 3 if i == j or j == k else 6
        w = perm*(1 if i == 0 else 2)*(1 if j == 0 else 2)*(1 if k == 0 else 2)
        weights += w
        q = add(add(k2[i], k2[j]), k2[k])
        e = sqrt(add(sc(sq(N), m*m), mul(sq(v), q)))
        invsum = add(invsum, sc(div(pt(1), e), w))
        hopsum = add(hopsum, sc(div(q, e), w))
    assert weights == 65**3
    em = sc(mul(N, invsum), -2*m*m/F(weights))
    et = sc(mul(v, hopsum), -F(2, weights))
    b = sub(add(sc(em, eps), sc(et, nu)), pt(eps*om))
    r = sqrt(add(pt(om*om), sq(b)))
    delta_hop = mul(sub(pt(1), sc(sqrt(pt(2)), F(1, 2))), iv.neg(et))
    delta_force = sc(delta_hop, nu)
    magnitude = div(sc(delta_force, 4*om*om), r)
    assert delta_hop[0] > F(267, 100)
    assert delta_force[0] > F(936, 1000)
    assert magnitude[0] > F(266, 1000)
    # q_K in [1/2,1) for K >= 2. This is an exact non-closure statement,
    # not a finite-time or measurement-resolution certificate.
    energy_K2 = mul(add(pt(1), sc(z, nu/F(2))), delta_hop)
    assert energy_K2[0] > F(313, 100)
    checks = [pauli_identity(K) for K in (1, 2, 3)]
    return {
      'schema': 'same_macro_state_source_response_1197_v1', 'status': 'PASS',
      'check_groups': 2,
      'source_interval': {'weights': weights, 'Ehop': iv.pair(et), 'b': iv.pair(b),
          'r': iv.pair(r), 'delta_Ehop': iv.pair(delta_hop),
          'delta_force': iv.pair(delta_force),
          'negative_acceleration_coefficient_magnitude': iv.pair(magnitude),
          'full_preparation_energy_per_site_copy_lower_at_K2': iv.pair(energy_K2),
          'all_K_at_least_2_acceleration_difference_abs_lower': '133/1000'},
      'exact_Pauli_checks': checks,
      'scope': 'Two equal coarse initial states have different exact finite-K relation curvature; only rejects this memoryless source summary. No certified finite-time instrument gap or GR exclusion.',
      'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'inherited_interval_code_sha256': hashlib.sha256((HERE.parent/'1196/certificate.py').read_bytes()).hexdigest()}


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--write', action='store_true'); a = p.parse_args()
    result = calculate()
    if a.write:
        with (HERE/'results.json').open('x', encoding='utf-8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2); f.write('\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
