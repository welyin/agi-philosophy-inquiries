"""Round 1065: finite calibrations, not a numerical proof of NSS or dimension."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

BASE = Path(__file__).resolve().parent

def compare(a, b, path='root'):
    if isinstance(a, dict):
        assert set(a) == set(b), path
        for key in a:
            compare(a[key], b[key], path + '.' + key)
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, f'{path}[{i}]')
    elif isinstance(a, float):
        assert math.isclose(a, b, abs_tol=1e-11, rel_tol=1e-9), (path, a, b)
    else:
        assert a == b, (path, a, b)

def run():
    ln2 = math.log(2)
    slow = []
    for L in (4., 8., 16., 64., 256., 1024., 4096.):
        c, half = 1/L, 1/(L+ln2)
        assert 0 < half < c
        slow.append({'L': L, 'cost': c, 'half_cost': half,
                     'ratio': half/c})
    ratio_witness = []
    for q in (.5, .9, .99, .999):
        L = max(4., 2*q*ln2/(1-q))
        ratio = L/(L+ln2)
        assert ratio > q
        ratio_witness.append({'proposed_q': q, 'L': L, 'larger_ratio': ratio})

    window = []
    for L0 in (4., 6., 8., 12., 20., 32.):
        a, M = 1/L0, 1/3
        delta = ln2/(L0*(L0-ln2))
        bound = math.floor((M-a)/delta)+1
        escape = math.floor((L0-3)/ln2)+1
        assert escape <= bound
        increments = []
        for k in range(escape-1):
            L = L0-k*ln2
            assert L-ln2 >= 3
            gain = 1/(L-ln2)-1/L
            assert gain >= delta-1e-14
            increments.append(gain)
        error_budget = delta/4
        robust_bound = math.floor((M-a)/(delta-error_budget))+1
        assert robust_bound >= bound
        window.append({'L0': L0, 'cost_threshold': a, 'delta': delta,
                       'exact_escape_doublings': escape, 'guaranteed_bound': bound,
                       'increment_loss_bound': error_budget,
                       'robust_bound': robust_bound,
                       'minimum_observed_in_window_increment': min(increments) if increments else None})

    # Matrix checks use a declared SU(2) endpoint model; no physical space is inferred.
    I = np.eye(2, dtype=complex)
    sigma = np.array([[[0,1],[1,0]], [[0,-1j],[1j,0]], [[1,0],[0,-1]]], complex)
    def element(x):
        t = np.linalg.norm(x)
        return I.copy() if t == 0 else math.cos(t)*I+1j*math.sin(t)*np.einsum('j,jab->ab', x/t, sigma)
    def principal_root(h):
        return (I+h)/math.sqrt(2+float(np.trace(h).real))
    rng = np.random.default_rng(1065)
    residuals = []
    for t in (.01, .03, .08, .13, .17):
        axis = rng.normal(size=3); axis /= np.linalg.norm(axis)
        h = element(t*axis); r = principal_root(h)
        A = element(np.array([.12,-.07,.09])); B = A @ h
        residuals.extend([np.linalg.norm(r@r-h), np.linalg.norm(A@r-B@r.conj().T),
                          np.linalg.norm(principal_root(h@h)-h),
                          np.linalg.norm(principal_root(h.conj().T)-r.conj().T)])
        assert 1/math.log(math.e/(t/2)) < 1/math.log(math.e/t)
    assert max(residuals) < 1e-12

    # Finite quotients check algebra only, not compactness/connectedness of a solenoid.
    finite = []
    for k in range(1, 7):
        modulus = 3**k
        inv2 = pow(2, -1, modulus)
        roots = [(z*inv2)%modulus for z in range(modulus)]
        assert sorted(roots) == list(range(modulus))
        def valuation(z):
            if z == 0: return k
            v = 0
            while z%3 == 0: z//=3; v+=1
            return v
        for z, root in enumerate(roots):
            assert (2*root)%modulus == z
            assert roots[(2*z)%modulus] == z
            assert valuation(root) == valuation(z)
        finite.append({'modulus': modulus, 'all_elements_checked': modulus,
                       'root_permutation': True, 'valuation_preserved': True})

    return {'round': 1065, 'scope': 'Finite checks of a proved strict-cost NSS bridge; no empirical or cognitive-axiom validation.',
            'slow_cost': slow, 'no_uniform_q_witnesses': ratio_witness,
            'finite_windows': window, 'su2_handshake_max_residual': float(max(residuals)),
            'odd_solenoid_finite_quotients': finite,
            'infinite_dimensional_and_topological_claims_proved_analytically': True}

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    path = BASE/'results.json'
    if args.write:
        with path.open('x', encoding='utf8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2); f.write('\n')
    else:
        compare(result, json.loads(path.read_text(encoding='utf8')))
    print(json.dumps({'round': 1065, 'passed': True, 'matrix_residual': result['su2_handshake_max_residual'],
                      'finite_windows': len(result['finite_windows']), 'empirical_groups': 0}, ensure_ascii=False))
