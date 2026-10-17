"""Exact fractions certify the finite screening; NumPy audits the paid unitary."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent

def compute():
    p = [F(v, 10) for v in (4, 3, 1, 2)]
    q = [F(v, 10) for v in (3, 4, 2, 1)]
    g = [F(v, 15) for v in (8, 4, 2, 1)]
    mt = lambda r, t: sum(max(F(0), a-t*b) for a, b in zip(r, g))
    assert sorted(p) == sorted(q)
    energy = lambda r: sum(i*a for i, a in enumerate(r))
    assert energy(p) == energy(q) == F(11, 10)
    tests = []
    for source, target, t in [(p, q, F(3, 4)), (q, p, F(3, 2))]:
        gap = mt(target, t)-mt(source, t)
        assert gap == F(1, 10)
        bound = gap-F(2, 1000)-t*F(1, 1000)
        assert bound > F(9, 100)
        tests.append({'direction': 'p_to_q' if source == p else 'q_to_p',
                      't': str(t), 'source_M': str(mt(source, t)),
                      'target_M': str(mt(target, t)), 'exact_gap': str(gap),
                      'conditional_error_lower': str(bound)})
    # Each transition is checked before any floating point work.
    pairs = [((0, 1), (1, 0)), ((1, 1), (0, 2)),
             ((2, 1), (3, 0)), ((3, 1), (2, 2))]
    perm = list(range(12))
    for a, b in pairs:
        assert sum(a) == sum(b)
        ia, ib = 3*a[0]+a[1], 3*b[0]+b[1]
        perm[ia], perm[ib] = ib, ia
    assert sorted(perm) == list(range(12))
    U = np.zeros((12, 12))
    U[perm, np.arange(12)] = 1
    H = np.diag([i+j for i in range(4) for j in range(3)])
    assert np.array_equal(U.T @ U, np.eye(12))
    assert np.array_equal(U @ H, H @ U)
    branches = []
    for source, target in [(p, q), (q, p)]:
        initial = [F(0)]*12
        for i, prob in enumerate(source): initial[3*i+1] = prob
        final = [F(0)]*12
        for i, j in enumerate(perm): final[j] = initial[i]
        reduced = [sum(final[3*i:3*i+3]) for i in range(4)]
        aux = [sum(final[3*i+j] for i in range(4)) for j in range(3)]
        assert reduced == target and aux == [F(1, 2), F(0), F(1, 2)]
        assert sum(j*x for j, x in enumerate(aux)) == 1
        assert sum((j-1)**2*x for j, x in enumerate(aux)) == 1
        assert aux[1] == 0
        # Exact distributions imply entropy/correlation identities analytically;
        # these numeric checks only audit that bookkeeping independently of labels.
        ent = lambda r: -sum(float(x)*math.log(float(x)) for x in r if x)
        mutual = ent(reduced)+ent(aux)-ent(final)
        assert abs(ent(final)-ent(initial)) < 2e-15
        assert abs(mutual-math.log(2)) < 2e-15
        branches.append({'source': [str(x) for x in source],
            'joint_final': [str(x) for x in final], 'material_final': [str(x) for x in reduced],
            'auxiliary_final': [str(x) for x in aux], 'auxiliary_mean': '1',
            'auxiliary_variance': '1', 'auxiliary_return_probability': '0',
            'auxiliary_entropy_increase': 'log(2)', 'final_mutual_information': 'log(2)'})
    return {'status': 'PASS_BOUNDED_ROLE_SCREEN_ONLY', 'mean_energy': str(energy(p)),
            'same_entropy_exact_reason': 'identical eigenvalue multisets',
            'same_relative_entropy_exact_reason': 'same H, mean energy, entropy and Gibbs reference',
            'monotone_tests': tests, 'energy_preserving_permutation': perm,
            'assisted_conversions': branches, 'geometry_claim': False,
            'actual_device_permission_claim': False, 'new_axiom': False}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--write', action='store_true')
    args = ap.parse_args()
    result = compute()
    target = HERE / 'results.json'
    if args.write:
        assert not target.exists(), 'Refuse to overwrite a saved scientific result'
        target.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    else:
        assert json.loads(target.read_text(encoding='utf-8')) == result
    print(json.dumps({'status': result['status'], 'mode': 'WRITE_ONCE' if args.write else 'READ_ONLY_MATCH',
                      'exact_conversion_gap': '1/10', 'assisted_return_probability': '0'}))

if __name__ == '__main__':
    main()
