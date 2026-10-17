"""Exact finite classical dictionary checks. Default replay never writes files."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json

D = Path(__file__).resolve().parent
ROOT = D.parents[3]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def radial(L, alpha, s):
    assert L > 0 and L + s > 0
    return L + s, alpha * L / (L + s)

def det(a, b, c, d):
    return a * d - b * c

def compute():
    baseline = json.loads((D / 'boundary_charge_sources.json').read_text(encoding='utf-8-sig'))
    for path, digest in baseline['source_sha256'].items():
        assert sha(ROOT / path) == digest, path
    goal = ROOT / 'research_cognition_physics/archive_1234_/_shared/goal_prompt.txt'
    assert sha(goal) == baseline['goal_sha256']
    rows = []
    for G in [F(1, 2), F(1), F(2)]:
        for L in [F(1), F(2), F(3)]:
            for s in [F(-1, 4), F(1, 4), F(1, 2)]:
                for alpha in [F(-2, 3), F(3, 7)]:
                    L1, a1 = radial(L, alpha, s)
                    assert radial(L1, a1, -s) == (L, alpha)
                    t = F(1, 7)
                    assert radial(*radial(L, alpha, s), t) == radial(L, alpha, s + t)
                    # Jacobian in (L, alpha), verified with the weighted 2-form.
                    J = det(F(1), F(0), alpha * s / L1**2, L / L1)
                    correct_ratio = L1 * J / L
                    assert correct_ratio == 1
                    bad_ratio = L1 / L
                    assert bad_ratio != 1
                    # q'=q+b is a different fixed family: match this one basepoint.
                    q = L**2 / (2 * G)
                    b = (L1**2 - L**2) / (2 * G)
                    assert 2 * G * (q + b) == L1**2
                    # Hamilton equations for H=alpha*L/G, with dq/dL=L/G.
                    q_dot = L / G
                    alpha_dot = -alpha / L
                    assert G * q_dot / L == 1
                    assert alpha_dot * L + alpha == 0
                    rows.append({'G': str(G), 'L': str(L), 's': str(s), 'alpha': str(alpha),
                                 'bad_pullback_ratio': str(bad_ratio),
                                 'correct_pullback_ratio': str(correct_ratio)})
    assert len(rows) == 54
    G, L, alpha, s = F(1), F(2), F(3, 7), F(1, 2)
    L1, a1 = radial(L, alpha, s)
    b = (L1**2 - L**2) / (2 * G)
    assert (L1, a1, b, alpha - a1) == (F(5, 2), F(12, 35), F(9, 8), F(3, 35))
    # Integrate the image of a finite rectangle, not just a single Jacobian.
    lo, hi, width = F(2), F(3), F(1)
    initial_volume = width * (hi**2 - lo**2) / (2 * G)
    bad_volume = width * ((hi + s)**2 - (lo + s)**2) / (2 * G)
    # At L'=x, transported alpha-width is width*(x-s)/x.
    image_lo, image_hi = lo + s, hi + s
    good_volume = width / G * (
        (image_hi**2 - image_lo**2) / 2 - s * (image_hi - image_lo))
    assert (initial_volume, bad_volume, good_volume) == (F(5, 2), F(3), F(5, 2))
    return {'status': 'PASS_CONDITIONAL_FINITE_BOUNDARY_DICTIONARY', 'cases': len(rows),
            'model': 'Classical reduced EH causal-diamond phase space; no native implementation',
            'witness': {'L_final': str(L1), 'fixed_area_parameter': str(b),
                        'alpha_fixed_area': str(alpha), 'alpha_fixed_radius': str(a1),
                        'alpha_difference': str(alpha - a1),
                        'wrong_local_symplectic_ratio': str(L1 / L),
                        'initial_rectangle_volume': str(initial_volume),
                        'wrong_rectangle_volume': str(bad_volume),
                        'transported_rectangle_volume': str(good_volume)},
            'rows': rows, 'sources_verified': len(baseline['source_sha256']),
            'new_physical_rounds': 0, 'new_axioms': 0}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    result = compute()
    path = D / 'boundary_charge_results.json'
    if args.write_results:
        assert not path.exists(), 'Preserve the saved result; use read-only replay.'
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        mode = 'initial_write'
    else:
        assert json.loads(path.read_text(encoding='utf-8-sig')) == result
        mode = 'read_only_match'
    print(json.dumps({'status': result['status'], 'mode': mode, 'cases': result['cases'],
                      'witness': result['witness']}, ensure_ascii=False))

if __name__ == '__main__':
    main()
