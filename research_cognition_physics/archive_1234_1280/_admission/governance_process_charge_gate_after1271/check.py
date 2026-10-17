"""Recompute a bounded candidate-qualification screen; default is read-only."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import math
import numpy as np


def exp_bounds(x, n):
    """Positive rational x; Taylor lower bound and a geometric upper tail."""
    term = F(1)
    total = term
    for k in range(1, n + 1):
        term *= x / k
        total += term
    next_term = term * x / (n + 1)
    ratio = x / (n + 2)
    assert 0 <= ratio < 1
    return total, total + next_term / (1 - ratio)


def rational_pair(x):
    return [x.numerator, x.denominator]


def run():
    lo_e, hi_e = exp_bounds(F(1), 20)
    lo_max, hi_max = exp_bounds(F(14, 9), 20)
    lo_9, _ = exp_bounds(F(9), 20)
    lo_07, _ = exp_bounds(F(7, 10), 8)
    assert lo_e > F(8, 3) and hi_e < 3 and hi_max < 5
    assert lo_9 > 3001 and lo_07 > 2
    support = F(3, 8) ** 2 * (1 - F(3, 8))
    gap = 1 - 10 * support
    assert support == F(45, 512) and gap == F(31, 256)
    assert F(1, 100) * 9 + F(1, 100) == F(1, 10)

    # Finite representation of the complete-instrument lemma, not the II1 model.
    I = np.eye(2)
    p = [np.diag([1.0, 0.0]), np.diag([0.0, 1.0])]
    C = np.kron(p[1], I)
    K = []
    for r in range(2):
        for m in range(2):
            reset = np.zeros((2, 2))
            reset[0, m] = 1
            K.append(np.kron(p[r], reset))
    effect = sum(k.T @ k for k in K)
    out_I = sum(k @ k.T for k in K)
    assert np.array_equal(effect, np.eye(4))
    assert not np.array_equal(out_I, np.eye(4))
    assert all(np.array_equal(k @ C, C @ k) for k in K)
    rho = np.diag([3/8, 1/8, 1/8, 3/8])
    out = sum(k @ rho @ k.T for k in K)
    branch = [sum(K[2*r+m] @ rho @ K[2*r+m].T for m in range(2))
              for r in range(2)]
    probs = [float(np.trace(b)) for b in branch]
    means = [float(np.trace(C @ b) / np.trace(b)) for b in branch]
    assert probs == [0.5, 0.5] and means == [0.0, 1.0]
    assert float(np.trace(C @ out)) == float(np.trace(C @ rho)) == 0.5
    for B in (C, np.eye(4)-C, 2*C+3*np.eye(4)):
        assert np.array_equal(sum(k.T @ B @ k for k in K), B)
    # A coherent full-rank input additionally checks that no classical-input
    # restriction has crept into this finite algebraic identity.
    v = np.array([1, 1j, -1, 1j]) / 2
    coherent = 0.5*np.outer(v, v.conj()) + np.eye(4)/8
    coherent_out = sum(k @ coherent @ k.T for k in K)
    assert abs(np.trace(C @ (coherent_out-coherent))) < 1e-14

    # Same-side complete energy: |q=0,a=1> <-> |q=1,a=0>.
    q = np.kron(p[1], I)
    a = np.kron(I, p[1])
    U = np.eye(4)[[0, 2, 1, 3]]
    before = np.diag([0.0, 1.0, 0.0, 0.0])
    after = U @ before @ U.T
    dq, da = [float(np.trace(H @ (after-before))) for H in (q, a)]
    assert np.array_equal(U.T @ U, np.eye(4))
    assert np.array_equal(U @ (q+a), (q+a) @ U)
    assert (dq, da) == (1.0, -1.0)

    # Diagnostics only: the certificates above do not rely on these floats.
    xstar = 1.5 + math.atan(1/(2*math.pi))/math.pi
    maximum = 2*math.exp(xstar)*math.sin(math.pi*(xstar-1))**2
    return {
        'schema': 'governance_process_charge_screen_results_v1',
        'status': 'PASS_WITH_SCOPE',
        'new_scientific_rounds': 0, 'new_formal_check_groups': 0,
        'new_adopted_axioms': 0,
        'rational_certificates': {
            'exp_1_lower': rational_pair(lo_e),
            'exp_1_upper': rational_pair(hi_e),
            'exp_14_over_9_upper': rational_pair(hi_max),
            'exp_9_lower': rational_pair(lo_9),
            'exp_07_lower': rational_pair(lo_07),
            'support_trace_strict_upper': rational_pair(support),
            'forward_half_L1_strict_lower': rational_pair(gap),
            'entropy_cost_at_delta_001_strict_lower': [9, 10],
        },
        'complete_instrument_diagnostic': {
            'TP': True, 'density_unital': False,
            'charge_distribution_preserved': True,
            'branch_probabilities': probs,
            'conditional_charge_means': means,
            'unconditional_charge_mean': 0.5,
            'coherent_input_charge_preserved': True,
        },
        'same_side_source_diagnostic': {
            'q_change': dq, 'auxiliary_change': da,
            'complete_source_change': dq+da,
            'conditional_left_area_change_beta_one': -(dq+da),
        },
        'floating_diagnostics_not_proof': {
            'd0_supremum': round(maximum, 12),
            'target_support_trace': round(math.exp(-2)-math.exp(-3), 12),
        },
        'scope': 'Exact rational certificates and finite lemma diagnostics only; '
                 'no actual gravitational instrument, no universal governance no-go, '
                 'no finite physical resource or field-equation certificate.'
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    destination = Path(__file__).with_name('results.json')
    if args.write:
        if destination.exists():
            raise SystemExit('Refusing to overwrite saved results')
        destination.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
        print('WROTE', destination.name)
    else:
        saved = json.loads(destination.read_text(encoding='utf-8'))
        assert saved == result, 'Recomputation differs from saved result'
        print('PASS: saved screen results reproduced; files unchanged')
