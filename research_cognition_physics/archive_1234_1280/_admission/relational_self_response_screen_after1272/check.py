"""Finite Z2 relational encoding diagnostic. Read-only unless --write is explicit."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import math
import sys
import numpy as np

D = Path(__file__).resolve().parent
ROOT = D.parents[3]
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
P0, P1 = (I2 + Z) / 2, (I2 - Z) / 2


def tensor(*xs):
    out = np.array([[1.0]], dtype=complex)
    for x in xs:
        out = np.kron(out, x)
    return out


def at(op, site):
    return tensor(*(op if i == site else I2 for i in range(3)))


def norm(x):
    return float(np.linalg.norm(x, 2))


def unitary(h, t):
    vals, vecs = np.linalg.eigh(h)
    return (vecs * np.exp(-1j * t * vals)) @ vecs.conj().T


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def check():
    src = json.loads((D / 'sources.json').read_text(encoding='utf-8'))
    for p, digest in src['local_sha256'].items():
        assert sha(ROOT / p) == digest, p
    assert sha(ROOT / 'research_cognition_physics/archive_1234_/_shared/goal_prompt.txt') == src['goal_sha256']
    coeff_sum = Q(1, 4) + Q(1, 5) + Q(1, 7) + Q(1, 4)
    assert coeff_sum == Q(59, 70)
    spectral_bound = 2 - coeff_sum
    sine_lower = Q(1, 4) - Q(1, 4)**3 / 6
    probability_lower = sine_lower**2
    assert probability_lower == Q(9025, 147456) > Q(3, 50)

    ids, idp = np.eye(8), np.eye(64)
    cs, flip = tensor(Z, Z, Z), tensor(X, X, X)
    gamma_op = tensor(flip, cs)
    e = np.zeros((64, 8), dtype=complex)
    e[:8, :] = ids / math.sqrt(2)
    e[56:64, :] = cs / math.sqrt(2)
    align_r = 0.25 * (2 * ids - at(Z, 0) @ at(Z, 1) - at(Z, 1) @ at(Z, 2))
    align = tensor(align_r, ids)
    terms0 = [(Q(1, 4), tensor(Z, Z, I2), (0, 1), 0),
              (Q(1, 5), tensor(I2, Z, Z), (1, 2), 1),
              (Q(1, 7), tensor(X, X, I2), (0, 1), 0)]
    odd = tensor(X, X, X)
    h0 = sum(float(c) * op for c, op, _, _ in terms0)
    f3 = (ids - at(Z, 2)) / 2
    fpar = (ids - cs) / 2
    uin = np.zeros(8, dtype=complex)
    uin[0] = 1
    psi0 = e @ uin
    residuals = {}

    def zero(name, value):
        residuals[name] = norm(value)
        assert residuals[name] < 2e-12, (name, residuals[name])

    zero('isometry', e.conj().T @ e - ids)
    zero('total_even_code', gamma_op @ e - e)
    zero('alignment_on_code', align @ e)
    zero('h0_odd_commutator', h0 @ odd - odd @ h0)
    zero('h0_z3_commutator', h0 @ at(Z, 2) - at(Z, 2) @ h0)
    zero('odd_system_parity', cs @ odd @ cs + odd)
    cases = []
    rng = np.random.default_rng(1272)
    passive = rng.normal(size=16) + 1j * rng.normal(size=16)
    passive /= np.linalg.norm(passive)
    for gam in (Q(0), Q(1, 4)):
        tag = str(gam)
        terms = terms0 + [(gam, odd, (0, 1, 2), 1)]
        h = np.zeros((8, 8), dtype=complex)
        hp = align.copy() + 2 * idp
        for j, (c, op, support, anchor) in enumerate(terms):
            assert anchor in support
            c_support = tensor(*(Z if i in support else I2 for i in range(3)))
            zero(tag + '_local_conjugation_' + str(j), c_support @ op @ c_support - cs @ op @ cs)
            lifted = float(c) * (tensor(at(P0, anchor), op) + tensor(at(P1, anchor), c_support @ op @ c_support))
            h += float(c) * op
            hp += lifted
            zero(tag + '_local_intertwiner_' + str(j), lifted @ e - e @ (float(c) * op))
        zero(tag + '_hermitian', hp - hp.conj().T)
        zero(tag + '_full_symmetry', hp @ gamma_op - gamma_op @ hp)
        zero(tag + '_h_intertwiner', hp @ e - e @ (h + 2 * ids))
        for i in range(3):
            zr = tensor(at(Z, i), ids)
            zero(tag + '_reference_z_conserved_' + str(i), hp @ zr - zr @ hp)
        vals = np.linalg.eigvalsh(hp)
        assert float(vals[0]) >= float(spectral_bound) - 2e-12
        up, ul = unitary(hp, 1), unitary(h + 2 * ids, 1)
        zero(tag + '_u_intertwiner', up @ e - e @ ul)
        zero(tag + '_complete_energy_conservation', up.conj().T @ hp @ up - hp)
        ep = tensor(e, I2)
        zero(tag + '_passive_partner', tensor(up, I2) @ ep @ passive - ep @ tensor(ul, I2) @ passive)
        # Test a non-invariant local effect as well as the local witness.
        effect = (ids + at(X, 0)) / 2
        lifted_effect = tensor(at(P0, 0), effect) + tensor(at(P1, 0), cs @ effect @ cs)
        zero(tag + '_effect_intertwiner', lifted_effect @ e - e @ effect)
        zero(tag + '_effect_symmetry', lifted_effect @ gamma_op - gamma_op @ lifted_effect)
        ev = np.linalg.eigvalsh(lifted_effect)
        assert ev[0] > -2e-12 and ev[-1] < 1 + 2e-12
        state = up @ psi0
        p3 = float(np.vdot(state, tensor(ids, f3) @ state).real)
        pp = float(np.vdot(state, tensor(ids, fpar) @ state).real)
        expected = math.sin(float(gam))**2
        assert abs(p3 - expected) < 2e-12 and abs(pp - expected) < 2e-12
        zero(tag + '_final_total_even', gamma_op @ state - state)
        # The reduced reference generally changes; no unchanged catalyst claim.
        ref0 = psi0.reshape(8, 8) @ psi0.reshape(8, 8).conj().T
        ref1 = state.reshape(8, 8) @ state.reshape(8, 8).conj().T
        ref_distance = float(np.sum(np.abs(np.linalg.eigvalsh(ref1 - ref0))) / 2)
        assert abs(ref_distance - expected) < 2e-12
        cases.append({'gamma': tag, 'local_p3': p3, 'system_odd_probability': pp,
                      'reference_marginal_trace_distance': ref_distance,
                      'full_h_min_eigenvalue': float(vals[0])})
    # Covariant encoder impossibility is analytic; this pair witnesses the failure.
    plus = np.zeros(8, dtype=complex)
    plus[0] = plus[4] = 1 / math.sqrt(2)
    a, b = e @ plus, e @ cs @ plus
    overlap = abs(np.vdot(a, b))
    assert overlap < 2e-12
    assert cases[1]['local_p3'] - cases[0]['local_p3'] > float(probability_lower)
    return {'status': 'PASS_RELATIONAL_SELF_RESPONSE_SCOPE',
            'kind': 'finite_mathematical_encoding_not_actual_gr_parent',
            'physical_dimension': 64, 'logical_dimension': 8,
            'exact_spectral_lower_bound': str(spectral_bound),
            'exact_probability_gap_lower_bound': str(probability_lower),
            'cases': cases, 'max_numerical_residual': max(residuals.values()),
            'residuals': residuals, 'noncovariant_encoder_witness_overlap': float(overlap),
            'new_scientific_rounds': 0, 'new_adopted_axioms': 0,
            'original_1200_parent_modified': False, 'update_law_selected': False}


def compare(a, b):
    if isinstance(a, dict):
        assert set(a) == set(b)
        for k in a:
            compare(a[k], b[k])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            compare(x, y)
    elif isinstance(a, float):
        assert abs(a - b) < 1e-10
    else:
        assert a == b


if __name__ == '__main__':
    result = check()
    out = D / 'results.json'
    if sys.argv[1:] == ['--write']:
        assert not out.exists(), 'Refusing to overwrite saved results'
        out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    else:
        assert not sys.argv[1:], 'Use no arguments or --write'
        compare(result, json.loads(out.read_text(encoding='utf-8')))
    print(json.dumps({k: result[k] for k in ['status', 'cases', 'max_numerical_residual']}, ensure_ascii=False))
