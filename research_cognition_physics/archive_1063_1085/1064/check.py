"""Finite calibrations of the literal-member interface classification.

Default: recompute and compare; --write exclusively creates results.json.
No historical files or libraries beyond NumPy are modified.
"""
from pathlib import Path
from itertools import permutations
import argparse
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent


def maximum(x):
    return float(np.max(np.abs(x)))


def row_swap(x, d, n, a, b):
    tail = x.shape[1:]
    return np.swapaxes(x.reshape((d,) * n + tail), a, b).reshape(x.shape)


def act(x, a, axis, d, n):
    y = x.reshape((d,) * n + x.shape[1:])
    y = np.moveaxis(np.tensordot(a, y, axes=(1, axis)), 0, axis)
    return y.reshape(x.shape)


def generators(d):
    out = []
    for i in range(d):
        for j in range(i + 1, d):
            x = np.zeros((d, d), complex)
            x[i, j] = x[j, i] = 1
            y = np.zeros((d, d), complex)
            y[i, j], y[j, i] = -1j, 1j
            out.extend((x, y))
    for i in range(d - 1):
        z = np.zeros((d, d), complex)
        z[i, i], z[i + 1, i + 1] = 1, -1
        out.append(z)
    return out


def epsilon(d):
    x = np.zeros(d ** d)
    for p in permutations(range(d)):
        parity = sum(p[i] > p[j] for i in range(d) for j in range(i + 1, d))
        x[np.ravel_multi_index(p, (d,) * d)] = (-1) ** parity / math.sqrt(math.factorial(d))
    return x


def localized_code(d):
    u = np.kron(epsilon(d), epsilon(d))
    w = row_swap(u, d, 2 * d, 0, d)
    v = np.column_stack((u, (w - u / d) / math.sqrt(1 - 1 / d ** 2)))
    return np.kron(np.eye(d), v), v, float(u @ w)


def local_calibration(d):
    n = 2 * d + 1
    w, v, overlap = localized_code(d)
    r = {'epsilon_overlap': abs(overlap - 1 / d),
         'private_isometry': maximum(v.conj().T @ v - np.eye(2)),
         'code_isometry': maximum(w.conj().T @ w - np.eye(2 * d))}
    r['singlet_generators'] = max(maximum(sum(act(v, t, k, d, 2*d)
                                                   for k in range(2*d))) for t in generators(d))
    r['full_collective_intertwiner'] = max(maximum(
        sum(act(w, t, k, d, n) for k in range(n)) - w @ np.kron(t, np.eye(2)))
        for t in generators(d))
    ext = np.kron(w, np.eye(d))
    logical_swap = np.arange(d * 2 * d).reshape(d, 2, d).swapaxes(0, 2).ravel()
    expected = ext[:, logical_swap]
    actual = row_swap(ext, d, n + 1, 0, n)
    r['literal_exchange_intertwiner'] = maximum(actual - expected)
    theta = .371
    r['full_time_exchange_intertwiner'] = maximum(
        np.cos(theta) * ext - 1j * np.sin(theta) * actual
        - (np.cos(theta) * ext - 1j * np.sin(theta) * expected))
    # Wrong physical member genuinely leaks out of the original fixed code.
    wrong = row_swap(ext, d, n + 1, 1, n)
    leakage = float(np.linalg.norm(wrong - ext @ (ext.conj().T @ wrong)))
    assert leakage > 1
    moved = row_swap(w, d, n, 0, 1)
    moved_ext = np.kron(moved, np.eye(d))
    r['moved_code_new_member_intertwiner'] = maximum(
        row_swap(moved_ext, d, n + 1, 1, n) - moved_ext[:, logical_swap])
    fixed_code_handoff_leak = float(np.linalg.norm(moved - w @ (w.conj().T @ moved)))
    assert fixed_code_handoff_leak > 1
    # Full original-space implementation of an invariant private update.
    hm = np.array([[.3, .4 + .2j], [.4 - .2j, -.3]])
    hr = v @ hm @ v.conj().T
    physical_internal_action = np.kron(np.eye(d), hr @ v)
    r['private_internal_intertwiner'] = maximum(
        physical_internal_action - w @ np.kron(np.eye(d), hm))
    heff = np.kron(np.kron(np.eye(d), hm), np.eye(d))
    seff = np.eye(d * 2 * d)[:, logical_swap]
    r['private_external_commutator'] = maximum(heff @ seff - seff @ heff)
    return {'d': d, 'members': n, 'private_dimension_used': 2,
            'complete_minimal_singlet_dimension': math.comb(2*d, d) // (d+1),
            'wrong_member_fixed_code_leak_frobenius': leakage,
            'same_code_handoff_leak_frobenius': fixed_code_handoff_leak,
            'residuals': r}


def swap_indices(d, n, a, b):
    return np.arange(d ** n).reshape((d,) * n).swapaxes(a, b).ravel()


def pipeline(d, gateway):
    # Rows of the complete unitary permutation, never a few selected inputs.
    p = np.arange(d ** 5)
    for a, b in ((gateway, 3), (gateway, 1), (gateway, 4)):
        p = p[swap_indices(d, 5, a, b)]
    return p


def dynamic_calibration(d):
    pa, pb = pipeline(d, 0), pipeline(d, 2)
    old = np.array(np.unravel_index(np.arange(d ** 5), (d,) * 5)).T
    new = old[:, [4, 3, 2, 0, 1]]
    mapping = np.ravel_multi_index(new.T, (d,) * 5)
    assert np.array_equal(pa[mapping], np.arange(d ** 5))
    k = swap_indices(d, 5, 0, 2)
    assert np.array_equal(k[pb], pa[k])  # T_B K = K T_A.
    # Direct input->output identities, including all spectator basis states.
    assert np.array_equal(new[:, 1], old[:, 3])
    assert np.array_equal(new[:, 4], old[:, 1])
    assert len(set(pa.tolist())) == d ** 5
    return {'d': d, 'operator_basis_columns_checked': d ** 5,
            'unknown_input_to_memory': True, 'old_memory_to_output': True,
            'handoff_intertwining': True, 'complete_permutation_unitary': True,
            'external_direct_contacts': ['A-I', 'A-O'],
            'contact_switch_is_an_input': True,
            'orthogonal_memory_output_probability_gap': 1}


def commutator_calibration(d):
    rng = np.random.default_rng(106400 + d)
    dimr = 2
    a = rng.normal(size=(d * dimr, d * dimr)) + 1j * rng.normal(size=(d * dimr, d * dimr))
    q, _ = np.linalg.qr(a)
    p = q[:, :2] @ q[:, :2].conj().T
    p_ext = np.kron(p, np.eye(d))
    ids = np.arange(d * dimr * d).reshape(d, dimr, d).swapaxes(0, 2).ravel()
    s = np.eye(d * dimr * d)[ids]
    lhs = np.linalg.norm(p_ext @ s - s @ p_ext) ** 2
    rhs = 0.
    for i in range(d):
        for j in range(d):
            e = np.zeros((d, d)); e[i, j] = 1
            e = np.kron(e, np.eye(dimr))
            rhs += np.linalg.norm(p @ e - e @ p) ** 2
    error = abs(float(lhs - rhs))
    assert lhs > .1 and error < 1e-10
    return {'d': d, 'matrix_unit_decomposition_error': error,
            'nonlocalized_projector_commutator_squared': float(lhs)}


def run():
    local = [local_calibration(d) for d in (2, 3)]
    dynamic = [dynamic_calibration(d) for d in (2, 3, 4, 5)]
    decomposition = [commutator_calibration(d) for d in (2, 3, 4)]
    residual = max([v for row in local for v in row['residuals'].values()] +
                   [row['matrix_unit_decomposition_error'] for row in decomposition])
    assert residual < 1e-10
    return {'round': 1064, 'passed': True, 'new_scientific_groups': 1,
            'empirical_groups': 0, 'accepted_cognitive_axioms': 0,
            'scope': 'Exact fixed-code literal-member classification; weak scheduled examples are reused exchange tools, not autonomous cognition or spacetime.',
            'localized_codes': local, 'weak_dynamic_examples': dynamic,
            'matrix_unit_checks': decomposition, 'maximum_residual': residual}


def compare(a, b, path='root'):
    if isinstance(a, dict):
        assert a.keys() == b.keys(), path
        for k in a: compare(a[k], b[k], path + '.' + k)
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)): compare(x, y, path + f'[{i}]')
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-11), (path, a, b)
    else:
        assert a == b, (path, a, b)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    target = HERE / 'results.json'
    if args.write:
        with target.open('x', encoding='utf8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2); f.write('\n')
    else:
        compare(result, json.loads(target.read_text(encoding='utf8')))
    print(json.dumps({'round': 1064, 'passed': True, 'maximum_residual': result['maximum_residual']}, ensure_ascii=False))
