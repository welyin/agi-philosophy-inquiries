"""Finite calibrations of proof.md; default execution is read-only."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

ATOL = 8e-11


def maxabs(x):
    return float(np.max(np.abs(x)))


def swap_matrix(d, n, a, b):
    rows = np.arange(d ** n).reshape((d,) * n)
    perm = np.swapaxes(rows, a, b).reshape(-1)
    return np.eye(d ** n)[perm]


def encoding(d):
    n = d + 1
    q = np.zeros((n, d))
    for k in range(1, n):
        q[:k, k - 1] = 1 / math.sqrt(k * (k + 1))
        q[k, k - 1] = -k / math.sqrt(k * (k + 1))
    e = np.zeros((d ** n, d, n))
    for r in range(n):
        sites = [j for j in range(n) if j != r]
        for p in itertools.permutations(range(d)):
            sign = (-1) ** sum(p[i] > p[j] for i in range(d) for j in range(i + 1, d))
            for g in range(d):
                x = [0] * n
                x[r] = g
                for site, value in zip(sites, p):
                    x[site] = value
                index = np.ravel_multi_index(tuple(x), (d,) * n)
                e[index, g, r] = (-1) ** r * sign / math.sqrt(math.factorial(d))
    w = math.sqrt(d / n) * np.einsum('xgr,rm->xgm', e, q)
    return e, q, w.reshape(d ** n, d * d)


def su_generators(d):
    out = []
    for i in range(d):
        for j in range(i + 1, d):
            x = np.zeros((d, d), complex)
            x[i, j] = x[j, i] = 0.5
            out.append(x)
            y = np.zeros((d, d), complex)
            y[i, j], y[j, i] = -0.5j, 0.5j
            out.append(y)
    for k in range(1, d):
        diag = np.zeros(d)
        diag[:k], diag[k] = 1, -k
        out.append(np.diag(diag) / math.sqrt(2 * k * (k + 1)))
    return out


def axis_apply(a, x, axis):
    return np.moveaxis(np.tensordot(a, x, axes=(1, axis)), 0, axis)


def evolution(h, time):
    values, vectors = np.linalg.eigh(h)
    return (vectors * np.exp(-1j * time * values)) @ vectors.conj().T


def run():
    rng = np.random.default_rng(20261009)
    rows, residuals = [], []
    for d in (2, 3, 4):
        n = d + 1
        e, q, w = encoding(d)
        gram = np.einsum('xgr,xhs->grhs', e, e)
        target = np.einsum('gh,rs->grhs', np.eye(d), (n * np.eye(n) - np.ones((n, n))) / d)
        errs = {'gram': maxabs(gram - target), 'isometry': maxabs(w.T @ w - np.eye(d * d))}
        wt = w.reshape((d,) * n + (d, d))
        all_internal = []
        for a, b in itertools.combinations(range(n), 2):
            perm = np.eye(n)
            perm[[a, b]] = perm[[b, a]]
            s = -q.T @ perm @ q
            all_internal.append(s)
            expected = np.einsum('xgm,mk->xgk', w.reshape(d ** n, d, d), s)
            residuals.append(maxabs(np.swapaxes(wt, a, b).reshape(d ** n, d, d) - expected))
        errs['internal_swap'] = max(residuals[-len(all_internal):])
        comm = float(np.linalg.norm(all_internal[0] @ all_internal[1] - all_internal[1] @ all_internal[0]))
        assert comm > 0.1
        gen_error = 0.0
        for t in su_generators(d):
            left = sum(axis_apply(t, wt, j) for j in range(n))
            right = np.einsum('xhm,hg->xgm', w.reshape(d ** n, d, d), t).reshape(wt.shape)
            gen_error = max(gen_error, maxabs(left - right))
        errs['collective_generator'] = gen_error
        wp = np.kron(w, np.eye(d))
        xp = wp.reshape((d,) * (n + 1) + (d, d, d))
        actual = sum(np.swapaxes(xp, j, n) for j in range(n))
        expected = xp + np.swapaxes(xp, n + 1, n + 3)
        errs['cross_primitive'] = maxabs(actual - expected)
        if d <= 3:
            ww = np.kron(w, w)
            xx = ww.reshape((d,) * (2 * n) + (d,) * 4)
            actual = sum(np.swapaxes(xx, i, j) for i in range(n) for j in range(n, 2 * n))
            expected = ((n * n - 1) / d) * xx + np.swapaxes(xx, 2 * n, 2 * n + 2)
            errs['cross_two_blocks'] = maxabs(actual - expected)
        residuals.extend(errs.values())
        rows.append({'d': d, 'minimum_members': n, 'private_dimension': d,
                     'physical_dimension': d ** n, 'code_dimension': d * d,
                     'private_commutator_norm': comm, 'residuals': errs})

    processes = []
    for d in (2, 3):
        n = d + 1
        _, q, w = encoding(d)
        wp = np.kron(w, np.eye(d))
        k = sum(swap_matrix(d, n + 1, i, n) for i in range(n))
        internal = swap_matrix(d, n + 1, 0, 1)
        perm = np.eye(n)
        perm[[0, 1]] = perm[[1, 0]]
        sm = -q.T @ perm @ q
        heff = np.eye(d ** 3) + swap_matrix(d, 3, 0, 2) + 0.37 * np.kron(np.eye(d), np.kron(sm, np.eye(d)))
        h = k + 0.37 * internal
        residuals.append(maxabs(h @ wp - wp @ heff))
        z = rng.normal(size=(d ** 3, 2)) + 1j * rng.normal(size=(d ** 3, 2))
        z /= np.linalg.norm(z)
        m = np.arange(1, d + 1) + 1j * np.arange(d, 0, -1)
        m /= np.linalg.norm(m)
        probs = []
        for t in (math.pi / 3, math.pi / 2, 2 * math.pi / 3):
            u = evolution(h, t)
            ue = evolution(heff, t)
            residuals.append(maxabs(u @ wp - wp @ ue))
            residuals.append(float(np.linalg.norm(u @ wp @ z - wp @ ue @ z)))
            pair = []
            for label in (0, 1):
                g = np.eye(d)[label]
                b = np.eye(d)[0]
                initial = wp @ np.kron(np.kron(g, m), b)
                final = (u @ initial).reshape(d ** n, d)
                # C stays |0>; expectation of (I-S_BC)/2 uses rho_B[0,0].
                p = (float(np.vdot(final, final).real) - float(np.vdot(final[:, 0], final[:, 0]).real)) / 2
                target = 0.0 if label == 0 else math.sin(t) ** 2 / 2
                residuals.append(abs(p - target))
                pair.append(p)
            probs.append({'time': t, 'p_zero': pair[0], 'p_one': pair[1]})
        processes.append({'d': d, 'microscopic_dimension': d ** (n + 1), 'reference_dimension': 2,
                          'internal_coefficient': 0.37, 'probabilities': probs})

    resources = []
    for d in (2, 3, 4, 5):
        s = swap_matrix(d, 2, 0, 1)
        rho = (np.eye(d * d) - s) / (d * (d - 1))
        marginal = np.einsum('ibjb->ij', rho.reshape(d, d, d, d))
        constraints = np.vstack([np.kron(t, np.eye(d)) + np.kron(np.eye(d), t) for t in su_generators(d)])
        singular = np.linalg.svd(constraints, compute_uv=False)
        nullity = int(np.sum(singular < 1e-10))
        assert nullity == (1 if d == 2 else 0)
        x = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
        u, _ = np.linalg.qr(x)
        u[:, 0] /= np.linalg.det(u)
        uu = np.kron(u, u)
        omega = np.eye(d).reshape(-1) / math.sqrt(d)
        residuals.extend([maxabs(marginal - np.eye(d) / d), maxabs(uu @ rho @ uu.conj().T - rho),
                          maxabs(np.kron(u, u.conj()) @ omega - omega), abs(float(np.trace(rho @ s).real) + 1)])
        # Independent full three-port comparison, including a common rotation.
        sab = swap_matrix(d, 3, 0, 1)
        effect = (np.eye(d ** 3) - swap_matrix(d, 3, 1, 2)) / 2
        unitary = math.cos(0.73) * np.eye(d ** 3) - 1j * math.sin(0.73) * sab
        for label in (0, 1):
            state = np.kron(np.kron(np.eye(d)[label], np.eye(d)[0]), np.eye(d)[0])
            rotated = np.kron(np.kron(u, u), u) @ state
            out = unitary @ state
            out_rotated = unitary @ rotated
            p = float(np.vdot(out, effect @ out).real)
            pr = float(np.vdot(out_rotated, effect @ out_rotated).real)
            residuals.extend([abs(pr - p), abs(p - label * math.sin(0.73) ** 2 / 2)])
        resources.append({'d': d, 'two_fundamental_singlet_dimension': nullity,
                          'mixed_antisymmetric_rank': d * (d - 1) // 2,
                          'mixed_antisymmetric_purity': float(np.trace(rho @ rho).real),
                          'swap_mean': float(np.trace(rho @ s).real)})
    tree = [{'d': d, 'depth': level, 'leaves': (d + 1) ** level,
             'private_factors': ((d + 1) ** level - 1) // d,
             'code_log_d_dimension': 1 + ((d + 1) ** level - 1) // d}
            for d in (2, 3, 4) for level in (1, 2, 3)]
    assert max(residuals) < ATOL, max(residuals)
    return {'round': 1063, 'scientific_calibration_groups': 1, 'new_empirical_groups': 0,
            'new_accepted_cognitive_axioms': 0,
            'scope': 'Finite rank calibration of the analytic minimal same-type SU(d) exchange hierarchy; not a full cognitive or spacetime countermodel.',
            'block_checks': rows, 'actual_processes': processes, 'pair_resource_checks': resources,
            'tree_dimension_formula': tree, 'maximum_residual': max(residuals),
            'all_assertions_passed': True}


def compare(a, b, path='result'):
    if isinstance(a, dict):
        assert a.keys() == b.keys(), path
        for key in a:
            compare(a[key], b[key], path + '/' + key)
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for index, (x, y) in enumerate(zip(a, b)):
            compare(x, y, path + '/' + str(index))
    elif isinstance(a, float):
        assert math.isclose(a, b, abs_tol=ATOL, rel_tol=1e-10), (path, a, b)
    else:
        assert a == b, (path, a, b)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--save-exclusive', action='store_true')
    args = parser.parse_args()
    result = run()
    output = Path(__file__).with_name('results.json')
    if args.save_exclusive:
        with output.open('x', encoding='utf8') as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
    else:
        compare(result, json.loads(output.read_text(encoding='utf8')))
    print(json.dumps({'passed': True, 'maximum_residual': result['maximum_residual'],
                      'mode': 'save-exclusive' if args.save_exclusive else 'read-only'}, ensure_ascii=False))
