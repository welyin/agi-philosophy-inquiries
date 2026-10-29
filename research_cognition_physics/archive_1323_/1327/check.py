"""Round 1327: complete CAR parent, reference dressing, sources and holdouts.

NumPy only. Default is read-only comparison with the saved deterministic result.
--write-new creates results.json exclusively; no historical file is overwritten.
"""
from __future__ import annotations
import argparse
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def masks(m, q):
    return [sum(1 << j for j in c) for c in itertools.combinations(range(m), q)]


def hop(s, i, j):
    if not (s >> j) & 1:
        return None
    sign = (-1) ** ((s & ((1 << j) - 1)).bit_count())
    t = s ^ (1 << j)
    if (t >> i) & 1:
        return None
    sign *= (-1) ** ((t & ((1 << i) - 1)).bit_count())
    return t | (1 << i), sign


def norm(a):
    return float(np.linalg.norm(a, ord=2))


def evolve(h, t):
    w, v = np.linalg.eigh(h)
    return (v * np.exp(-1j * t * w)) @ v.conj().T


def parent(n, q, N, M, rule):
    basis = []
    for s in masks(3 * n, q):
        k = (s >> n).bit_count()
        if k <= N:
            basis += [(s, N - k - h, h) for h in range(N - k + 1)]
    lookup = {x: j for j, x in enumerate(basis)}
    d = len(basis)
    t = np.zeros((d, d), complex)
    jm = np.zeros_like(t)
    nr, nl, nb, jz = [], [], [], []
    for c, (s, a, b) in enumerate(basis):
        k = (s >> n).bit_count()
        nr.append(q - k); nl.append(k); nb.append(a + b)
        n0 = ((s >> n) & ((1 << n) - 1)).bit_count()
        n1 = (s >> (2 * n)).bit_count()
        jz.append((a - b + n0 - n1) / 2)
        if a:
            jm[lookup[(s, a - 1, b + 1)], c] += math.sqrt(a * (b + 1))
        for i in range(n):
            flipped = hop(s, 2 * n + i, n + i)
            if flipped:
                ss, sign = flipped
                jm[lookup[(ss, a, b)], c] += sign
        for alpha, occ in enumerate((a, b)):
            if not occ:
                continue
            factor = (1 / math.sqrt(a + b) if rule == 'balance'
                      else 1 / math.sqrt(N) if rule == 'total' else 1)
            for i in range(n):
                for j in range(n):
                    moved = hop(s, (alpha + 1) * n + i, j)
                    if moved:
                        ss, sign = moved
                        state = (ss, a - (alpha == 0), b - (alpha == 1))
                        t[lookup[state], c] += sign * M[i, j] * math.sqrt(occ) * factor
    return dict(basis=basis, lookup=lookup, T=t, Jm=jm,
                nr=np.diag(nr), nl=np.diag(nl), nb=np.diag(nb), Jz=np.diag(jz))


def logical(n, q, M):
    bs = masks(2 * n, q)
    ix = {s: j for j, s in enumerate(bs)}
    t = np.zeros((len(bs), len(bs)), complex)
    k = np.array([(s >> n).bit_count() for s in bs])
    for c, s in enumerate(bs):
        for i in range(n):
            for j in range(n):
                out = hop(s, n + i, j)
                if out:
                    ss, sign = out
                    t[ix[ss], c] += sign * M[i, j]
    return bs, t, np.diag(k)


def dressing(p, bs, N):
    top = np.zeros((len(p['basis']), len(bs)), complex)
    for j, s in enumerate(bs):
        # Physical highest-weight mask has exactly the same R,L0 ordering.
        candidates = [z for z in p['lookup'] if z[0] == s and z[2] == 0]
        assert len(candidates) == 1
        top[p['lookup'][candidates[0]], j] = 1
    blocks = [top]
    for h in range(1, N + 1):
        blocks.append(p['Jm'] @ blocks[-1] / math.sqrt(h * (N - h + 1)))
    return np.concatenate(blocks, axis=1)


def run():
    checks = []
    cases = []
    def eq(name, a, b, tol=2e-10):
        error = norm(np.asarray(a) - np.asarray(b))
        checks.append(dict(name=name, residual=error, tolerance=tol))
        if error > tol:
            raise AssertionError((name, error, tol))

    # q>n cases include logical R_i and L_i simultaneously occupied.
    configurations = [(2, 1, 1), (2, 1, 3), (2, 2, 2), (2, 2, 3),
                      (2, 3, 2), (2, 4, 2), (3, 2, 2), (3, 3, 3)]
    for n, q, N in configurations:
        rng = np.random.default_rng(1327000 + 100 * n + 10 * q + N)
        M = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / 7
        p = parent(n, q, N, M, 'balance')
        old = parent(n, q, N, M, 'total')
        bs, t, nl = logical(n, q, M)
        D = dressing(p, bs, N)
        I = np.eye(len(bs)); IR = np.eye(N + 1)
        lift = lambda x: np.kron(IR, x)
        label = f'n{n}.q{q}.N{N}'
        eq(label + '.isometry', D.conj().T @ D, np.eye(D.shape[1]))
        eq(label + '.forward', p['T'] @ D, D @ lift(t))
        eq(label + '.backward', p['T'].conj().T @ D, D @ lift(t.conj().T))
        eq(label + '.right_source', p['nr'] @ D, D @ lift(q * I - nl))
        eq(label + '.left_source', p['nl'] @ D, D @ lift(nl))
        eq(label + '.reference_source', p['nb'] @ D, D @ lift(N * I - nl))
        eq(label + '.total_number', (p['nb'] + p['nl']) @ D, N * D)
        Jp = p['Jm'].conj().T
        casimir = p['Jz'] @ p['Jz'] + (Jp @ p['Jm'] + p['Jm'] @ Jp) / 2
        eq(label + '.highest_spin', casimir @ D, N * (N + 2) / 4 * D)
        multiplicity = int(np.count_nonzero(abs(np.linalg.eigvalsh(casimir) - N * (N + 2) / 4) < 1e-9))
        assert multiplicity == D.shape[1]
        eq(label + '.scalar_SU2', p['T'] @ p['Jm'], p['Jm'] @ p['T'])
        oldt = t @ np.diag(np.sqrt((N - np.diag(nl)) / N))
        eq(label + '.old_extension', old['T'] @ D, D @ lift(oldt))
        if q == 1:
            eq(label + '.one_fermion_identical_full_parent', p['T'], old['T'])
        E0, omega, momentum = 3.0, 0.7, 0.19
        H0 = (E0 + momentum) * p['nr'] + (E0 + omega - momentum) * p['nl'] + omega * p['nb']
        H = H0 + p['T'] + p['T'].conj().T
        logical_H = (E0 * q + omega * N) * I + momentum * (q * I - 2 * nl) + t + t.conj().T
        eq(label + '.complete_generator', H @ D, D @ lift(logical_H))
        U = evolve(H, 0.83)
        eq(label + '.all_code_dynamics', U @ D, D @ lift(evolve(logical_H, 0.83)))
        eq(label + '.momentum_source', (p['nr'] - p['nl']) @ D, D @ lift(q * I - 2 * nl))
        B = np.zeros((n, n), complex); B[0, -1] = 0.6 + 0.8j
        dparent = parent(n, q, N, B, 'balance')['T']
        _, dlogical, _ = logical(n, q, B)
        eq(label + '.offdiagonal_source', (dparent + dparent.conj().T) @ D,
           D @ lift(dlogical + dlogical.conj().T))
        records = []
        target_records = []
        for k in range(min(q, n) + 1):
            P = np.diag((np.diag(p['nl']) == k).astype(float))
            LP = np.diag((np.diag(nl) == k).astype(float))
            eq(label + f'.record_{k}', P @ D, D @ lift(LP))
            records.append(P @ U @ D)
            target_records.append(D @ lift(LP @ evolve(logical_H, 0.83)))
        record_map = np.concatenate(records, axis=0)
        eq(label + '.complete_record_and_partners', record_map, np.concatenate(target_records, axis=0))
        eq(label + '.record_isometry', record_map.conj().T @ record_map, np.eye(D.shape[1]))
        # These are physical full-input map identities, so arbitrary partners are covered.
        readout_work = None
        if q <= n:
            right_cols = [j for j, s in enumerate(bs) if (s >> n) == 0]
            columns = [h * len(bs) + j for h in range(N + 1) for j in right_cols]
            E = D[:, columns]
            rho = (U @ E) @ (U @ E).conj().T / len(columns)
            after = np.zeros_like(rho)
            for k in range(q + 1):
                P = np.diag((np.diag(p['nl']) == k).astype(float))
                after += P @ rho @ P
            W = float(np.trace(H @ (after - rho)).real)
            interaction = float(np.trace((p['T'] + p['T'].conj().T) @ rho).real)
            eq(label + '.record_work', np.array([[W]]), np.array([[-interaction]]))
            eq(label + '.energy_conservation', U.conj().T @ H @ U, H)
            readout_work = W
        cases.append(dict(n=n, q=q, N=N, full_parent_dimension=len(p['basis']),
                          dressed_dimension=D.shape[1], dark_dimension=len(p['basis']) - D.shape[1],
                          occupancy_readout_work=readout_work))

    holdouts = []
    for N in (2, 3):
        n = q = 2; g = 1.0
        M = np.eye(n) * g
        p = parent(n, q, N, M, 'balance')
        old = parent(n, q, N, M, 'total')
        v = np.zeros(len(p['basis']), complex)
        v[p['lookup'][((1 << n) - 1, N, 0)]] = 1
        end = np.diag(p['nl']) == q
        H = p['T'] + p['T'].conj().T
        H_old = old['T'] + old['T'].conj().T
        T = np.pi / (2 * g)
        Omega = g * math.sqrt(2 * (2 * N - 1) / N)
        t_opt = np.pi / Omega
        cap = 4 * N * (N - 1) / (2 * N - 1) ** 2
        row = dict(N=N, g=g, calibration_time=T, old_optimal_time=t_opt,
                   old_analytic_global_cap=cap, unavoidable_failure=1-cap)
        for key, hh, tt in [('balance_same_clock', H, T), ('old_same_clock', H_old, T), ('old_optimal', H_old, t_opt)]:
            out = evolve(hh, tt) @ v
            prob = float(np.vdot(out[end], out[end]).real)
            mean_left = float(np.vdot(out, p['nl'] @ out).real)
            row[key] = dict(all_delivered=prob, mean_left=mean_left, reference_consumed=mean_left)
            predicted = (1.0 if key == 'balance_same_clock' else
                         cap * math.sin(Omega * tt / 2) ** 4)
            eq(f'holdout.N{N}.{key}', np.array([[prob]]), np.array([[predicted]]))
            eq(f'holdout.N{N}.{key}.resource', np.array([[float(np.vdot(out,p['nb']@out).real)+mean_left]]), np.array([[N]]))
        holdouts.append(row)

    return dict(round=1327, status='PASS', new_axioms=0, stage_complete=False,
                check_count=len(checks), max_residual=max(x['residual'] for x in checks),
                cases=cases, holdouts=holdouts, checks=checks,
                scope='Finite complete CAR parent; arbitrary unknown reference and partners on declared dressed code. No local gauge-field or empirical parameter identification.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-new', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    target = HERE / 'results.json'
    if args.write_new:
        with target.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    elif target.exists():
        saved = json.loads(target.read_text(encoding='utf-8'))
        assert saved['status'] == 'PASS' and saved['check_count'] == result['check_count']
        assert [(x['n'], x['q'], x['N'], x['full_parent_dimension']) for x in saved['cases']] == [(x['n'], x['q'], x['N'], x['full_parent_dimension']) for x in result['cases']]
        for a, b in zip(saved['holdouts'], result['holdouts']):
            for key in ('balance_same_clock', 'old_same_clock', 'old_optimal'):
                assert abs(a[key]['all_delivered'] - b[key]['all_delivered']) < 1e-10
    else:
        raise FileNotFoundError('Create results once with --write-new, then use read-only --check')
    print(json.dumps({k: result[k] for k in ('round','status','check_count','max_residual','cases','holdouts')}, ensure_ascii=False))


if __name__ == '__main__':
    main()
