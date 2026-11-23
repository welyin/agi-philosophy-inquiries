"""Independent 1327 audit, using combinatorial dressing (no author imports).

Default/--check is read-only. --write-new exclusively creates the result.
The physical parent keeps all CAR occupations at fixed Q=q and K=N.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def basis(modes, particles):
    return list(itertools.combinations(range(modes), particles))


def fermion_move(occupied, target, source):
    if source not in occupied:
        return None
    old_position = occupied.index(source)
    remaining = occupied[:old_position] + occupied[old_position + 1:]
    if target in remaining:
        return None
    new_position = sum(x < target for x in remaining)
    out = remaining[:new_position] + (target,) + remaining[new_position:]
    return out, (-1) ** (old_position + new_position)


def build_parent(n, q, N, matrix, rule):
    states = []
    for occupied in basis(3 * n, q):
        k = sum(x >= n for x in occupied)
        if k <= N:
            for blue in range(N - k + 1):
                states.append((occupied, N - k - blue, blue))
    index = {state: j for j, state in enumerate(states)}
    dim = len(states)
    forward = np.zeros((dim, dim), dtype=complex)
    raising = np.zeros_like(forward)
    left, right, photons, magnetic = [], [], [], []
    for col, (occupied, red, blue) in enumerate(states):
        k = sum(x >= n for x in occupied)
        left.append(k)
        right.append(q - k)
        photons.append(red + blue)
        n_red = sum(n <= x < 2 * n for x in occupied)
        n_blue = sum(x >= 2 * n for x in occupied)
        magnetic.append((red - blue + n_red - n_blue) / 2)
        if blue:
            row = index[(occupied, red + 1, blue - 1)]
            raising[row, col] += math.sqrt(blue * (red + 1))
        for flavor in range(n):
            moved = fermion_move(occupied, n + flavor, 2 * n + flavor)
            if moved:
                out, sign = moved
                raising[index[(out, red, blue)], col] += sign
        if red + blue == 0:
            continue
        denominator = math.sqrt(red + blue if rule == 'balance' else N)
        for color, number in enumerate((red, blue)):
            if not number:
                continue
            for i in range(n):
                for j in range(n):
                    moved = fermion_move(occupied, n * (color + 1) + i, j)
                    if moved:
                        out, sign = moved
                        row = index[(out, red - (color == 0), blue - (color == 1))]
                        forward[row, col] += sign * matrix[i, j] * math.sqrt(number) / denominator
    return dict(states=states, index=index, T=forward, Jplus=raising,
                nr=np.diag(right), nl=np.diag(left), nb=np.diag(photons),
                Jz=np.diag(magnetic))


def build_logical(n, q, matrix):
    states = basis(2 * n, q)
    index = {state: j for j, state in enumerate(states)}
    T = np.zeros((len(states), len(states)), complex)
    for col, occupied in enumerate(states):
        for i in range(n):
            for j in range(n):
                moved = fermion_move(occupied, n + i, j)
                if moved:
                    out, sign = moved
                    T[index[out], col] += sign * matrix[i, j]
    k = np.diag([sum(x >= n for x in occupied) for occupied in states])
    return states, T, k


def combinatorial_dressing(n, N, physical, logical):
    """Explicit binomial splitting; no SU(2) ladder is used to build D."""
    D = np.zeros((len(physical['states']), (N + 1) * len(logical)), complex)
    for h in range(N + 1):
        for col, occupied in enumerate(logical):
            r = tuple(x for x in occupied if x < n)
            flavors = tuple(x - n for x in occupied if x >= n)
            k = len(flavors)
            for colors in itertools.product((0, 1), repeat=k):
                blue = h - sum(colors)
                if not 0 <= blue <= N - k:
                    continue
                sequence = r + tuple(n * (color + 1) + f for color, f in zip(colors, flavors))
                inversions = sum(sequence[i] > sequence[j]
                                 for i in range(len(sequence)) for j in range(i + 1, len(sequence)))
                state = (tuple(sorted(sequence)), N - k - blue, blue)
                amplitude = (-1) ** inversions * math.sqrt(math.comb(N - k, blue) / math.comb(N, h))
                D[physical['index'][state], h * len(logical) + col] += amplitude
    return D


def unitary(H, time):
    values, vectors = np.linalg.eigh(H)
    return (vectors * np.exp(-1j * time * values)) @ vectors.conj().T


def compute():
    checks = []
    cases = []
    holdouts = []

    def close(name, left, right, tolerance=3e-10):
        # Frobenius residual dominates the operator norm for these finite matrices.
        residual = float(np.linalg.norm(np.asarray(left) - np.asarray(right)))
        checks.append(dict(name=name, residual=residual, tolerance=tolerance))
        if residual > tolerance:
            raise AssertionError((name, residual, tolerance))

    configs = [(2, 1, 1), (2, 1, 3), (2, 2, 2), (2, 2, 3),
               (2, 3, 2), (2, 4, 2), (3, 2, 3), (3, 3, 3)]
    for n, q, N in configs:
        rng = np.random.default_rng(711327 + 100 * n + 10 * q + N)
        M = (rng.uniform(-1, 1, (n, n)) + 1j * rng.uniform(-1, 1, (n, n))) / 5
        parent = build_parent(n, q, N, M, 'balance')
        alternative = build_parent(n, q, N, M, 'total')
        logical, T, k = build_logical(n, q, M)
        D = combinatorial_dressing(n, N, parent, logical)
        I = np.eye(len(logical))
        lift = lambda a: np.kron(np.eye(N + 1), a)
        tag = f'n{n}.q{q}.N{N}'
        expected_dimension = sum(math.comb(2 * n, j) * math.comb(n, q - j) * (N - j + 1)
                                 for j in range(min(N, q, 2 * n) + 1) if 0 <= q - j <= n)
        assert len(parent['states']) == expected_dimension
        close(tag + '.combinatorial_isometry', D.conj().T @ D, np.eye(D.shape[1]))
        close(tag + '.forward_full_map', parent['T'] @ D, D @ lift(T))
        close(tag + '.backward_full_map', parent['T'].conj().T @ D, D @ lift(T.conj().T))
        close(tag + '.R_source', parent['nr'] @ D, D @ lift(q * I - k))
        close(tag + '.L_source', parent['nl'] @ D, D @ lift(k))
        close(tag + '.reference_source', parent['nb'] @ D, D @ lift(N * I - k))
        close(tag + '.whole_K', parent['nb'] + parent['nl'], N * np.eye(expected_dimension))
        close(tag + '.whole_Q', parent['nr'] + parent['nl'], q * np.eye(expected_dimension))
        close(tag + '.old_occupation_dependence', alternative['T'] @ D,
              D @ lift(T @ np.diag(np.sqrt((N - np.diag(k)) / N))))
        if q == 1:
            close(tag + '.one_particle_full_parent_identity', parent['T'], alternative['T'])
        Jp = parent['Jplus']
        Jm = Jp.conj().T
        Jz = parent['Jz']
        Casimir = Jz @ Jz + (Jp @ Jm + Jm @ Jp) / 2
        close(tag + '.maximum_spin', Casimir @ D, N * (N + 2) / 4 * D)
        maximum_spin_multiplicity = int(np.count_nonzero(np.abs(np.linalg.eigvalsh(Casimir) - N * (N + 2) / 4) < 1e-9))
        assert maximum_spin_multiplicity == D.shape[1]
        close(tag + '.SU2_scalar_forward', Jp @ parent['T'], parent['T'] @ Jp)
        for i in range(n):
            for j in range(n):
                for component, coefficient in [('real', 1), ('imaginary', 1j)]:
                    dM = np.zeros((n, n), complex)
                    dM[i, j] = coefficient
                    pT = build_parent(n, q, N, dM, 'balance')['T']
                    lT = build_logical(n, q, dM)[1]
                    close(tag + f'.source.M_{i}{j}.{component}', (pT + pT.conj().T) @ D,
                          D @ lift(lT + lT.conj().T))
        ER, EL, omega, time = 3.2, 4.1, 0.6, 0.73
        V = parent['T'] + parent['T'].conj().T
        H = ER * parent['nr'] + EL * parent['nl'] + omega * parent['nb'] + V
        h = (ER * q + omega * N) * I + (EL - ER - omega) * k + T + T.conj().T
        close(tag + '.complete_H', H @ D, D @ lift(h))
        U = unitary(H, time)
        u = unitary(h, time)
        close(tag + '.entire_code_time_evolution', U @ D, D @ lift(u))
        P = D @ D.conj().T
        close(tag + '.dark_code_no_leakage', (np.eye(expected_dimension) - P) @ H @ P, np.zeros_like(H))
        close(tag + '.full_parent_energy', U.conj().T @ H @ U, H)
        instrument = []
        for kk in range(q + 1):
            selector = np.diag((np.diag(parent['nl']) == kk).astype(float))
            logical_selector = np.diag((np.diag(k) == kk).astype(float))
            branch = selector @ U @ D
            close(tag + f'.record_branch_{kk}', branch, D @ lift(logical_selector @ u))
            instrument.append(branch)
        instrument = np.concatenate(instrument, axis=0)
        close(tag + '.record_complete_isometry', instrument.conj().T @ instrument, np.eye(D.shape[1]))
        reference_record_work = None
        if q <= n:
            pure_right = [j for j, occupied in enumerate(logical) if all(x < n for x in occupied)]
            columns = [hh * len(logical) + j for hh in range(N + 1) for j in pure_right]
            E = D[:, columns]
            # General pure state of all unknown right/reference input and a 3-dimensional partner.
            input_coefficients = rng.normal(size=(len(columns), 3)) + 1j * rng.normal(size=(len(columns), 3))
            input_coefficients /= np.linalg.norm(input_coefficients)
            initial = E @ input_coefficients
            evolved = U @ initial
            reference_record_work = 0.0
            recorded_norm = 0.0
            for kk in range(q + 1):
                branch = (np.diag(parent['nl']) == kk)[:, None] * evolved
                reference_record_work += np.vdot(branch, H @ branch).real
                recorded_norm += np.vdot(branch, branch).real
            reference_record_work -= np.vdot(evolved, H @ evolved).real
            close(tag + '.explicit_unknown_partner_branch_normalization', recorded_norm, 1.0)
            close(tag + '.explicit_partner_energy', np.vdot(initial, H @ initial), np.vdot(evolved, H @ evolved))
            close(tag + '.record_work_not_deleted', reference_record_work, -np.vdot(evolved, V @ evolved).real)
            delta_b = np.vdot(evolved, parent['nb'] @ evolved) - np.vdot(initial, parent['nb'] @ initial)
            delta_l = np.vdot(evolved, parent['nl'] @ evolved) - np.vdot(initial, parent['nl'] @ initial)
            close(tag + '.partner_reference_consumption', delta_b + delta_l, 0)
        cases.append(dict(n=n, q=q, N=N, physical_dimension=expected_dimension,
                          dressed_dimension=D.shape[1], dark_dimension=expected_dimension-D.shape[1],
                          maximum_spin_multiplicity=maximum_spin_multiplicity,
                          sample_record_work=reference_record_work))

    for N in (2, 3):
        n = q = 2
        M = np.eye(n)
        parent = build_parent(n, q, N, M, 'balance')
        alternative = build_parent(n, q, N, M, 'total')
        logical, _, k = build_logical(n, q, M)
        D = combinatorial_dressing(n, N, parent, logical)
        input_index = logical.index((0, 1))
        E = D[:, [h * len(logical) + input_index for h in range(N + 1)]]
        output_selector = np.diag((np.diag(parent['nl']) == 2).astype(float))
        frequency = math.sqrt(2 * (2 * N - 1) / N)
        cap = 4 * N * (N - 1) / (2 * N - 1) ** 2
        row = dict(N=N, predicted_old_cap=cap, minimum_old_failure=1-cap, observations=[])
        for rule, time in [('balance', math.pi/2), ('total', math.pi/2), ('total', math.pi/frequency)]:
            p = parent if rule == 'balance' else alternative
            out = unitary(p['T'] + p['T'].conj().T, time) @ E
            success = out.conj().T @ output_selector @ out
            expected = 1.0 if rule == 'balance' else cap * math.sin(frequency*time/2)**4
            close(f'holdout.N{N}.{rule}.t{time:.8f}.all_unknown_reference_success', success, expected*np.eye(N+1))
            close(f'holdout.N{N}.{rule}.t{time:.8f}.whole_reference', out.conj().T @ (parent['nb'] + parent['nl']) @ out, N*np.eye(N+1))
            mean_left = float(np.trace(out.conj().T @ parent['nl'] @ out).real / (N+1))
            row['observations'].append(dict(rule=rule, time=time, success=float(np.trace(success).real/(N+1)),
                                            mean_left=mean_left, reference_consumed=mean_left))
        holdouts.append(row)
    return dict(round=1327, status='PASS', stage_complete=False, new_axioms=0,
                construction='Explicit binomial reference splitting with independent tuple-CAR signs; no author imports or highest-weight ladder in dressing.',
                residual_norm='Frobenius', check_count=len(checks), max_residual=max(x['residual'] for x in checks),
                cases=cases, holdouts=holdouts, checks=checks)


def compare(saved, current, path='root'):
    if isinstance(current, dict):
        assert saved.keys() == current.keys(), path
        for key in current:
            compare(saved[key], current[key], path + '.' + key)
    elif isinstance(current, list):
        assert len(saved) == len(current), path
        for i, (a, b) in enumerate(zip(saved, current)):
            compare(a, b, path + f'[{i}]')
    elif isinstance(current, float):
        assert math.isfinite(saved) and abs(saved-current) <= 1e-9 * (1 + abs(current)), path
    else:
        assert saved == current, path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--write-new', action='store_true')
    arguments = parser.parse_args()
    result = compute()
    destination = HERE / 'independent_results.json'
    if arguments.write_new:
        with destination.open('x', encoding='utf-8') as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
    else:
        compare(json.loads(destination.read_text(encoding='utf-8')), result)
    print(json.dumps({key: result[key] for key in ('round', 'status', 'check_count', 'max_residual', 'cases', 'holdouts')}, ensure_ascii=False))


if __name__ == '__main__':
    main()
