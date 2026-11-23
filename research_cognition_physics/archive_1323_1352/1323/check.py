"""Round 1323: finite full-sector extraction checks; read-only by default."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = STAGE.parent.parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def same_result(actual, saved, path='result'):
    """Check all saved evidence, including secondary diagnostics and bindings."""
    if isinstance(actual, dict):
        assert isinstance(saved, dict) and actual.keys() == saved.keys(), path
        for key in actual:
            same_result(actual[key], saved[key], path + '.' + key)
    elif isinstance(actual, list):
        assert isinstance(saved, list) and len(actual) == len(saved), path
        for i, (a, b) in enumerate(zip(actual, saved)):
            same_result(a, b, path + '[' + str(i) + ']')
    elif isinstance(actual, float):
        assert isinstance(saved, (float, int)) and not isinstance(saved, bool), path
        assert np.isfinite(saved) and np.isclose(actual, saved, rtol=1e-10, atol=1e-10), path
    else:
        assert type(actual) is type(saved) and actual == saved, path


def annihilators(N):
    mats = [np.zeros((N, N + 1), complex) for _ in range(2)]
    for k in range(N + 1):
        if k:
            mats[0][k - 1, k] = np.sqrt(k)
        if k < N:
            mats[1][k, k] = np.sqrt(N - k)
    return mats


def extraction(N):
    if N == 0:
        return np.zeros((0, 1), complex)
    return np.vstack(annihilators(N)) / np.sqrt(N)


def spin(two_j):
    dim = two_j + 1
    up = np.zeros((dim, dim), complex)
    for k in range(two_j):
        up[k + 1, k] = np.sqrt((k + 1) * (two_j - k))
    return [(up + up.conj().T) / 2,
            (up - up.conj().T) / (2j),
            np.diag(np.arange(dim) - two_j / 2)]


def spin_left():
    # alpha=0 means up, while the reference occupation basis increases n_up.
    return [np.array([[0, 1], [1, 0]], complex) / 2,
            np.array([[0, -1j], [1j, 0]], complex) / 2,
            np.diag([.5, -.5])]


def evolve(H, t):
    vals, vecs = np.linalg.eigh(H)
    return (vecs * np.exp(-1j * t * vals)) @ vecs.conj().T


def parent(N, M, normalized=True, er=1., el=2., omega=1.):
    n = M.shape[0]
    nr = (N + 1) * n
    nl = 2 * N * n
    W = extraction(N)
    T = np.kron(W, M) * (1 if normalized else np.sqrt(N))
    V = np.zeros((nr + nl, nr + nl), complex)
    V[nr:, :nr], V[:nr, nr:] = T, T.conj().T
    PR = np.diag([1.] * nr + [0.] * nl)
    PL = np.eye(nr + nl) - PR
    NB = N * PR + (N - 1) * PL
    H0 = er * PR + el * PL + omega * NB
    return H0 + V, V, PR, PL, NB, W, nr


def run():
    checks = []

    def close(name, value, expected=0., tol=3e-10):
        residual = float(np.linalg.norm(np.asarray(value) - np.asarray(expected)))
        if not residual <= tol:
            raise AssertionError((name, residual, tol))
        checks.append(dict(name=name, residual=residual, tolerance=tol, passed=True))

    rng = np.random.default_rng(1323)
    capacity = []
    for N in range(1, 5):
        W = extraction(N)
        close(f'N{N}.isometry', W.conj().T @ W, np.eye(N + 1))
        close(f'N{N}.dark_dimension', 2 * N - np.linalg.matrix_rank(W), N - 1)
        for i, (a, b, c) in enumerate(zip(spin_left(), spin(N - 1), spin(N))):
            Jout = np.kron(a, np.eye(N)) + np.kron(np.eye(2), b)
            close(f'N{N}.intertwiner{i}', Jout @ W - W @ c)
        for dr in (1, 2, 3):
            capacity.append(dict(N=N, right_dimension=dr,
                                 input_per_flavor=dr * (N + 1), output_per_flavor=2 * N,
                                 recovery_error_lower_bound=max(0., 1 - 2 * N / (dr * (N + 1)))))
        for n in (3, 4, 5):
            tag = f'N{N}.n{n}'
            a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
            b = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
            u, _ = np.linalg.qr(a)
            v, _ = np.linalg.qr(b)
            M = u @ np.diag(np.linspace(.025, .1, n)) @ v.conj().T
            H, V, PR, PL, NB, W, nr = parent(N, M)
            d = H.shape[0]
            close(tag + '.hermitian', H - H.conj().T)
            close(tag + '.full_number', NB + PL, N * np.eye(d))
            close(tag + '.commutator', H @ (NB + PL) - (NB + PL) @ H)
            close(tag + '.energy_resonance', (H - V) @ V - V @ (H - V))
            assert np.linalg.eigvalsh(H)[0] > 0
            t = 1.37
            U = evolve(V, t)
            close(tag + '.unitarity', U.conj().T @ U, np.eye(d))
            u1, s, vh = np.linalg.svd(M)
            A = np.kron(np.eye(N + 1), vh.conj().T @ np.diag(np.cos(t * s)) @ vh)
            B = np.kron(W, -1j * u1 @ np.diag(np.sin(t * s)) @ vh)
            model = np.vstack([A, B])
            close(tag + '.all_input_formula', U[:, :nr], model)
            close(tag + '.success_failure', A.conj().T @ A + B.conj().T @ B, np.eye(nr))
            Pdark = np.zeros((d, d), complex)
            Pdark[nr:, nr:] = np.kron(np.eye(2 * N) - W @ W.conj().T, np.eye(n))
            close(tag + '.dark_retained', U @ Pdark, Pdark)
            psi = rng.normal(size=nr * 2) + 1j * rng.normal(size=nr * 2)
            psi /= np.linalg.norm(psi)
            initial = np.zeros(d * 2, complex)
            initial[:nr * 2] = psi
            final = np.kron(U, np.eye(2)) @ initial
            close(tag + '.unknown_partner', final, np.kron(model, np.eye(2)) @ psi)

            def mean(O, state):
                return float(np.vdot(state, np.kron(O, np.eye(2)) @ state).real)

            p = mean(PL, final)
            close(tag + '.reference_balance', mean(NB, final) - mean(NB, initial), -p)
            close(tag + '.right_balance', mean(PR, final) - mean(PR, initial), -p)
            close(tag + '.total_energy', mean(H, final), mean(H, initial))
            for param, deriv in [('er', PR), ('el', PL), ('omega', NB)]:
                delta = 1e-5
                defaults = dict(er=1., el=2., omega=1.)
                plus, minus = dict(defaults), dict(defaults)
                plus[param] += delta
                minus[param] -= delta
                finite = (parent(N, M, **plus)[0] - parent(N, M, **minus)[0]) / (2 * delta)
                close(tag + '.source_' + param, finite, deriv, 3e-8)

    # Independent common clock calibrated at N=1, without refitting at N=4.
    clock_results = []
    gamma = .2
    time = np.pi / (2 * gamma)
    for n in (3, 4, 5):
        for N in range(1, 5):
            for normalized in (True, False):
                _, V, _, PL, _, W, nr = parent(N, gamma * np.eye(n), normalized)
                U = evolve(V, time)
                success = U[nr:, :nr]
                rate = 1. if normalized else np.sqrt(N)
                expected = np.sin(rate * gamma * time) ** 2
                close(f'clock.N{N}.n{n}.{normalized}', success.conj().T @ success,
                      expected * np.eye(nr))
                if normalized:
                    close(f'clock.N{N}.n{n}.full_transfer', success, -1j * np.kron(W, np.eye(n)))
                clock_results.append(dict(N=N, n=n, normalized=normalized,
                                          probability=float(np.trace(success.conj().T @ success).real / nr)))

    # The whole endpoint amplitude, not just its modulus, admits periodic aliases.
    for N in (1, 2, 4):
        _, V, _, _, _, W, nr = parent(N, (4 * (N - 1) + 1) * gamma * np.eye(1))
        close(f'alias.N{N}.same_endpoint', evolve(V, time)[nr:, :nr], -1j * W)
    alias_rates = [4 * (N - 1) + 1 for N in (1, 2, 4)]
    assert len(set(alias_rates)) > 1

    # Sequential extraction has a fixed coefficient and keeps joint information.
    for N in range(2, 5):
        seq = np.kron(np.eye(2), extraction(N - 1)) @ extraction(N)
        direct = np.vstack([b @ a for a in annihilators(N)
                            for b in annihilators(N - 1)]) / np.sqrt(N * (N - 1))
        close(f'sequence.N{N}.coefficient', seq, direct)
        close(f'sequence.N{N}.isometry', seq.conj().T @ seq, np.eye(N + 1))
    cat = np.array([1., 0., 1.]) / np.sqrt(2)
    bell = (extraction(2) @ cat).reshape(2, 2)
    close('not_cloning.left_mixed', bell @ bell.conj().T, np.eye(2) / 2)
    close('not_cloning.joint_pure', np.linalg.norm(bell), 1.)

    # Hom-space checks use full generator constraints rather than only W_N.
    hom = []
    for N in (1, 2, 3):
        for two_j in (0, 1, 2):
            dr = two_j + 1
            m, r = (N + 1) * dr, 2 * N
            constraints = []
            for a, b, c, d0 in zip(spin_left(), spin(N - 1), spin(N), spin(two_j)):
                Jo = np.kron(a, np.eye(N)) + np.kron(np.eye(2), b)
                Ji = np.kron(c, np.eye(dr)) + np.kron(np.eye(N + 1), d0)
                constraints.append(np.kron(np.eye(m), Jo) - np.kron(Ji.T, np.eye(r)))
            nullity = m * r - np.linalg.matrix_rank(np.vstack(constraints), tol=1e-9)
            expected = 1 if two_j == 0 else 0 if two_j == 1 else (1 if N == 1 else 2)
            close(f'hom.N{N}.dr{dr}', nullity, expected)
            hom.append(dict(N=N, right_dimension=dr, complex_intertwiners=int(nullity)))

    xs = np.arange(1., 5.)
    coefficients = np.polynomial.polynomial.polyfit(xs, xs ** -.5, 3)
    close('finite_window.polynomial_interpolation',
          np.polynomial.polynomial.polyval(xs, coefficients), xs ** -.5)
    # Construct the right-only sector through the same full parent routine.
    Hvac, Vvac, PRvac, PLvac, NBvac, _, _ = parent(0, .2 * np.eye(5))
    close('N0.complete_dark_block', Hvac, np.eye(5))
    close('N0.no_conversion', Vvac)
    close('N0.no_reference_or_left', NBvac + PLvac)

    binding_paths = [HERE / 'proof.md', HERE / 'check.py',
                     STAGE / '_admission/round1323_admission.md',
                     STAGE.parent / 'archive_1313_1322/1313/proof.md',
                     STAGE.parent / 'archive_1313_1322/1319/proof.md',
                     STAGE.parent / 'archive_1201_1218/1215/proof.md']
    return dict(round=1323, status='PASS', scope='declared finite full-unknown-reference extraction contract',
                scientific_new_axioms=0, checks=checks, count=len(checks),
                max_residual=max(c['residual'] for c in checks), capacity=capacity,
                fixed_clock=clock_results, intertwiner_dimensions=hom,
                endpoint_alias_rate_ratios=alias_rates,
                finite_window_polynomial_coefficients=coefficients.tolist(),
                source_sha256={p.relative_to(ROOT).as_posix(): digest(p) for p in binding_paths},
                limitations=['not a complete local Higgs field construction',
                             'finite balance window; no continuum or all-UV claim',
                             'general flavor matrices retain failures unless old singular-value conditions hold'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check', action='store_true')
    group.add_argument('--write-new', action='store_true')
    args = parser.parse_args()
    result = run()
    target = HERE / 'results.json'
    if args.write_new:
        with target.open('x', encoding='utf8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    else:
        saved = json.loads(target.read_text(encoding='utf8'))
        same_result(result, saved)
    print(json.dumps(dict(round=1323, mode='write-new' if args.write_new else 'read-only',
                          count=result['count'], max_residual=result['max_residual'], passed=True)))


if __name__ == '__main__':
    main()
