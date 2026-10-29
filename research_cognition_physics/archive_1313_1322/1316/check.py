"""1316 full Dirac chains, shared response and exact light-pole bounds.

Only NumPy and the standard library. No empirical fit or imported author model.
Default/--check is read-only; --write-new creates results exclusively.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
TOL = 3e-8


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def unitary(n, rng):
    return np.linalg.qr(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))[0]


def invsqrt(a):
    s, v = np.linalg.eigh(a)
    assert min(s) > 0
    return (v / np.sqrt(s)) @ v.conj().T


def geneig(a, b):
    w = invsqrt(b)
    return np.linalg.eigvalsh(w @ a @ w)


def chain(lengths, q, m):
    n, count = len(lengths), sum(lengths)
    t = np.zeros((count, count + n))
    dtq = np.zeros_like(t)
    p = np.zeros((n, count + n))
    v0 = np.zeros((count + n, n))
    s, formula_s, f, heavy = [], [], [], []
    row = col = 0
    for a, size in enumerate(lengths):
        b = np.zeros((size, size + 1))
        for j in range(size):
            b[j, j], b[j, j + 1] = 1, -q
            dtq[row + j, col + j + 1] = -m
        t[row:row + size, col:col + size + 1] = m * b
        p[a, col + size] = 1
        zero = q ** (-np.arange(size + 1, dtype=float))
        zero /= np.linalg.norm(zero)
        v0[col:col + size + 1, a] = zero
        f.append(zero[-1])
        k = b @ b.T
        x = np.linalg.solve(k, np.eye(size)[:, -1])
        s.append(q * q * float(x @ x) / m ** 2)
        if q == 1:
            formula_s.append(size * (2 * size + 1) / (6 * m ** 2 * (size + 1)))
        else:
            numerator = sum((q ** j - q ** (-j)) ** 2 for j in range(1, size + 1))
            formula_s.append(numerator / (q ** (size + 1) - q ** (-size - 1)) ** 2 / m ** 2)
        heavy += list(m ** 2 * (1 + q * q - 2 * q * np.cos(np.arange(1, size + 1) * np.pi / (size + 1))))
        row += size
        col += size + 1
    return dict(t=t, p=p, v0=v0, f=np.array(f), s=np.array(s), formula_s=np.array(formula_s),
                g=t @ t.T, heavy=np.sort(heavy), dtq=dtq)


def response(x, y0, z=0., portal=None):
    zmap = y0 @ x['p'] if portal is None else portal
    g = x['g']
    return zmap @ x['t'].T @ np.linalg.solve(g, np.linalg.solve(g - z * np.eye(len(g)), x['t'])) @ zmap.conj().T


def full_mass(x, y0, rho, portal=None):
    zmap = y0 @ x['p'] if portal is None else portal
    return np.vstack([rho * zmap, x['t']])


def run():
    rng = np.random.default_rng(1316)
    errors, cases, sources, deletion = {}, [], [], []

    def eq(name, a, b, scale=None, tol=TOL):
        scale = max(1., float(np.linalg.norm(b))) if scale is None else scale
        error = float(np.linalg.norm(np.asarray(a) - np.asarray(b))) / scale
        errors[name] = max(errors.get(name, 0.), error)
        assert error < tol, (name, error)

    def psd(name, a, scale=1.):
        violation = max(0., -float(np.linalg.eigvalsh((a + a.conj().T) / 2).min())) / max(1., scale)
        errors[name] = max(errors.get(name, 0.), violation)
        assert violation < TOL, (name, violation)

    for n in (3, 4, 5):
        lengths = list(range(1, n + 1))
        for q in (1., 1.25, 1.75):
            m = 3.
            x = chain(lengths, q, m)
            y0 = unitary(n, rng) @ np.diag(np.linspace(.3, .8, n)) @ unitary(n, rng)
            y = y0 @ np.diag(x['f'])
            d = y @ y.conj().T
            c0 = response(x, y0)
            gap = float(np.linalg.eigvalsh(x['g']).min())
            eq('null_mode', x['t'] @ x['v0'], np.zeros((sum(lengths), n)))
            eq('zero_orthonormal', x['v0'].T @ x['v0'], np.eye(n))
            eq('endpoint_zero', x['p'] @ x['v0'], np.diag(x['f']))
            eq('heavy_spectrum', np.linalg.eigvalsh(x['g']), x['heavy'])
            eq('closed_response', x['s'], x['formula_s'])
            eq('matched_C6', c0, y0 @ np.diag(x['s']) @ y0.conj().T)
            predicted = np.sort(x['s'] / x['f'] ** 2)
            eq('portal_free_generalized_spectrum', geneig(c0, d), predicted)
            # Change all portal singular values and both mixing dictionaries.
            y02 = unitary(n, rng) @ np.diag(np.linspace(.6, 1.9, n)) @ unitary(n, rng)
            y2 = y02 @ np.diag(x['f'])
            eq('independent_portal_same_spectrum', geneig(response(x, y02), y2 @ y2.conj().T), predicted)
            for zeta in (.08, .65):
                rho = np.sqrt(zeta * gap / np.linalg.norm(d, 2))
                mat = full_mass(x, y0, rho)
                hright = mat.conj().T @ mat
                allpoles = np.linalg.eigvalsh(hright)
                poles = allpoles[:n]
                assert np.max(poles) < gap and allpoles[n] >= gap - TOL
                assert np.max(poles) <= zeta * gap + TOL
                low = geneig(rho ** 2 * d, np.eye(n) + rho ** 2 * c0 / (1 - zeta))
                high = geneig(rho ** 2 * d, np.eye(n) + rho ** 2 * c0)
                psd('physical_pole_lower_bound', np.diag(poles - low), np.linalg.norm(hright))
                psd('physical_pole_upper_bound', np.diag(high - poles), np.linalg.norm(hright))
                for pole in poles:
                    cz = response(x, y0, float(pole))
                    effective = rho ** 2 * d - pole * (np.eye(n) + rho ** 2 * cz)
                    eq('exact_light_pole', np.min(np.abs(np.linalg.eigvalsh(effective))), 0.,
                       scale=max(1., rho ** 2 * np.linalg.norm(d)))
                    psd('response_lower_bound', cz - c0)
                    psd('response_upper_bound', c0 / (1 - zeta) - cz)
                z = zeta * gap * .43
                hleft = mat @ mat.conj().T
                schur = hleft[:n, :n] - z * np.eye(n) - hleft[:n, n:] @ np.linalg.solve(hleft[n:, n:] - z * np.eye(len(x['g'])), hleft[n:, :n])
                eq('full_Schur_identity', schur, rho ** 2 * d - z * (np.eye(n) + rho ** 2 * response(x, y0, z)))
                cases.append(dict(n=n, lengths=lengths, q=q, m=m, rho=float(rho), zeta=zeta,
                                  gap_squared=gap, light_poles=poles.tolist(), lower=low.tolist(), upper=high.tolist(),
                                  predicted_generalized_spectrum=predicted.tolist(),
                                  low_D_condition=float(np.linalg.cond(d))))
            x2 = chain(lengths, q, 2 * m)
            eq('independent_scale_same_Y', y0 @ np.diag(x2['f']), y)
            eq('independent_scale_C6', response(x2, y0), c0 / 4)
            # A genuinely additional endpoint coupling orthogonal to the light code.
            k = np.zeros((n, len(x['g'])), complex)
            k[0, -1] = .4 / m
            portal = y0 @ x['p'] + k @ x['t']
            eq('multiport_same_Y', portal @ x['v0'], y)
            newc = response(x, y0, portal=portal)
            diff = float(np.linalg.norm(geneig(newc, d) - predicted))
            assert diff > 1e-3
            deletion.append(dict(n=n, q=q, same_Y=True, generalized_spectrum_change=diff))

    # Full-source differentiation, before eliminating any parameter or heavy field.
    n, lengths, q, m, rho = 3, [1, 2, 3], 1.4, 2.7, .4
    y0 = unitary(n, rng) @ np.diag([.31, .53, .77]) @ unitary(n, rng)
    x = chain(lengths, q, m)
    mat = full_mass(x, y0, rho)
    h = mat.conj().T @ mat
    eig, vec = np.linalg.eigh(h)
    step = 1e-5
    for name in ('m', 'q', 'rho', 'ReY01', 'ImY01'):
        if name == 'm':
            plus = full_mass(chain(lengths, q, m + step), y0, rho)
            minus = full_mass(chain(lengths, q, m - step), y0, rho)
            deriv = np.vstack([np.zeros((n, mat.shape[1])), x['t'] / m])
        elif name == 'q':
            plus = full_mass(chain(lengths, q + step, m), y0, rho)
            minus = full_mass(chain(lengths, q - step, m), y0, rho)
            deriv = np.vstack([np.zeros((n, mat.shape[1])), x['dtq']])
        elif name == 'rho':
            plus, minus = full_mass(x, y0, rho + step), full_mass(x, y0, rho - step)
            deriv = np.vstack([y0 @ x['p'], np.zeros_like(x['t'])])
        else:
            dy = np.zeros_like(y0)
            dy[0, 1] = 1 if name == 'ReY01' else 1j
            plus, minus = full_mass(x, y0 + step * dy, rho), full_mass(x, y0 - step * dy, rho)
            deriv = np.vstack([rho * dy @ x['p'], np.zeros_like(x['t'])])
        eq('mass_source_' + name, (plus - minus) / (2 * step), deriv)
        dh = deriv.conj().T @ mat + mat.conj().T @ deriv
        dhnum = (plus.conj().T @ plus - minus.conj().T @ minus) / (2 * step)
        eq('complete_squared_source_' + name, dhnum, dh)
        analytic = np.real(np.diag(vec.conj().T @ dh @ vec))
        numeric = (np.linalg.eigvalsh(plus.conj().T @ plus) - np.linalg.eigvalsh(minus.conj().T @ minus)) / (2 * step)
        eq('full_pole_HF_' + name, numeric, analytic, tol=3e-7)
        sources.append(name)

    # Same nontrivial light target can be admitted for every finite chain count.
    # This is an explicitly calibrated construction, NOT an empirical prediction.
    calibrated = []
    for n in (3, 4, 5):
        x = chain(list(range(1, n + 1)), 1.3, 5.)
        u = np.exp(2j * np.pi * np.outer(np.arange(n), np.arange(n)) / n) / np.sqrt(n)
        target = u @ np.diag(np.arange(1, n + 1) * .01)
        y0 = target @ np.diag(1 / x['f'])
        eq('arbitrary_target_Y', y0 @ np.diag(x['f']), target)
        inv = float(np.imag(u[0, 0] * u[1, 1] * u[0, 1].conjugate() * u[1, 0].conjugate()))
        assert abs(inv) > 1e-3
        calibrated.append(dict(n=n, portal_norm=float(np.linalg.norm(y0, 2)),
                               mixing_quartet=inv, phase_instrument_certified=False))

    rcharges = dict(original_left=1, Higgs=2, N_conjugate=1, added_S_left=3)
    assert (rcharges['original_left'] + rcharges['Higgs'] + rcharges['N_conjugate']) % 4 == 0
    assert (rcharges['N_conjugate'] + rcharges['added_S_left']) % 4 == 0
    assert 2 * rcharges['N_conjugate'] % 4 == 2 * rcharges['added_S_left'] % 4 == 2
    assert sum((rcharges['N_conjugate'], rcharges['added_S_left'])) % 4 == 0
    return {'schema': 'round1316_common_clockwork_response_v1', 'status': 'PASS', 'seed': 1316,
            'cases': cases, 'source_parameters': sources, 'multiport_deletions': deletion,
            'calibrated_light_targets': calibrated, 'r_charges_mod4_all_left': rcharges,
            'max_normalized_errors': errors, 'ordinary_tolerance': TOL, 'HF_tolerance': 3e-7,
            'new_axioms': 0, 'stage_complete': False,
            'scope': 'tree-level declared finite-chain parent and exact quadratic mass poles; no empirical parameter fit',
            'bindings': {str(p.relative_to(HERE.parent)): sha(p) for p in
                         (HERE/'proof.md', HERE/'check.py', HERE.parent/'_admission/round1316_admission.md')}}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--check', action='store_true')
    p.add_argument('--write-new', action='store_true')
    args = p.parse_args()
    result = run()
    out = HERE/'results.json'
    if args.write_new:
        with out.open('x', encoding='utf8', newline='\n') as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write('\n')
    else:
        saved = json.loads(out.read_text(encoding='utf8'))
        assert saved['status'] == result['status'] == 'PASS'
        assert saved['bindings'] == result['bindings']
        assert len(saved['cases']) == len(result['cases']) == 18
        for old, new in zip(saved['cases'], result['cases']):
            for key in ('n', 'lengths', 'q', 'm', 'zeta'):
                assert old[key] == new[key]
            for key in ('light_poles', 'lower', 'upper', 'predicted_generalized_spectrum'):
                assert np.allclose(old[key], new[key], atol=1e-9, rtol=1e-8)
        assert saved['source_parameters'] == result['source_parameters']
    print(json.dumps({'status': 'PASS', 'cases': len(result['cases']),
                      'source_parameters': len(result['source_parameters']),
                      'multiport_deletions': len(result['multiport_deletions']),
                      'maximum_normalized_error': max(result['max_normalized_errors'].values()),
                      'mode': 'write-new' if args.write_new else 'read-only'}))


if __name__ == '__main__':
    main()
