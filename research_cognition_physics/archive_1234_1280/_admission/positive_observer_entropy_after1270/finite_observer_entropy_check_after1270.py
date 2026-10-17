"""Finite positive-observer entropy diagnostic; not a gravitational model.

Only NumPy and the standard library. Default: compare saved results read-only.
--write-once creates results next to this file and refuses to overwrite.
"""
from pathlib import Path
import json
import math
import sys
import numpy as np

RESULT = Path(__file__).with_suffix('.json')


def entropy(a):
    ev = np.linalg.eigvalsh((a + a.conj().T) / 2)
    assert ev.min() > -2e-13
    pos = ev[ev > 0]
    return float(-np.sum(pos * np.log(pos)))


def sqrtm(a):
    ev, v = np.linalg.eigh(a)
    assert ev.min() > -2e-13
    return (v * np.sqrt(np.maximum(ev, 0))) @ v.conj().T


def gamma_pdf(x, k):
    out = np.zeros_like(x, dtype=float)
    mask = x > 0
    out[mask] = np.exp((k - 1) * np.log(x[mask]) - x[mask] - math.lgamma(k))
    return out


def calculate(order):
    beta = math.log(2)
    energies = np.array([0., .5, 1.5])
    z = float(np.exp(-beta * energies).sum())
    psi = np.array([1., 1j, 2.], dtype=complex) / math.sqrt(6)
    pure = np.outer(psi, psi.conj())
    diagonal = np.diag([.2, .3, .5]).astype(complex)
    mixed = .65 * pure + .35 * diagonal
    states = {'commuting': diagonal, 'pure_coherent': pure, 'mixed_coherent': mixed}
    nodes, weights = np.polynomial.legendre.leggauss(order)
    rows = []
    for k in [4, 6]:
        # Integer Gamma entropy, Euler constant in binary64 precision.
        digamma = sum(1 / n for n in range(1, k)) - 0.5772156649015329
        h_unit = k + math.lgamma(k) + (1 - k) * digamma
        for width in [1., 2., 4., 8., 16.]:
            shifts = energies / width
            # Integrate in u=y/width; split at every support threshold.
            cuts = sorted(set([0., *shifts, *np.arange(1., 81., 1.)]))
            xs = np.concatenate([(lo + hi) / 2 + (hi - lo) * nodes / 2
                                 for lo, hi in zip(cuts[:-1], cuts[1:])])
            ws = np.concatenate([(hi - lo) * weights / 2
                                 for lo, hi in zip(cuts[:-1], cuts[1:])])
            gp = np.stack([gamma_pdf(xs - shift, k) for shift in shifts], axis=1)
            g0 = gamma_pdf(xs, k)
            for name, rho in states.items():
                pi = rho.diagonal().real
                srho = entropy(rho)
                root = sqrtm(rho)
                p_y = gp @ pi
                mass = float(np.dot(ws, p_y))
                mean_y = float(np.dot(ws, xs * p_y)) * width
                h_y = float(-np.dot(ws, p_y * np.log(p_y))) + math.log(width)
                h_q = h_unit + math.log(width)
                s_fiber, avg_r, iyj, kl_shift, fiber_match = 0., 0., 0., 0., 0.
                for u, weight, g, prob in zip(xs, ws, gp, p_y):
                    if prob < 1e-280:
                        continue
                    amp = np.sqrt(g)
                    d_u = amp[:, None] * rho * amp[None, :]
                    s_fiber += weight * entropy(d_u)
                    # Independent reference purification calculation, without dephasing rho.
                    r_u = root.T @ np.diag(g) @ root.conj()
                    avg_r += weight * prob * entropy(r_u / prob)
                    fiber_match = max(fiber_match, abs(entropy(d_u / prob) - entropy(r_u / prob)))
                    mask = g > 0
                    iyj += weight * float(np.dot(pi[mask] * g[mask], np.log(g[mask] / prob)))
                    base = math.exp((k - 1) * math.log(u) - u - math.lgamma(k))
                    kl_shift += weight * float(np.dot(pi[mask] * g[mask], np.log(g[mask] / base)))
                iyr = srho - avg_r
                exact = s_fiber + math.log(width) + math.log(beta) - math.log(z) - beta * mean_y
                drel = -srho + beta * float(pi @ energies) + math.log(z)
                predicted = -drel - beta * k * width + h_q + math.log(beta)
                gap = iyj - iyr
                bound = float(pi @ energies**2) / (2 * (k - 2) * width**2)
                assert abs(mass - 1) < 2e-10
                assert abs(mean_y - (k * width + float(pi @ energies))) < 2e-9
                assert abs(h_y - h_q - iyj) < 2e-9
                assert abs(exact - predicted - gap) < 2e-9
                assert -2e-10 <= gap <= iyj + 2e-10
                assert iyj <= kl_shift + 2e-10 and kl_shift <= bound + 2e-10
                assert fiber_match < 2e-12
                if name == 'commuting':
                    assert abs(gap) < 2e-9
                if name == 'pure_coherent':
                    assert abs(iyr) < 2e-9
                rows.append(dict(k=k, width=width, state=name, observer_mean=k*width,
                    entropy_exact=exact, entropy_leading=predicted, gap=gap,
                    I_Y_energy=iyj, I_Y_reference=iyr, weighted_shift_KL=kl_shift,
                    bound=bound, mass_error=abs(mass-1),
                    identity_error=abs(exact-predicted-gap), fiber_spectrum_error=fiber_match))
    return rows


def main():
    low, high = calculate(32), calculate(64)
    keys = ['entropy_exact', 'entropy_leading', 'gap', 'I_Y_energy',
            'I_Y_reference', 'weighted_shift_KL', 'bound']
    convergence = max(abs(a[key]-b[key]) for a,b in zip(low,high) for key in keys)
    assert convergence < 2e-9, convergence
    result = dict(status='PASS_FINITE_POSITIVE_OBSERVER_ENTROPY', scope='Finite inner-action algebra only; fixed trace; no spacetime or QFT entropy claim',
        cases=len(high), quadrature_comparison_error=convergence,
        max_identity_error=max(row['identity_error'] for row in high),
        proof_not_replaced_by_numerics=True, rows=high)
    if sys.argv[1:] == ['--write-once']:
        with RESULT.open('x', encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write('\n')
        action='WRITE_ONCE'
    elif not sys.argv[1:]:
        saved=json.loads(RESULT.read_text('utf-8'))
        assert saved['status']==result['status'] and saved['cases']==result['cases']
        for a,b in zip(saved['rows'],high):
            for key in keys:
                assert abs(a[key]-b[key])<2e-9, (a['state'],key)
        action='READ_ONLY_MATCH'
    else:
        raise SystemExit('Expected no arguments or --write-once')
    print(json.dumps({key:result[key] for key in ['status','cases','quadrature_comparison_error','max_identity_error']} | {'action':action},ensure_ascii=False))


if __name__=='__main__':
    main()
