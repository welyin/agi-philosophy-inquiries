"""Finite source-record/recovery diagnostics. Default is read-only."""
from pathlib import Path
import argparse
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent


def trace_norm(a):
    return float(np.linalg.svd(a, compute_uv=False).sum())


def choi(kraus):
    return sum(np.outer(k.reshape(-1, order='F'),
                        k.reshape(-1, order='F').conj()) for k in kraus) / 2


def isometry(rng, rows, cols):
    z = rng.normal(size=(rows, cols)) + 1j*rng.normal(size=(rows, cols))
    return np.linalg.qr(z)[0]


def run():
    bell = np.array([1., 0., 0., 1.], dtype=complex) / math.sqrt(2)
    phi = np.outer(bell, bell.conj())
    z = np.diag([1., -1.])
    cases = []
    a = 1.
    for s in [0., .1, .5, 1., 2., 10.]:
        r = math.hypot(a, s)
        y = np.array([r, -r])
        p = np.array([(1+a/r)/2, (1-a/r)/2])
        q = p[::-1]
        ks = [np.diag(np.sqrt([p[j], q[j]])) for j in range(2)]
        b = float(np.sqrt(p*q).sum())
        bound = (1-s/r)/2
        c = choi(ks)
        fe = float((bell.conj()@c@bell).real)
        dist = trace_norm(c-phi)/2
        assert abs(y@p-a) < 1e-12 and abs(y@q+a) < 1e-12
        assert abs(((y-a)**2)@p-s*s) < 1e-12
        assert abs(((y+a)**2)@q-s*s) < 1e-12
        assert np.linalg.norm(sum(k.conj().T@k for k in ks)-np.eye(2)) < 1e-12
        assert abs(b-s/r) < 1e-12 and abs(fe-(1+b)/2) < 1e-12
        assert abs(dist-bound) < 1e-12
        # Channel identity certifies the analytic half-diamond value; no SDP claim.
        for i in range(2):
            for j in range(2):
                e = np.zeros((2, 2)); e[i,j] = 1
                image = sum(k@e@k.conj().T for k in ks)
                target = (1+b)/2*e+(1-b)/2*z@e@z
                assert np.linalg.norm(image-target) < 1e-12
        cases.append({'s_over_a': s, 'B': b, 'Bell_infidelity': 1-fe,
                      'Bell_half_trace_distance': dist,
                      'analytic_half_diamond': bound})

    # Non-QND instruments: 4 outcomes, 2 hidden Kraus labels, qout dimension 3.
    # Each outcome gets an arbitrary CPTP recovery 3 -> 2 via environment dim 3.
    rng = np.random.default_rng(105910)
    worst_moment_excess = worst_cp_excess = worst_fe_excess = 0.
    for _ in range(64):
        v = isometry(rng, 24, 2).reshape(4, 2, 3, 2)
        vals = rng.normal(size=4)
        probs = np.array([[sum(np.vdot(k[:,i], k[:,i]).real for k in branch)
                           for i in range(2)] for branch in v])
        assert np.max(abs(probs.sum(axis=0)-1)) < 1e-12
        means = vals@probs
        std = np.sqrt(((vals[:,None]-means)**2*probs).sum(axis=0))
        delta = abs(means[0]-means[1])
        B = float(np.sqrt(probs[:,0]*probs[:,1]).sum())
        cap = float(std.sum()/math.hypot(delta, std.sum()))
        worst_moment_excess = max(worst_moment_excess, B-cap)
        all_k = []
        for branch, ps in zip(v, probs):
            rec = isometry(rng, 6, 3).reshape(3, 2, 3)
            composed = [r@k for r in rec for k in branch]
            cross = sum(np.outer(k[:,0], k[:,1].conj()) for k in composed)
            worst_cp_excess = max(worst_cp_excess,
                                 trace_norm(cross)-math.sqrt(ps[0]*ps[1]))
            all_k.extend(composed)
        c = choi(all_k)
        fe = float((bell.conj()@c@bell).real)
        worst_fe_excess = max(worst_fe_excess, fe-(1+B)/2)
        assert np.linalg.eigvalsh(c).min() > -1e-12
        assert abs(np.trace(c)-1) < 1e-12
        assert trace_norm(c-phi)/2 >= 1-fe-1e-12
    assert max(worst_moment_excess, worst_cp_excess, worst_fe_excess) < 1e-12

    eta = 4156812479/38839500000
    d = 1-1/math.sqrt(5)
    source = d-eta
    assert source > 0
    return {'passed': True, 'new_science_groups': 0, 'new_empirical_groups': 0,
            'saturating_abstract_instruments': cases,
            'non_QND_recovery_diagnostic': {'cases': 64,
                'maximum_positive_moment_excess': worst_moment_excess,
                'maximum_positive_CP_block_excess': worst_cp_excess,
                'maximum_positive_Bell_fidelity_excess': worst_fe_excess},
            'round1053_transport': {'d': d, 'eta_upper': eta,
                'half_mean_separation_lower': source,
                'mean_separation_lower': 2*source,
                'illustrative_s_equals_a_recovery_lower': (1-1/math.sqrt(2))/2},
            'actual_gravitational_instrument_certified': False,
            'universal_statement_justified_by': 'analytic proof, not sampling'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--save', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'results.json'
    if args.save:
        with path.open('x', encoding='utf-8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write('\n')
    else:
        saved = json.loads(path.read_text(encoding='utf-8'))
        assert saved == result
    print(json.dumps(result, ensure_ascii=False))
