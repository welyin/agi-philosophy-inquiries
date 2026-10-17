"""Small existing-material calibration diagnostic; default is read-only.

Exact integer/Fraction witnesses carry the two-endpoint calculation.
Finite matrices diagnose the analytic full-input identities, not experiments.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def mm(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def tr(a):
    return list(map(list, zip(*a)))


def exact_witness(z):
    assert z in (-1, 1)
    u = [[0]*4 for _ in range(4)]
    for a in range(2):
        for b in range(2):
            u[2*b+a][2*a+b] = z if a == b == 1 else 1
    identity = [[int(i == j) for j in range(4)] for i in range(4)]
    assert mm(tr(u), u) == identity
    j0 = [[int(i == 2*j) for j in range(2)] for i in range(4)]
    moved = [[int(i == j) for j in range(2)] for i in range(4)]
    assert mm(u, j0) == moved
    # Final C state for input +, with its 1/sqrt(2) normalization included.
    rho = [[Q(1, 2), Q(z, 2)], [Q(z, 2), Q(1, 2)]]
    effect = [[Q(1, 2), Q(3, 10)], [Q(3, 10), Q(1, 2)]]
    return sum(mm(effect, rho)[j][j] for j in range(2))


def unitary(z):
    u = np.zeros((4, 4), dtype=complex)
    for a in range(2):
        for b in range(2):
            u[2*b+a, 2*a+b] = z if a == b == 1 else 1
    return u


def calculate():
    for rel, digest in json.loads((HERE/'sources.json').read_text(encoding='utf-8'))['local_sha256'].items():
        assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest() == digest, rel
    pplus, pminus = exact_witness(1), exact_witness(-1)
    assert (pplus, pminus, (pplus-pminus)/2) == (Q(4, 5), Q(1, 5), Q(3, 10))
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    i2 = np.eye(2, dtype=complex)
    ep = (i2 + 0.6*x)/2
    em = i2-ep
    projplus, projminus = (i2+x)/2, (i2-x)/2
    kp = np.sqrt(0.8)*projplus + np.sqrt(0.2)*projminus
    km = np.sqrt(0.2)*projplus + np.sqrt(0.8)*projminus
    # Actual inherited complete instrument: C tensor two record factors.
    v = np.zeros((8, 2), dtype=complex)
    for out in range(2):
        v[4*out] = kp[out]
        v[4*out+3] = km[out]
    r = max(np.linalg.norm(v.conj().T@v-i2), np.linalg.norm(kp.conj().T@kp-ep),
            np.linalg.norm(km.conj().T@km-em))
    j0 = np.zeros((4, 2), dtype=complex)
    j0[0, 0], j0[2, 1] = 1, 1
    jseq = np.zeros((8, 2), dtype=complex)
    jseq[1, 0], jseq[5, 1] = 1, 1  # A unknown, B=0, C=1.
    plus = np.array([1, 1], dtype=complex)/np.sqrt(2)
    phase = np.diag([1, 1j])
    r = max(r, np.linalg.norm(phase@x@phase.conj().T-y))
    swap = unitary(1)
    n = np.diag([0, 1])
    nt = np.kron(n, i2)+np.kron(i2, n)
    rng = np.random.default_rng(1272003)
    worst_bound_slack = 1.
    phases = [0., np.pi, np.pi/2, -np.pi/2] + list(rng.uniform(-np.pi, np.pi, 20))
    for phi in phases:
        z = np.exp(1j*phi)
        u = unitary(z)
        r = max(r, np.linalg.norm(u.conj().T@u-np.eye(4)),
                np.linalg.norm(u@swap-swap@u), np.linalg.norm(u@nt-nt@u),
                np.linalg.norm(u@j0-swap@j0))
        seq = np.kron(i2, u)@np.kron(u, i2)@jseq
        expected = np.zeros((8, 2), dtype=complex)
        expected[2, 0], expected[3, 1] = 1, z
        r = max(r, np.linalg.norm(seq-expected))
        target = np.trace(np.kron(np.eye(4), ep)@np.outer(seq@plus, (seq@plus).conj())).real
        state = np.array([1, z], dtype=complex)/np.sqrt(2)
        px = np.vdot(state, ep@state).real
        py = np.vdot(state, (i2+0.6*y)/2@state).real
        r = max(r, abs(target-px), abs(px-(0.5+0.3*z.real)), abs(py-(0.5+0.3*z.imag)),
                abs((10/3)*complex(px-0.5, py-0.5)-z))
        # Reference test is additional diagnosis. The isometry identity above is universal.
        psi = rng.normal(size=6)+1j*rng.normal(size=6)
        psi /= np.linalg.norm(psi)
        out = np.kron(seq, np.eye(3))@psi
        exact = np.kron(expected, np.eye(3))@psi
        r = max(r, np.linalg.norm(out-exact), abs(np.vdot(out, out)-1))
        eta = 0.0001
        delta = 10*np.sqrt(2)*eta/3
        for sx in (-1, 1):
            for sy in (-1, 1):
                w = (10/3)*complex(px+sx*eta-0.5, py+sy*eta-0.5)
                zhat = w/abs(w)
                uh = unitary(zhat)
                full = np.kron(i2, u)@np.kron(u, i2)
                fullhat = np.kron(i2, uh)@np.kron(uh, i2)
                difference = np.linalg.norm(fullhat-full, ord=2)
                assert abs(zhat-z) <= 2*delta+1e-13
                assert difference <= 4*delta+1e-13
                worst_bound_slack = min(worst_bound_slack, 4*delta-difference)
    assert r < 1e-12
    eta, xi = Q(1, 10000), Q(1, 1000)
    budget = 40*np.sqrt(2)*float(eta)/3+float(xi)
    assert eta+xi == Q(11, 10000) and budget < 0.003
    return {'status': 'PASS_SCOPED_CALIBRATION_TRANSFER',
        'exact': {'z_plus_target': str(pplus), 'z_minus_target': str(pminus),
                  'target_gap': str(pplus-pminus), 'minimax_error_lower_bound': str((pplus-pminus)/2),
                  'target_only_budget': str(eta+xi)},
        'numerical_diagnostics': {'phase_points': len(phases), 'max_identity_residual': float(r),
            'full_process_budget': float(budget), 'smallest_sample_bound_slack': float(worst_bound_slack)},
        'scope': {'same_B_no_reset': True, 'same_pair_rule': True, 'all_reference_identity_analytic': True,
            'independent_geometry_test': False, 'physical_experiment': False, 'new_scientific_rounds': 0,
            'source_and_clock_certificate': 'Inherited 1112/1115 finite logical boundary permissions; not newly simulated.'}}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--write', action='store_true')
    args = p.parse_args()
    got = calculate()
    target = HERE/'results.json'
    if args.write:
        with target.open('x', encoding='utf-8') as f:
            json.dump(got, f, ensure_ascii=False, indent=2)
            f.write('\n')
    else:
        old = json.loads(target.read_text(encoding='utf-8'))
        assert got['status'] == old['status'] and got['exact'] == old['exact'] and got['scope'] == old['scope']
        assert got['numerical_diagnostics']['phase_points'] == old['numerical_diagnostics']['phase_points']
        assert abs(got['numerical_diagnostics']['full_process_budget']-old['numerical_diagnostics']['full_process_budget']) < 1e-14
    print(json.dumps(got, ensure_ascii=False))


if __name__ == '__main__':
    main()
