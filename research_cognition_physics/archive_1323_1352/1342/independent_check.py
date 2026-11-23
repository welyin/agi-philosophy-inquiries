"""Independent audit of round 1342 (no imports from the author implementation)."""
import json
import math
from pathlib import Path

import numpy as np


def adj(a):
    return a.conj().T


def hermitian_basis(n):
    basis = []
    for i in range(n):
        for j in range(i + 1):
            m = np.zeros((n, n), dtype=complex)
            if i == j:
                m[i, i] = 1
                basis.append(m)
            else:
                m[i, j] = 1
                basis.extend([m + adj(m), 1j * (m - adj(m))])
    return basis


def zero_ray(forms):
    """Start with the mixed barycentre; reduce rank inside its affine section."""
    d = len(forms[0])
    state = np.eye(d, dtype=complex) / d
    for rank in range(d, 1, -1):
        p, u = np.linalg.eigh(state)
        u = u[:, -rank:]
        p = p[-rank:]
        assert p.min() > 1e-12
        vars = hermitian_basis(rank)
        rows = [[np.trace(adj(u) @ f @ u @ m).real for m in vars] for f in forms]
        _, _, vh = np.linalg.svd(rows, full_matrices=True)
        delta = sum(a * m for a, m in zip(vh[-1], vars))
        assert np.linalg.norm(np.asarray(rows) @ vh[-1]) < 1e-11
        invroot = 1 / np.sqrt(p)
        k = invroot[:, None] * delta * invroot[None, :]
        bound = np.linalg.eigvalsh(k)[-1]
        assert bound > 1e-8
        state = u @ (np.diag(p) - delta / bound) @ adj(u)
    val, v = np.linalg.eigh(state)
    assert abs(val[-1] - 1) < 5e-10
    return v[:, -1]


def populations(x):
    n = x.shape[0]
    ps = [np.outer(c, c.conj()) for c in x.T]
    forms = [np.eye(n), ps[1] - ps[0], ps[2] - ps[0]]
    sub = np.eye(n, dtype=complex)
    vs = []
    for _ in range(n - 1):
        v = zero_ray([adj(sub) @ f @ sub for f in forms])
        vs.append(sub @ v)
        _, _, vh = np.linalg.svd(v.conj()[None, :], full_matrices=True)
        sub = sub @ adj(vh)[:, 1:]
    vs.append(sub[:, 0])
    return np.column_stack(vs)


def check_triple(x, p, q):
    ratio = 1 + np.sqrt(2)
    v = populations(x)
    co = adj(v) @ x
    pop_error = np.max(np.abs(np.abs(co) - np.abs(co[:, :1])))
    assert pop_error < 5e-9
    t1 = np.angle(co[:, 1] * co[:, 0].conj())
    t2 = np.angle(co[:, 2] * co[:, 0].conj())
    # Nearest rational-grid point and modular inverse. There is no search.
    cells = np.rint(q * (t2 - ratio * t1) / (2 * np.pi)).astype(np.int64) % q
    branch = (pow(p, -1, q) * cells) % q
    omega = t1 + 2 * np.pi * branch
    u1 = v @ (np.exp(1j * omega) * co[:, 0])
    u2 = v @ (np.exp(1j * ratio * omega) * co[:, 0])
    first = np.linalg.norm(u1 - x[:, 1])
    last = np.linalg.norm(u2 - x[:, 2])
    span = np.ptp(omega)
    assert first < 2e-8
    assert last < 3 * np.pi / q + 1e-8
    assert span <= 2 * np.pi * q + 1e-8
    return float(pop_error), float(first), float(last), float(span)


def run():
    rng = np.random.default_rng(22041342)
    triples = []
    for n in (1, 2, 3, 5, 8):
        for _ in range(3):
            x = rng.normal(size=(n, 3)) + 1j * rng.normal(size=(n, 3))
            triples.append(x / np.linalg.norm(x, axis=0))
    triples.extend([np.eye(2, dtype=complex)[:, [0, 1, 1]],
                    np.eye(3, dtype=complex),
                    np.eye(5, dtype=complex)[:, [1, 1, 1]]])
    p, q = 985, 408
    assert abs((1 + np.sqrt(2)) - p / q) < 1 / q ** 2
    obs = np.max([check_triple(x, p, q) for x in triples], axis=0)
    endpoint_bound = 3 * np.pi / q
    angle = lambda v, w: np.arccos(np.clip(abs(np.vdot(v, w)), 0, 1))
    worst = 0.
    for _ in range(120):
        e = rng.normal(size=4)
        a = rng.normal(size=4) + 1j * rng.normal(size=4)
        a /= np.linalg.norm(a)
        t, s = rng.uniform(-3, 3, size=2)
        b, c = np.exp(-1j * e * t) * a, np.exp(-1j * e * (t + s)) * a
        excess = abs(angle(a, b) - angle(b, c)) - np.ptp(e) * abs(t - s) / 2
        worst = max(worst, excess)
    assert worst < 1e-12
    # Common diagonal centre shifts leave ray motion unchanged and cannot reduce width.
    # An orthogonal intermediate qubit sets gap = k*pi with k positive and odd.
    odds = np.arange(1, 32, 2)
    qubit_floor = np.min(np.arccos(np.clip(abs(np.cos(np.sqrt(2) * odds * np.pi / 2)), 0, 1)))
    threshold = 2 * (math.acos(.87) - math.acos(.96))
    return dict(triples=len(triples), denominator=q, uniform_last_error=endpoint_bound,
                last_error_measured=float(obs[2]), first_error_measured=float(obs[1]),
                population_error=float(obs[0]), guaranteed_dimensionless_width=2*np.pi*q,
                max_width_used=float(obs[3]), bandwidth_checks=120, worst_violation=worst,
                finite_qubit_floor_angle_Btau_at_most_31pi=float(qubit_floor),
                data1341_threshold=threshold)


if __name__ == "__main__":
    data = run()
    path = Path(__file__).with_name("independent_results.json")
    if path.exists():
        old = json.loads(path.read_text(encoding="utf-8"))
        assert old.keys() == data.keys()
        for k, v in data.items():
            if k in ("first_error_measured", "population_error", "last_error_measured", "max_width_used"):
                continue  # SVD null choices depend on LAPACK; theorem bounds are asserted above.
            assert np.isclose(old[k], v)
    print(json.dumps(data, indent=2))

