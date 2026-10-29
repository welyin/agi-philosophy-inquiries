"""Count-zero finite-source feasibility probe; stdout only, no saved results.

This is not a round-1326 implementation or a physical Higgs identification.
The finite boson boundary is part of this probe's exact model.
"""
import json
import numpy as np


def build(B):
    states = [(a, b) for N in range(B + 1) for a in range(N + 1) for b in [N - a]]
    index = {x: k for k, x in enumerate(states)}
    d = len(states)
    a0 = np.zeros((d, d))
    a1 = np.zeros_like(a0)
    for j, (a, b) in enumerate(states):
        if a: a0[index[a - 1, b], j] = np.sqrt(a)
        if b: a1[index[a, b - 1], j] = np.sqrt(b)
    number = np.array([sum(x) for x in states], float)
    return states, (a0, a1), number


def interaction(ops, number, gd, gu, normalized):
    d = len(number)
    V = np.zeros((4*d, 4*d))
    down = np.zeros(d)
    down[number > 0] = 1 / np.sqrt(number[number > 0])
    up = 1 / np.sqrt(number + 2)
    if not normalized: down = np.ones(d); up = np.ones(d)
    for loc, a in [(1, ops[0]), (2, ops[1])]:
        V[loc*d:(loc+1)*d, :d] = gd * a * down
    V[d:2*d, 3*d:] = gu * ops[1].T * up
    V[2*d:3*d, 3*d:] = -gu * ops[0].T * up
    return V + V.T


def main():
    B = 18
    _, ops, number = build(B)
    d = len(number)
    p = .005
    D = np.diag([p, -p, -p, p])
    source = -(ops[0] + ops[0].T)
    family = {}
    # One calibration of leading low blocks, then fixed coefficients.
    R0 = np.diag(number) + source
    e0, u0 = np.linalg.eigh(R0)
    chi0 = u0[:, 0]
    E0 = np.kron(np.eye(4), chi0[:, None])
    Vnorm = interaction(ops, number, .013, .011, True)
    norm0 = E0.T @ Vnorm @ E0
    bareunit = interaction(ops, number, 1., 1., False)
    bare0 = E0.T @ bareunit @ E0
    gd = float(norm0[1, 0] / bare0[1, 0])
    gu = float(norm0[2, 3] / bare0[2, 3])
    family['normalized'] = (Vnorm, .013)
    family['bare_linear'] = (interaction(ops, number, gd, gu, False), np.sqrt(B+2)*max(gd, gu))
    rows = []
    for f in [1., 2.]:
        R = np.diag(number) + f*source
        er, ur = np.linalg.eigh(R)
        chi = ur[:, 0]
        E = np.kron(np.eye(4), chi[:, None])
        for label, (V, g_bound) in family.items():
            H = np.kron(np.eye(4), R) + np.kron(D, np.eye(d)) + V
            eig, vectors = np.linalg.eigh(H)
            A = D + E.T @ V @ E
            approx = np.linalg.eigvalsh(A)
            full = eig[:4] - er[0]
            gap = er[1] - er[0]
            denominator = gap - 2*(p+g_bound)
            bound = g_bound*g_bound / denominator
            source_full = np.kron(np.eye(4), source)
            reference_source = float(chi @ source @ chi)
            forces = np.array([vectors[:, j] @ source_full @ vectors[:, j] - reference_source for j in range(4)])
            assert denominator > 0
            assert np.all(full <= approx + 1e-10)
            assert np.all(full >= approx - bound - 1e-10)
            assert eig[4] - eig[3] >= denominator - 1e-10
            rows.append(dict(menu=label, f=f, reference_ground=er[0], reference_gap=gap,
                             exact_low_energies=full.tolist(), first_order_low_energies=approx.tolist(),
                             energy_bound=bound, max_energy_residual=float(max(abs(full-approx))),
                             exact_excitation_source=forces.tolist(),
                             mean_reference=float(chi @ np.diag(number) @ chi)))
    print(json.dumps(dict(status='COUNT_ZERO_FEASIBILITY_ONLY', boson_cutoff=B, dimension=4*d,
                          p=p, normalized_couplings=[.013,.011], bare_once_calibrated=[gd,gu],
                          calibration='first low block at f=1; fixed for f=2',
                          Higgs_identity_claimed=False, rows=rows), indent=2))


if __name__ == '__main__':
    main()
