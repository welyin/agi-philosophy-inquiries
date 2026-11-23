"""Round 1344 independent charged multiplicity check. No imports from the main analysis."""
import argparse
import json
from pathlib import Path

import numpy as np


# Siebel, WUB-DIS 2003-11, table B.7: every entry tabulated at theta3=148 degrees.
# Table B.7 quotes statistical errors. The shared two-jet control is DELPHI paper table 3.
BINS = [
    (105, .777, .048), (107, .804, .035), (110, .797, .031),
    (113, .800, .030), (117, .854, .031), (121, .762, .026),
    (126, .816, .026), (131, .759, .024), (136, .754, .023), (141, .701, .022),
]
N2, N2_STAT = .570, .001


def casimir_weights(kind):
    # Work on |q> x |qbar> x |g>. The amplitude tensor is T[a][i,j].
    matrices = []
    if kind == "SU3":
        for i, j in [(0, 1), (0, 2), (1, 2)]:
            E = np.zeros((3, 3), complex)
            E[i, j] = .5
            matrices.extend([E + E.T, 1j * (E.T - E)])
        matrices.extend([np.diag([1., -1., 0.]) / 2,
                         np.diag([1., 1., -2.]) / np.sqrt(12)])
    else:
        for j, k in [(1, 2), (2, 0), (0, 1)]:
            A = np.zeros((3, 3), complex)
            A[j, k], A[k, j] = -1j, 1j
            matrices.append(A)
    T = np.array(matrices)
    C = np.trace(sum(M @ M for M in T)).real / 3
    index = np.trace(T[0] @ T[0]).real
    state = T.transpose(1, 2, 0).copy()
    state /= np.linalg.norm(state)
    pair_color = np.zeros((3, 3))
    for M in T:
        adj = np.einsum("aij,bji->ab", T, [M @ X - X @ M for X in T]) / index
        assert np.allclose(adj, adj.conj().T)
        parts = [
            np.einsum("ij,jka->ika", M, state),
            -np.einsum("ija,jk->ika", state, M),
            np.einsum("ab,ijb->ija", adj, state),
        ]
        assert np.linalg.norm(sum(parts)) < 1e-12
        for i in range(3):
            for j in range(3):
                pair_color[i,j] += np.vdot(parts[i], parts[j]).real
    return -np.array([pair_color[1,2], pair_color[0,2], pair_color[0,1]]) / C


def cone(b, beta):
    if beta == 0:
        return b / 2
    tangentsq = np.tan(np.deg2rad(beta))**2
    return np.log1p(.5 * b * tangentsq) / np.log1p(tangentsq)


def quadrature(theta, beta):
    # Direct integration over the polar cosine and the azimuth of the soft gluon.
    nodes, weights = np.polynomial.legendre.leggauss(32)
    low = np.cos(np.deg2rad(beta))
    mu = (low + 1 + (1 - low) * nodes)[:, None] / 2
    sinpolar = np.sqrt(1 - mu**2)
    phi = (np.arange(512) + .5) * 2 * np.pi / 512
    b = 1 - np.cos(theta)
    reference = 2 / mu
    antenna = b / ((1 - sinpolar * np.cos(phi)) *
                   (1 - sinpolar * np.cos(phi - theta)))
    return np.dot(weights, np.mean(antenna, axis=1)) / np.dot(weights, reference[:,0])


def geometry(angles, color, beta):
    angles = np.deg2rad(angles)
    # Leg labels ordered by energy. Angle i is opposite leg i.
    x = 2 * np.sin(angles) / sum(np.sin(angles))
    assignment = []
    values = []
    dipoles = cone(1 - np.cos(angles), beta)
    for g in range(3):
        q, qb = [i for i in range(3) if i != g]
        assignment.append((x[q]**2 + x[qb]**2) / ((1 - x[q]) * (1 - x[qb])))
        values.append(color[0] * dipoles[q] + color[1] * dipoles[qb] + color[2] * dipoles[g])
    return float(np.average(values, weights=assignment))


def run():
    colors = {name:casimir_weights(name) for name in ("SU3", "SO3")}
    assert np.allclose(colors["SU3"], [9/8, 9/8, -1/8])
    assert np.allclose(colors["SO3"], [1/2]*3)
    tests=[]
    for degree in [0., 23., 61., 117., 151., 180.]:
        theta = np.deg2rad(degree)
        for opening in [10., 20., 30., 40.]:
            tests.append(abs(cone(1 - np.cos(theta), opening) - quadrature(theta, opening)))
    assert max(tests)<2e-12
    cap = float(1.5 * cone(1.5, 30))
    assert np.isclose(geometry([120]*3,colors["SO3"],30),cap)
    records=[]
    for theta2, n, error in BINS:
        angles = [360 - 148 - theta2, theta2, 148]
        outputs = {name: {str(beta):geometry(angles, color, beta) for beta in (0,30)}
                   for name, color in colors.items()}
        R = n / N2
        # The single-bin statistical propagation does not create a global independent-bin test.
        stat = np.hypot(error / N2, n * N2_STAT / N2**2)
        assert outputs["SO3"]["30"] <= cap + 1e-12
        records.append(dict(theta2=theta2,N3=n,N3_stat=error,N3_over_N2=R,
                            ratio_stat_only=stat, candidate_ratio=outputs,
                            min_hadron_bridge_factor_vs_SO3_LO_cap=R/cap))
    return dict(source_table="WUB-DIS 2003-11, Table B.7; every tabulated theta3=148 deg",
                N2=N2,N2_stat=N2_STAT, theta3=148,
                cone_half_angle_deg=30,jet_resolution=.015,
                normalized_anti_color_dots={n:v.tolist() for n,v in colors.items()},
                max_angular_quadrature_error=max(tests),
                SO3_30deg_max_ratio=cap,records=records)


def compare(a,b):
    if isinstance(a, dict):
        assert a.keys()==b.keys()
        for k in a:
            compare(a[k], b[k])
    elif isinstance(a, list):
        assert len(a)==len(b)
        for p,q in zip(a,b):
            compare(p,q)
    elif isinstance(a,str):
        assert a==b
    else:
        assert np.isclose(a,b,atol=2e-12,rtol=2e-12), (a,b)


if __name__ == "__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    result=run()
    out=Path(__file__).with_name("independent_results.json")
    if args.check:
        compare(result,json.loads(out.read_text(encoding="utf-8")))
    else:
        with out.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write("\n")
    print("PASS", len(result["records"]), result["SO3_30deg_max_ratio"])
