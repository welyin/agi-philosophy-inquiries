"""Finite graph, sharp dephasing bound and matched spurion check; NumPy only."""
import argparse
import json
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT / "inputs.json").read_text(encoding="utf-8-sig"))


def mixing(s12, s13, s23, delta):
    c12, c13, c23 = np.sqrt(1.0 - np.square([s12, s13, s23]))
    z = np.exp(1j * delta)
    return np.array([
        [c12*c13, s12*c13, s13/z],
        [-s12*c23-c12*s23*s13*z, c12*c23-s12*s23*s13*z, s23*c13],
        [s12*s23-c12*c23*s13*z, -c12*s23-s12*c23*s13*z, c23*c13],
    ])


def gap(unitary):
    return 1.0 - np.linalg.svd(abs(unitary)**2, compute_uv=False)[1]


def raw_graph():
    box = DATA["edge_work_box"]
    lower = {label: (center - box["widths"]*sum(errors))/(1+box["amplitude_transfer"])
             for label, center, errors in DATA["edges"]}
    vertices = "uctdsb"
    dist = np.full((6, 6), np.inf)
    np.fill_diagonal(dist, 0)
    for edge, x in lower.items():
        i, j = map(vertices.index, edge)
        dist[i, j] = dist[j, i] = 1.0/x**2
    for k in range(6):
        dist = np.minimum(dist, dist[:, k, None]+dist[None, k, :])
    i, j = np.unravel_index(np.argmax(dist), dist.shape)
    diameter = np.sqrt(dist[i, j])
    norm = box["norm_ceiling_example"]
    return {"min_amplitudes": lower, "worst_pair": vertices[i]+vertices[j],
            "sqrt_resistance_diameter": diameter,
            "required_Esum_for_diagonal_spread": box["desired_strength_spread"]/(diameter*norm)}


def entry_lowers(low, high):
    """Triangle bounds for ALL phase values; rounded down with arithmetic slack."""
    a, b, c = low
    A, B, C = high
    x, y, z = np.sqrt(1-np.square(high))
    X, Y, Z = np.sqrt(1-np.square(low))

    def diff(p0, p1, q0, q1):
        return max(0., p0-q1, q0-p1)

    return np.array([
        [x*y, a*y, b],
        [diff(a*z, A*Z, x*c*b, X*C*B), diff(x*z, X*Z, a*c*b, A*C*B), c*y],
        [diff(a*c, A*C, x*z*b, X*Z*B), diff(x*c, X*C, a*z*b, A*Z*B), z*y],
    ])


def positive_definite(matrix):
    """Exact rational Schur complement certificate."""
    active = matrix
    while active:
        pivot = active[0][0]
        if pivot <= 0:
            return False
        active = [[active[i][j]-active[i][0]*active[0][j]/pivot
                   for j in range(1, len(active))] for i in range(1, len(active))]
    return True


def certify_box(low, high):
    floor = np.maximum(0, np.floor(10**7 * (entry_lowers(low, high)**2-1e-12))).astype(int)
    lap = [[Fraction(1, 6) for _ in range(6)] for _ in range(6)]
    for i in range(3):
        for j in range(3):
            a, b = i, j+3
            w = Fraction(int(floor[i, j]), 10**7)
            lap[a][a] += w
            lap[b][b] += w
            lap[a][b] -= w
            lap[b][a] -= w
    lo, hi = 0, 10**7
    while hi-lo > 1:
        mid = (lo+hi)//2
        possible = [[lap[i][j]-(Fraction(mid, 10**7) if i == j else 0)
                     for j in range(6)] for i in range(6)]
        if positive_definite(possible):
            lo = mid
        else:
            hi = mid
    return lo/10**7


def angular_checks():
    angles = DATA["angles"]
    mult = angles["widths"]
    q = np.array(angles["quark"])
    u = mixing(*q[:, 0])
    out = {"quark": {"central_gap": gap(u), "box_gap_at_least": certify_box(
        q[:3, 0]-mult*q[:3, 1], q[:3, 0]+mult*q[:3, 1])}}
    def reactor(t):
        return np.sqrt((1-np.sqrt(1-t))/2)
    re, re_err = angles["reactor_amplitude"]
    for name in ("NO", "IO"):
        solar, se = angles[name]["solar"]
        atm, *ae = angles[name]["atmospheric"]
        local = [np.sqrt(solar), reactor(re), np.sqrt(atm)]
        low = [np.sqrt(solar-mult*se), reactor(re-mult*re_err), np.sqrt(atm-mult*max(ae))]
        high = [np.sqrt(solar+mult*se), reactor(re+mult*re_err), np.sqrt(atm+mult*max(ae))]
        ext = [gap(mixing(*local, delta)) for delta in [0., np.pi]]
        out[name] = {"central_min_gap_over_phase": min(ext), "endpoint_gaps": ext,
                     "box_gap_at_least": certify_box(low, high)}
    eig = np.linalg.svd(abs(u)**2)
    diag1, diag2 = np.diag(eig[0][:, 1]), np.diag(eig[2][1])
    g = diag1+u@diag2@u.conj().T
    rotated = u.conj().T@g@u
    loss = sum(np.linalg.norm(h-np.diag(np.diag(h)))**2 for h in [g, rotated])
    ratio = loss/np.linalg.norm(g-np.trace(g)*np.eye(3)/3)**2
    assert abs(ratio-gap(u)) < 1e-12
    return out, u


def example(unitary):
    y = np.array(DATA["up_yukawa"])
    coeff = DATA["epsilon"]
    a = np.diag(1+coeff*y**2)@unitary
    up, down = a@a.conj().T, a.conj().T@a
    spread = np.ptp(np.linalg.eigvalsh(up))
    sb = abs(down[1, 2])
    exact = np.diag(up).real
    sb_check = (unitary[2, 1].conjugate()*unitary[2, 2]+(exact[1]-exact[0])/
                (exact[2]-exact[0])*unitary[1, 1].conjugate()*unitary[1, 2])
    assert abs(sb/spread-abs(sb_check)) < 1e-12
    return {"max_minus_min_strength": spread, "sb_return_abs": sb,
            "sb_abs_over_spread": sb/spread,
            "up_offdiag_norm": np.linalg.norm(up-np.diag(np.diag(up))),
            "down_offdiag_norm": np.linalg.norm(down-np.diag(np.diag(down)))}


def run():
    two_bases, unitary = angular_checks()
    return {"direct_edges": raw_graph(), "two_bases": two_bases, "spurion": example(unitary)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = json.loads(json.dumps(run()))
    path = ROOT/"results.json"
    if args.write_new:
        with path.open("x", encoding="utf-8") as outfile:
            json.dump(result, outfile, ensure_ascii=False, indent=2)
            outfile.write("\n")
    elif args.check:
        saved = json.loads(path.read_text(encoding="utf-8-sig"))
        def same(left, right):
            if isinstance(left, dict):
                return left.keys() == right.keys() and all(same(left[k], right[k]) for k in left)
            if isinstance(left, list):
                return len(left) == len(right) and all(same(a,b) for a,b in zip(left,right))
            if isinstance(left, (float, int)):
                return np.isclose(left, right, rtol=2e-10, atol=3e-14)
            return left == right
        assert same(result, saved)
        print("1476 verified")
    else:
        print(json.dumps(result, indent=2))
