"""Dimension-independent variational bounds specialized to three observed families."""
import argparse
import json
from math import sin, sqrt
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
LOW = json.loads((HERE/"../1386/inputs.json").read_text(encoding="utf-8"))
HIGH = {"up": [[2.87e-6, .07e-6], [.00145, .00003], [.4454, .0048]],
        "down": [[.65e-5, .01e-5], [1.29e-4, .02e-4], [.00606, .00006]],
        "angles": [[.22708, .00082], [.00419, .00009], [.04743, .00047], [1.139, .023]]}


def ceiling_loss(up, down, angles, angular, width):
    """Upper bound on the second-smallest visible eigenvalue for any allowed weights."""
    c_over_t = (up[1][0]+width*up[1][1])/(up[2][0]-width*up[2][1])
    s_over_b = (down[1][0]+width*down[1][1])/(down[2][0]-width*down[2][1])
    s13, s23 = (angles[i][0]+width*angles[i][1] for i in (1, 2))
    if angular:
        s13, s23 = sin(s13), sin(s23)
    miss_tb = 1-(1-s13*s13)*(1-s23*s23)
    return 1-c_over_t**2, 1-s_over_b**2, miss_tb


def matrix(s, delta):
    x, z, y = s
    cx, cz, cy = np.sqrt(1-np.array([x, z, y])**2)
    phase = np.exp(1j*delta)
    return np.array([[cx*cz, x*cz, z/phase],
                     [-x*cy-cx*y*z*phase, cx*cy-x*y*z*phase, y*cz],
                     [x*y-cx*cy*z*phase, -cx*y-x*cy*z*phase, cy*cz]])


def compute():
    results = {}
    for label, up, down, angles, angular in (
            ("MZ", LOW["yukawa"]["up"], LOW["yukawa"]["down"],
             LOW["quark_sines_and_phase"], False),
            ("high", HIGH["up"], HIGH["down"], HIGH["angles"], True)):
        bound_rows = {}
        for width in (0, 1, 2):
            good_u, good_d, miss = ceiling_loss(up, down, angles, angular, width)
            d1 = good_u+good_d-1
            d2 = d1-good_d*miss
            bound_rows[str(width)] = {"largest_remainder_at_least": d1,
                                     "second_remainder_at_least": d2}
        s = [sin(a[0]) if angular else a[0] for a in angles[:3]]
        v = matrix(s, angles[3][0])
        assert np.allclose(v.conj().T@v, np.eye(3))
        hu = np.diag(np.square(np.array(up)[:,0]))
        hd = v @ np.diag(np.square(np.array(down)[:,0])) @ v.conj().T
        # Explicit balanced rank-two completion, after rescaling to the spectral edge.
        x, y = 1/hu[-1,-1], 1/down[2][0]**2
        a = x*hu+y*hd
        radius = np.linalg.eigvalsh(a)[-1]
        d = np.eye(3)-a/radius
        eig, rot = np.linalg.eigh(d)
        assert abs(eig[0]) < 1e-13 and eig[1] >= bound_rows["0"]["second_remainder_at_least"]
        delta = (eig[1]-eig[0])*(eig[2]-eig[1])*(eig[2]-eig[0])
        j = abs((rot[0,0]*rot[1,1]*rot[0,1].conj()*rot[1,0].conj()).imag)
        basis = np.eye(3, dtype=complex)
        projectors = [np.outer(col, col.conj()) for col in basis]
        cycle = (basis[0] @ d @ projectors[1] @ d @ projectors[2] @ d @ basis[0])
        reverse = (basis[0] @ d @ projectors[2] @ d @ projectors[1] @ d @ basis[0])
        odd = abs((d[0,1]*d[1,2]*d[2,0]).imag)
        difference = abs(abs(1+1j*cycle)**2-abs(1-1j*cycle)**2)/4
        assert abs(reverse-cycle.conj()) < 1e-15 and abs(difference-odd) < 2e-16
        assert j <= 1/(6*sqrt(3))
        assert np.isclose(j*delta, odd, rtol=2e-7, atol=1e-22)
        results[label] = {"one_or_two_marginal_width_bounds": bound_rows,
                          "balanced_rank_two_eigenvalues": eig.tolist(),
                          "CP_odd_product_absolute": float(odd),
                          "J_of_remainder_absolute": float(j)}
    return results


def check(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            check(a[key], b[key])
    else:
        assert np.allclose(a, b, rtol=2e-7, atol=2e-13)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = compute()
    if args.check:
        check(result, json.loads((HERE/"independent_results.json").read_text(encoding="utf-8")))
        print("1397 independent capacity and spectral-CP certificate verified")
    else:
        print(json.dumps(result, indent=2, allow_nan=False))
