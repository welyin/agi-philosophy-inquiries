"""Review of the oriented cofactor. No use of round 1562 author modules."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


BASE = Path(__file__).resolve().parent


def read(relative):
    return json.loads((BASE / relative).read_text(encoding="utf-8-sig"))


def cofactor(matrix):
    # C_ij = conjugated signed minor with row i and column j deleted.
    mask = np.arange(3)
    return np.array([
        [(-1)**(i+j)*np.linalg.det(matrix[np.ix_(mask != i, mask != j)]).conjugate()
         for j in range(3)] for i in range(3)
    ])


def masses(tensors):
    h = np.diag(tensors["H_diag"])
    entries = np.array(tensors["F_real_upper"]) + 1j*np.array(tensors["F_imag_upper"])
    f = np.zeros((3, 3), dtype=complex)
    iu, ju = np.triu_indices(3)
    f[iu, ju] = f[ju, iu] = entries
    r = tensors["r"]
    s = complex(tensors["s_real"], tensors["s_imag"])
    return r*(h+s*f), h+f


def phase_review(y, r):
    k = (y[0]/y)**2
    eps = np.sum(k[1:]*r/(1-k[1:]*(1+r)))
    assert eps < .002
    for phi in np.linspace(0, math.pi, 641):
        z = k*(r*np.exp(1j*phi)-1)
        factors = 1+z
        incremental_phase = sum(np.angle(factors))
        w = (1+3*z)/factors
        derivative = np.sum(np.imag(w))
        step = 1e-6
        # Strip the positive radial h^3; evaluate the invariant logarithmic derivative.
        plus, minus = [np.prod(1+np.exp(sign*2*step)*z) for sign in (1, -1)]
        fd = np.angle(plus/minus)/(2*step)
        assert abs(derivative-fd) < 3e-9
        assert phi-3e-14 <= incremental_phase <= phi+eps*math.sin(phi)+3e-14
        assert abs(abs(w[0]-3)-2/r) < 4e-15
        if .01 <= phi <= 3.13:
            assert np.min(np.sign(w.imag)) > 0
    return float(eps)


def audit():
    previous = read("../1556/results.json")["candidates"]
    hi = read("../1397/inputs.json")["high_scale"]
    table_down = np.array(hi["down"])
    table_up = np.array(hi["up"])
    edm = read("../1389/inputs.json")["edm"]
    slope, stat, syst = edm["abs_slope_e_fm"]
    theta = edm["limit_90_ecm"]/((slope-2*(stat+syst))*edm["fm_in_cm"])
    out = {"illustrative_theta_limit_near_zero_rad": theta, "candidates": {}}
    for name in ("release_down_A", "release_down_B"):
        up, down = masses(previous[name]["tensors"])
        yu, yc, yt = np.linalg.svd(up, compute_uv=False)[::-1]
        y = np.linalg.svd(down, compute_uv=False)[::-1]
        d, s, b = y
        r0 = table_down[0, 0]/d
        rmin, rmax = (table_down[0, 0]+np.array([-2, 2])*table_down[0, 1])/d
        det = np.linalg.det(down)
        # Random independent U(3) bases, including nontrivial determinant phases.
        gen = np.random.default_rng(3562)
        def basis():
            a = gen.normal(size=(3,3))+1j*gen.normal(size=(3,3))
            return np.linalg.qr(a)[0]
        left, right = basis(), basis()
        rotated = left@down@right.conj().T
        orient = np.linalg.det(left)*np.linalg.det(right).conjugate()
        expected = left@cofactor(down)@right.conj().T/orient
        assert np.linalg.norm(cofactor(rotated)-expected) < 5e-19
        for phi in (0., .7, 2.4, math.pi):
            c = d*d*(r0*np.exp(1j*phi)-1)/det.conjugate()
            effective = down+c*cofactor(down)
            changed = rotated+c*orient*cofactor(rotated)
            assert np.linalg.norm(changed-left@effective@right.conj().T) < 2e-17
            z = d*d/y**2*(r0*np.exp(1j*phi)-1)
            assert np.allclose(y*abs(1+z), np.linalg.svd(effective, compute_uv=False)[::-1],
                               rtol=1e-11, atol=0.)
            assert abs(np.linalg.det(effective)/det-np.prod(1+z)) < 2e-12
        eps = phase_review(y, rmax)
        phase_review(y, r0)
        phase_review(y, rmin)
        umax, _, tmax = table_up @ np.array([1, 2])
        upper_from_top = (tmax+yu)/(yc*yt)
        other_lower = np.array([yc-upper_from_top*yu*yt, yt-upper_from_top*yu*yc])
        assert min(other_lower) > umax
        c_up = (umax+yu)/(yc*yt)
        c_down_min = (rmin-1)*d/(s*b)
        out["candidates"][name] = {
            "up_yukawas": [yu, yc, yt], "down_yukawas": y.tolist(),
            "r_central": r0, "r_two_printed_widths": [rmin, rmax],
            "phase_margin_coefficient_upper": eps,
            "up_reassignment_other_two_lower": other_lower.tolist(),
            "coefficient_abs_ratio_upper_two_widths": c_up/c_down_min,
            "kappa_tilde_d_upper_near_zero": 2*theta/rmin,
            "kappa_d_zero_phase_central": 3-2/r0,
            "kappa_d_pi_phase_central": 3+2/r0,
        }
    return out


def check_equal(now, saved):
    if isinstance(now, dict):
        assert now.keys() == saved.keys()
        for key in now:
            check_equal(now[key], saved[key])
    else:
        assert np.allclose(now, saved, rtol=2e-11, atol=1e-20), (now, saved)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = audit()
    result = BASE / "independent_results.json"
    if args.write:
        with result.open("x", encoding="utf-8") as file:
            json.dump(data, file, indent=2)
            file.write("\n")
    else:
        check_equal(data, read(result))
    print("Independent round 1562 checks passed.")
