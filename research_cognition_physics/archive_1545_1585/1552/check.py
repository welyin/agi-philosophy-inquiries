"""Exact rational Bernstein cover of one shared SO(10) charged family."""
import argparse
import json
import math
import runpy
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


@dataclass(frozen=True)
class G:
    """Gaussian integer for the global polynomial certificates."""
    r: int
    i: int = 0

    def __add__(self, other):
        if isinstance(other, int):
            return G(self.r + other, self.i)
        return G(self.r + other.r, self.i + other.i)

    __radd__ = __add__

    def __neg__(self):
        return G(-self.r, -self.i)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if isinstance(other, int):
            return G(self.r * other, self.i * other)
        if not isinstance(other, G):
            return NotImplemented
        return G(self.r * other.r - self.i * other.i, self.r * other.i + self.i * other.r)

    __rmul__ = __mul__

    def conjugate(self):
        return G(self.r, -self.i)

    def abs2(self):
        return self.r ** 2 + self.i ** 2


def adjugate(m):
    out = np.empty((3, 3), dtype=object)
    for i in range(3):
        for j in range(3):
            a, b = (j + 1) % 3, (j + 2) % 3
            c, d = (i + 1) % 3, (i + 2) % 3
            out[i, j] = m[a, c] * m[b, d] - m[a, d] * m[b, c]
    return out


def norm2(m):
    return sum(v.abs2() for v in m.flat)


def from_pairs(rows):
    return np.array([[G(*v) for v in row] for row in rows], dtype=object)


def bernstein(samples):
    """Interpolate the exact degree-n polynomial at j/n into Bernstein form."""
    degree = len(samples) - 1
    coeff = [F(0)] * (degree + 1)
    product = [1]
    diffs = list(samples)
    for k in range(degree + 1):
        weight = F(diffs[0], math.factorial(k))
        for j, val in enumerate(product):
            coeff[j] += weight * val
        if k < degree:
            product = [-k * x + degree * y for x, y in zip(product + [0], [0] + product)]
            diffs = [b - a for a, b in zip(diffs, diffs[1:])]
    return [sum(coeff[k] * F(math.comb(j, k), math.comb(degree, k)) for k in range(j + 1))
            for j in range(degree + 1)]


def certify(ctrl, remaining, depth=0):
    if min(ctrl) > 0:
        return (1, depth)
    if remaining == 0:
        raise ArithmeticError("Could not certify positivity")
    row = ctrl
    left, right = [row[0]], [row[-1]]
    while len(row) > 1:
        row = [(a + b) / 2 for a, b in zip(row, row[1:])]
        left.append(row[0])
        right.append(row[-1])
    l, dl = certify(left, remaining - 1, depth + 1)
    r, dr = certify(right[::-1], remaining - 1, depth + 1)
    return l + r, max(dl, dr)


def charged_reference(cfg, point):
    parent = read(HERE / cfg["historical"])
    h = read(HERE / cfg["quarks"])["high_scale"]
    up, down = [np.array(h[name])[:, 0] for name in ("up", "down")]
    lep = np.array(read(HERE / cfg["charged_leptons"])["scales"]["1e16"])[3:, 0]
    old = runpy.run_path(str((HERE / cfg["historical"]).with_name("check.py")))
    u0, d0 = old["matrices"](np.array(h["angles_theta12_theta13_theta23_delta"])[:, 0],
                               up, down, parent["charged_only_completion_relative_phases"])
    z, err = old["charged_solution"](u0, d0, lep, parent["charged_only_start_x_y"])
    c = lep[-1] / np.linalg.norm(d0 + z * u0, 2)
    def convert(v):
        a = np.array(v, dtype=float)
        return (a[..., 0] + 1j * a[..., 1]) / point["scale"]
    ur, dr = convert(point["U"]), convert(point["D"])
    cr, jr = convert(point["c"]), convert(point["d"])
    differences = [np.linalg.norm(ur - u0 * up[-1]), np.linalg.norm(dr - d0 * down[-1]),
                   abs(cr - c / down[-1]), abs(jr - c * z / up[-1])]
    assert max(differences) < 1e-11 and err < 2e-11
    er = cr * dr + jr * ur
    rel = np.max(np.abs(np.linalg.svd(er, compute_uv=False)[::-1] / lep - 1))
    assert rel < 2e-10
    return ur, dr, er, cr, jr


def evaluate():
    cfg = read(HERE / "inputs.json")
    pt = cfg["point"]
    scale = pt["scale"]
    u, d = from_pairs(pt["U"]), from_pairs(pt["D"])
    c, j = G(*pt["c"]), G(*pt["d"])
    e = d * c + u * j
    direction = np.array([G(*v) for v in pt["electron"]], dtype=object)
    norm_e = norm2(direction)
    h = e @ e.conj().T
    he = h @ direction
    ray = sum(a.conjugate() * b for a, b in zip(direction, he))
    assert ray.i == 0
    residual_sq = F(norm2(he * norm_e - direction * ray.r), norm_e ** 3)
    separation = F(norm2(adjugate(e)), 3 * norm2(e)) - F(ray.r, norm_e)
    assert separation > 0 and residual_sq / separation**2 < F(1, 10**18)
    bnum, bden = cfg["shape_fraction"]
    anum, aden = cfg["electron_fraction"]
    degree = 24
    certificate = {}
    for sign in [-1, 1]:
        data = [[], [], []]
        for i in range(degree + 1):
            t = 2 * i - degree
            nz = G(sign * degree, sign * t)
            dz = G(degree, -t)
            x = d * (scale * dz) - e * nz
            y = u * (scale * nz * dz * j) + x * (3 * scale * dz + c * nz)
            adj_x = adjugate(x)
            p = y @ adj_x @ y.T
            a, b = norm2(p), norm2(adjugate(p))
            data[0].append(((x @ adj_x)[0, 0]).abs2())
            data[1].append(bnum * a**2 - bden * b)
            data[2].append(anum * a * norm_e - aden * norm2(p.conj().T @ direction))
        certs = [certify(bernstein(vals), cfg["bernstein_max_bisections"]) for vals in data]
        certificate[str(sign)] = dict(zip(("det_X", "spectral_shape", "electron_response"), certs))

    ur, dr, er, cr, jr = charged_reference(cfg, pt)
    phases = np.exp(2j * np.pi * np.arange(4096) / 4096)[:, None, None]
    x = dr - phases * er
    y = phases * jr * ur + (3 + phases * cr) * x
    pmat = y @ np.linalg.solve(x, y.swapaxes(-1, -2))
    vals, axes = np.linalg.eigh(pmat @ pmat.conj().swapaxes(-1, -2))
    electron = np.linalg.svd(er)[0][:, -1]
    pe = np.abs(electron.conj() @ axes[:, :, -1].T)**2
    rho = (vals[:, 1] - vals[:, 0]) / (vals[:, 2] - vals[:, 0])
    beta = bnum / bden
    ratio = math.sqrt(beta) / (1 - 2 * beta)
    pe_bound = (math.sqrt((anum/aden) / (1 - 2 * beta)) + 1e-9)**2
    oscillation = read(HERE / cfg["oscillation"])
    held = {}
    for ordering in ("NO", "IO"):
        obs = oscillation[ordering]
        if ordering == "NO":
            lo, hi = obs["solar"][1] / obs["atm_abs"][2], obs["solar"][2] / obs["atm_abs"][1]
            pe_min = obs["s13"][1]
        else:
            lo = 1 - obs["solar"][2] / obs["atm_abs"][1]
            hi = 1 - obs["solar"][1] / obs["atm_abs"][2]
            pe_min = obs["s12"][1] * (1 - obs["s13"][2])
        root = math.sqrt(lo)
        held[ordering] = {"lower_gap_fraction_box": [lo, hi], "heaviest_electron_fraction_floor": pe_min,
                          "relative_shape_correction_floor": (root - ratio) / (1 + root)}
        assert lo > ratio**2 and pe_min > pe_bound
    db = read(HERE / cfg["reactor"])["daya_bay"]
    assert db["sin_sq_2theta13"] - 2*db["one_sigma"] > 4 * pe_bound * (1-pe_bound)
    return {"exact_chart_certificates_leaves_depth": certificate, "charged_match_scale_GeV": 1e16,
            "singular_ratio_mid_over_high_upper": ratio, "lower_gap_fraction_upper": ratio**2,
            "heaviest_electron_fraction_upper": pe_bound, "NO_reactor_amplitude_upper":4*pe_bound*(1-pe_bound),
            "DayaBay_two_printed_width_floor":db["sin_sq_2theta13"] - 2*db["one_sigma"],
            "oscillation_holdouts": held,
            "grid_diagnostic_only":{"lower_gap_fraction":[float(rho.min()),float(rho.max())],
                                    "heaviest_electron_fraction":[float(pe.min()),float(pe.max())]}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    check = parser.parse_args().check
    output = evaluate()
    path = HERE / "results.json"
    if check:
        expected = read(path)
        diagnostic = output.pop("grid_diagnostic_only")
        reference = expected.pop("grid_diagnostic_only")
        assert all(np.allclose(value, reference[key], atol=1e-11, rtol=1e-8) for key,value in diagnostic.items())
        assert json.loads(json.dumps(output)) == expected
        print("1552 exact phase certificate and charged/oscillation replay passed")
    else:
        path.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
