"""Independent finite-box and positive-rank audit of the oscillation-tensor test."""

import argparse
import json
from dataclasses import dataclass
from decimal import Decimal
from itertools import product
from math import isqrt
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SCALE = 10**12
ONE_RESULT = HERE / "independent_results.json"


@dataclass(frozen=True)
class Box:
    """Inclusive interval with integer endpoints in units of 1/SCALE."""
    lo: int
    hi: int

    def __add__(self, other):
        return Box(self.lo + other.lo, self.hi + other.hi)

    def __neg__(self):
        return Box(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        values = [x*y for x, y in product((self.lo, self.hi), (other.lo, other.hi))]
        return Box(min(values)//SCALE, -(-max(values)//SCALE))

    def __truediv__(self, other):
        assert other.lo * other.hi > 0
        endpoints = list(product((self.lo, self.hi), (other.lo, other.hi)))
        return Box(min(SCALE*x//y for x, y in endpoints),
                   max(-(-SCALE*x//y) for x, y in endpoints))

    def sqrt(self):
        assert self.lo >= 0
        upper = isqrt(SCALE*self.hi)
        return Box(isqrt(SCALE*self.lo), upper + int(upper*upper != SCALE*self.hi))

ONE = Box(SCALE, SCALE)


def scaled(value):
    number = Decimal(str(value)) * SCALE
    assert number == int(number)
    return int(number)


def diagonal_enclosures(params, order, cp_sign):
    x, y, z, r = params
    s12, c12 = x.sqrt(), (ONE - x).sqrt()
    s23, c23 = z.sqrt(), (ONE - z).sqrt()
    p = y.sqrt() if cp_sign == 1 else -y.sqrt()
    v = (ONE - y).sqrt()
    # Mass-squared shift sets the first eigenvalue to zero; |Delta31| = 1.
    first = [c12*v, -s12*c23 - c12*s23*p, s12*s23 - c12*c23*p]
    second = [s12*v, c12*c23 - s12*s23*p, -c12*s23 - s12*c23*p]
    third = [p, s23*v, c23*v]
    delta = ONE if order == 1 else -ONE
    ds = []
    for row in range(3):
        j, k = (row+1) % 3, (row+2) % 3
        off = r*second[j]*second[k] + delta*third[j]*third[k]
        assert off.lo * off.hi > 0
        # adj(H) = r Delta u_1 u_1^T and the off-diagonal minor fixes d_row.
        ds.append(-r*delta*first[j]*first[k]/off)
    return ds


def audit_continuum(settings):
    names = ["sin2_theta12", "sin2_theta13", "sin2_theta23", "solar_over_abs31"]
    outer = tuple(Box(*map(scaled, settings[k])) for k in names)
    record = {}
    for order, sign in product((-1, 1), repeat=2):
        ds = diagonal_enclosures(outer, order, sign)
        left, right = ds[0] - ds[1], ds[0] - ds[2]
        min_gap = max(min(left.lo, -right.hi), min(right.lo, -left.hi))
        assert min_gap > scaled(0.0119016)
        record[f"ordering_{order}_cos_phase_{sign}"] = {
            "certified_cells": 1, "max_bisection_depth": 0,
            "min_separation_over_abs31": min_gap/SCALE,
        }
    return record


def pmns(x, y, z, phase):
    def rotation(i, j, sine2, complex_phase=0.):
        mat = np.eye(3, dtype=complex)
        c, s = np.sqrt(1-sine2), np.sqrt(sine2)*np.exp(-1j*complex_phase)
        mat[i, i] = mat[j, j] = c
        mat[i, j], mat[j, i] = s, -s.conjugate()
        return mat
    return rotation(1, 2, z) @ rotation(0, 2, y, phase) @ rotation(0, 1, x)


def audit_rank(box, weights):
    rng = np.random.default_rng(1625)
    samples = 0
    for order, phase in product((-1, 1), (0., np.pi, 0.37, 2.1, -0.8)):
        for _ in range(4):
            x, y, z, r = [rng.uniform(*box[key]) for key in
                          ("sin2_theta12", "sin2_theta13", "sin2_theta23", "solar_over_abs31")]
            U = pmns(x, y, z, phase)
            ev = np.array([0., r, float(order)])
            H = (U*ev) @ U.conj().T
            cp = (H[0, 1]*H[1, 2]*H[2, 0]).imag
            j = (1-y)*np.sqrt(x*(1-x)*z*(1-z)*y)*np.sin(phase)
            assert abs(cp - r*order*(order-r)*j) < 1e-14
            hzero = H + (-min(ev))*np.eye(3)
            null = U[:, np.argmin(ev)]
            assert np.linalg.norm(hzero @ null) < 1e-14
            assert (null.conj() @ (weights * null)).real > 0
            H = hzero + rng.uniform(0.01, 0.2)*np.eye(3)
            # Find the largest admissible b by PSD feasibility in the original basis.
            low, high = 0., min(H.diagonal().real / weights)
            for _ in range(75):
                b = (low+high)/2
                if np.linalg.eigvalsh(H-np.diag(b*weights))[0] >= 0:
                    low = b
                else:
                    high = b
            assert low > 0
            R = H-np.diag(low*weights)
            # Build a rectangular factor from the best principal pair; no eigensystem of R.
            pairs = ((0, 1), (0, 2), (1, 2))
            idx = max(pairs, key=lambda pair: np.linalg.det(R[np.ix_(pair, pair)]).real)
            chol = np.linalg.cholesky(R[np.ix_(idx, idx)])
            columns = np.linalg.solve(chol, R[list(idx), :]).conj().T
            assert np.linalg.norm(columns @ columns.conj().T - R, 2) < 2e-12
            assert np.linalg.matrix_rank(columns, tol=1e-8) == 2
            samples += 1
    return samples


def run():
    cfg = json.loads((HERE / "inputs.json").read_text(encoding="utf-8-sig"))
    box = cfg["box"]
    branches = audit_continuum(box)
    minimum = min(v["min_separation_over_abs31"] for v in branches.values())
    family = json.loads((HERE / cfg["family_inputs"]).read_text(encoding="utf-8-sig"))
    fit = json.loads((HERE / cfg["family_output"]).read_text(encoding="utf-8-sig"))
    weights = np.array([family["pole_MeV"]["electron"][0], family["pole_MeV"]["muon"][0],
                        fit["roots_MeV"][1]])
    x, y, z = [box[k] for k in ("sin2_theta12", "sin2_theta13", "sin2_theta23")]
    radicand = np.prod([x[0], 1-x[1], z[0], 1-z[1], y[0]])
    return {
        "continuous_endpoints_certified": branches,
        "box_separation_eV2_rounded_down": int(1e12 * minimum*box["abs31_eV2"][0])/1e12,
        "positive_bound_abs_J_per_abs_sin_phase_rounded_down": int(1e8*(1-y[1])*np.sqrt(radicand))/1e8,
        "positive_rank_two_witnesses_and_zero_mass_null_checks": audit_rank(box, weights/max(weights)),
        "interval_endpoint_unit": "1e-12 (outward integer arithmetic)",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.check:
        assert result == json.loads(ONE_RESULT.read_text(encoding="utf-8"))
        print("1625 independent PASS: continuous box, CP, both orderings and rank-two construction")
    else:
        with ONE_RESULT.open("x", encoding="utf-8") as file:
            file.write(json.dumps(result, indent=2) + "\n")
