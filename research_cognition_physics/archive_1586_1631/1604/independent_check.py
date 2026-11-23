"""Conditional R2 dipole test, with independent quadrature for the finite-mass kernel."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent
OUTPUT = HERE / "independent_results.json"


def old(number):
    return json.loads((ARCHIVE / str(number) / "inputs.json").read_text(encoding="utf-8-sig"))


def kernel_closed(x):
    if x == 0:
        return 0.25
    if x == 1:
        return 7 / 72
    fs = (x + 1) / (4 * (1 - x)**2) + x * math.log(x) / (2 * (1 - x)**3)
    ff = (x*x - 5*x - 2) / (12 * (x - 1)**3) + x * math.log(x) / (2 * (x - 1)**4)
    return (5 / 3) * fs - ff


def kernel_quadrature(x):
    nodes, weights = leggauss(160)
    t = (nodes + 1) / 2
    values = t * (1 - t) * (5 / 3 - t) / (1 - t + x*t)
    return float(weights @ values) / 4


def evaluate():
    angles = old(1370)
    atmospheric = old(1580)["theta23_sin2"]
    phase_space = (1 - (0.00051099895 / 0.1056584)**2)**3
    normalization = (9 / 137.035999 * 0.1056584**5 * 2.1969803e-6 * phase_space
                     / (1024 * math.pi**4 * 6.582119569e-25))
    kstar = kernel_closed((180 / 1000)**2)
    assert abs(kstar - kernel_quadrature((180 / 1000)**2)) < 1e-13
    for x in (0, 0.015, 0.5, 1, 2, 4):
        assert abs(kernel_closed(x) - kernel_quadrature(x)) < 1e-12
    out = {}
    for order in ("NO", "IO"):
        x = angles[order]["s12"][1]
        z = angles[order]["s13"][2 if order == "NO" else 1]
        v = atmospheric[order][2 if order == "NO" else 1]
        if order == "NO":
            prod = math.sqrt(1-z) * (math.sqrt(x*(1-x)*(1-v)) - (1-x)*math.sqrt(v*z))
        else:
            prod = math.sqrt(z*(1-z)*v)
        mass_ratio = (normalization * (prod*kstar)**2 / 1.5e-13)**0.25
        width_ratio = 1 / (16 * math.pi * mass_ratio**2)
        out[order] = {
            "abs_b_e_b_mu_lower": prod,
            "M_over_Y_TeV_lower": mass_ratio / 1000,
            "M_over_Y_TeV_lower_at_amplitude_pollution_0p2": mass_ratio * math.sqrt(0.8) / 1000,
            "Gamma_visible_over_M_cubed_GeV_minus2_upper": width_ratio,
            "Gamma_visible_MeV_upper_at_1_TeV": width_ratio * 1000**3 * 1000
        }
    zmax = angles["IO"]["s13"][2]
    pmax = math.sqrt(zmax * (1 - zmax) * atmospheric["IO"][2])
    helicities = 1 + (0.00051099895 / 0.1056584)**2
    out["IO"]["abs_b_e_b_mu_upper"] = pmax
    out["IO"]["M_over_Y_TeV_sufficient_without_pollution"] = (
        normalization * (pmax / 4)**2 * helicities / 1.5e-13)**0.25 / 1000
    # Preserve charged-lepton marginals with another left-handed Yukawa texture.
    b = np.sqrt([0.65, 0.20, 0.15])
    ranked = np.outer(np.array([1, 1, 0]) / math.sqrt(2), b)
    separated = np.array([[b[0], 0, b[2]], [0, b[1], 0], [0, 0, 0]])
    kernels = np.diag([kernel_closed((m/1000)**2) for m in (0.002, 1.3, 180)])
    assert np.allclose(np.sum(abs(ranked)**2, axis=0), np.sum(abs(separated)**2, axis=0))
    assert np.linalg.matrix_rank(ranked) == 1 and np.linalg.matrix_rank(separated) == 2
    return {"kernel_floor": kstar, "K0": 0.25, "K1": 7/72,
            "normalization_GeV4": normalization, "ordering": out,
            "equal_column_strength_example": {
                "squared_norms": (b*b).tolist(),
                "rank_one_weighted_e_mu": float((ranked.T @ kernels @ ranked)[0, 1]),
                "rank_two_weighted_e_mu": float((separated.T @ kernels @ separated)[0, 1])}}


def matches(a, b):
    if isinstance(a, dict):
        return isinstance(b, dict) and a.keys() == b.keys() and all(matches(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return isinstance(b, list) and len(a) == len(b) and all(matches(x, y) for x, y in zip(a, b))
    return math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-15)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--check", action="store_true")
    action.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = evaluate()
    if args.write_new:
        with OUTPUT.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2, ensure_ascii=False)
            stream.write("\n")
    elif not matches(result, json.loads(OUTPUT.read_text(encoding="utf-8"))):
        raise SystemExit("Stored values differ from recalculation.")
    print("1604 independent check passed")
