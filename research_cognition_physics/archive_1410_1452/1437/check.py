"""Verify the Z4 selection, neutral directions, and single-pole width ratios."""

import argparse
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def enumerate_terms(max_degree, residue):
    terms = set()
    for counts in itertools.product(range(max_degree + 1), repeat=4):
        h1, conjugate_h1, h2, conjugate_h2 = counts
        degree = sum(counts)
        net = h2 - conjugate_h2
        if degree <= max_degree and h1 - conjugate_h1 + net == 1 and net % 4 == residue:
            terms.add((degree, net))
    return [list(item) for item in sorted(terms)]


def verify():
    p = json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))
    eta = p["bound_on_absolute_net_winding_fraction"]
    width_floors = {
        "integer_lifts": 4.0,
        "mu_winding_only": (2 - 4 * eta) ** 2,
        "both_windings": ((2 - 4 * eta) / (1 + 4 * eta)) ** 2,
    }
    mass_a = 1000 * p["design_A_mass_GeV"]
    mu, tau = p["mu_mass_MeV"], p["tau_mass_MeV"]
    phase_ratio = math.sqrt((1 - (2 * mu / mass_a) ** 2) / (1 - (2 * tau / mass_a) ** 2))
    unit_ratio = (mu / tau) ** 2 * phase_ratio
    terms = p["counterexample_mu_terms"]
    beta = p["counterexample_tan_beta"]
    fractions = [item["mass_fraction"] for item in terms]
    mass = sum(fractions)
    radial = sum(item["degree"] * item["mass_fraction"] for item in terms)
    axial = sum(item["mass_fraction"] * ((item["net_H2"] - 1) * beta + item["net_H2"] / beta)
                for item in terms)
    mu_terms = enumerate_terms(5, p["required_H2_charges"]["mu"])
    tau_terms = enumerate_terms(5, p["required_H2_charges"]["tau"])
    assert mu_terms == [[3, 2], [5, -2], [5, 2]]
    assert tau_terms == [[1, 1], [3, 1], [5, 1]]
    assert (mass, radial, axial) == (1.0, 3.0, 0.0)
    # Also test the exact Schur complement for the three actual heavy Yukawa links.
    import numpy as np
    h1, h2 = 0.18 + 0.05j, 0.23 - 0.02j
    y1, y2, y3, y4 = 0.7, 0.6, 0.5, -0.4
    heavy = np.array([[2.1, y2 * h1.conjugate()], [y4 * h1, 1.8]])
    matched = -np.array([y1 * h2, 0]) @ np.linalg.solve(heavy, np.array([0, y3 * h2]))
    expected = y1 * y2 * y3 * h2**2 * h1.conjugate() / (2.1 * 1.8 - y2 * y4 * abs(h1)**2)
    assert abs(matched - expected) < 1e-14
    return {
        "allowed_degree_and_netH2_through_five": {"mu": mu_terms, "tau": tau_terms},
        "phase_space_ratio_at_design_mass": phase_ratio,
        "unit_strength_width_ratio": unit_ratio,
        "phase_and_mass_corrected_width_floors": width_floors,
        "illustrative_width_ratio_floors": {k: v * unit_ratio for k, v in width_floors.items()},
        "common_relative_error_for_disjoint_unit_and_floor": {
            k: (v - 1) / (v + 1) for k, v in width_floors.items()},
        "counterexample": {"mass": mass, "aligned_h": radial, "odd_A": axial},
        "heavy_messenger_schur_checked": True,
    }


def comparable(x):
    if isinstance(x, float):
        return float(format(x, ".12g"))
    if isinstance(x, dict):
        return {k: comparable(v) for k, v in x.items()}
    return x


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--check", action="store_true")
    modes.add_argument("--write", action="store_true")
    options = parser.parse_args()
    result = comparable(verify())
    target = HERE / "results.json"
    if options.write:
        with target.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")
    elif result != json.loads(target.read_text(encoding="utf-8")):
        raise ValueError("Saved results differ")
    print(json.dumps(result))
