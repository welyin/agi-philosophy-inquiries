"""Inventory Weyl charges, enumerate singlet completions, and match a two-source seesaw."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def anomaly(charges):
    return sum(charges), sum(q**3 for q in charges)


def cancelled(linear, cubic):
    return 2*linear % 12 == 0 and 182*cubic % 72 == 0


def compute():
    cfg = read("inputs.json")
    old = read(cfg["flavor_mediator_inputs"])
    q = old["flavor_charge"]
    tens, fives = cfg["su5_generation_charges_10"], cfg["su5_generation_charges_5bar"]
    pairs = [q[name] for name in ["P", "Pbar", "T", "Tbar"]]
    sectors = {
        "families": [n for n in tens for _ in range(10)]+[n for n in fives for _ in range(5)],
        "two_complete_heavy_pairs": [n for n in pairs for _ in range(5)],
        "flavor_superpartners": [q[name] for name in cfg["neutral_charged_under_flavor"]],
        "neutral_stabilization_and_higgsinos": cfg["optional_stabilization_multipliers"] +
                                                [n for n in cfg["weak_doublet_higgsinos"] for _ in range(2)],
    }
    charges = sum(sectors.values(), [])
    linear, cubic = anomaly(charges)
    mixed = 3*sum(tens) + sum(fives) + sum(pairs)
    assert (linear, cubic, mixed) == (58, 388, 12)
    assert all((182*(n+12)**3 - 182*n**3) % 72 == 0 for n in range(12))
    solutions = {}
    for count in range(3):
        sol = []
        for extra in itertools.combinations_with_replacement(range(12), count):
            a, b = anomaly(extra)
            assert cancelled(linear+a, cubic+b) == ((linear+a) % 6 == 0 and (cubic+b) % 36 == 0)
            if cancelled(linear+a, cubic+b):
                sol.append(list(extra))
        solutions[str(count)] = sol
    assert solutions["1"] == [[2], [8]] and len(solutions["2"]) == 6
    assert all(cancelled(*anomaly([n, (-n) % 12])) for n in range(12))
    assert all(cancelled(*anomaly([n])) for n in [0, 6])

    eps, x = old["scalar_magnitude_over_cutoff"], old["real_X_over_cutoff"]
    s = eps*np.exp(1j*math.pi*old["scalar_angle_over_pi"])
    f = s**4 / x**2
    assert all(abs((s*np.exp(-2j*math.pi*g/12))**4 /
                   (x*np.exp(2j*math.pi*4*g/12))**2 - f) < 1e-14 for g in range(12))
    r_mass = read(cfg["target_neutrino_witness"])["joint_witness"]["R"]
    coeff = 6*(r_mass+1/r_mass) + 16
    r = 8/(coeff+math.sqrt(coeff**2-1936))
    a, b = np.array(cfg["source_a"]), np.array(cfg["source_b"])
    k = cfg["majorana_cross_coefficient"]
    # The original factorized row phase is removed in the common family basis.
    row = np.diag(s**np.array(fives))
    phase = np.diag((s/eps)**np.array(fives))
    inv_abs_rows = eps**(-np.array(fives, dtype=float))
    y0 = inv_abs_rows * math.sqrt(r) * b
    y2 = inv_abs_rows * a / math.sqrt(abs(f))
    yd = row @ np.column_stack([y0, s**2 * (y2 + k*y0)])
    def mass(shift):
        m = np.array([[1, k*s**2], [k*s**2, (k*k+shift)*s**4 + x*x]])
        reduced = phase.conj() @ (yd @ np.linalg.solve(m, yd.T)) @ phase.conj()
        z = s**4/(x*x + shift*s**4)
        reduced_expected = r*np.outer(b,b) + z/abs(f)*np.outer(a,a)
        np.testing.assert_allclose(reduced, reduced_expected, atol=8e-15)
        masses = np.linalg.svd(reduced, compute_uv=False)
        return z, reduced, masses
    z, reduced, masses = mass(0)
    eta = -float(np.angle(z))
    target = np.outer(a,a) + r*np.exp(1j*eta)*np.outer(b,b)
    np.testing.assert_allclose(reduced*np.exp(1j*eta), target, atol=8e-15)
    assert masses[2] < 2e-15 and abs(masses[0]/masses[1]-r_mass) < 1e-12
    z_bad, _, m_bad = mass(cfg["allowed_S4_contamination"])
    harmful_terms = []
    for u,v in itertools.combinations_with_replacement(["S", "X", "A2", "A3"], 2):
        if (2*cfg["neutral_pair_example"][1]+q[u]+q[v])%12 == 0:
            harmful_terms.append([u,v])
    assert len(harmful_terms) == 3
    return {
        "sector_linear_and_cubic": {key: list(anomaly(val)) for key,val in sectors.items()},
        "total_linear_and_cubic": [linear,cubic],
        "residues_mod_6_and_36": [linear % 6, cubic % 36],
        "unit_color_and_weak_complete_pairs": [mixed, mixed],
        "unit_color_and_weak_triplets_only": [mixed, mixed-sum(pairs)],
        "singlet_solutions_by_multiplicity": solutions,
        "two_scalar_insertions_for_charge_2_Majorana_pair": harmful_terms,
        "rank_two_witness": {
            "charge_assignment": {"column_for_b": 0, "column_for_a": 2},
            "relative_light_source_phase_deg": round(math.degrees(eta), 10),
            "mass_ratio": round(masses[0]/masses[1], 10),
            "small_root_relative_amplitude": round(r, 12),
            "max_diff_to_registered_two_source_matrix": round(float(np.max(abs(reduced*np.exp(1j*eta)-target))), 14),
        },
        "same_charges_allowed_counterexample": {
            "relative_light_source_phase_deg": round(-np.angle(z_bad, deg=True), 8),
            "mass_ratio_for_same_real_Yukawas": round(m_bad[0]/m_bad[1], 8)
        }
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    output = compute()
    if args.check:
        assert output == read("results.json")
        print("1688 checks passed")
    else:
        print(json.dumps(output, ensure_ascii=False, indent=2, allow_nan=False))
