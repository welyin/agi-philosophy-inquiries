"""Check color invariants, delivered partners and finite charge-period records."""
from itertools import combinations, product
from pathlib import Path
import argparse
import json
import math

import numpy as np


def generators(group):
    rotations, symmetric = [], []
    for i, j in combinations(range(3), 2):
        u = np.zeros((3, 3), complex)
        u[i, j], u[j, i] = 1, -1
        rotations.append(u)
        v = np.zeros_like(u)
        v[i, j] = v[j, i] = 1j
        symmetric.append(v)
    if group == "SU3":
        return rotations + symmetric + [1j * np.diag([1, -1, 0]), 1j * np.diag([0, 1, -1])]
    return rotations


def tensor(parts):
    out = np.array([[1.]])
    for part in parts:
        out = np.kron(out, part)
    return out


def projector(group, m, n):
    dimension = 3 ** (m+n)
    penalty = np.zeros((dimension, dimension), complex)
    for gen in generators(group):
        action = np.zeros_like(penalty)
        for k in range(m+n):
            legs = [np.eye(3)] * (m+n)
            legs[k] = gen if k < m else gen.conj()
            action += tensor(legs)
        penalty += action.conj().T @ action
    if group == "mu3_SO3":
        phase = np.exp(2j * np.pi * (m-n) / 3)
        penalty += abs(phase-1)**2 * np.eye(dimension)
    vals, vecs = np.linalg.eigh(penalty)
    kept = vecs[:, vals < 1e-9]
    return kept @ kept.conj().T


def retained_partner_prob(projection):
    # Each row labels active color; columns represent distinct undiscarded antipartners.
    original = np.eye(len(projection)) / np.sqrt(len(projection))
    passed = projection @ original
    failed = original - passed
    assert np.isclose(np.linalg.norm(passed)**2 + np.linalg.norm(failed)**2, 1)
    return float(np.linalg.norm(passed)**2)


def calculate():
    groups = {"SU3": (0, 2), "SO3": (1, 3), "mu3_SO3": (0, 3)}
    out = {}
    for group, (pair, four) in groups.items():
        multiplicities = {}
        for m, n in ((1, 1), (3, 0), (2, 0), (2, 2)):
            proj = projector(group, m, n)
            rank = int(round(np.trace(proj).real))
            assert rank == {(1, 1): 1, (3, 0): 1, (2, 0): pair, (2, 2): four}[(m, n)]
            rate = retained_partner_prob(proj)
            assert np.isclose(rate, rank / (3 ** (m+n)))
            multiplicities[f"{m}R_{n}dual"] = rank
        passed = pair / 9
        dark = passed * np.sin(np.pi/3)**2
        out[group] = {"singlet_multiplicities": multiplicities,
                      "pair_pass_and_minus": float(dark),
                      "pair_pass_and_plus": float(passed - dark),
                      "pair_rejected": 1-passed, "four_leg_pass": four / 81}

    unit_charge_error = .08 / 3
    phase_error = .02
    visibility_floor = .95
    angles = [2 * math.pi * (1/3 + 2*i*unit_charge_error) + j*phase_error
              for i,j in product((-1,1), repeat=2)]
    assert 0 < min(angles) <= max(angles) < math.pi and min(angles) > math.pi/2
    min_minus = (1 - visibility_floor * math.cos(min(angles))) / 2
    return {"round": 1352, "contracts": out, "finite_charge_error": unit_charge_error,
            "path_phase_error_radians": phase_error, "visibility_floor": visibility_floor,
            "true_SO3_joint_floor": min_minus/9,
            "with_efficiency_0.90_floor": .90 * min_minus/9,
            "full_four_leg_gap": 1/81}


def check_saved(actual, expected):
    if isinstance(actual, dict):
        assert actual.keys() == expected.keys()
        for k in actual:
            check_saved(actual[k], expected[k])
    elif isinstance(actual, float):
        assert np.isclose(actual, expected,
                          rtol=1e-12, atol=1e-12)
    else:
        assert actual == expected


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = calculate()
    path = Path(__file__).with_name("results.json")
    if args.write:
        with path.open("x", encoding="utf8") as out:
            out.write(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    else:
        check_saved(result, json.loads(path.read_text(encoding="utf8")))
    print(json.dumps(result))
