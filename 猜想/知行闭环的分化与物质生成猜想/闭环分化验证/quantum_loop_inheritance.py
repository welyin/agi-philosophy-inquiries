"""Compare a fixed two-register quantum feedback cycle and its split records."""

import argparse
import json
from pathlib import Path

import numpy as np


RESULT = Path(__file__).with_name("quantum_loop_inheritance_results.json")
I = np.eye(2)
X = np.array([[0, 1], [1, 0]])
Z = np.diag([1, -1])
PLUS = np.array([1, 1]) / np.sqrt(2)
MINUS = np.array([1, -1]) / np.sqrt(2)
ZERO = np.array([1, 0])


def controlled_flip(control):
    gate = np.zeros((4, 4))
    for s in range(2):
        for r in range(2):
            target_s = s ^ r if control == "R" else s
            target_r = s ^ r if control == "S" else r
            gate[2 * target_s + target_r, 2 * s + r] = 1
    return gate


def density(vector):
    return np.outer(vector, vector.conj())


def run_checks():
    feedback = controlled_flip("R")
    record = controlled_flip("S")
    plus_readout = np.kron(density(PLUS), I)
    seed = np.kron(ZERO, ZERO)
    dephase_ops = (np.kron(Z, I), np.kron(I, Z))
    angles = np.linspace(0, np.pi, 33)
    for theta in angles:
        c, s = np.cos(theta), np.sin(theta)
        rotate = np.array([[c, -s], [s, c]])
        record_after_rotation = record @ np.kron(rotate, I)
        cycle = record_after_rotation @ feedback
        assert np.allclose(cycle.T @ cycle, np.eye(4))
        first = cycle @ seed
        second = cycle @ first
        intact = second @ plus_readout @ second
        assert np.isclose(intact, (1 + np.sin(2 * theta) ** 2) / 2)
        for sign in dephase_ops:
            rho = density(first)
            erased = (rho + sign @ rho @ sign) / 2
            assert np.isclose(np.trace(plus_readout @ cycle @ erased @ cycle.T), 0.5)
        cut = record_after_rotation @ first
        assert np.isclose(cut @ plus_readout @ cut,
                          (1 + np.sin(2 * theta) * np.cos(2 * theta)) / 2)

    # Right before the copying gate the logical input may be any qubit.
    bell_plus, bell_minus = [record @ np.kron(v, ZERO) for v in (PLUS, MINUS)]
    for state in (bell_plus, bell_minus):
        rho = density(state).reshape(2, 2, 2, 2)
        for child in (np.trace(rho, axis1=0, axis2=2),
                      np.trace(rho, axis1=1, axis2=3)):
            assert np.allclose(child, I / 2)
    assert np.isclose(bell_plus @ np.kron(X, X) @ bell_plus, 1)
    assert np.isclose(bell_minus @ np.kron(X, X) @ bell_minus, -1)

    # Recover the complex logical basis operators, including both coherences.
    embed = record[:, [0, 2]]
    for i in range(2):
        for j in range(2):
            unit = np.zeros((2, 2))
            unit[i, j] = 1
            assert np.allclose(record.T @ (embed @ unit @ embed.T) @ record,
                               np.kron(unit, density(ZERO)))
    return {"angles_checked": len(angles), "pi_over_four_plus_probabilities":
            {"intact": 1, "dephased_S_or_R": 0.5, "second_feedback_removed": 0.5},
            "local_bell_marginals_equal": True, "joint_xx_expectations": [1, -1],
            "logical_matrix_units_recovered": 4}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    results = run_checks()
    if args.check:
        if results != json.loads(RESULT.read_text(encoding="utf-8")):
            raise SystemExit("Saved result differs from recomputation.")
        print("Verified", RESULT.name)
    else:
        RESULT.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
