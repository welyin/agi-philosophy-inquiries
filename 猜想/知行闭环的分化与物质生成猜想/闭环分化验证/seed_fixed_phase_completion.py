"""Test coherent reuse of one fixed phase from note 46, and the addressable-gate extension."""

import argparse
import json
from pathlib import Path

import numpy as np

from public_choice_quantum_return import CNOT, H, I, X, Y, Z, equal


RESULT = Path(__file__).with_name("seed_fixed_phase_completion_results.json")
PAULIS = (X, Y, Z)


def rotation(axis, angle):
    generator = sum(coefficient * p for coefficient, p in zip(axis, PAULIS))
    return np.cos(angle / 2) * I - 1j * np.sin(angle / 2) * generator


def closest_power(half_angle, target, budget):
    # Physical distance permits an overall sign for rotations in SU(2).
    residuals = half_angle * np.arange(budget + 1) - target / 2
    return int(np.argmax(np.abs(np.cos(residuals))))


def run():
    rz = np.exp(-1j * np.pi / 6) * np.diag([1, np.exp(1j * np.pi / 3)])
    rx = H @ rz @ H
    word = rz @ rx
    equal(np.linalg.det(word), 1)
    equal(np.trace(word), 1.5)
    axis = np.array([np.sqrt(3), 1, np.sqrt(3)]) / np.sqrt(7)
    conjugate = rz @ word @ rz.conj().T
    angle = 2 * np.arccos(0.75)
    equal(word, rotation(axis, angle))
    other = np.array([
        axis[0] / 2 - np.sqrt(3) * axis[1] / 2,
        np.sqrt(3) * axis[0] / 2 + axis[1] / 2,
        axis[2],
    ])
    equal(conjugate, rotation(other, angle))
    directions = np.array([axis, other, np.cross(axis, other)])
    assert np.linalg.matrix_rank(directions) == 3
    target = rotation(axis, np.pi / 2) @ rotation(other, np.pi / 5)
    trials = []
    for budget in (32, 256, 2048):
        a = closest_power(angle / 2, np.pi / 2, budget)
        b = closest_power(angle / 2, np.pi / 5, budget)
        actual = np.linalg.matrix_power(word, a) @ np.linalg.matrix_power(conjugate, b)
        error = min(np.linalg.norm(actual - sign * target, ord=2) for sign in (1, -1))
        trials.append({"per_axis_search_limit": budget, "word_powers": [a, b],
                       "operator_norm_error_up_to_sign": round(float(error), 9)})
    assert trials[-1]["operator_norm_error_up_to_sign"] < 0.005
    # Is a full local calibration possible in this particular qubit tensor model?
    # Test linear span without assuming exact preparation of the Y eigenstate.
    zero = I[:, 0]
    one = I[:, 1]
    plus = H @ zero
    phased_plus = rz @ plus
    projectors = [np.outer(v, v.conj()) for v in (zero, one, plus, phased_plus)]
    design = np.array([[np.trace(p @ rho).real for p in (I, X, Y, Z)]
                       for rho in projectors])
    equal(abs(np.linalg.det(design)), np.sqrt(3))
    ranks = [int(np.linalg.matrix_rank(design)), int(np.linalg.matrix_rank(np.kron(design, design)))]
    assert ranks == [4, 16]
    # Physical effect for phased_plus is available by reversing its preparation.
    equal(np.linalg.matrix_power(rz, 5) @ rz, -I)
    equal((H @ rz.conj().T) @ phased_plus, zero)
    # The CNOT has to be independently routable for the multi-qubit conclusion.
    equal(CNOT @ np.kron(plus, zero), np.array([1, 0, 0, 1]) / np.sqrt(2))
    return {
        "one_reusable_nonreal_phase_in_radians": "pi/3",
        "trace_of_determinant_one_word": "3/2",
        "axis_and_conjugate_and_bracket_rank": 3,
        "one_and_two_wire_exact_tomography_ranks": ranks,
        "actual_words_for_one_two_axis_target": trials,
        "target_uses_arbitrary_angles_only_in_comparison": True,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.check:
        if json.loads(RESULT.read_text(encoding="utf-8")) != result:
            raise SystemExit("Saved fixed-phase results differ.")
        print("Verified", RESULT.name)
    else:
        RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
