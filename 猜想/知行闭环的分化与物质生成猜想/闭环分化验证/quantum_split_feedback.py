"""Local action and inherited-relation checks on the two active child loops."""

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from quantum_loop_inheritance import I, X, Z, ZERO, controlled_flip, density


RESULT = Path(__file__).with_name("quantum_split_feedback_results.json")
Y = np.array([[0, -1j], [1j, 0]])
BITS = tuple(product(range(2), repeat=2))
LOCAL_ID = np.eye(4)
TOL = 1e-12


def check(actual, target):
    np.testing.assert_allclose(actual, target, atol=TOL, rtol=0)


def marginals(state):
    rho = density(state).reshape(4, 4, 4, 4)
    return np.trace(rho, axis1=1, axis2=3), np.trace(rho, axis1=0, axis2=2)


def run():
    # Basis order: S_A, R_A, S_B, R_B. Both return memories start blank.
    local_input = np.kron(I, ZERO[:, None])
    blank_returns = np.kron(local_input, local_input)
    phi_plus = np.array([1, 0, 0, 1]) / np.sqrt(2)
    psi_minus = np.array([0, 1, -1, 0]) / np.sqrt(2)
    seeds = np.column_stack([phi_plus, psi_minus])
    rotate = np.array([[1, -1], [1, 1]]) / np.sqrt(2)
    without_feedback = controlled_flip("S") @ np.kron(rotate, I)
    cycle = without_feedback @ controlled_flip("R")
    first = np.kron(cycle, cycle) @ blank_returns
    action = np.kron(X, X)
    inherited = np.kron(X, I)
    return_bit = np.kron(I, Z)
    relation = np.kron(inherited, inherited)
    check(cycle @ np.kron(Z, I) @ cycle.T, -np.kron(Y, Y))
    check(cycle @ (-np.kron(Y, Y)) @ cycle.T, inherited)
    check(np.linalg.matrix_power(cycle, 6), -LOCAL_ID)
    check(inherited @ action, action @ inherited)
    check(inherited @ return_bit, return_bit @ inherited)

    induced = {}
    cuts = {}
    partials = {}
    memory_pulse = np.kron(I, Y)
    for pulses in BITS:
        pa, pb = pulses
        pulse = np.kron(memory_pulse if pa else LOCAL_ID,
                        memory_pulse if pb else LOCAL_ID)
        for removed in BITS:
            ca, cb = removed
            last = np.kron(without_feedback if ca else cycle,
                           without_feedback if cb else cycle)
            frame = last @ pulse @ first
            distinguished = frame @ seeds
            left = distinguished.conj().T @ np.kron(action, LOCAL_ID) @ distinguished
            right = distinguished.conj().T @ np.kron(LOCAL_ID, action) @ distinguished
            expected = ((-1) ** pa, (-1) ** pb)
            for col in range(2):
                check(left[col, col], 0 if ca else expected[0])
                check(right[col, col], 0 if cb else expected[1])
            for state in distinguished.T:
                reductions = marginals(state)
                for site in range(2):
                    key = (site, pulses[site], removed[site])
                    if key in partials:
                        check(reductions[site], partials[key])
                    else:
                        partials[key] = reductions[site]
            if pulses == (0, 0):
                cuts[f"cut_{ca}{cb}"] = [0 if ca else 1, 0 if cb else 1]
            if removed != (0, 0):
                continue
            induced[f"return_{pa}{pb}"] = list(expected)
            # Simultaneous eigenvalue tests on the entire parent logical subspace.
            check(np.kron(action, LOCAL_ID) @ frame, expected[0] * frame)
            check(np.kron(LOCAL_ID, action) @ frame, expected[1] * frame)
            check(relation @ frame, frame @ np.kron(Z, Z))
            check(relation @ distinguished, distinguished @ Z)

            # The action readout gives the correct pulse syndromes without harming
            # unknown inputs, including a reference purifying an arbitrary input.
            pa_projector = (LOCAL_ID + expected[0] * action) / 2
            pb_projector = (LOCAL_ID + expected[1] * action) / 2
            diagnosed = np.kron(pa_projector, pb_projector) @ frame
            undo = first.conj().T @ pulse.conj().T @ last.conj().T
            check(undo @ diagnosed, np.eye(4))

    # Natural return distribution untouched by pulse, measurements or reset.
    second = np.kron(cycle, cycle) @ first
    for a, b in BITS:
        projector = np.kron((LOCAL_ID + (-1)**a * return_bit) / 2,
                            (LOCAL_ID + (-1)**b * return_bit) / 2)
        check(second.conj().T @ projector @ second, np.eye(4) / 4)
    return {
        "local_action_expectations_for_perturbed_returns": induced,
        "local_action_expectations_after_feedback_cuts": cuts,
        "old_joint_parity_on_phi_plus_and_psi_minus": [1, -1],
        "natural_returns_for_next_cycle_00_01_10_11": [0.25] * 4,
        "exact_recovery_on_parent_action_input_dimension": 4,
        "local_cycle_sixth_power_is_minus_identity": True,
        "pulse_and_cut_settings_tested": len(BITS) ** 2,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    answer = run()
    if args.check:
        if json.loads(RESULT.read_text(encoding="utf-8")) != answer:
            raise SystemExit("Saved answer differs from recomputation.")
        print("Verified", RESULT.name)
    else:
        RESULT.write_text(json.dumps(answer, indent=2) + "\n", encoding="utf-8")
