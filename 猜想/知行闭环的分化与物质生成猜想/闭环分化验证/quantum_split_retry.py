"""Coherent early/late retry using the unchanged return test and swap from note 38."""

import argparse
import json
from pathlib import Path

import numpy as np

from quantum_loop_inheritance import I, X, Z, ZERO, controlled_flip
from quantum_split_router import check_equal


RESULT = Path(__file__).with_name("quantum_split_retry_results.json")
Y = np.array([[0, -1j], [1j, 0]])


def clock_blocks(loop):
    # Flag sectors 00,01,10,11; data order S_A,R_A,S_C,R_C.
    eye = np.eye(16)
    projector = np.kron((np.eye(4) - np.kron(Y, Y).real) / 2, np.diag([1, 0, 0, 0]))
    swap = eye.reshape(2, 2, 2, 2, 16).swapaxes(1, 2).reshape(16, 16)
    pair = np.kron(loop, loop)
    evolution = np.kron(np.diag([1, 0, 0, 0]), np.kron(loop, np.eye(4)))
    evolution += np.kron(np.diag([0, 1, 1, 1]), pair)
    blocks = []
    for target in (2, 1):
        scatter = np.eye(64)
        src = slice(0, 16)
        dst = slice(16 * target, 16 * (target + 1))
        scatter[src, src] = eye - projector
        scatter[dst, dst] = eye - swap @ projector @ swap.T
        scatter[dst, src] = swap @ projector
        scatter[src, dst] = projector @ swap.T
        check_equal(scatter @ scatter.T, np.eye(64))
        blocks.append(scatter @ evolution)
    return blocks[0], blocks[1], evolution


def run():
    rotate = np.array([[1, -1], [1, 1]]) / np.sqrt(2)
    loop = controlled_flip("S") @ np.kron(rotate, I) @ controlled_flip("R")
    source = np.kron(I, ZERO[:, None])
    pass_test = (np.eye(4) - np.kron(Y, Y).real) / 2
    ancestor = np.column_stack(([1, 0, 0, 1], [0, -1, 1, 0])) / np.sqrt(2)
    single = []
    expected = [Z, -X, np.zeros((2, 2)), -Z, X, np.zeros((2, 2))]
    for t, sign in enumerate(expected, 1):
        encoded = np.linalg.matrix_power(loop, t) @ source
        effect = encoded.T @ pass_test @ encoded
        check_equal(effect, (I + sign) / 2)
        check_equal(ancestor.T @ np.kron(effect, I) @ ancestor, I / 2)
        single.append(np.round(np.linalg.eigvalsh(effect), 8).tolist())

    early, late, idle = clock_blocks(loop)
    # A fixed unitary is the cyclic pointer shift composed with these blocks.
    ticks = (early, idle, idle, late, idle, idle, idle, idle)
    for tick in ticks:
        check_equal(tick.T @ tick, np.eye(64))
    seeds = np.eye(64)[:, [0, 8]]
    path = seeds.copy()
    propagator = np.eye(64)
    for t in range(4):
        path = ticks[t] @ path
        propagator = ticks[t] @ propagator
        if t == 0:
            blocks = path.reshape(4, 16, 2)
            check_equal(blocks[2].T @ blocks[2], (I + Z) / 2)
            check_equal(blocks[0].T @ blocks[0], (I - Z) / 2)
    branches = path.reshape(4, 16, 2)
    early_effect = branches[2].T @ branches[2]
    late_effect = branches[1].T @ branches[1]
    check_equal(early_effect + late_effect, I)
    check_equal(early_effect - late_effect, Z)
    check_equal(branches[[0, 3]], np.zeros((2, 16, 2)))
    check_equal(propagator.T @ path, seeds)
    # This identity lets the inverse recover arbitrary inputs with unknown partners.
    transfer_x = path @ X @ path.T
    check_equal(path.T @ transfer_x @ path, X)
    signs = np.column_stack(([1, 1], [1, -1])) / np.sqrt(2)
    late_early_states = (path @ signs).reshape(4, 16, 2)
    reduced = [s.T @ s for s in np.moveaxis(late_early_states, 2, 0)]
    check_equal(reduced[0], reduced[1])

    # B goes through its own V once per tick, including both attempts.
    b_map = np.linalg.matrix_power(loop, 4) @ source
    joint = np.kron(path, b_map)
    flag_sign = np.kron(np.diag([0, -1, 1, 0]), np.eye(16))
    b_read = np.linalg.matrix_power(loop, 4) @ np.kron(Z, I) @ np.linalg.matrix_power(loop.T, 4)
    check_equal(joint.T @ np.kron(flag_sign, b_read) @ joint, np.kron(Z, Z))
    check_equal(joint.T @ np.kron(transfer_x, np.eye(4)) @ joint, np.kron(X, I))
    for basis in np.eye(16).reshape(16, 4, 4):
        check_equal(joint.T @ np.kron(np.eye(64), basis) @ joint,
                    np.kron(I, b_map.T @ basis @ b_map))
    # The original root used a special encoding. Reading the birth time on this
    # code leaves B holding the root qubit up to a known correction, even though
    # tracing time cannot recover a freely supplied local |+> or |-> at A.
    routed_ancestor = (joint @ ancestor).reshape(4, 16, 4, 2)
    for sector, bit, correction in ((2, 0, Z), (1, 1, X)):
        child_state = branches[sector, :, bit]
        conditional_b = np.tensordot(child_state.conj(), routed_ancestor[sector],
                                     axes=([0], [0]))
        check_equal(correction @ b_map.T @ conditional_b, I / np.sqrt(2))
    for i in range(4, len(ticks)):
        path = ticks[i] @ path
        branches = path.reshape(4, 16, 2)
        check_equal(branches[[0, 3]], np.zeros((2, 16, 2)))
        pair = np.linalg.matrix_power(np.kron(loop, loop), i - 3)
        check_equal(branches[[1, 2]], pair @ (propagator @ seeds).reshape(4, 16, 2)[[1, 2]])
    return {"single_window_effect_spectra": single,
            "earliest_trial": 1, "retry_trial": 4, "certain_birth_by_trial": 4,
            "ancestor_code_probability_by_birth_time": {"early": 0.5, "late": 0.5},
            "three_local_updates_checked_after_retry": len(ticks) - 4,
            "unperturbed_B_operator_basis_checks": 16,
            "coherent_recovery_input_dimension": 2,
            "root_qubit_recoverable_at_B_given_birth_time": True,
            "timing_trace_identical_for_plus_and_minus": True}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.check:
        if result != json.loads(RESULT.read_text(encoding="utf-8")):
            raise SystemExit("Saved answers differ from recomputation.")
        print("Verified", RESULT.name)
    else:
        RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
