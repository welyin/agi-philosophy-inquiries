"""Three actual return samples control a fourth birth through two local relays."""

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from quantum_recursive_ports import (
    ENCODE, READY, ROUTER, V, advance_to_window, apply_axes, tick,
)


RESULT = Path(__file__).with_name("quantum_consensus_dispatch_results.json")
BITS = tuple(product(range(2), repeat=3))
IDENTITY = np.eye(2)
X = IDENTITY[::-1]
CNOT = np.eye(4)[[0, 1, 3, 2]]
SWAP = np.eye(4)[[0, 2, 1, 3]]


def equal(value, expected):
    np.testing.assert_allclose(value, expected, rtol=0, atol=1e-12)


def local_sample():
    """Physical isometry on S,R,M with M blank. Return its columns on K."""
    # Control iY on S by the physical R immediately after V.
    iy = np.array([[0, 1], [-1, 0]])
    control_return = np.kron(IDENTITY, np.diag([1, 0])) + np.kron(iy, np.diag([0, 1]))
    full = np.kron(control_return @ V, IDENTITY).reshape(2, 2, 2, 8)
    full = apply_axes(full, SWAP, (1, 2))
    full = apply_axes(full, SWAP, (0, 1)).reshape(8, 8)
    equal(full.T @ full, np.eye(8))
    columns = full[:, [0, 2]]
    equal(columns, np.kron(READY, np.ones((2, 1)) / np.sqrt(2)))
    return columns.reshape(4, 2, 2)


def inherited_state():
    entered, clock = advance_to_window(0, 2, 7)
    actual, _ = tick(entered, clock, 7)
    # Axes now A,C,B,root. The occupied spectator retained its old record.
    inherited = actual[1, 1, 0, :, :, :, 2, :].transpose(0, 2, 1, 3)
    inherited = apply_axes(inherited, V.T, (2,))
    code = np.einsum("act,tbx->acbx", ENCODE.reshape(2, 2, 2), ENCODE.reshape(2, 2, 2))
    equal(inherited, np.einsum("ia,jc,kb,acbx->ijkx", READY, READY, READY, code))
    return inherited[:2, :2, :2, :], code


def add_message(state, site, courier):
    # The stored logical Z is now R. The sample occupies a separate local M.
    # A two-bit cell is the control for the first gate, using its low (R) bit.
    control_r = np.kron(IDENTITY, CNOT)
    state = apply_axes(state, control_r, (site, courier))
    return apply_axes(state, CNOT, (site + 3, courier))


def feedback_interval(state):
    # Six actual ticks on every cell; the separate evidence and cables persist.
    for _ in range(6):
        for site in range(3):
            state = apply_axes(state, V, (site,))
    return state


def route():
    inherited, code = inherited_state()
    sample = local_sample()
    prepared = np.einsum("ima,jnc,kob,acbx->ijkmnox", sample, sample, sample, inherited)
    state = np.zeros((4, 4, 4, 2, 2, 2, 2, 2, 2))
    state[:, :, :, :, :, :, 0, 0, :] = prepared
    # Courier 6 travels A->C->B and 7 travels B->C->A.
    state = add_message(state, 0, 6)
    state = add_message(state, 2, 7)
    state = feedback_interval(state)
    state = add_message(state, 1, 6)
    state = add_message(state, 1, 7)
    state = feedback_interval(state)
    state = add_message(state, 2, 6)
    state = add_message(state, 0, 7)
    expected = np.zeros_like(state)
    for z in BITS:
        for e in BITS:
            g = (sum(z) + sum(e)) % 2
            expected[z + e + (g, g)] = code[z] / np.sqrt(8)
    equal(state, expected)
    return state, code


def birth_map(open_port, spare=0):
    # Rows are activity flag, parent S/R, spare S/R. Include the real router.
    before = np.zeros((2, 4, 4, 4))
    before[0, :, spare, :] = V
    gate = ROUTER if open_port else np.eye(32)
    return gate @ before.reshape(32, 4)


def verify():
    state, code = route()
    ports = [birth_map(False), birth_map(True)]
    # Opening the contact cannot erase a prior occupied receiver.
    for spare in (1, 2, 3):
        entered = birth_map(True, spare) @ READY
        equal(entered, birth_map(False, spare) @ READY)
    expected_split = np.kron(READY, READY) @ ENCODE
    equal(ports[1][16:] @ READY, expected_split)
    equal(ports[1][:16] @ READY, np.zeros_like(expected_split))
    final = np.zeros((8, 2, 32, 4, 32, 2))
    for index, e in enumerate(BITS):
        for g in range(2):
            amplitudes = state[(slice(None),) * 3 + e + (g, g)]
            block = np.einsum("ai,cj,bk,ijkx->acbx", ports[g == 0], V,
                              ports[g == 1], amplitudes, optimize=True)
            final[index, g] = block
            flags = block.reshape(2, 16, 4, 2, 16, 2)
            for fa, fb in product(range(2), repeat=2):
                rows = flags[fa, :, :, fb].reshape(-1, 2)
                expected = np.zeros((2, 2))
                x = (g + sum(e)) % 2
                if (fa, fb) == (1 - g, g):
                    expected[x, x] = 1 / 8
                equal(rows.T @ rows, expected)
    matrix = final.reshape(-1, 2)
    equal(matrix.T @ matrix, IDENTITY)
    for g in range(2):
        sub = final[:, g].reshape(-1, 2)
        equal(sub.T @ sub, IDENTITY / 2)
    # Explicit external reference with a non-real relative phase.
    partner = np.diag([1, 1j]) / np.sqrt(2)
    recovered = matrix.T @ (matrix @ partner)
    equal(recovered, partner)
    phase_dephased = sum(
        np.abs(matrix[:, x]) @ np.abs(matrix[:, 1 - x]) for x in range(2)
    )
    equal(phase_dephased, 0)
    return {
        "leaf_order": ["A", "C", "B"],
        "initial_history": "note 40 j=7 A admits blank C; old D is occupied",
        "B_frame_alignment_after_j7": "one local V dagger",
        "fresh_return_source": "physical R after the next local V; extract with controlled iY and two swaps",
        "sample_bits": 3, "zero_initialized_sample_qubits": 3,
        "zero_initialized_courier_qubits": 2,
        "communication_path": ["A-C-B", "B-C-A"],
        "directed_single_qubit_transmissions": 4, "parallel_edge_rounds": 2,
        "local_V_ticks_per_communication_round": 6,
        "birth_window_after_extract_local_V_tick": 13,
        "candidate_spares": 2, "blank_qubits_per_candidate_spare": 2,
        "candidate_activity_flags": 2,
        "checked_root_Z_and_evidence_configurations": 16,
        "each_sample_triple_probability_for_any_root": 0.125,
        "ancestor_and_evidence_rule": "g=x xor e_A xor e_C xor e_B",
        "flags_A_B_at_g0_and_g1": [[1, 0], [0, 1]],
        "marginal_gate_probabilities_for_any_root": [0.5, 0.5],
        "occupied_spare_basis_cases_checked": 3,
        "both_root_columns_relayed_with_phases_verified": True,
        "complete_output_isometry_and_inverse_on_code_with_partner_verified": True,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = verify()
    if args.check:
        if json.loads(RESULT.read_text(encoding="utf-8")) != result:
            raise SystemExit("Saved results differ from recomputation.")
        print("Verified", RESULT.name)
    else:
        RESULT.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
