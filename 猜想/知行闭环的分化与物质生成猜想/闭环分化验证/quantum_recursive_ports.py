"""Reversible four-cell test of birth-port alignment and a second split."""

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from quantum_loop_inheritance import I, Z, controlled_flip
from quantum_split_router import check_equal


RESULT = Path(__file__).with_name("quantum_recursive_ports_results.json")
# Three activity bits f_B,f_C,f_D, four cells A,B,C,D, and logical input columns.
SHAPE = (2, 2, 2, 4, 4, 4, 4, 2)
CLOCK = 16
A, B, C, D = range(3, 7)
BASIS = np.eye(4)
READY = BASIS[:, :2]  # |0,x> for logical x
ROTATION = np.array([[1, -1], [1, 1]]) / np.sqrt(2)
V = controlled_flip("S") @ np.kron(ROTATION, I) @ controlled_flip("R")
Z_RETURN = np.kron(I, Z)
ENCODE = np.array([[1, 0], [0, -1], [0, 1], [1, 0]]) / np.sqrt(2)


def apply_axes(state, gate, axes):
    """Apply a matrix on the named factors and keep all other indices."""
    new_axes = tuple(range(len(axes)))
    moved = np.moveaxis(state, axes, new_axes)
    flat = moved.reshape(gate.shape[1], -1)
    updated = (gate @ flat).reshape(moved.shape)
    return np.moveaxis(updated, new_axes, axes)


def controlled(gate):
    identity = np.eye(len(gate))
    zeros = np.zeros_like(gate)
    return np.block([[identity, zeros], [zeros, gate]])


def reversible_router():
    # W from note 38, followed by an exchange of action/return in EACH child.
    identity = np.eye(16)
    w = identity.reshape(2, 2, 2, 2, 16).swapaxes(1, 2).reshape(16, 16)
    swap_sr = BASIS[[0, 2, 1, 3]]
    j = np.kron(swap_sr, swap_sr) @ w
    y = np.array([[0, -1j], [1j, 0]])
    pred = (BASIS - np.kron(y, y).real) / 2
    source = np.kron(pred, np.diag([1, 0, 0, 0]))
    arrival = j @ source @ j.T
    router = np.block([[identity - source, source @ j.T],
                       [j @ source, identity - arrival]])
    check_equal(router.T @ router, np.eye(32))
    return router, pred


ROUTER, ACCEPT = reversible_router()
ACTIVE_UPDATE = controlled(V)
ACTIVE_B_ROUTER = controlled(ROUTER)


def tick(state, position, second_at, alter_a_return=False):
    """A clock sector maps to its successor; no external activity decisions."""
    state = apply_axes(state, V, (A,))
    for flag, cell in enumerate((B, C, D)):
        state = apply_axes(state, ACTIVE_UPDATE, (flag, cell))
    if position == 0:
        state = apply_axes(state, ROUTER, (0, A, B))
    if position == second_at:
        if alter_a_return:  # Only for the counterfactual protocol.
            state = apply_axes(state, Z_RETURN, (A,))
        state = apply_axes(state, ROUTER, (1, A, C))
        state = apply_axes(state, ACTIVE_B_ROUTER, (0, 2, B, D))
    return state, (position + 1) % CLOCK


def seed(c_value, d_value):
    state = np.zeros(SHAPE)
    for x in (0, 1):
        state[0, 0, 0, x, 0, c_value, d_value, x] = 1
    return state


def pair_code(c_value, d_value, steps):
    """Independent prediction from the logical isometry, before the second port."""
    local = np.linalg.matrix_power(V, steps) @ READY
    output = np.zeros(SHAPE)
    for a, b in product(range(2), repeat=2):
        output[1, 0, 0, :, :, c_value, d_value, :] += np.einsum(
            "i,j,x->ijx", local[:, a], local[:, b], ENCODE[2*a + b])
    return output


def local_pair_output(bit, spare, altered):
    # V^7 = -V on each parent. Both minus signs cancel below.
    if spare != 0 or altered:
        active = V @ READY[:, bit]
        if altered:
            active = Z_RETURN @ active
        return np.outer(active, BASIS[:, spare])
    return (np.kron(READY, READY) @ ENCODE[:, bit]).reshape(4, 4)


def expected_seventh(c_value, d_value, alter_a=False):
    """Product of the two local branch maps on an entangled ancestor."""
    output = np.zeros(SHAPE)
    flags = (1, int(c_value == 0 and not alter_a), int(d_value == 0))
    for a, b in product(range(2), repeat=2):
        ac = local_pair_output(a, c_value, alter_a)
        bd = local_pair_output(b, d_value, False)
        output[flags] += np.einsum("ik,jl,x->ijklx", ac, bd, ENCODE[2*a+b])
    return output


def advance_to_window(c_value, d_value, second_at):
    out, pos = tick(seed(c_value, d_value), 0, second_at)
    check_equal(out, pair_code(c_value, d_value, 0))
    for t in range(1, second_at):
        out, pos = tick(out, pos, second_at)
        check_equal(out, pair_code(c_value, d_value, t))
    assert pos == second_at
    return out, pos


def verify():
    for gate in (V, ACTIVE_UPDATE, ROUTER, ACTIVE_B_ROUTER):
        check_equal(gate.T @ gate, np.eye(len(gate)))
    check_equal(np.linalg.matrix_power(V, 6), -BASIS)

    # Checking every spare basis column with its phase also covers coherent
    # superpositions and external references by linearity.
    for c_value, d_value in product(range(4), repeat=2):
        entered, pos = advance_to_window(c_value, d_value, 7)
        complete, _ = tick(entered, pos, 7)
        check_equal(complete, expected_seventh(c_value, d_value))
        matrix = complete.reshape(-1, 2)
        check_equal(matrix.T @ matrix, I)
    entered, pos = advance_to_window(0, 2, 7)
    altered, _ = tick(entered, pos, 7, alter_a_return=True)
    check_equal(altered, expected_seventh(0, 2, alter_a=True))

    # Both spares blank. At child tick two a different return subspace
    # meets the SAME gate. Propagate all four orthogonal growth flags coherently.
    twice = np.linalg.matrix_power(V, 2) @ READY
    check_equal(twice.T @ ACCEPT @ twice, I/2)
    entered, pos = advance_to_window(0, 0, 2)
    coherent, _ = tick(entered, pos, 2)
    check_equal(coherent[0], np.zeros_like(coherent[0]))
    for c_on, d_on in product(range(2), repeat=2):
        branch = coherent[1, c_on, d_on].reshape(-1, 2)
        check_equal(branch.T @ branch, I/4)
    check_equal(coherent.reshape(-1, 2).T @ coherent.reshape(-1, 2), I)
    return {
        "clock_positions": CLOCK,
        "full_unitary_dimension_by_factorization": CLOCK * np.prod(SHAPE[:-1]).item(),
        "spare_basis_configurations_checked": 16,
        "logical_code_dimension": 2,
        "child_ticks_until_second_window": 7,
        "active_ABCD_for_C_blank_D_occupied": [1, 1, 1, 0],
        "active_ABCD_for_C_occupied_D_blank": [1, 1, 0, 1],
        "active_ABCD_if_both_blank": [1, 1, 1, 1],
        "active_ABCD_after_A_return_phase_flip_with_D_occupied": [1, 1, 0, 0],
        "window_at_child_tick_two_flag_probabilities_00_01_10_11": [0.25]*4,
        "all_unitary_factors_and_logical_branch_grams_within_1e-12": True,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = verify()
    if args.check:
        if json.loads(RESULT.read_text(encoding="utf-8")) != result:
            raise SystemExit("Saved answers differ from recomputation.")
        print("Verified", RESULT.name)
    else:
        RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
