"""Join two independently supplied real orientation carriers with the birth router."""

import argparse
import json
from pathlib import Path

import numpy as np

from quantum_recursive_ports import ROUTER, apply_axes, controlled
from public_choice_quantum_return import CNOT, I, Z, equal


RESULT = Path(__file__).with_name("seed_orientation_join_results.json")
J = np.array([[0, 1], [-1, 0]])
JA, JB = np.kron(J, I), np.kron(I, J)
SAME = (np.eye(4) - JA @ JB) / 2
OPPOSITE = np.eye(4) - SAME
FLIP_B = np.kron(I, Z)


def scatter_write_unscatter(state, inverse=False):
    if inverse:
        state = apply_axes(state, controlled(FLIP_B), (3, 1))
    state = apply_axes(state, ROUTER, (0, 1, 2))
    state = apply_axes(state, CNOT, (0, 3))
    state = apply_axes(state, ROUTER.T, (0, 1, 2))
    if not inverse:
        state = apply_axes(state, controlled(FLIP_B), (3, 1))
    return state


def verify():
    # Axes: activity flag, both input carriers, spare cell, new memory, input columns.
    seed = np.zeros((2, 4, 4, 2, 4))
    seed[0, :, 0, 0] = np.eye(4)
    output = scatter_write_unscatter(seed)
    expected = np.zeros_like(seed)
    expected[0, :, 0, 0] = SAME
    expected[0, :, 0, 1] = FLIP_B @ OPPOSITE
    equal(output, expected)
    equal(scatter_write_unscatter(output, inverse=True), seed)
    encoding = output[0, :, 0].reshape(8, 4)
    equal(encoding.T @ encoding, np.eye(4))
    equal(encoding @ encoding.T, np.kron(SAME, I))
    equal(encoding @ JA, np.kron(JA, I) @ encoding)
    equal(encoding @ JB, np.kron(JB, Z) @ encoding)
    equal(np.kron(JA-JB, I) @ encoding, 0)

    # Take a genuine inaccessible orthogonal copy of m, remove that copy,
    # and undo the original gates on the accessible output, branch by branch.
    omega = np.eye(4).ravel() / 2
    post = np.zeros((16, 16))
    for bit, sector in enumerate((SAME, OPPOSITE)):
        branch = np.zeros_like(output)
        branch[:, :, :, bit] = output[:, :, :, bit]
        decoded = scatter_write_unscatter(branch, inverse=True)
        target = np.zeros_like(seed)
        target[0, :, 0, 0] = sector
        equal(decoded, target)
        branch_bell = sector.ravel() / 2
        post += np.outer(branch_bell, branch_bell)
        equal(branch_bell @ branch_bell, 0.5)
        equal(omega @ np.kron(np.eye(4), sector.T) @ omega, 0.5)
    equal(omega @ post @ omega, 0.5)
    return {
        "input_real_dimension": 4,
        "orientation_ranks": [int(np.trace(SAME)), int(np.trace(OPPOSITE))],
        "router_applications_in_forward_protocol": 2,
        "blank_two_bit_spare_and_activity_flag_returned": True,
        "new_coherent_memory_bits": 1,
        "operator_columns_reversed_exactly": 4,
        "A_orientation_sign_for_m_0_1": [1, 1],
        "B_orientation_sign_for_m_0_1": [1, -1],
        "four_dimensional_bell_squared_fidelity_after_orientation_bit_leak": 0.5,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = verify()
    if args.check:
        if result != json.loads(RESULT.read_text(encoding="utf-8")):
            raise SystemExit("Saved results differ from recomputation.")
        print("Verified", RESULT.name)
    else:
        RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
