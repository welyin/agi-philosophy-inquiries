"""Check hereditary input wiring, local feedback and lineage observables."""

import argparse
import json
from functools import reduce
from itertools import product
from pathlib import Path

import numpy as np

from quantum_loop_inheritance import I, Z, ZERO, controlled_flip


RESULT = Path(__file__).with_name("quantum_recursive_lineage_results.json")
Y = np.array([[0, -1j], [1j, 0]])
Y_PLUS = np.array([1, 1j]) / np.sqrt(2)
Y_MINUS = Y_PLUS.conj()
K = np.kron(ZERO[:, None], I)  # logical input on the return port


def check(observed, expected):
    np.testing.assert_allclose(observed, expected, atol=1e-12, rtol=0)


def tensor(factors):
    return reduce(np.kron, factors)


def all_shapes(leaves):
    if leaves == 1:
        yield None
        return
    for left_size in range(1, leaves):
        for left in all_shapes(left_size):
            for right in all_shapes(leaves - left_size):
                yield (left, right)


def tree_code(shape, birth):
    if shape is None:
        return I, (1,)
    left, left_signs = tree_code(shape[0], birth)
    right, right_signs = tree_code(shape[1], birth)
    return np.kron(left, right) @ birth, left_signs + tuple(-s for s in right_signs)


def after_losing(state, leaves, lost):
    amplitudes = np.moveaxis(state.reshape([2] * leaves), lost, -1).reshape(-1, 2)
    return amplitudes @ amplitudes.conj().T


def run():
    rotate = np.array([[1, -1], [1, 1]]) / np.sqrt(2)
    cut = controlled_flip("S") @ np.kron(rotate, I)
    cycle = cut @ controlled_flip("R")
    trigger = -np.kron(Y, Y)
    birth = cycle @ K
    check(trigger @ birth, birth)
    check(birth @ Y_PLUS, np.kron(Y_PLUS, Y_MINUS))
    check(birth @ Y_MINUS, np.kron(Y_MINUS, Y_PLUS))

    # Source order S_A,R_A,S_C,R_C. Transfer R_A to S_C and then
    # exchange action and return locally at each newborn, as in the proposal.
    blank = np.kron(ZERO, ZERO)[:, None]
    original = np.kron(birth, blank).reshape(2, 2, 2, 2, 2)
    routed = original.swapaxes(1, 2)
    rewired = routed.swapaxes(0, 1).swapaxes(2, 3).reshape(16, 2)
    check(rewired, np.kron(K, K) @ birth)

    shape_counts = {}
    loss_tests = 0
    for leaves in range(2, 7):
        count = 0
        for shape in all_shapes(leaves):
            count += 1
            code, signs = tree_code(shape, birth)
            check(code.conj().T @ code, I)
            for root_sign, root_state in ((1, Y_PLUS), (-1, Y_MINUS)):
                expected = tensor([Y_PLUS if root_sign * s == 1 else Y_MINUS for s in signs])
                check(code @ root_state, expected)
            check(tensor([Z] * leaves) @ code, code @ Z)
            for lost in range(leaves):
                check(after_losing(code[:, 0], leaves, lost),
                      after_losing(code[:, 1], leaves, lost))
                overlap = np.moveaxis((code @ Y_PLUS).reshape([2] * leaves), lost, -1).reshape(-1, 2)
                other = np.moveaxis((code @ Y_MINUS).reshape([2] * leaves), lost, -1).reshape(-1, 2)
                check(overlap @ other.conj().T, 0)
                loss_tests += 1
            # Joint decoding of a logical qubit still entangled with a reference.
            reference_input = np.array([1, 0, 0, 1]) / np.sqrt(2)
            encoded = np.kron(I, code) @ reference_input
            check(np.kron(I, code.conj().T) @ encoded, reference_input)
        shape_counts[str(leaves)] = count

    # During A's second birth, B is a spectator (in its freely evolving frame).
    split_a = np.kron(birth, I)
    for row, col in product(range(4), repeat=2):
        unit = np.zeros((4, 4))
        unit[row, col] = 1
        after = (split_a @ unit @ split_a.conj().T).reshape(4, 2, 4, 2)
        check(after.trace(axis1=0, axis2=2), unit.reshape(2, 2, 2, 2).trace(axis1=0, axis2=2))

    # At local birth frames, all three loops have the same known action input.
    # Each logical return has zero Z expectation on the whole ancestral code.
    code3 = split_a @ birth
    cut_results = {}
    for cuts in product(range(2), repeat=3):
        frame = tensor([(cut if removed else cycle) @ K for removed in cuts]) @ code3
        for i, removed in enumerate(cuts):
            observable = tensor([trigger if j == i else np.eye(4) for j in range(3)])
            check(frame.conj().T @ observable @ frame, 0 if removed else I)
        cut_results["".join(map(str, cuts))] = [1 - bit for bit in cuts]

    return {
        "hereditary_input_frame_verified": True,
        "ordered_binary_tree_counts_by_leaves": shape_counts,
        "one_leaf_loss_tests": loss_tests,
        "joint_logical_Z_is_leaf_Z_parity": True,
        "sibling_matrix_units_tested": 16,
        "local_birth_frame_trigger_expectations_with_cuts": cut_results,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    results = run()
    if args.check:
        if json.loads(RESULT.read_text(encoding="utf-8")) != results:
            raise SystemExit("Saved results differ from recomputation.")
        print("Verified", RESULT.name)
    else:
        RESULT.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
