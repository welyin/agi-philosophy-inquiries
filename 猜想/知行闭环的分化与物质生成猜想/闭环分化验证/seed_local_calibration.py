"""Use the loop's own return and its actual birth flag to test local calibration."""

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from quantum_consensus_dispatch import equal, local_sample
from quantum_recursive_ports import I, Z, READY, ROUTER, V


RESULT = Path(__file__).with_name("seed_local_calibration_results.json")
X = I[::-1]
Y = np.array([[0, -1j], [1j, 0]])
J = (1j * Y).real


def coin_branches():
    # The sampling loop has its own blank and one blank witness. With logical
    # input zero it ends in |00> at the parent and a coherent fair return.
    actual = local_sample()[:, :, 0]
    equal(actual[1:], 0)
    return actual[0]


def actions(coin):
    first = [coin[0] * I, coin[1] * J]
    return [first, [Z @ gate for gate in first]]


def channel(operators, state):
    return sum(gate @ state @ gate.conj().T for gate in operators)


def run():
    coin = coin_branches()
    private_actions = actions(coin)
    bell = (V @ READY)[:, 0]
    input_pair = np.outer(bell, bell)
    blank = np.eye(4)[:, 0]
    basis = [I, X, Z]
    for operators in private_actions:
        equal(sum(g.T @ g for g in operators), I)
        for observable in basis:
            equal(channel(operators, observable), np.trace(observable) * I / 2)

    pair_outputs = []
    for index, operators in enumerate(private_actions):
        pair = channel([np.kron(g, I) for g in operators], input_pair)
        equal(pair, (np.eye(4) + (-1)**(index + 1) * np.kron(Y, Y)) / 4)
        pair_outputs.append(pair)
        before = np.kron(np.diag([1, 0]), np.kron(pair, np.outer(blank, blank)))
        after = ROUTER @ before @ ROUTER.T
        flag = [np.trace(block).real for block in (after[:16, :16], after[16:, 16:])]
        equal(flag, [index, 1 - index])

    difference = pair_outputs[0] - pair_outputs[1]
    for left, right in product(basis, repeat=2):
        equal(np.trace(difference @ np.kron(left, right)), 0)
    equal(np.linalg.svd(difference, compute_uv=False).sum() / 2, 1)
    y_plus = (I + Y) / 2
    y_positive = [np.trace(channel(gates, y_plus) @ y_plus).real for gates in private_actions]
    equal(y_positive, [1, 0])
    return {
        "loop_return_coin_probabilities": [0.5, 0.5],
        "real_isolated_operator_basis_checked": ["I", "X", "Z"],
        "real_product_directions_with_identical_statistics": 9,
        "input_pair_source": "first V from the parent's blank |00>",
        "original_router_birth_probabilities": [1, 0],
        "trace_distance_of_pair_outputs": 1,
        "local_Y_plus_readings_with_added_Y_preparation_and_readout": [1, 0],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.check:
        if json.loads(RESULT.read_text(encoding="utf-8")) != result:
            raise SystemExit("Saved calibration results differ.")
        print("Verified", RESULT.name)
    else:
        RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
