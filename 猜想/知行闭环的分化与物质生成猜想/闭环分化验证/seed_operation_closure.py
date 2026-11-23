"""Audit the real working gates and test phase orientation at the actual split port."""

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from quantum_consensus_dispatch import CNOT, SWAP, local_sample
from quantum_living_tree_recovery import leaf_correction
from quantum_loop_inheritance import I, X, Z, PLUS, ZERO
from quantum_recursive_ports import ACCEPT, ENCODE, READY, ROUTER, V, controlled


RESULT = Path(__file__).with_name("seed_operation_closure_results.json")
Y = np.array([[0, -1j], [1j, 0]])
S = np.diag([1, 1j])


def check(actual, expected):
    np.testing.assert_allclose(actual, expected, atol=1e-12, rtol=0)


def density(psi):
    return np.outer(psi, psi.conj())


def trace_distance(a, b):
    return float(np.abs(np.linalg.eigvalsh(a - b)).sum() / 2)


def transcripts(oracle, rounds):
    """Keep unnormalized amplitudes so every branch retains its probability."""
    tree = {"": np.kron(ZERO, ZERO)}
    readouts = tuple(np.kron(I, np.diag(column)) for column in I)
    for turn in range(rounds):
        children = {}
        for history, state in tree.items():
            phase = np.kron(I, oracle) if turn % 2 else np.kron(oracle, I)
            feedback = np.kron(X, I) if history.endswith("1") else np.eye(4)
            before_readout = V @ phase @ feedback @ V @ state
            for outcome, projection in enumerate(readouts):
                children[history + str(outcome)] = projection @ before_readout
        tree = children
    return tree


def split_flag(first_sign, second_sign):
    gates = {1: S, -1: S.conj()}
    witnesses = np.kron(gates[first_sign] @ PLUS, gates[second_sign] @ PLUS)
    # Five ticks before the normal birth tick; the six give a global minus sign.
    source = np.linalg.matrix_power(V, 6) @ witnesses
    blank_spare = np.kron(ZERO, ZERO)
    arrival = ROUTER @ np.kron(ZERO, np.kron(source, blank_spare))
    return np.sum(np.abs(arrival.reshape(2, -1))**2, axis=1)


def run():
    iy = 1j * Y
    check(iy.imag, 0)
    check(ACCEPT, (np.eye(4) - np.kron(Y, Y)) / 2)
    checked = {
        "V": V, "V_inverse": V.conj().T, "split": ROUTER, "split_accept": ACCEPT,
        "controlled_split": controlled(ROUTER), "swap": SWAP, "CNOT": CNOT,
        "conditional_return_iY": iy, "aligned_return_sampler": local_sample(),
        "birth_isometry": ENCODE, "ready_port": READY,
        "living_tree_correction_0": leaf_correction(0),
        "living_tree_correction_1": leaf_correction(1),
    }
    for matrix in checked.values():
        check(np.asarray(matrix).imag, 0)
    check(np.linalg.matrix_power(V, 6), -np.eye(4))
    bell = ENCODE[:, 0].reshape(2, 2)
    marginal = bell @ bell.conj().T
    check(marginal, I / 2)
    y_plus = S @ PLUS
    target = density(y_plus)
    check(trace_distance(marginal, target), 0.5)
    for real_state in (marginal, density(ZERO), density(PLUS)):
        check(np.trace(real_state @ target), 0.5)
    check(np.trace(density(S.conj() @ PLUS) @ target), 0)

    positive = transcripts(S, 5)
    negative = transcripts(S.conj(), 5)
    for history, state in positive.items():
        check(state.conj(), negative[history])
    probs = np.array([np.vdot(v, v).real for v in positive.values()])
    check(probs.sum(), 1)
    assert np.count_nonzero(probs > 1e-12) > 1
    flags = {}
    for left, right in product((1, -1), repeat=2):
        expected_flag = int(left != right)
        check(split_flag(left, right), I[expected_flag])
        flags[f"{left:+d},{right:+d}"] = expected_flag
    return {
        "audited_real_matrices": list(checked),
        "bare_feedback_ticks_returning_minus_identity": 6,
        "minimum_native_qubit_trace_distance_to_y_plus": 0.5,
        "max_squared_fidelity_from_real_native_states": 0.5,
        "equal_prior_single_oracle_orientation_success_with_real_resources": 0.5,
        "audit_oracle_queries_per_trial": 5,
        "tested_adaptive_record_histories": len(positive),
        "split_activity_flag_for_two_supplied_phase_signs": flags,
        "equal_prior_single_oracle_success_given_calibrated_Y_readout": 1,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.check:
        if json.loads(RESULT.read_text(encoding="utf-8")) != result:
            raise SystemExit("Saved results differ from recomputation.")
        print("Verified seed operation closure and four relative-orientation inputs.")
    else:
        RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
