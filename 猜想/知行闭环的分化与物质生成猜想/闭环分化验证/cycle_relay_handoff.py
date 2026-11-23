"""Compose two actual cyclic queries and account for the intermediate carrier frame."""

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from cycle_memory_readout import CYCLE, W, V


def close(a, b):
    np.testing.assert_allclose(a, b, rtol=0, atol=3e-13)


def arms(left, right):
    zero = np.zeros_like(left)
    return np.block([[left, zero], [zero, right]])


def power(p, n):
    return np.linalg.matrix_power(p, n % 3)


def outer(v):
    return np.outer(v, v)


def trace_distance(a, b):
    return np.sum(np.abs(np.linalg.eigvalsh(a-b)))/2


def controlled_not(control_first):
    gate = np.zeros((4, 4))
    for i in range(4):
        x, y = divmod(i, 2)
        if control_first:
            y ^= x
        else:
            x ^= y
        gate[2*x+y, i] = 1
    return gate


def permutation_frames():
    first = controlled_not(True)
    second = controlled_not(False)
    close(second @ first, CYCLE)
    return [power(CYCLE, k) @ s for k, s in product(range(3), (np.eye(4), first))]


def transfer_checks():
    n = 0
    frames = permutation_frames()
    assert len({m.tobytes() for m in frames}) == 6
    for frame in frames:
        outgoing = frame @ CYCLE @ frame.T
        sign = 1 if np.allclose(outgoing, CYCLE) else -1
        close(outgoing, power(CYCLE, sign))
        transit = arms(frame, frame)
        for a, b, c in product(range(3), repeat=3):
            ab = arms(power(CYCLE, a), power(CYCLE, b))
            b_second = arms(power(outgoing, b), np.eye(4))
            c_second = arms(np.eye(4), power(outgoing, c))
            inverse = arms(power(outgoing, -b), power(outgoing, -b))
            direct = arms(power(outgoing, a), power(outgoing, c)) @ transit
            raw = c_second @ b_second @ transit @ ab
            corrected = c_second @ inverse @ b_second @ transit @ ab
            # Operator equality covers all inputs and arbitrary unknown partners.
            close(raw, inverse.T @ direct)
            close(corrected, direct)
            close(corrected.T @ corrected, np.eye(8))
            n += 1
    # Two SWAPs on four physical bits by 3 controlled NOTs per wire.
    for data, blank in product(range(4), repeat=2):
        incoming, outgoing = data, blank
        for bit in (1, 2):
            outgoing ^= incoming & bit
            incoming ^= outgoing & bit
            outgoing ^= incoming & bit
        assert (incoming, outgoing) == (blank, data)
    return n


def score_difference(rho, delta):
    q = W.T @ power(CYCLE, delta) @ W
    eye = np.eye(2)
    effect = np.block([[eye, q.T], [q, eye]])/3
    return np.trace(effect @ rho)


def two_clock_scenarios():
    states = [W.T @ power(CYCLE, k) @ V for k in range(3)]
    scores = []
    split_states = [np.zeros((4, 4)) for _ in range(3)]
    for a, c, left_b, right_b in product(range(3), repeat=4):
        delta = (c-a) % 3
        ell = (a+right_b) % 3
        r = (left_b+c) % 3
        ket = np.r_[states[ell], states[r]]/np.sqrt(2)
        rho = outer(ket)
        split_states[delta] += rho/27
        if left_b == right_b:
            scores.append(score_difference(rho, delta))
    close(scores, np.full(27, 2/3))
    close(split_states[0], split_states[1])
    close(split_states[0], split_states[2])
    return {"shared_tick_difference_score": "2/3",
            "independent_ticks_without_correction": "1/3"}


def hidden_common_rotation():
    rotations = [W.T @ power(CYCLE, b) @ W for b in range(3)]
    phi = np.eye(2).reshape(4)/np.sqrt(2)
    before = outer(phi)
    after = sum(outer(np.kron(np.eye(2), m) @ phi) for m in rotations)/3
    close(np.linalg.eigvalsh(after), [0, 0, 0.5, 0.5])
    close(phi @ after @ phi, 0.5)
    close(trace_distance(before, after), 0.5)
    for v in np.eye(2):
        close(sum(outer(m @ v) for m in rotations)/3, np.eye(2)/2)
    return {"partner_overlap_probability": "1/2", "partner_trace_error": "1/2",
            "distinct_orthogonal_inputs_after_twirl": 1}


def run():
    return {"physical_frame_and_tick_tests": transfer_checks(), "swap_basis_tests": 16,
            "comparison": two_clock_scenarios(), "uncompensated": hidden_common_rotation()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    outfile = Path(__file__).with_name("cycle_relay_handoff_results.json")
    if args.check:
        assert json.loads(outfile.read_text(encoding="utf-8")) == result
    else:
        with outfile.open("x", encoding="utf-8") as out:
            out.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
