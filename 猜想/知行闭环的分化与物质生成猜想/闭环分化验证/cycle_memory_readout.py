"""Two-route comparison of three-cycle counts with a rewound readout."""
import argparse
import json
from pathlib import Path

import numpy as np

H = np.array([[1., 1.], [1., -1.]]) / np.sqrt(2)
I4 = np.eye(4)
CYCLE = np.zeros((4, 4))
for x in range(2):
    for y in range(2):
        CYCLE[2*y+(x ^ y), 2*x+y] = 1
POWERS = [I4, CYCLE, CYCLE.T]
V = (I4[:, 1] - I4[:, 3]) / np.sqrt(2)
W = np.column_stack((V, (I4[:, 1] + I4[:, 3] - 2*I4[:, 2]) / np.sqrt(6)))
CODE = np.kron(np.eye(2), W)


def circuit():
    """Gate on the work address, path and two-bit memory; all columns retained."""
    gates = np.eye(32)
    for k in range(4):
        controlled = np.eye(8)
        if k < 3:
            controlled[4:, 4:] = POWERS[k].T
        gates[8*k:8*k+8, 8*k:8*k+8] = np.kron(H, I4) @ controlled
    return gates @ np.kron(np.kron(H, H), np.eye(8))


def code_outputs():
    states = {}
    for a in range(3):
        for b in range(3):
            states[a, b] = np.concatenate((POWERS[a] @ V, POWERS[b] @ V)) / np.sqrt(2)
    return states


def matrix_checks():
    u = circuit()
    states = code_outputs()
    start = np.zeros((32, 8))
    start[:8] = np.eye(8)
    encoded = u @ start @ CODE
    initial = u @ start @ start.T @ u.T
    outcomes = []
    for k in range(3):
        projector = np.zeros((32, 32))
        projector[8*k:8*k+4, 8*k:8*k+4] = I4
        outcomes.append(projector)
    good = sum(outcomes)
    bad = np.eye(32) - good
    s = 3/8
    np.testing.assert_allclose(encoded.T @ good @ encoded, s*np.eye(4), atol=1e-14)
    np.testing.assert_allclose(initial @ encoded, encoded, atol=1e-14)
    # A good component after a miss and either rewind answer is a scalar
    # multiple of the original one, simultaneously on the entire code.
    for project, factor in ((initial, 1-s), (np.eye(32)-initial, -(1-s))):
        np.testing.assert_allclose(good @ project @ bad @ encoded,
                                   factor*(good @ encoded), atol=1e-14)
    for project, factor in ((initial, 1-s), (np.eye(32)-initial, s)):
        np.testing.assert_allclose(bad @ project @ bad @ encoded,
                                   factor*(bad @ encoded), atol=1e-14)
    retry = 2*s*(1-s)
    recovered_good = sum(encoded.T @ bad @ proj @ good @ proj @ bad @ encoded
                         for proj in (initial, np.eye(32)-initial))
    np.testing.assert_allclose(recovered_good, (1-s)*retry*np.eye(4), atol=1e-14)
    for (a, b), state in states.items():
        delta = (b-a) % 3
        prepared = u @ start @ state
        hits = [np.linalg.norm(proj @ prepared)**2 for proj in outcomes]
        np.testing.assert_allclose(np.array(hits)/s,
                                   [2/3 if k == delta else 1/6 for k in range(3)])
        x_prob = 4*np.array(hits)
        x_only = (x_prob[delta] + (1-x_prob).sum()/2) / 3
        assert np.isclose(x_only, 7/12)
    for delta in range(3):
        density = sum(np.outer(states[a, (a+delta) % 3], states[a, (a+delta) % 3])
                      for a in range(3))/3
        eig = np.linalg.eigvalsh(CODE.T @ density @ CODE)
        np.testing.assert_allclose(eig, [0, 0, 0.5, 0.5], atol=1e-14)
    return s, retry


def run():
    np.testing.assert_array_equal(CYCLE @ CYCLE @ CYCLE, I4)
    np.testing.assert_allclose(W.T @ (I4 + CYCLE + CYCLE.T) @ W, 0, atol=1e-14)
    j_full = (CYCLE-CYCLE.T) / np.sqrt(3)
    j = W.T @ j_full @ W
    np.testing.assert_allclose(j @ j, -np.eye(2), atol=1e-14)
    np.testing.assert_allclose(j.T, -j, atol=1e-14)
    np.testing.assert_allclose(W.T @ CYCLE @ W, (-np.eye(2)+np.sqrt(3)*j)/2, atol=1e-14)
    center = (I4[:, 1] + I4[:, 2] + I4[:, 3]) / 3
    quarter_turn = center + j_full @ (I4[:, 1] - center)
    assert np.isclose(min(quarter_turn), 1/3-1/np.sqrt(3))
    s, retry = matrix_checks()
    finite = []
    for rounds in (1, 2, 4, 8, 16):
        unfinished = (1-s)*(1-retry)**(rounds-1)
        score = (2-unfinished)/3  # on a timeout guess a uniform delta
        finite.append({"tests": rounds, "timeout": round(unfinished, 10),
                       "guaranteed_score": round(score, 10)})
    return {"response_settings": 9, "real_effective_dimension": 4,
            "classical_triangle_quarter_turn_minimum": round(float(min(quarter_turn)), 10),
            "one_binary_comparison": "7/12", "full_readout": "2/3",
            "first_good": "3/8", "next_good_after_bad": "15/32",
            "mean_good_tests": "7/3", "finite_caps": finite}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    output = Path(__file__).with_name("cycle_memory_readout_results.json")
    if args.check:
        assert json.loads(output.read_text(encoding="utf-8")) == result
    else:
        with output.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
