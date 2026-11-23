"""One returned route carrying two independent cyclic marks; scheduled resource injection."""

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np


OMEGA = np.exp(2j * np.pi / 3)
H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
X = np.array([[0, 1], [1, 0]])
Z = np.diag([1, -1])
S = np.diag([1, 1j])
PLUS = np.ones(2) / np.sqrt(2)


def project(ket):
    return np.outer(ket, ket.conj())


def real_operator(matrix):
    return np.block([[matrix.real, -matrix.imag], [matrix.imag, matrix.real]])


def close(left, right):
    np.testing.assert_allclose(left, right, atol=2e-12, rtol=0)


def cyclic_counter(state):
    a, b = state
    b ^= a
    a ^= b
    return a, b


def score_characters():
    strata = {}
    for left, right in product(range(3), repeat=2):
        states = [np.zeros((2, 2), dtype=complex) for _ in range(3)]
        for a, b in product(range(3), repeat=2):
            angle = OMEGA ** (left*a + right*b)
            ket = np.array([1, angle]) / np.sqrt(2)
            states[(b-a) % 3] += project(ket) / 3
        if (left + right) % 3 == 0 and right != 0:
            effects = [2*state/3 for state in states]
            close(sum(effects), np.eye(2))
            # Equal-prior dual certificate: I/3 dominates every prior-weighted state.
            for density in states:
                assert np.linalg.eigvalsh((np.eye(2)-density)/3).min() >= -2e-12
        else:
            close(states[0], states[1])
            close(states[0], states[2])
            effects = [np.eye(2)/3 for _ in states]
        score = sum(np.trace(e @ state).real for e, state in zip(effects, states)) / 3
        expected = 2/3 if right and (left+right) % 3 == 0 else 1/3
        close(score, expected)
        strata[f"{left},{right}"] = round(score, 10)
    return strata


def four_real_dimensions():
    template = [project(np.array([1, OMEGA**i]) / np.sqrt(2)) for i in range(3)]
    effects = [real_operator(2*state/3) for state in template]
    close(sum(effects), np.eye(4))
    probabilities = []
    for a, b in product(range(3), repeat=2):
        matrix = np.diag([OMEGA**a, OMEGA**b])
        real_state = real_operator(matrix) @ np.r_[PLUS, 0*PLUS]
        probs = [real_state @ effect @ real_state for effect in effects]
        expected = np.full(3, 1/6)
        expected[(b-a) % 3] = 2/3
        close(probs, expected)
        probabilities.append(probs[(b-a) % 3])
    return round(float(np.mean(probabilities)), 10)


def fixed_readout():
    # Randomly cycle the names so each true difference has the same success.
    y = np.array([[0, -1j], [1j, 0]])
    plus_y = (np.eye(2) + y)/2
    scores = []
    for a, b in product(range(3), repeat=2):
        diff = (b - a) % 3
        success = 0.0
        for shift in range(3):
            shifted = (diff + shift) % 3
            ket = np.array([1, OMEGA**shifted])/np.sqrt(2)
            prob_plus = float(np.vdot(ket, plus_y @ ket).real)
            if shifted == 1:
                success += prob_plus / 3
            if shifted == 2:
                success += (1 - prob_plus) / 3
        scores.append(success)
    close(scores, np.full(9, (1+np.sqrt(3)/2)/3))
    return round(scores[0], 10)


def real_coefficients(r, s):
    return np.array([sum(r**a * s**((a+delta) % 3) for a in range(3))/3
                     for delta in range(3)])


def robust_real_bound():
    rows = []
    for error in (0.0, 0.01, 0.05, 0.25, 0.5):
        low = (1-2*error)**(1/3)
        certified = 1/3 + error*(1-low)/9
        best = 1/3
        for r, s in product(np.linspace(low, 1, 19), repeat=2):
            coeff = real_coefficients(r, s)
            optimum = 1/3 + np.ptp(coeff)/6
            closed = 1/3+(1-r)*(1-s)*(1+max(r,s)+r*s)/18
            close(optimum, closed)
            best = max(best, optimum)
        close(best, certified)
        rows.append({"cycle_trace_error": error, "sharp_real_score": round(certified, 10)})
    return rows


def injected_programs():
    # Alice's reference share is measured in X, then Bob uses its reported sign.
    pair = np.array([[np.sqrt(3)/2, 0], [0, 1/2]])
    program = None
    for sign in (1, -1):
        steered = np.array([1, sign]) @ pair / np.sqrt(2)
        close(np.vdot(steered, steered), 0.5)
        corrected = (np.eye(2) if sign == 1 else Z) @ (np.sqrt(2)*steered)
        candidate = Z @ H @ S @ corrected
        close(candidate[1]/candidate[0], OMEGA)
        program = candidate
    gate = np.diag([1, OMEGA])
    conjugate = X @ program
    close(conjugate[1]/conjugate[0], OMEGA.conjugate())
    # CNOT from unknown route occupancy onto program; read the program in Z.
    forward = [np.diag(program), np.diag(program[::-1])]
    reverse = [np.diag(conjugate), np.diag(conjugate[::-1])]
    branches = [(forward[0], gate, 0.5),
                (reverse[0] @ forward[1], gate, 0.25),
                (reverse[1] @ forward[1], np.eye(2), 0.25)]
    count = 0
    for action, expected, weight in branches:
        close(action.conj().T @ action, weight*np.eye(2))
        for i, j in product(range(2), repeat=2):
            unit = np.eye(2)[:, [i]] @ np.eye(2)[[j], :]
            close(action @ unit @ action.conj().T,
                  weight*(expected @ unit @ expected.conj().T))
            count += 1
    return {"matrix_units_verified": count, "window_success": 0.75,
            "abort_returns_identity": 0.25, "mean_reference_pairs_until_success": 2}


def run():
    orbit = [(0, 1)]
    for _ in range(3):
        orbit.append(cyclic_counter(orbit[-1]))
    assert orbit == [(0, 1), (1, 1), (1, 0), (0, 1)]
    assert cyclic_counter((0, 0)) == (0, 0)
    return {"classical_cycle": orbit, "two_complex_path_scores": score_characters(),
            "four_real_dimensions_score": four_real_dimensions(),
            "fixed_detector_balanced_score": fixed_readout(),
            "approximate_closure": robust_real_bound(), "gate_supply": injected_programs()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    results = run()
    out = Path(__file__).with_name("cyclic_port_comparison_results.json")
    if args.check:
        assert json.loads(out.read_text(encoding="utf-8")) == json.loads(json.dumps(results))
    else:
        with out.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results))
