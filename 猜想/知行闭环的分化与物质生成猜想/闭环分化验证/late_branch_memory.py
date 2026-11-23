"""A delayed choice with zero or one frozen branch record returned on time."""

import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import hypot, sqrt
from pathlib import Path

import numpy as np

from growing_frontier import build_frontiers, reports


INPUTS = tuple(product(range(2), repeat=2))
PAULI = np.array([[[0, 1], [1, 0]], [[1, 0], [0, -1]], [[0, -1j], [1j, 0]]])


def tree_archives():
    depth = 3
    generation = build_frontiers(depth)[depth]
    archives = {}
    for wanted in INPUTS:
        found = []
        for side in range(2):
            half = [path[1:] for path in generation if path[0] == side]
            preparation = dict.fromkeys(half, 1)
            preparation[half[0]] = wanted[side]
            result, _ = reports({depth - 1: preparation}, 1, 1)
            # Subtree birth in the growing tree is one time unit later.
            local_time, value = result[depth - 1]
            assert 1 + local_time == 5 and 1 + local_time + 1 == 6
            found.append(value)
        assert tuple(found) == wanted
        archives[wanted] = tuple(found)
    return archives


def arrival_patterns():
    # The reader sends both requests. Slow replies are valid but miss the deadline.
    question_time, deadline = 7, 9
    slow = 4
    fast = 2
    flags = []
    for delays in ((slow, slow), (fast, slow), (slow, fast)):
        arrivals = [i for i, delay in enumerate(delays) if question_time + delay <= deadline]
        flags.append(arrivals[0] if arrivals else None)
    assert flags == [None, 0, 1]
    return flags


def classical_scores(archives, flags):
    scores = Counter()
    for encoded in product(range(2), repeat=4):
        successes = []
        for flag in flags:
            histories = defaultdict(Counter)
            for word, memory in zip(INPUTS, encoded):
                reply = None if flag is None else archives[word][flag]
                for question in range(2):
                    histories[memory, question, reply][word[question]] += 1
            successes.append(Fraction(sum(max(tally.values()) for tally in histories.values()), 8))
        scores[successes[0], (successes[1] + successes[2]) / 2] += 1
    assert scores == {
        (Fraction(1, 2), Fraction(3, 4)): 2,
        (Fraction(3, 4), Fraction(7, 8)): 12,
        (Fraction(1, 2), Fraction(1)): 2,
    }
    return scores


def optimal_qubit(p):
    if p <= 2 / 3:
        a = 1 / sqrt(2)
        c = 0.
        score = 0.5 + p / 4 + (2 - p) / (4 * sqrt(2))
    else:
        norm = hypot(2 * (1 - p), p)
        a = 2 * (1 - p) / norm
        c = sqrt(max(0., 1 - 2 * a * a))
        score = 0.5 + p / 4 + norm / 4
    return a, c, score


def quantum_score(p, flags, archives):
    a, c, optimum = optimal_qubit(p)
    states = {}
    for word in INPUTS:
        x, y = (1 - 2 * bit for bit in word)
        bloch = np.array([x * a, y * a, x * y * c])
        density = (np.eye(2) + np.tensordot(bloch, PAULI, axes=1)) / 2
        assert np.allclose(np.linalg.eigvalsh(density), [0, 1])
        states[word] = density
    total = 0.
    for flag, weight in zip(flags, (1 - p, p / 2, p / 2)):
        for question in range(2):
            for reply in ((None,) if flag is None else (0, 1)):
                # Positive-eigenvalue effect optimally separates the two
                # unnormalized hypotheses in this observable history.
                hypotheses = [np.zeros((2, 2), dtype=complex) for _ in range(2)]
                for word, density in states.items():
                    if flag is None or archives[word][flag] == reply:
                        hypotheses[word[question]] += density / 4
                values, vectors = np.linalg.eigh(hypotheses[0] - hypotheses[1])
                keep = vectors[:, values > 0]
                effect = keep @ keep.conj().T
                correct = np.trace(effect @ hypotheses[0] + (np.eye(2) - effect) @ hypotheses[1])
                assert abs(correct.imag) < 1e-12
                total += weight * correct.real / 2
    assert abs(total - optimum) < 1e-12
    return float(total)


def run():
    archives = tree_archives()
    flags = arrival_patterns()
    classes = classical_scores(archives, flags)
    samples = {}
    for p in (Fraction(0), Fraction(1, 2), Fraction(2, 3), Fraction(4, 5), Fraction(1)):
        classical = max((1 - p) * none + p * one for none, one in classes)
        quantum = quantum_score(float(p), flags, archives)
        samples[str(p)] = {"classical": str(classical), "qubit": round(quantum, 12)}
    assert samples["4/5"]["classical"] == "9/10"
    assert abs(samples["4/5"]["qubit"] - (7 + sqrt(5)) / 10) < 1e-12
    return {"branch_preparations": len(archives), "return_patterns": len(flags),
            "one_bit_encoders": sum(classes.values()), "verified_profiles": samples}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    path = Path(__file__).with_name("late_branch_memory_results.json")
    if args.check:
        assert json.loads(path.read_text(encoding="utf-8")) == result
    else:
        with path.open("x", encoding="utf-8") as output:
            output.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
