"""Compare two actual real cyclic carriers across a shared or separately offset middle clock."""

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
import json
from math import sqrt
from pathlib import Path

import numpy as np

from cycle_memory_readout import CYCLE, V, W


def pair_carriers():
    orbit = [(W.T @ np.linalg.matrix_power(CYCLE, k) @ V) for k in range(3)]
    return {(a, b): np.concatenate((orbit[a], orbit[b])) / sqrt(2)
            for a, b in product(range(3), repeat=2)}


def direct_ensembles(pairs, separate):
    states = np.zeros((3, 16, 16))
    offsets = range(3) if separate else (0,)
    for a, b, c in product(range(3), repeat=3):
        for u, v in product(offsets, repeat=2):
            left = pairs[(a+u) % 3, (b+u) % 3]
            right = pairs[(b+v) % 3, (c+v) % 3]
            psi = np.kron(left, right)
            states[(c-a) % 3] += np.outer(psi, psi) / (9*len(offsets)**2)
    for state in states:
        np.testing.assert_allclose(np.trace(state), 1, atol=1e-13)
    return states


def measurement_certificate(states):
    mean = sum(states) / 3
    eig, vectors = np.linalg.eigh(mean)
    support = vectors[:, eig > 1e-11]
    inverse = (support / np.sqrt(eig[eig > 1e-11])) @ support.T
    complement = np.eye(mean.shape[0]) - support @ support.T
    measures = [(inverse @ state @ inverse + complement)/3 for state in states]
    np.testing.assert_allclose(sum(measures), np.eye(mean.shape[0]), atol=2e-12)
    primal = sum(np.trace(e @ r) for e, r in zip(measures, states)) / 3
    dual = sum(state @ e for e, state in zip(measures, states))/3
    np.testing.assert_allclose(dual, dual.T, atol=2e-12)
    for effect, state in zip(measures, states):
        assert min(np.linalg.eigvalsh(effect)) > -2e-12
        assert min(np.linalg.eigvalsh(dual-state/3)) > -2e-12
    np.testing.assert_allclose(np.trace(dual), primal, atol=2e-12)
    scores = [float(np.trace(state @ e)) for state, e in zip(states, measures)]
    np.testing.assert_allclose(scores, primal, atol=2e-12)
    return {"dimension": mean.shape[0], "ensemble_support_rank": support.shape[1],
            "optimal_score": round(float(primal), 12),
            "conditional_scores": [round(x, 12) for x in scores]}


def character_certificate(separate):
    blocks = defaultdict(lambda: [0, 0, 0])
    for first, second, h, g in product(range(2), range(2), (-1, 1), (-1, 1)):
        a = h*(1-first)
        b = h*first + g*(1-second)
        c = g*second
        key = ((a+c) % 3, b % 3)
        if separate:
            key += (h, g)
        blocks[key][c % 3] += 1
    shape_counts = Counter(tuple(sorted(row, reverse=True)) for row in blocks.values())
    assert sum(sum(row) for row in blocks.values()) == 16
    score = sum(sum(sqrt(n) for n in row)**2 for row in blocks.values()) / 48
    return {"block_shapes": {str(k): v for k, v in sorted(shape_counts.items())},
            "score": round(score, 12)}


def report_only(pairs):
    conditional = np.zeros((3, 4, 4))
    for (a, b), state in pairs.items():
        conditional[(b-a) % 3] += np.outer(state, state)/3
    effects = (4/3) * conditional
    np.testing.assert_allclose(sum(effects), np.eye(4), atol=2e-12)
    joint = np.zeros((3, 3, 3))
    for a, b, c in product(range(3), repeat=3):
        l = pairs[a, b]
        r = pairs[b, c]
        lp = np.array([l @ e @ l for e in effects])
        rp = np.array([r @ e @ r for e in effects])
        np.testing.assert_allclose(lp[(b-a) % 3], 2/3)
        np.testing.assert_allclose(rp[(c-b) % 3], 2/3)
        joint[(c-a) % 3] += np.outer(lp, rp)/27
    posterior = joint / joint.sum(axis=0)
    for i, j in product(range(3), repeat=2):
        expected = np.full(3, 1/4)
        expected[(i+j) % 3] = 1/2
        np.testing.assert_allclose(posterior[:, i, j], expected)
    return float(np.round(joint.max(axis=0).sum(), 12))


def six_outcomes():
    # All amplitudes are integers times 1/sqrt(8) after H on the route bit.
    initial = np.rint(np.sqrt(2)*V).astype(int)
    turns = [np.linalg.matrix_power(CYCLE, j).astype(int) @ initial for j in range(3)]
    weights = {}
    for a, b in product(range(3), repeat=2):
        weights[a, b] = np.concatenate([(turns[a][1:] + sign*turns[b][1:])**2
                                         for sign in (1, -1)])
        assert weights[a, b].sum() == 8
    local = np.zeros((3, 6), dtype=int)
    for (a, b), row in weights.items():
        local[(b-a) % 3] += row
    local_score = Fraction(int(local.max(axis=0).sum()), 9*8)
    assert local_score == Fraction(7, 12)
    scores = {"each_local_difference": str(local_score)}
    for label, offsets in (("shared_middle", (0,)), ("independent_frame_shifts", range(3))):
        counts = np.zeros((3, 6, 6), dtype=int)
        for a, b, c, u, v in product(range(3), range(3), range(3), offsets, offsets):
            left = weights[(a+u) % 3, (b+u) % 3]
            right = weights[(b+v) % 3, (c+v) % 3]
            counts[(c-a) % 3] += np.outer(left, right)
        assert counts.sum() == 27*64*len(offsets)**2
        scores[label] = str(Fraction(int(counts.max(axis=0).sum()), int(counts.sum())))
    assert scores["shared_middle"] == "49/96"
    assert scores["independent_frame_shifts"] == "11/24"
    return scores


def run():
    pairs = pair_carriers()
    result = {"ordinary_two_answers": report_only(pairs), "six_local_outcomes": six_outcomes()}
    for separate, identifier, predicted in (
            (False, "shared_middle", (26+4*sqrt(2))/48),
            (True, "independent_frame_shifts", 1/2)):
        certificate = measurement_certificate(direct_ensembles(pairs, separate))
        exact = character_certificate(separate)
        np.testing.assert_allclose([certificate["optimal_score"], exact["score"]],
                                   predicted, atol=2e-12, rtol=0)
        result[identifier] = {"real_matrix": certificate, "characters": exact,
                              "raw_input_triples_and_shifts": 243 if separate else 27}
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    path = Path(__file__).with_name("shared_middle_cycle_results.json")
    if args.check:
        assert json.loads(path.read_text(encoding="utf-8")) == result
    else:
        with path.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
