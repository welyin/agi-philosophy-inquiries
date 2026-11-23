"""Late pair parity after a single addressed XOR. Check the Helstrom ensembles directly."""

import argparse
import itertools
import json
from pathlib import Path

import numpy as np


def uniform_edges(n):
    return {(i, j): 2 / (n * (n - 1)) for i, j in itertools.combinations(range(n), 2)}


def spectral_strategy(markers, edges):
    gram = markers @ markers.conj().T
    weights = np.zeros(gram.shape)
    for (i, j), probability in edges.items():
        weights[i, j] = weights[j, i] = probability * abs(gram[i, j])
    values, vectors = np.linalg.eigh(weights)
    amplitude = np.abs(vectors[:, -1])
    return 0.5 + float(values[-1]) / 2, amplitude


def query_ensemble(initial, markers):
    # Computational-basis XOR on the target, with any retained auxiliary system.
    n, _, aux = initial.shape
    gram = markers @ markers.conj().T
    tables = list(itertools.product((0, 1), repeat=n))
    answers = []
    for table in tables:
        after = initial.copy()
        for site, bit in enumerate(table):
            if bit:
                after[site] = initial[site, ::-1]
        ket = after.reshape(-1)
        density = np.outer(ket, ket.conj()).reshape(n, 2*aux, n, 2*aux)
        answers.append((density * gram[:, None, :, None]).reshape(ket.size, ket.size))
    return tables, answers


def helstrom_by_edge(initial, markers, edges):
    tables, answers = query_ensemble(initial, markers)
    n = len(markers)
    negative_norms = np.linalg.norm((initial[:, 0] - initial[:, 1]) / np.sqrt(2), axis=1)
    score = 0
    for (i, j), probability in edges.items():
        difference = sum((-1)**(bits[i] ^ bits[j]) * state
                         for bits, state in zip(tables, answers)) / (2**(n-1))
        direct = 0.5 + np.abs(np.linalg.eigvalsh(difference)).sum() / 4
        expected = 0.5 + abs(np.vdot(markers[j], markers[i])) * negative_norms[i] * negative_norms[j]
        np.testing.assert_allclose(direct, expected, atol=1e-12)
        score += probability * direct
    return float(score)


def score_winner(markers, edges):
    theory, amplitude = spectral_strategy(markers, edges)
    initial = amplitude[:, None, None] * np.array([1, -1])[None, :, None] / np.sqrt(2)
    actual = helstrom_by_edge(initial, markers, edges)
    np.testing.assert_allclose(actual, theory, atol=1e-12)
    return round(actual, 12)


def extremal_markers(n, r):
    first = np.array([1., 0.])
    second = np.array([r, np.sqrt(1 - r*r)])
    bisector = (first + second) / np.sqrt(2 + 2*r)
    return np.array([first, second] + [bisector] * (n-2))


def run():
    rng = np.random.default_rng(24033)
    random_trials = 0
    for n in (3, 4):
        for _ in range(16):
            markers = rng.normal(size=(n, 3)) + 1j * rng.normal(size=(n, 3))
            markers /= np.linalg.norm(markers, axis=1)[:, None]
            initial = rng.normal(size=(n, 2, 2)) + 1j * rng.normal(size=(n, 2, 2))
            initial /= np.linalg.norm(initial)
            edges = dict.fromkeys(uniform_edges(n))
            for e, weight in zip(edges, rng.dirichlet(np.ones(len(edges)))):
                edges[e] = weight
            actual = helstrom_by_edge(initial, markers, edges)
            optimal, _ = spectral_strategy(markers, edges)
            assert actual <= optimal + 1e-12
            score_winner(markers, edges)
            random_trials += 1

    equality_cases = 0
    for n in (2, 3, 4, 7, 12):
        for r in np.linspace(0, 1, 9):
            markers = extremal_markers(n, r)
            bound = 0.5 + 1/n - (1-r)/(n*(n-1))
            predicted, _ = spectral_strategy(markers, uniform_edges(n))
            np.testing.assert_allclose(predicted, bound, atol=1e-12)
            if n <= 4:
                score_winner(markers, uniform_edges(n))
            equality_cases += 1

    uniform = uniform_edges(3)
    identical = np.ones((3, 1))
    phases = np.exp(1j * np.array([0, 1.02, -2.23]))[:, None]
    ideal = score_winner(identical, uniform)
    assert abs(ideal - 5/6) < 1e-12
    assert score_winner(phases, uniform) == ideal
    disconnected = {(0, 1): 0.8, (2, 3): 0.2}
    markers = np.eye(3)[[0, 0, 1, 2]]
    assert score_winner(markers, disconnected) == score_winner(np.ones((4, 1)), disconnected)
    return {"three_site_ideal": ideal, "three_site_distinct_records": score_winner(np.eye(3), uniform),
            "three_site_one_distinguishable_pair_maximum": score_winner(extremal_markers(3, 0), uniform),
            "three_site_worst_overlap_0p8_maximum": score_winner(extremal_markers(3, 0.8), uniform),
            "disconnected_unseen_edge": score_winner(markers, disconnected),
            "random_complex_ancillary_trials": random_trials, "sharp_bound_instances": equality_cases}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    output = Path(__file__).with_name("late_pair_comparison_results.json")
    if args.check:
        assert json.loads(output.read_text(encoding="utf-8")) == result
    else:
        with output.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
