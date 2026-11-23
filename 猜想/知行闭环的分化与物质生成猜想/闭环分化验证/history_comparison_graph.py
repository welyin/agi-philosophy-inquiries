"""Check two-history tests of a single history-preserving channel."""

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


def pair_tests(records):
    # Row i contains the environment ket for the definite history i.
    gram = records @ records.conj().T
    assert np.allclose(gram.diagonal(), 1)
    return (1 + gram.real) / 2


def coherent_return(records):
    mean = np.mean(records, axis=0)
    return float(np.vdot(mean, mean).real)


def components(n, edges):
    groups = list(range(n))
    for a, b in edges:
        old, new = groups[b], groups[a]
        groups = [new if value == old else value for value in groups]
    return groups


def check_all_small_graphs():
    n = 5
    possible = list(itertools.combinations(range(n), 2))
    disconnected = 0
    for bits in itertools.product((0, 1), repeat=len(possible)):
        edges = [e for e, bit in zip(possible, bits) if bit]
        groups = components(n, edges)
        records = np.eye(n)[groups]
        scores = pair_tests(records)
        assert all(scores[i, j] == 1 for i, j in edges)
        if len(set(groups)) > 1:
            assert scores.min() == 0.5
            disconnected += 1
        else:
            assert np.all(scores == 1)
    return {"graphs": 2**len(possible), "disconnected_examples": disconnected}


def path_test(h):
    angles = np.linspace(0, math.pi / 2, h + 1)
    records = np.stack((np.cos(angles), np.sin(angles)), axis=1)
    scores = pair_tests(records)
    edge_score = (1 + math.cos(math.pi / (2 * h))) / 2
    assert np.allclose(np.diag(scores, 1), edge_score)
    assert np.isclose(scores[0, -1], 0.5)
    assert np.isclose(h * math.acos(2 * edge_score - 1), math.pi / 2)
    return {"edges": h, "each_neighbor_success": round(edge_score, 12),
            "endpoint_success": 0.5, "global_uniform_success": round(coherent_return(records), 12)}


def average_and_partner_checks():
    rng = np.random.default_rng(2307)
    trials = 0
    for n in (3, 5, 8, 16):
        for jitter in (0.01, 0.3, 1.0):
            for _ in range(5):
                records = jitter * (rng.normal(size=(n, 4)) + 1j * rng.normal(size=(n, 4)))
                records[:, 0] += 1
                records /= np.linalg.norm(records, axis=1, keepdims=True)
                edges = {(i, i+1) for i in range(n-1)}
                edges.update((a, b) for a, b in itertools.combinations(range(n), 2)
                             if rng.random() < 0.25)
                laplacian = np.zeros((n, n))
                for i, j in edges:
                    laplacian[i, i] += 1
                    laplacian[j, j] += 1
                    laplacian[i, j] -= 1
                    laplacian[j, i] -= 1
                scores = pair_tests(records)
                gap = np.linalg.eigvalsh(laplacian)[1]
                uniform = coherent_return(records)
                assert 1 - uniform <= 4 * sum(1-scores[i, j] for i, j in edges)/(n*gap) + 1e-12
                gram = records @ records.conj().T
                assert np.isclose(gram.sum().real / n**2, uniform)
                # A generic joint pure input with a three-dimensional retained partner.
                data = rng.normal(size=(n, 3)) + 1j * rng.normal(size=(n, 3))
                data /= np.linalg.norm(data)
                blocks = np.einsum('ia,jb->iajb', data, data.conj())
                difference = (blocks * (gram - 1)[:, None, :, None]).reshape(3*n, 3*n)
                distance = np.abs(np.linalg.eigvalsh(difference)).sum()/2
                mean_bound = math.sqrt(max(0, 2*n*(1-math.sqrt(uniform))))
                anchor_bound = min(max(np.linalg.norm(records - e, axis=1)) for e in records)
                assert distance <= min(1, mean_bound, anchor_bound) + 1e-12
                trials += 1
    return trials


def run():
    d = 64
    records = np.ones((d, 1), dtype=complex)
    records[-1] = -1
    uniform = coherent_return(records)
    fixed_score = pair_tests(records)[0, -1]
    assert uniform == (1 - 2/d)**2 and fixed_score == 0
    return {"component_checks": check_all_small_graphs(),
            "path_examples": [path_test(n) for n in (4, 16, 64)],
            "random_partner_and_spectral_checks": average_and_partner_checks(),
            "one_hidden_sign": {"histories": d, "uniform_success": uniform,
                                "fixed_phase_pair_success": float(fixed_score),
                                "optimal_pair_visibility": 1,
                                "worst_trace_error": 1}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    path = Path(__file__).with_name("history_comparison_graph_results.json")
    if args.check:
        assert json.loads(path.read_text(encoding="utf-8")) == result
    else:
        with path.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
