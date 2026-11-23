"""Orbit spans, reversible storage, and finite tests for descendants of a cyclic source."""
import argparse
import json
from fractions import Fraction
from functools import reduce
from itertools import product
from pathlib import Path

import numpy as np

from cycle_memory_readout import CYCLE, V, W


P = W.T @ CYCLE @ W
START = W.T @ V
STATES = [np.linalg.matrix_power(P, i) @ START for i in range(3)]


def tensor(items):
    return reduce(np.kron, items)


def partitions(n, minimum=1):
    if n == 0:
        yield ()
        return
    for k in range(minimum, n+1):
        for tail in partitions(n-k, k):
            yield (k,) + tail


def check_partition(parts):
    # Each block has one repeatedly applied count, distinct blocks choose independently.
    grouped = [np.column_stack([tensor([ket]*k) for ket in STATES]) for k in parts]
    labels = list(product(range(3), repeat=len(parts)))
    columns = np.column_stack([tensor([group[:, j] for group, j in zip(grouped, label)])
                               for label in labels])
    predicted = int(np.prod([2 if k == 1 else 3 for k in parts]))
    analytic_gram = tensor([np.eye(3) + ((-0.5)**k)*(np.ones((3, 3))-np.eye(3))
                            for k in parts])
    np.testing.assert_allclose(columns.T @ columns, analytic_gram, atol=1e-11)

    # Columns in physical memory. Coordinates in an arbitrary real orthonormal
    # code basis model a compressor when a general joint isometry is available.
    vals, vecs = np.linalg.eigh(columns @ columns.T)
    basis = vecs[:, vals > 1e-9]
    assert basis.shape[1] == predicted
    encoded = np.zeros((2*((predicted+1)//2), len(labels)))
    encoded[:predicted] = basis.T @ columns
    np.testing.assert_allclose(basis @ encoded[:predicted], columns, atol=1e-11)
    for group_index in range(len(parts)):
        gate = tensor([tensor([P if j == group_index else np.eye(2)]*k)
                       for j, k in enumerate(parts)])
        compressed_gate = basis.T @ gate @ basis
        np.testing.assert_allclose(gate @ basis, basis @ compressed_gate, atol=1e-11)
        np.testing.assert_allclose(compressed_gate.T @ compressed_gate,
                                   np.eye(predicted), atol=1e-11)
    return predicted, len(labels)


def synchronous_readout(n):
    tests = n//3
    # Check using the actual trine preparations and three native yes/no projections.
    q = Fraction(1, 4**tests)
    hit = 1-q+q*q/3
    returned = hit+(1-hit)*Fraction(1, 4**n)
    for truth in range(3):
        live = [abs(STATES[truth] @ ket)**(2*tests) for ket in STATES]
        success = 0
        for hits in product((0, 1), repeat=3):
            if not hits[truth]:
                continue
            prob = np.prod([p if yes else 1-p for p, yes in zip(live, hits)])
            success += prob/sum(hits)
        np.testing.assert_allclose(success, float(hit), atol=1e-12)
    return {"n": n, "guess_success": round(float(hit), 10),
            "full_projection_success": round(float(returned), 10)}


def independent_tests():
    unit_vectors = [np.kron(ket, ket) for ket in STATES]
    verify = 2*sum(np.outer(u, u) for u in unit_vectors)/3
    ev = np.linalg.eigvalsh(verify)
    np.testing.assert_allclose(ev, [0, .5, .5, 1], atol=1e-12)
    sectors = []
    for n in range(1, 7):
        # Independent real cyclic shifts act on this one real plane only by
        # the signed total phase. The complex vector is a calculation basis.
        axis = np.array([1, 1j])/np.sqrt(2)
        signs = [1 if j % 2 == 0 else -1 for j in range(n)]
        vector = tensor([axis if s == 1 else axis.conj() for s in signs])
        eigenvalue = np.vdot(axis, P @ axis)
        for exponents in product(range(3), repeat=n):
            # Apply the gates in the tensor factors without materializing 729 matrices.
            actual = tensor([np.linalg.matrix_power(P, k) @ (axis if s == 1 else axis.conj())
                             for s, k in zip(signs, exponents)])
            phase = eigenvalue ** (sum(s*k for s, k in zip(signs, exponents)) % 3)
            np.testing.assert_allclose(actual, phase*vector, atol=1e-12)
        sectors.append(3**n)
    return {"single_site_verification_eigenvalues": [0, 0.5, 0.5, 1],
            "signed_sector_settings_checked": sum(sectors)}


def run():
    by_size = {}
    total = 0
    for n in range(1, 7):
        rows = []
        for parts in partitions(n):
            dimension, count = check_partition(parts)
            rows.append({"group_sizes": list(parts), "orbit_labels": count,
                         "code_dimension": dimension, "extra_memory_dimension": (dimension+1)//2})
            total += 1
        by_size[str(n)] = rows
    return {"partitions_checked": total, "cases_by_size": by_size,
            "independent_projection": independent_tests(),
            "synchronous_classical_log": [synchronous_readout(n) for n in (3, 6, 12, 18)]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    file = Path(__file__).with_name("ancestral_cycle_span_results.json")
    if args.check:
        assert json.loads(file.read_text(encoding="utf-8")) == result
    else:
        with file.open("x", encoding="utf-8") as out:
            out.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in result if k != "cases_by_size"}))
