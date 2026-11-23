"""Audit the actual seed, arbitrary additive generators and finite-time readout."""

import argparse
import hashlib
import json
import sys
from functools import reduce
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "猜想" / "知行闭环的分化与物质生成猜想" / "闭环分化验证"))
from quantum_recursive_lineage import K, all_shapes, controlled_flip
from quantum_recursive_ports import READY, V

I = np.eye(2)
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1, -1])
YZERO = np.array([1, 1j]) / np.sqrt(2)
ATOL = 3e-12


def equal(left, right):
    np.testing.assert_allclose(left, right, atol=ATOL, rtol=0)


def product(factors):
    return reduce(np.kron, factors)


def build(shape, split):
    if shape is None:
        return I, [1]
    left, ls = build(shape[0], split)
    right, rs = build(shape[1], split)
    return np.kron(left, right) @ split, ls + [-s for s in rs]


def pauli(vec):
    x, q, z = vec
    return (x * X + q * Y + z * Z) / 2


def defect(code, parent, leaves):
    count = len(leaves)
    total = np.zeros((2 ** count, 2 ** count), dtype=complex)
    for site, generator in enumerate(leaves):
        factors = [I] * count
        factors[site] = pauli(generator)
        total += product(factors)
    return total @ code - code @ pauli(parent)


def reduced_after_losing(state, count, lost):
    columns = np.moveaxis(state.reshape([2] * count), lost, -1).reshape(-1, 2)
    return columns @ columns.conj().T


def evaluate():
    split = V @ READY
    rotation = np.array([[1, -1], [1, 1]]) / np.sqrt(2)
    equal(split, controlled_flip("S") @ np.kron(rotation, I) @ controlled_flip("R") @ K)
    rng = np.random.default_rng(1819)
    cases = []
    reductions = 0
    for count in range(2, 7):
        shapes = list(all_shapes(count))
        for shape in shapes:
            code, signs = build(shape, split)
            plus, minus = (code @ vector for vector in (YZERO, YZERO.conj()))
            equal(plus, product([YZERO if s == 1 else YZERO.conj() for s in signs]))
            equal(minus, plus.conj())
            parent = rng.normal(size=3)
            leaves = rng.normal(size=(count, 3))
            d = defect(code, parent, leaves)
            delta = np.dot(signs, leaves[:, 1]) - parent[1]
            transverse = np.hypot(leaves[:, 0], leaves[:, 2])
            base = delta**2 + parent[0]**2 + parent[2]**2 + transverse @ transverse
            offset = 2 * np.prod(transverse) if count == 2 else 0.0
            equal(np.linalg.eigvalsh(d.conj().T @ d), (base + np.array([-offset, offset])) / 4)
            # Only the diagonal part is used for the readout experiment.
            phase = 0.61 * float(np.dot(signs, leaves[:, 1]))
            a = np.exp(-0.5j * phase)
            output = (a * plus + a.conjugate() * minus) / np.sqrt(2)
            for site in range(count):
                factors = [Z] * count
                factors[site] = X
                equal(np.vdot(output, product(factors) @ output), signs[site] * np.sin(phase))
            equal(np.vdot(output, product([Z] * count) @ output), np.cos(phase))
            # Deliberately nonorthogonal and distinct retained reference vectors.
            p = np.kron(plus, np.array([1, 0]))
            m = np.kron(minus, np.array([1, 2j]) / np.sqrt(5))
            before = (p + m) / np.sqrt(2)
            after = (a * p + a.conjugate() * m) / np.sqrt(2)
            for site in range(count):
                equal(reduced_after_losing(after, count + 1, site),
                      reduced_after_losing(before, count + 1, site))
                reductions += 1
        cases.append(len(shapes))
    axis = lambda q: np.array([0.0, q, 0.0])
    for parent, children, size in ((1, [2, 0], .5), (2, [1, 0], .5)):
        equal(np.linalg.norm(defect(split, axis(parent), list(map(axis, children))), 2), size)
    code, _ = build(((None, None), None), split)
    equal(defect(code, axis(1), list(map(axis, [1, 0, 0]))), 0)
    interaction = split @ (X / 2) @ split.conj().T
    equal(interaction @ split, split @ (X / 2))
    for first, second in ((0, 2), (1, 3)):
        equal(interaction.reshape(2, 2, 2, 2).trace(axis1=first, axis2=second), 0)
    # A single interrogation interval can miss an entire relative revolution.
    for time, target in ((1, 0), (.5, 1)):
        reference = code[:, 0]
        acted = code @ (np.cos(np.pi*time) * I - 1j * np.sin(np.pi*time) * Y)[:, 0]
        singular = np.linalg.eigvalsh(np.outer(reference, reference.conj()) -
                                      np.outer(acted, acted.conj()))
        equal(sum(abs(singular)) / 2, target)
    # Same spectrum at the R1818 swap endpoint and at the R1819 birth.
    for axis_spectrum in rng.normal(size=(24, 3)):
        axis_spectrum /= np.linalg.norm(axis_spectrum)
        traceless_parent = -axis_spectrum
        left = axis(traceless_parent[1])
        dynamical_error = np.linalg.norm(defect(split, traceless_parent, [left, axis(0)]), 2)
        energy = (I - 2 * pauli(axis_spectrum)) / 2
        plus_energy = np.vdot(YZERO, energy @ YZERO).real
        minus_energy = np.vdot(YZERO.conj(), energy @ YZERO.conj()).real
        record_error = abs(plus_energy - minus_energy) / 2
        equal(record_error**2 + dynamical_error**2, .25)
    return {
        "original_topic40_and_topic49_ready_maps_equal": True,
        "leaf_counts": list(range(2, 7)),
        "all_ordered_binary_shapes_by_leaf_count": cases,
        "dense_generator_and_two_parity_certificates": sum(cases),
        "single_leaf_losses_with_retained_reference": reductions,
        "absolute_test_tolerance": ATOL,
        "compensating_internal_vertex_operator_norms": [.5, .5],
        "corresponding_root_to_all_leaves_operator_norm": 0,
        "channel_distance_at_1_and_half_for_2pi_alias": [0, 1],
        "transverse_joint_generator_example_has_zero_marginal_traces": True,
        "same_port_spectral_axes_checked_at_unit_gap": 24,
        "squared_sum_of_two_optimal_endpoint_errors_per_gap_squared": .25,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = evaluate()
    path = HERE / "independent_results.json"
    payload = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.check:
        if json.loads(path.read_bytes()) != result:
            raise SystemExit("1819 independent results differ")
    else:
        path.write_bytes(payload.encode("utf-8"))
    print("1819 independent PASS", hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12])
