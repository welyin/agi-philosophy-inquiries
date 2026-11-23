"""Verify optimal logical-qubit recovery after copying the actual split location."""

import argparse
import json
from pathlib import Path

import numpy as np


RESULT = Path(__file__).with_name("lineage_location_recovery_results.json")
TOL = 1e-10


def density(vector):
    vector = np.asarray(vector, dtype=complex)
    return np.outer(vector, vector.conj())


def verify(tau, parities):
    size = len(tau)
    identity = np.eye(size)
    zero = np.zeros_like(identity)
    projectors = [np.diag(np.array(parities) == k).astype(float) for k in (0, 1)]
    routes = [np.block([[projectors[k], zero], [zero, projectors[1-k]]])
              for k in (0, 1)]
    # Store an actual orthogonal location copy as the last output factor.
    location_copy = sum(np.kron(q, np.eye(2)[:, [k]]) for k, q in enumerate(routes))
    assert np.allclose(location_copy.conj().T @ location_copy, np.eye(2*size), atol=TOL)

    corrections = []
    block_norms = []
    for k in (0, 1):
        cross = projectors[k] @ tau @ projectors[1-k]
        left, singular, right_h = np.linalg.svd(cross)
        unitary = left @ right_h
        corrections.append(np.block([[identity, zero], [zero, unitary]]))
        block_norms.append(float(singular.sum()))
    visibility = sum(block_norms)
    assert abs(block_norms[0] - block_norms[1]) < TOL

    def stored(operator):
        before = np.kron(operator, tau)
        with_location = location_copy @ before @ location_copy.conj().T
        without_location = np.trace(with_location.reshape(2*size, 2, 2*size, 2),
                                    axis1=1, axis2=3)
        block_sum = sum(q @ before @ q for q in routes)
        assert np.allclose(without_location, block_sum, atol=TOL, rtol=0)
        return without_location

    def recover(operator):
        state = stored(operator)
        corrected = sum(u @ q @ state @ q @ u.conj().T
                        for q, u in zip(routes, corrections))
        return np.trace(corrected.reshape(2, size, 2, size), axis1=1, axis2=3)

    # Four matrix units test the channel with any untouched partner by linearity.
    bell_recovered = np.zeros((4, 4), dtype=complex)
    for x in (0, 1):
        for y in (0, 1):
            unit = np.outer(np.eye(2)[x], np.eye(2)[y])
            recovered = recover(unit)
            assert np.allclose(recovered, unit*(1 if x == y else visibility),
                               atol=TOL, rtol=0)
            bell_recovered += np.kron(unit, recovered) / 2
    bell = np.array([1, 0, 0, 1]) / np.sqrt(2)
    bell_score = float(np.vdot(bell, bell_recovered @ bell).real)

    difference = stored(np.array([[0, 1], [1, 0]], dtype=complex))
    distance = np.linalg.svd(difference, compute_uv=False).sum() / 2
    assert abs(distance - visibility) < TOL
    assert abs(bell_score - (1+distance)/2) < TOL
    probs = [float(np.trace(p @ tau).real) for p in projectors]
    assert visibility <= 2*np.sqrt(probs[0]*probs[1]) + TOL
    return {"location_0_probability_for_input_0_1": [round(v, 10) for v in probs],
            "trace_distance_between_encoded_plus_minus": round(float(distance), 10),
            "optimal_worst_case_squared_fidelity": round(bell_score, 10)}


def calculate():
    rng = np.random.default_rng(4201)
    matrix = rng.normal(size=(8, 5)) + 1j*rng.normal(size=(8, 5))
    mixed = matrix @ matrix.conj().T
    mixed /= np.trace(mixed)
    bell_key = density(np.array([1, 0, 0, 1])/np.sqrt(2))
    reduced_key = np.trace(bell_key.reshape(2, 2, 2, 2), axis1=1, axis2=3)
    cases = {
        "fixed_evidence": (density([1, 0]), [0, 1], 0.5),
        "fair_classical_evidence": (np.eye(2)/2, [0, 1], 0.5),
        "fair_coherent_evidence": (density(np.array([1, 1])/np.sqrt(2)), [0, 1], 1),
        "partly_coherent_evidence": (np.array([[0.5, 0.25j], [-0.25j, 0.5]]),
                                    [0, 1], 0.75),
        "biased_pure_evidence": (density([np.sqrt(0.1), 1j*np.sqrt(0.9)]),
                                 [0, 1], 0.8),
        "bell_evidence_both_parts_accessible": (bell_key, [0, 0, 1, 1], 1),
        "bell_evidence_second_part_inaccessible": (reduced_key, [0, 1], 0.5),
        "three_endpoints_classical_independent": (np.eye(8)/8,
                                                 [i.bit_count()%2 for i in range(8)], 0.5),
        "three_endpoints_coherent_independent": (density(np.ones(8)/np.sqrt(8)),
                                                 [i.bit_count()%2 for i in range(8)], 1),
        "three_endpoints_mixed": (mixed, [i.bit_count()%2 for i in range(8)], None),
    }
    results = {}
    for name, (tau, parity, expected) in cases.items():
        row = verify(tau, parity)
        if expected is not None:
            assert abs(row["optimal_worst_case_squared_fidelity"] - expected) < TOL
        results[name] = row
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    calculated = calculate()
    if args.check:
        if json.loads(RESULT.read_text(encoding="utf-8")) != calculated:
            raise SystemExit("Saved results differ from the recalculated results.")
        print(f"Verified {len(calculated)} cases, including all four logical matrix units.")
    else:
        RESULT.write_text(json.dumps(calculated, indent=2) + "\n", encoding="utf-8")
