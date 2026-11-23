"""Check finite classical recodings and the same two-record task on a real plane."""

import argparse
import json
from collections import deque
from math import factorial
from pathlib import Path

import numpy as np


RESULT = Path(__file__).with_name("recoverable_record_paths_results.json")


def reachable_encodings(records, labels):
    start = tuple(range(records))
    queue = deque([start])
    visited = {start}
    while queue:
        current = queue.popleft()
        for empty in set(range(labels)).difference(current):
            for record in range(records):
                updated = list(current)
                updated[record] = empty
                target = tuple(updated)
                if target not in visited:
                    visited.add(target)
                    queue.append(target)
    return len(visited)


def check_transfer(before, after, labels, fraction):
    records = len(before)
    identity = np.eye(records)
    basis = np.eye(labels)
    encoding = (1 - fraction) * basis[:, before] + fraction * basis[:, after]
    # Both endpoints used by the moving record must have the same reading.
    decoder = np.zeros((records, labels))
    for record, (origin, target) in enumerate(zip(before, after)):
        decoder[record, [origin, target]] = 1
    assert np.array_equal(decoder.sum(axis=0), np.ones(labels))
    assert np.array_equal(decoder @ encoding, identity)
    partner = np.eye(3)
    assert np.array_equal(np.kron(decoder, partner) @ np.kron(encoding, partner),
                          np.eye(3 * records))


def run_checks():
    component_sizes = {}
    segment_samples = 0
    fractions = np.linspace(0, 1, 17)
    for records in range(2, 6):
        full = reachable_encodings(records, records)
        spare = reachable_encodings(records, records + 1)
        assert full == 1 and spare == factorial(records + 1)
        component_sizes[str(records)] = {"full": full, "one_spare": spare}
        assignment = list(range(records))
        for record, new_label in ((0, records), (1, 0), (0, 1)):
            updated = assignment.copy()
            updated[record] = new_label
            for fraction in fractions:
                check_transfer(assignment, updated, records + 1, fraction)
                segment_samples += 1
            assignment = updated
        assert assignment == [1, 0, *range(2, records)]

    swap = np.array([[0, 1], [1, 0]])
    identity = np.eye(2)
    worst_error = 0.0
    max_real_error = 0.0
    for fraction in fractions:
        encoding = (1 - fraction) * identity + fraction * swap
        recovered = encoding @ encoding
        worst_error = max(worst_error, float(1 - np.min(np.diag(recovered))))

        angle = np.pi * fraction / 2
        rotation = np.array([[np.cos(angle), -np.sin(angle)],
                             [np.sin(angle), np.cos(angle)]])
        joint = np.kron(rotation, np.eye(3))
        max_real_error = max(max_real_error,
                             float(np.max(np.abs(joint.T @ joint - np.eye(6)))))
    assert worst_error == 0.5 and max_real_error < 1e-12
    assert np.allclose(rotation ** 2, swap, rtol=0, atol=1e-12)
    return {
        "records_to_reachable_deterministic_encodings": component_sizes,
        "exact_segment_samples_including_classical_partner": segment_samples,
        "two_label_swap_peak_error": worst_error,
        "real_plane_inverse_and_endpoint_within_1e-12": True,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    results = run_checks()
    if args.check:
        if json.loads(RESULT.read_text(encoding="utf-8")) != results:
            raise SystemExit("Saved result differs from recomputation.")
        print("Verified", RESULT.name)
    else:
        RESULT.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
