"""Construct finite packets from the actual four-leg source and verify transfer."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


BASE = Path(__file__).resolve().parents[1]
EARLIER = runpy.run_path(str(BASE / "1355" / "check.py"))
NU = 0.01


def matchings(ports):
    if not ports:
        yield ()
        return
    first, *tail = ports
    for index, partner in enumerate(tail):
        for rest in matchings(tail[:index] + tail[index + 1:]):
            yield ((first, partner), *rest)


def source_columns(k, local_packet):
    """The same wiring for each candidate, with the physical, unnormalized C."""
    indices = np.indices((3,) * (2 * k)).reshape(2 * k, -1)
    cup = np.eye(3) / math.sqrt(3)
    cols = []
    for matching in matchings(list(range(2 * k))):
        vector = np.ones(indices.shape[1], dtype=complex)
        left, right = [], []
        for i, j in matching:
            if j < k:
                left.append((i, j))
            elif i >= k:
                right.append((i, j))
            else:
                vector *= cup[indices[i], indices[j]]
        assert len(left) == len(right)
        for (a, b), (c, d) in zip(sorted(left), sorted(right)):
            vector *= local_packet[indices[a], indices[b], indices[c], indices[d]]
        cols.append(vector)
    return np.column_stack(cols)


def test_case(k, group):
    device = EARLIER["Device"](group)
    packet = device.matrix(EARLIER["THETA"])[:, 2].reshape((3,) * 4)
    matrix = source_columns(k, packet)
    routes = matrix.shape[1]
    gram = matrix.conj().T @ matrix
    eig, modes = np.linalg.eigh(gram)
    expected = {("SU3", 3): 6, ("SO3", 3): 15,
                ("SU3", 4): 23, ("SO3", 4): 91}[(group, k)]
    assert max(abs(eig[:-expected]), default=0.0) < 1e-11
    gap = eig[-expected]
    assert gap > NU
    filt = (modes * np.sqrt(NU / np.maximum(eig, NU))) @ modes.conj().T
    encode = (matrix @ filt) / math.sqrt(routes)
    record = encode.conj().T @ encode
    mixed_success = np.trace(record).real / routes
    assert np.allclose(record @ record, NU / routes * record, atol=2e-13)
    assert abs(mixed_success - expected * NU / routes**2) < 1e-13
    # Keep an external Bell partner across both heralded operations.
    logical = modes[:, -2:] / math.sqrt(2)
    assert np.linalg.norm(record @ logical - (NU / routes) * logical) < 1e-12
    assert np.allclose(logical.conj().T @ logical, np.eye(2)/2)
    return {"routes": routes, "rank": expected, "minimum_positive_gram": float(gap),
            "one_way_on_code": NU / routes, "mixed_trial_success": float(mixed_success),
            "round_trip_on_code": (NU / routes)**2}


def compare(current, saved):
    if isinstance(current, dict):
        assert current.keys() == saved.keys()
        for key in current:
            compare(current[key], saved[key])
    else:
        np.testing.assert_allclose(current, saved, rtol=2e-10, atol=1.e-13)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--write-new", action="store_true")
    args = cli.parse_args()
    results = {str(k): {g: test_case(k, g) for g in ("SU3", "SO3")} for k in (3, 4)}
    for row in results.values():
        diff = row["SO3"]["mixed_trial_success"] - row["SU3"]["mixed_trial_success"]
        row["maximum_equal_halfwidth_for_separation_exclusive"] = diff / 2
    path = Path(__file__).with_name("results.json")
    if args.write_new:
        with path.open("x", encoding="utf-8") as out:
            json.dump(results, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        compare(results, json.loads(path.read_bytes()))
    print(json.dumps(results))
