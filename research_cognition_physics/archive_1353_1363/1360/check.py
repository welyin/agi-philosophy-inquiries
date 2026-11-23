"""Contract the actual balanced preparation with a supplied trivalent matter packet."""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np


ARCHIVE = Path(__file__).resolve().parents[1]
BALANCED = runpy.run_path(str(ARCHIVE / "1358" / "check.py"))


def volume():
    tensor = np.zeros((3, 3, 3))
    for word in itertools.permutations(range(3)):
        inversions = sum(a > b for a, b in itertools.combinations(word, 2))
        tensor[word] = (-1)**inversions / math.sqrt(6)
    return tensor.ravel()


def test(group, p):
    device = BALANCED["EARLIER"]["Device"](group)
    theta = BALANCED["EARLIER"]["THETA"]
    packet = device.matrix(theta)[:, 2].reshape((3,) * 4)
    matrix = BALANCED["source_columns"](p, packet)
    routes = matrix.shape[1]
    eig, modes = np.linalg.eigh(matrix.conj().T @ matrix)
    nu = BALANCED["NU"]
    weights = np.sqrt(nu / np.maximum(eig, nu))
    encoded = matrix @ ((modes * weights) @ modes.conj().T) / math.sqrt(routes)
    # Each physical cap is the inverse of the normalized particle-antiparticle pair.
    target_size = 3**(2*p - 3)
    blocks = encoded.reshape(target_size, 27, routes)
    delivered = np.einsum("amr,m->ar", blocks, volume()) / math.sqrt(27)
    gram = delivered.conj().T @ delivered
    spec, basis = np.linalg.eigh(gram)
    ideal = nu / (routes * 27)
    rank = sum(spec > ideal/2)
    expected = 1 if p == 3 else {"SU3": 3, "SO3": 6}[group]
    assert rank == expected and np.allclose(spec[-rank:], ideal, rtol=2e-9, atol=1e-14)
    assert np.max(np.abs(spec[:-rank]), initial=0) < 1e-13
    # A Bell pair with two distinct code directions checks coherence on a real multiplicity.
    if rank >= 2:
        bell = basis[:, -2:] / math.sqrt(2)
        assert np.linalg.norm(gram @ bell - ideal*bell) < 1e-12
    report = {"routes": routes, "physical_dimension": target_size, "neutral_rank": int(rank),
              "success_on_code": ideal, "raw_mixed_success": float(np.trace(gram).real / target_size)}
    if p == 4:
        swapped = delivered.reshape((3,)*5 + (routes,)).swapaxes(0, 1).reshape(delivered.shape)
        for sign, label in [(1, "sym"), (-1, "antisym")]:
            component = (delivered + sign * swapped)/2
            joint = np.linalg.norm(component)**2 / target_size
            projected_rank = joint * target_size / ideal
            assert abs(projected_rank - round(projected_rank)) < 1e-9
            weight = (3 + sign) / 6
            report[label] = {"invariant_rank": round(projected_rank), "swap_flag": weight,
                             "joint": float(joint), "given_swap_flag": float(joint / weight)}
    return report


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            compare(a[key], b[key])
    else:
        np.testing.assert_allclose(a, b, atol=1e-15, rtol=1e-9)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    report = {group: {str(p): test(group, p) for p in (3, 4)} for group in ("SU3", "SO3")}
    path = Path(__file__).with_name("results.json")
    if args.write_new:
        with path.open("x", encoding="utf8") as f:
            json.dump(report, f, indent=2, allow_nan=False)
            f.write("\n")
    else:
        compare(report, json.loads(path.read_bytes()))
    print(json.dumps(report))
