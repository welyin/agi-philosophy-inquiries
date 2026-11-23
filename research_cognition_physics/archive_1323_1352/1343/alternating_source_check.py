"""Sequential one-arm pulses; no point in the sequence opens both arms."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np


def generators(group):
    eye = np.eye(3)
    outer = lambda i, j: np.outer(eye[i], eye[j])
    if group == "rotation":
        return [1j*(outer(i, j)-outer(j, i)) for i, j in ((0, 1), (1, 2), (2, 0))]
    basis = []
    for i, j in itertools.combinations(range(3), 2):
        basis.extend([(outer(i, j)+outer(j, i))/2, 1j*(outer(i, j)-outer(j, i))/2])
    return basis+[np.diag([1., -1., 0.])/2, np.diag([1., 1., -2.])/(2*np.sqrt(3))]


def evolve(h, t):
    eigenvalues, columns = np.linalg.eigh(h)
    return (columns*np.exp(-1j*t*eigenvalues)) @ columns.conj().T


def read(group):
    basis = generators(group)
    casimir = np.trace(sum(t@t for t in basis)).real/3
    arms = []
    for left in (True, False):
        block = np.vstack([np.kron(t, np.eye(3)) if left else np.kron(np.eye(3), t)
                           for t in basis])/np.sqrt(casimir)
        h = np.zeros((9+len(block), 9+len(block)), dtype=complex)
        h[9:, :9], h[:9, 9:] = block, block.conj().T
        arms.append(h)
        assert abs(np.linalg.norm(h, ord=2)-1) < 1e-12
    theta, m = np.pi, 256
    target = evolve(arms[0]+arms[1], theta)[:, :9]
    pulse = evolve(arms[1], theta/m) @ evolve(arms[0], theta/m)
    seq = np.linalg.matrix_power(pulse, m)[:, :9]
    error = np.linalg.norm(seq-target, ord=2)
    assert error <= theta**2/m
    prob = lambda x: np.linalg.norm(x[9:])**2/9  # purifying partners retained by identity
    return dict(group=group, ideal=float(prob(target)), alternate=float(prob(seq)),
                worst_input_column_error=float(error), full_norm_upper=float(theta**2/m))


def run():
    a, b = read("unitary"), read("rotation")
    robust = .255-2*((np.pi+.005)**2/256 + 512*.00005)
    assert robust > .126 and a["alternate"]-b["alternate"] > robust
    return dict(arms=[a, b], pulse_count=512, guaranteed_gap=float(robust))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = Path(__file__).with_name("alternating_source_results.json")
    if args.write:
        with path.open("x", encoding="utf8") as f:
            json.dump(result, f, indent=2)
            f.write("\n")
    else:
        old = json.loads(path.read_text(encoding="utf8"))
        assert old.keys() == result.keys() and old["pulse_count"] == 512
        assert np.isclose(old["guaranteed_gap"], result["guaranteed_gap"])
        for row, ref in zip(result["arms"], old["arms"]):
            for k, v in row.items():
                assert v == ref[k] if k == "group" else np.isclose(v, ref[k])
    print(json.dumps(result))
