"""Inverse test for two coherent projectors at one common pulse angle."""
import argparse
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def rotation(i, j, sine, phase=0.):
    matrix = np.eye(3, dtype=complex)
    matrix[i, i] = matrix[j, j] = np.sqrt(1-sine**2)
    matrix[i, j] = sine*np.exp(-1j*phase)
    matrix[j, i] = -matrix[i, j].conjugate()
    return matrix


def compute():
    source = json.loads((HERE.parent/"1386"/"inputs.json").read_text(encoding="utf-8"))
    a, b, c, delta = [v[0] for v in source["quark_sines_and_phase"]]
    v = rotation(1, 2, c) @ rotation(0, 2, b, delta) @ rotation(0, 1, a)
    t = np.diag(np.exp(2j*np.pi*np.arange(3)/3))
    w = t.conj().T @ v @ t @ v.conj().T
    eigen = np.linalg.eigvals(w)
    cases = []
    for ref in np.argsort(np.angle(eigen)):
        pair = np.delete(eigen, ref)/eigen[ref]
        for sign in [1, -1]:
            s = sign*np.exp(.5j*np.angle(np.prod(pair)))
            overlap = (sum(pair)-2*s)/(s-1)**2
            assert abs(overlap.imag) < 1e-8
            if not 0 <= overlap.real <= 1:
                continue
            z = float(overlap.real)
            direction = np.sqrt([z, 1-z])
            gate1 = np.diag([s, 1.])
            gate2 = np.eye(2) + (s-1)*np.outer(direction, direction)
            product = gate1 @ gate2
            assert np.allclose(np.poly(product), np.poly(pair))
            comm_norm = np.linalg.norm(product-gate2@gate1, 2)
            assert np.isclose(comm_norm, abs(s-1)**2*np.sqrt(z*(1-z)))
            cases.append({"reference_phase": float(np.angle(eigen[ref])), "pulse_phase": float(np.angle(s)),
                          "direction_overlap_squared": z, "order_difference_norm": float(comm_norm)})
    assert cases and all(v["order_difference_norm"] > 1e-7 for v in cases)
    return {"round": 1413, "ordered_cubic_clock_common_phase_free": True,
            "smallest_absolute_equal_angle": min(abs(v["pulse_phase"]) for v in cases),
            "equal_pulse_solutions": cases, "free_unequal_spectral_gates_commute": True}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = compute()
    if args.check:
        saved = json.loads((HERE/"independent_results.json").read_text(encoding="utf-8"))
        for x, y in zip(result["equal_pulse_solutions"], saved["equal_pulse_solutions"], strict=True):
            assert x.keys() == y.keys() and np.allclose(list(x.values()), list(y.values()), rtol=1e-8, atol=1e-10)
        assert np.isclose(result["smallest_absolute_equal_angle"], saved["smallest_absolute_equal_angle"])
        print("1413 independent phase budgets and gate order verified")
    else:
        print(json.dumps(result, indent=2))
