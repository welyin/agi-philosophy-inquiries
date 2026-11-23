"""Independent finite matrices for partial port sharing, norm budgets and weak CP."""
import argparse
import cmath
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "independent_results.json"
I = np.eye(3, dtype=complex)
W = cmath.exp(2j * math.pi / 3)
T = I[:, [1, 2, 0]]
FOURIER = np.array([[W**(-j*k) for k in range(3)] for j in range(3)]) / math.sqrt(3)
assert np.allclose(T @ FOURIER, FOURIER @ np.diag([1, W, W**2]))


def read(relative):
    return json.loads((HERE / relative).read_text(encoding="utf-8"))


def close(x, y):
    return np.allclose(x, y, rtol=1e-10, atol=1e-13)


def phase_gate(vector, phi):
    u = np.array(vector, dtype=complex)
    u /= np.linalg.norm(u)
    return I + (cmath.exp(1j*phi) - 1)*np.outer(u, u.conj())


def relative_right(s):
    """Mode zero takes the smallest singular value; modes 1 and 2 remain untouched."""
    light, mid, heavy = s
    if heavy == 0:
        return {"path_gain": 0., "right_phase_rad": 0.}
    ratio = mid/heavy
    t = math.atan2(ratio*math.sin(math.pi/3), 1-ratio/2)
    gain = heavy / (2*math.sin(t+math.pi/3))
    phase = 2*math.asin(min(1., light/(2*gain))) - 2*t + 2*math.pi/3
    gate = phase_gate(FOURIER[:, 0], phase)
    a, b = gain, -gain*cmath.exp(2j*t)/W
    y = a*I + b*T @ gate
    assert np.linalg.matrix_rank(gate-I, tol=1e-12) <= 1
    assert close(np.linalg.norm(y @ FOURIER, axis=0), s)
    assert close(np.sort(np.linalg.svd(y, compute_uv=False)), s)
    mismatch = 2*math.sin(abs(cmath.phase(cmath.exp(1j*phase)))/4)
    assert gain*mismatch + 1e-14 >= max(0., (heavy-mid-light)/3)
    return {"path_gain": gain, "right_phase_rad": phase, "phase_insensitive_mismatch": mismatch,
            "misalignment_budget_relative_to_heavy": gain*mismatch/heavy}


def filter_gate(j, t):
    x = math.sqrt((4-t*t)/3)
    a, b = (x+t)/2, (x-t)/2
    return a*I - b*W**(-j)*T


def two_steps(s):
    scale = math.sqrt(s[-1])
    if not scale:
        return 0.
    first = scale * filter_gate(0, s[0]/s[-1])
    second = scale * filter_gate(1, s[1]/s[-1])
    assert close([np.linalg.norm(x, 2) for x in (first, second)], [scale]*2)
    assert close(np.linalg.norm(second @ first @ FOURIER, axis=0), s)
    return scale


def dataset(pairs):
    s = np.array([row[0] for row in pairs])
    e = np.array([row[1] for row in pairs])
    delta = max(0., s[2]-s[1]-s[0])/3
    lower = max(0., s[2]-s[1]-s[0]-2*sum(e))/3
    nearest = s + np.array([delta, delta, -delta])
    assert close(max(abs(nearest-s)), delta)
    assert nearest[-1] <= nearest[0]+nearest[1]+1e-13
    return {"one_step_op_norm_spectrum_distance": delta, "two_width_lower": lower,
            "two_step_minimax_norm": two_steps(s), "partial_right": relative_right(s)}


def weak_read(u, d):
    bases = []
    for gate, a, b in ((u,.87,.41+.32j), (d,.56,.65-.19j)):
        y = a*I + b*gate
        levels, vectors = np.linalg.eigh(y @ y.conj().T)
        assert np.min(np.diff(levels)) > .001
        bases.append(vectors)
    v = bases[0].conj().T @ bases[1]
    j = (v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag
    return {"abs_J": float(abs(j)), "probabilities": (abs(v)**2).tolist()}


def run():
    low = read("../1386/inputs.json")["yukawa"]
    high = read("../1397/inputs.json")["high_scale"]
    high_l = read("../1402/inputs.json")["scales"]["1e16"][-3:]
    samples = {"MZ": {s: low[s] for s in ("up", "down", "charged_leptons")},
               "1e16": {"up": high["up"], "down": high["down"], "charged_leptons": high_l}}
    result = {scale: {f: dataset(pairs) for f, pairs in family.items()}
              for scale, family in samples.items()}
    for edge in ([0,0,0],[0,0,1],[0,1,1],[1,1,1],[.5,1,2]):
        two_steps(edge)
        relative_right(edge)
    # Choose instead the eigenbasis of the actual up unitary for the CP witnesses.
    t = np.diag([1, W, W**2])
    one = phase_gate([1, 2j, 1+1j], .91)
    two = phase_gate([2-1j, 1+2j, 1], -1.33)
    rank1 = weak_read(t, t @ one)
    rank2 = weak_read(t, t @ one @ two)
    assert np.linalg.matrix_rank(one-I, tol=1e-11) == 1
    assert np.linalg.matrix_rank(one @ two-I, tol=1e-11) == 2
    assert rank1["abs_J"] < 1e-12 and rank2["abs_J"] > .0001
    assert np.max(np.min(rank1["probabilities"], axis=1)) > .001
    result["one_relative_direction"] = rank1
    result["two_relative_directions"] = rank2
    result["positive_lower_on_fixed_left_example"] = 2 - math.sqrt(3)
    return result


def equal(a, b):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(equal(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return len(a) == len(b) and all(equal(i, j) for i, j in zip(a, b))
    return math.isclose(a, b, rel_tol=2e-12, abs_tol=1e-13)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    flags = parser.add_mutually_exclusive_group()
    flags.add_argument("--write", action="store_true")
    flags.add_argument("--check", action="store_true")
    result = run()
    if parser.parse_args().write:
        OUTPUT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    else:
        assert equal(result, json.loads(OUTPUT.read_text(encoding="utf-8")))
        print("1411 partial right sharing, optimal two-stage norm and relative-rank CP verified")
