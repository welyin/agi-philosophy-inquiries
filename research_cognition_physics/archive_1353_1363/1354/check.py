"""Charge-conserving delivery with a finite reference. No postselection in totals."""
import argparse
import json
from pathlib import Path

import numpy as np


def selected_weights(p, k, rate):
    """Maximize the target overlap at fixed success, including zero supports."""
    a = np.roll(p, k)
    weights = np.zeros_like(p)
    live = p > 0
    remaining = rate
    while np.any(live):
        scale = remaining / p[live].sum()
        clipped = live & (scale * p >= a)
        if not np.any(clipped):
            weights[live] = scale * p[live]
            remaining = 0
            break
        weights[clipped] = a[clipped]
        remaining -= weights[clipped].sum()
        live &= ~clipped
    for j in np.flatnonzero(p == 0):
        weights[j] = min(max(remaining, 0), a[j])
        remaining -= weights[j]
    np.testing.assert_allclose(weights.sum(), rate, atol=1e-13)
    return weights


def trial(p, k, rate):
    ref = np.sqrt(p)
    caps = np.roll(p, k)
    w = selected_weights(p, k, rate)
    gain = np.sqrt(np.divide(w, caps, out=np.zeros_like(p), where=caps > 0))
    shift = np.roll(np.eye(len(p)), k, axis=0)
    success = np.diag(gain) @ shift
    # A failure leaves the input payload in place and keeps its charge.
    failure = np.diag(np.sqrt(np.maximum(0, 1 - np.diag(success.T @ success))))
    np.testing.assert_allclose(success.T @ success + failure.T @ failure, np.eye(len(p)), atol=1e-12)
    phase = np.exp(2j * np.pi / len(p))
    group = np.diag(phase ** np.arange(len(p)))
    np.testing.assert_allclose(group @ success, success @ (phase**k * group), atol=1e-12)
    out = success @ ref
    np.testing.assert_allclose(out @ out, rate, atol=1e-12)
    fidelity = float((ref @ out)**2 / rate)
    assert -1e-12 <= fidelity <= 1 + 1e-12
    fidelity = min(1.0, max(0.0, fidelity))
    return {"success": rate, "conditional_fidelity": fidelity, "kept_sector_weights": w.tolist(),
            "success_and_swap_minus": rate*(1-fidelity)/2, "failed": 1-rate}


def calculate():
    p = np.array([1/2, 1/3, 1/6])
    shift = np.roll(np.eye(3), 1, axis=0)
    phi = np.sqrt(p)
    repeated = [float(abs(phi @ np.linalg.matrix_power(shift, n) @ phi)**2) for n in (1, 2, 3)]
    matched = np.diag(np.sqrt(p)).reshape(-1)
    extended = np.kron(shift, np.eye(3))
    entangled = [float(abs(matched @ np.linalg.matrix_power(extended, n) @ matched)**2)
                 for n in (1, 2, 3)]
    np.testing.assert_allclose(entangled, [0, 0, 1], atol=1e-14)
    flat_purification = np.eye(3)/np.sqrt(3)
    moved = shift @ flat_purification
    np.testing.assert_allclose(moved @ moved.T, flat_purification @ flat_purification.T)
    np.testing.assert_allclose(np.vdot(moved, flat_purification), 0, atol=1e-14)
    # Test disjoint supports and noncoprime shifts for a separate Z4 register.
    disjoint = trial(np.array([.5,0,.5,0]), 1, .75)
    islands = trial(np.array([.1,.4,.1,.4]), 2, 1)
    np.testing.assert_allclose([disjoint["conditional_fidelity"], islands["conditional_fidelity"]], [0,1])
    return {"round":1354, "p":p.tolist(), "shift":1, "max_exact_return_success":min(np.roll(p,1)/p),
            "controlled_reference": [trial(p,1,r) for r in (1/3,.5,.75,1)],
            "same_unfiltered_gate_return_fidelities":repeated, "charge_matched_purification":entangled,
            "readout_abs_error":.005}


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    opts=parser.parse_args()
    data=calculate()
    path=Path(__file__).with_name("results.json")
    if opts.write:
        with path.open("x", encoding="utf8") as stream:
            json.dump(data,stream,ensure_ascii=False,indent=2)
            stream.write("\n")
    else:
        archived=json.loads(path.read_text(encoding="utf8"))
        assert data.keys()==archived.keys()
        for k in data.keys()-{"controlled_reference"}:
            np.testing.assert_allclose(data[k],archived[k],atol=1e-12,rtol=1e-12)
        for actual,stored in zip(data["controlled_reference"],archived["controlled_reference"],strict=True):
            assert actual.keys()==stored.keys()
            for k in actual:
                np.testing.assert_allclose(actual[k],stored[k],atol=1e-12,rtol=1e-12)
    print(json.dumps(data))
