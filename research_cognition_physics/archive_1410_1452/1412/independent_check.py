"""Independent joint bound on the symmetry and texture-zero repairs of a named mass block."""
import argparse
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def joint_bound(q, left, cap):
    target = q*q*(1-left)/(left + q*q*(1-left))
    radius = cap/(1-q*q)
    right = float(np.clip(target, max(0., left-radius), min(1., left+radius)))
    a, b = np.sqrt([1-left, left])
    c, d = np.sqrt([1-right, right])
    block = np.array([[a, -b], [b, a]]) @ np.diag([q, -1.]) @ np.array([[c, d], [-d, c]])
    zero_error = abs(q*a*c - b*d)
    assert abs(abs(block[0, 0])-zero_error) < 1e-14
    assert abs(abs(block[0, 1]**2-block[1, 0]**2) - (1-q*q)*abs(left-right)) < 1e-14
    assert np.allclose(np.linalg.svd(block, compute_uv=False), [1., q])
    return {"right_sine_squared": right, "smallest_diagonal": float(zero_error)}


def compute():
    data = json.loads((HERE.parent/"1386"/"inputs.json").read_text(encoding="utf-8"))
    u, c = data["yukawa"]["up"][:2]
    q = (u[0]+2*u[1])/(c[0]-2*c[1])
    t = .083-2*.008
    left = t*t/(1+t*t)
    cases = {str(cap): joint_bound(q, left, cap) for cap in [0., .0005, .001, .002, .004, .006]}
    assert cases["0.001"]["smallest_diagonal"] > .001
    return {"round": 1412, "q_upper": q, "left_sine_squared_lower": left, "joint_frontier": cases,
            "simultaneous_per_mille_repairs_excluded": True}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = compute()
    if args.check:
        assert result == json.loads((HERE/"independent_results.json").read_text(encoding="utf-8"))
        print("1412 independent mass block and joint frontier verified")
    else:
        print(json.dumps(result, indent=2))
