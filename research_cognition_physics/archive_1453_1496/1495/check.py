"""Finite real flavor rotations and the reactor overlap of the third mass ray."""

import argparse
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


def amplitude(p):
    return 4.0 * p * (1.0 - p)


def ray_angle(p):
    return math.asin(math.sqrt(p))


def same_report(computed, stored):
    if isinstance(computed, dict):
        return computed.keys() == stored.keys() and all(
            same_report(value, stored[key]) for key, value in computed.items()
        )
    if isinstance(computed, list):
        return len(computed) == len(stored) and all(
            same_report(a, b) for a, b in zip(computed, stored)
        )
    return math.isclose(computed, stored, rel_tol=1e-11, abs_tol=1e-12)


def certificate():
    data = json.loads((ROOT / "inputs.json").read_text(encoding="utf-8"))["reactor"]
    delta = data["working_half_widths"] * data["published_width"]
    lo, hi = data["sin2_2theta13"] - delta, data["sin2_2theta13"] + delta
    small = [(1.0 - math.sqrt(1.0 - a)) / 2.0 for a in [lo, hi]]
    windows = [small, [1.0 - small[1], 1.0 - small[0]]]
    # Order one for TS is impossible because ord(T) != ord(S).
    axial = [0.0, 1 / 3, 2 / 3, (3 - math.sqrt(5)) / 6,
             (3 + math.sqrt(5)) / 6, 1.0]
    all_rows = sorted({round(x, 15) for v in axial for x in (v, (1.0 - v) / 2.0)})
    candidates = []
    for p in all_rows:
        for left, right in windows:
            closest = min(max(left, p), right)
            candidates.append((abs(ray_angle(closest) - ray_angle(p)), p, closest))
    gap, nearest, target = min(candidates)
    positive_amplitudes = [amplitude(p) for p in all_rows if 0.0 < p < 1.0]
    # Explicit nearest real branch: icosahedron edge midpoint and cyclic coordinate rotation.
    phi = (1 + math.sqrt(5)) / 2
    w = np.array([1.0, 1.0 + phi, phi])
    w /= np.linalg.norm(w)
    omega = np.exp(2j * np.pi / 3)
    e = np.array([1, omega, omega**2]) / np.sqrt(3)
    real_p = abs(np.vdot(e, w))**2
    assert math.isclose(real_p, nearest, abs_tol=1e-13)
    s = 2 * np.outer(w, w) - np.eye(3)
    t = np.roll(np.eye(3), 1, axis=0)
    assert np.allclose(np.linalg.matrix_power(t @ s, 5), np.eye(3))
    # Reach the nearest allowed target on the projective geodesic, moving only the mass ray.
    c = np.vdot(e, w)
    perpendicular = (w - c * e) / math.sqrt(1.0 - abs(c)**2)
    adjusted = (c / abs(c)) * math.sqrt(target) * e + math.sqrt(1.0 - target) * perpendicular
    assert math.isclose(abs(np.vdot(e, adjusted))**2, target, abs_tol=1e-13)
    overlap_cost = math.acos(min(1.0, abs(np.vdot(w, adjusted))))
    assert math.isclose(overlap_cost, gap, abs_tol=1e-13)
    # Finite complex counterexample. The denominator is selected for this example only.
    n, k = 43, 19
    phase = np.exp(2j * math.pi * k / n)
    unitary = np.array([[0, phase, 0], [phase.conjugate(), 0, 0], [0, 0, -1]])
    axis = np.array([phase, 1, 0]) / np.sqrt(2)
    unit_row = np.ones(3) / np.sqrt(3)
    complex_p = abs(np.vdot(unit_row, axis))**2
    assert np.allclose(unitary @ unitary, np.eye(3)) and np.allclose(unitary @ axis, axis)
    assert np.isclose(np.linalg.det(unitary), 1) and lo <= amplitude(complex_p) <= hi
    return {
        "reactor_amplitude_work_range": [lo, hi], "p3e_windows": windows,
        "real_axial_probabilities": sorted(axial), "all_real_electron_probabilities": all_rows,
        "smallest_nonzero_real_amplitude": min(positive_amplitudes),
        "nearest_real_p": float(real_p), "nearest_target_p": target,
        "joint_ray_angle_floor_degrees": math.degrees(gap),
        "one_moving_projector_norm_floor": math.sin(gap),
        "either_ray_if_equal_budget_projector_norm_floor": math.sin(gap / 2),
        "minimum_additive_amplitude_change": min(
            abs(amplitude(p) - min(hi, max(lo, amplitude(p)))) for p in all_rows
        ),
        "complex_example": {"n": n, "k": k, "electron_probability": float(complex_p),
                            "oscillation_amplitude": float(amplitude(complex_p))},
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = certificate()
    rendered = json.dumps(result, ensure_ascii=False, indent=2)
    output = ROOT / "results.json"
    if args.write:
        with output.open("x", encoding="utf-8") as stream:
            stream.write(rendered + "\n")
    elif not same_report(result, json.loads(output.read_text(encoding="utf-8"))):
        raise ValueError("Saved results differ from the recomputed certificate")
    print(rendered)
