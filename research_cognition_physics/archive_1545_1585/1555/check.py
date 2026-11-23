"""Reproduce round 1555 using the numerical public likelihood, without plots."""
import argparse
import hashlib
import json
from math import pi
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def circular_margin(delta, interval):
    """Positive inside a closed arc; magnitude is angular distance to its ends."""
    lo, hi = np.asarray(interval) * pi
    mid = (lo + hi) / 2
    return float((hi - lo) / 2 - abs((delta - mid + pi) % (2 * pi) - pi))


def verify_original(path, grid):
    # This specific release stores a pickled NumPy dictionary. Verify the exact
    # pinned public bytes before permitting NumPy to deserialize that container.
    if hashlib.sha256(path.read_bytes()).hexdigest() != grid["sha256"]:
        raise ValueError("Unrecognized public file; refusing to deserialize")
    data = np.load(path, allow_pickle=True).item()
    x, m, q = (data[key] for key in ("sin2_theta23", "deltam232", "deltallh"))
    axes = [np.unique(axis) for axis in (x, m)]
    assert [len(axis) for axis in axes] == grid["source_shape"]
    idx = np.lexsort((m, x))
    np.testing.assert_array_equal(axes[0][:12], grid["grid_x"])
    np.testing.assert_array_equal(axes[1], grid["mass_squared_eV2"])
    np.testing.assert_array_equal(q[idx].reshape(grid["source_shape"])[:12], grid["q_rows"])


def compute(data):
    g = data["icecube"]
    x = np.asarray(g["grid_x"])
    q = np.asarray(g["q_rows"])
    m = np.asarray(g["mass_squared_eV2"])
    target = g["working_level"]
    assert q.shape == (len(x), len(m))
    assert np.all(np.diff(x) > 0) and np.all(np.diff(m) > 0)

    def profile(value):
        if value < x[0] or value > x[-1]:
            return None
        row = np.array([np.interp(value, x, column) for column in q.T])
        j = int(np.argmin(row))
        return float(row[j]), float(m[j])

    crossings = []
    for j, column in enumerate(q.T):
        for i, (left, right) in enumerate(zip(column, column[1:])):
            if left > target >= right:
                first = x[i] + (x[i + 1] - x[i]) * (left - target) / (left - right)
                crossings.append((float(first), float(m[j])))
    assert np.min(q[0]) > target and crossings
    first, crossing_mass = min(crossings)
    base, down, up = data["t2k_disappearance"]["NO"]
    old_range = [base - 2 * down, base + 2 * up]
    rows = []
    for name, delta, dp, theta, m1, mbb, old_chi in data["points"]:
        margin = [circular_margin(delta, arc) for arc in data["joint"]["delta_intervals_pi"]]
        width = data["theta_rounding"]
        low, high = theta - width, theta + width
        scan = [low, high] + [v for v in x if low < v < high]
        rounded = [profile(v) for v in scan]
        best = profile(theta)
        rounded_range = all(v is not None for v in rounded)
        rounded_min = min(v[0] for v in rounded) if rounded_range else None
        rows.append({
            "name": name, "delta_rad": delta, "sin2theta23": theta,
            "circular_margins_1_3_sigma_rad": margin,
            "all_rounded_phases_in_3sigma_outside_1sigma": margin[0] < -dp and margin[1] > dp,
            "inside_joint_theta23_1sigma": any(lo <= theta <= hi for lo, hi in data["joint"]["theta23_1sigma"]),
            "inside_old_pure_disappearance_two_widths": old_range[0] <= theta <= old_range[1],
            "ice_q_and_profile_mass": best, "ice_min_q_with_theta_rounding": rounded_min,
            "min_shift_theta_to_reach_gridded_q4_from_printed_value": first - theta,
            "min_q_reduction_for_rounded_point_when_covered": None if rounded_min is None else rounded_min - target
        })
    assert all(v["all_rounded_phases_in_3sigma_outside_1sigma"] for v in rows)
    assert not any(v["inside_joint_theta23_1sigma"] for v in rows)
    assert [v["inside_old_pure_disappearance_two_widths"] for v in rows] == [True, False, False, True]
    assert rows[1]["ice_q_and_profile_mass"] is None and rows[1]["sin2theta23"] + data["theta_rounding"] < x[0]
    assert all(v["ice_min_q_with_theta_rounding"] > target for v in (rows[0], rows[2], rows[3]))
    return {
        "rows": rows, "old_pure_disappearance_two_widths": old_range,
        "ice_first_published_grid_edge_and_q": [float(x[0]), float(q[0].min())],
        "ice_first_x_with_q_at_most_4": first, "ice_mass_at_first_crossing_eV2": crossing_mass,
        "input_sha256": hashlib.sha256((HERE / "inputs.json").read_bytes()).hexdigest()
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    parser.add_argument("--verify-public-map", type=Path)
    options = parser.parse_args()
    inputs = json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))
    if options.verify_public_map:
        verify_original(options.verify_public_map, inputs["icecube"])
    results = json.dumps(compute(inputs), ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    dest = HERE / "results.json"
    if options.write_new:
        with dest.open("x", encoding="utf-8", newline="\n") as out:
            out.write(results)
    else:
        assert json.loads(results) == json.loads(dest.read_text(encoding="utf-8")), "Result mismatch"
    print("1555 verified")
