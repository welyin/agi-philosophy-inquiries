"""Rate-only inversion and reciprocal charged-current transfer bounds."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def ratio_windows(inputs, widths):
    bounds = {}
    for label in ("PIENU", "NA62"):
        row = inputs[label]
        error = row.get("total_width")
        if error is None:
            error = math.hypot(row["stat"], row["syst"])
        value, theory = row["measurement"], row["standard"]
        theory_error = widths * row["theory_width"]
        bounds[label] = [
            math.sqrt((value - widths * error) / (theory + theory_error)),
            math.sqrt((value + widths * error) / (theory - theory_error)),
        ]
    tau = inputs["Belle_II"]
    value, error = tau["leptonic_tau_gmu_over_ge"], widths * tau["published_width"]
    bounds["Belle_II"] = [1 / (value + error), 1 / (value - error)]
    return bounds


def intersection(bounds):
    return [max(row[0] for row in bounds.values()), min(row[1] for row in bounds.values())]


def first_width(inputs, fixed=None):
    low, high = 0.0, 2.0
    for _ in range(60):
        mid = (low + high) / 2
        left, right = intersection(ratio_windows(inputs, mid))
        feasible = left <= right if fixed is None else left <= fixed <= right
        if feasible:
            high = mid
        else:
            low = mid
    return high


def strength_bounds(interval, projection, initial_capacity, target):
    gains = [1 + projection * (r * r - 1) for r in interval]
    capacity = max(gains[1], 1 / gains[0])
    assert capacity > 1
    visits = math.ceil(math.log(target / initial_capacity) / (2 * math.log(capacity)))
    return {
        "largest_one_return_condition": capacity,
        "mass_ratio_upper_after_four": initial_capacity * capacity ** 8,
        "necessary_return_count_per_leg": visits,
    }


def calculate():
    inputs = json.loads((HERE / "inputs.json").read_bytes())
    geometry = json.loads((HERE / inputs["active_geometry"]).read_bytes())
    target = geometry["required_mass_ratio_lower"]
    tables = {}
    for width in (1, 2):
        channels = ratio_windows(inputs, width)
        common = intersection(channels)
        assert common[0] <= common[1]
        endpoints = {}
        for name in inputs["endpoints"]:
            endpoint = geometry[name]
            projection = endpoint["geometry"]["electron_plane_probability"]
            base = endpoint["positive_compressed_propagation"]["unfiltered_capacity"]
            endpoints[name] = {
                "any_plane": strength_bounds(common, 1, base, target),
                "registered_plane": strength_bounds(common, projection, base, target),
            }
            if common[0] < 1:
                endpoints[name]["registered_plane_electron_suppressed"] = strength_bounds(
                    [common[0], 1], projection, base, target
                )
        mesons = intersection({key: channels[key] for key in ("PIENU", "NA62")})
        tables[str(width)] = {
            "channels": channels,
            "common": common,
            "tau_gmu_over_ge_using_mesons_only": [1/mesons[1], 1/mesons[0]],
            "endpoints": endpoints,
        }
    row_overlap = 0.25
    w = np.array([[1.0, 0, 0], [row_overlap, math.sqrt(1-row_overlap**2), 0], [0, 0, 1]])
    gram = w.T @ w
    assert np.allclose(np.diag(w @ w.T), 1)
    assert np.isclose(np.linalg.cond(gram), (1+row_overlap)/(1-row_overlap))
    result = {
        "required_mass_ratio": target,
        "common_r_least_width": first_width(inputs),
        "r_equals_one_least_width": first_width(inputs, 1),
        "windows": tables,
        "equal_rows_unseen_overlap_example": {
            "overlap": row_overlap,
            "relative_inclusive_rates": np.diag(w @ w.T).tolist(),
            "return_eigenvalues": np.linalg.eigvalsh(gram).round(12).tolist(),
        },
    }
    assert tables["1"]["common"][0] > 1 and tables["2"]["common"][0] < 1
    for item in tables["2"]["endpoints"].values():
        assert item["any_plane"]["mass_ratio_upper_after_four"] < target
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = calculate()
    output = HERE / "results.json"
    if args.check:
        if json.loads(output.read_bytes()) != result:
            raise SystemExit("saved output differs")
        print("1764 reciprocal transfer checks passed")
    else:
        output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
