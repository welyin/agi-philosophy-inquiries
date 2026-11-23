"""Quark mass shape versus the b-tag fraction. No CKM fit input."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def bounds(up, down):
    g = math.prod(up * down) ** (1 / 3)
    u2, u3 = up[1:]
    d2, d3 = down[1:]
    simple = (u2 * d2 / g) ** 2
    x = (g / u2) ** 2 - d2**2
    y = d3**2 - d2**2
    z = (g / u2) ** 2 - (g / u3) ** 2
    assert 0 < x < min(y, z)
    sharp = (1 - x / y) * (1 - x / z)
    cosine, sine = math.sqrt(sharp), math.sqrt(1 - sharp)
    v = np.array([[1, 0, 0], [0, cosine, -sine], [0, sine, cosine]])
    singular = np.linalg.svd(up[:, None] * v * down, compute_uv=False)
    assert abs(singular[1] / g - 1) < 1e-10
    assert 0 < sharp < simple
    return simple, sharp, singular[1] / g


def run():
    cfg = json.loads((HERE / "inputs.json").read_text(encoding="utf8"))
    masses = json.loads((HERE / cfg["mass_input"]).read_text(encoding="utf8"))["yukawa"]
    up, down = (np.array(masses[label]) for label in ("up", "down"))
    lower = cfg["cms_lower"]
    signs = np.array([-1, 1, -1])
    cases = []
    for width in cfg["box_widths"]:
        u, d = (a[:, 0] + a[:, 1] * width * signs for a in (up, down))
        loose, sharp_at_corner, witness = bounds(u, d)
        cases.append({
            "width": width, "uniform_upper_p": loose,
            "sharp_at_loose_corner_only": sharp_at_corner,
            "attainment_middle_over_g": round(witness, 10),
            "minimum_relative_weight_using_uniform_bound": lower * (1-loose) / (loose * (1-lower)),
            "upper_R_by_relative_weight": [w * loose / (1-loose+w*loose)
                                           for w in cfg["relative_bottom_weights"]],
        })
    box_bound = cases[2]["uniform_upper_p"]
    # Relative drift is defined as the ratio of the normalized middle singular
    # value between MZ and a reference scale with an exactly geometric spectrum.
    p_min = lower / (lower + cfg["relative_bottom_weights"][1] * (1 - lower))
    ratio = math.sqrt(box_bound / p_min)
    assert cases[2]["upper_R_by_relative_weight"][-1] < lower < cases[3]["uniform_upper_p"]
    return {"round": 1443, "cases": cases, "mass_only_inputs": 6,
            "CKM_fit_inputs": 0, "CMS_lower": lower,
            "p_min_for_1.05": p_min, "max_normalized_middle_for_1.05": ratio,
            "min_normalized_drift_from_exact_UV_for_1.05": 1 - ratio}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write-new", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    dest = HERE / "results.json"
    if args.write_new:
        with dest.open("x", encoding="utf8") as handle:
            json.dump(result, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
    else:
        assert result == json.loads(dest.read_text(encoding="utf8"))
    print("1443 OK", result["cases"][2]["uniform_upper_p"])
