"""Relate a W width count to the generator traces measured with four jets."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def generators():
    matrices = []
    for i, j in ((0, 1), (0, 2), (1, 2)):
        x = np.zeros((3, 3), complex)
        x[i, j] = x[j, i] = 0.5
        y = np.zeros_like(x)
        y[i, j], y[j, i] = -0.5j, 0.5j
        matrices.extend([x, y])
    t3 = np.diag([0.5, -0.5, 0])
    t8 = np.diag([1, 1, -2]) / math.sqrt(12)
    return matrices + [t3, t8]


def traces(ts):
    c = sum(t @ t for t in ts)
    denominator = np.trace(c @ c).real
    gram = np.array([[np.trace(a @ b).real for b in ts] for a in ts])
    comm = [a@b-b@a for a in ts for b in ts]
    return {"y": float(np.trace(gram@gram) / denominator),
            "x": float(sum(np.vdot(d, d).real for d in comm) / denominator),
            "effective_matter_dimension": float(np.trace(c).real**2 / denominator)}


def solve():
    cfg = json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))
    b0, err = cfg["hadronic_B"]
    multiple = cfg["widths"]
    b_lo, b_hi = b0 - multiple*err, b0 + multiple*err
    def ratio(b):
        return cfg["lepton_channels"]*b/(1-b)
    s_lo, s_hi = cfg["accessible_quark_row_norm"]
    k_lo, k_hi = cfg["inclusive_relative_transfer"]
    d_interval = [ratio(b_lo)/(s_hi*k_hi), ratio(b_hi)/(s_lo*k_lo)]
    allowed = list(range(math.ceil(d_interval[0]), math.floor(d_interval[1])+1))
    assert allowed == [3]

    opal, stat, systematic = cfg["OPAL_y"]
    y_max = opal + multiple*(stat+systematic)
    sharp = (3 - math.sqrt(3))/2
    t_min = 2*math.sqrt(3)-3
    ts = generators()
    u2 = [ts[0], ts[1], ts[-2]]
    blocks = {}
    for name, t in (("no_abelian", 0), ("minimum", t_min), ("equal_metric", 1), ("large_abelian", 12)):
        value = traces(u2 + [math.sqrt(t)*ts[-1]])
        assert abs(value["y"] - (6+2*t**2)/(9+2*t+t**2)) < 1e-12
        blocks[name] = {"coupling_weight_ratio": t, **value}
    assert abs(blocks["minimum"]["y"] - sharp) < 1e-12
    abelian = traces(ts[-2:])
    full = traces(ts)
    so3 = traces([ts[i] for i in (1,3,5)])
    assert abs(full["y"] - 3/8) < 1e-12 and abs(full["x"] - 9/4) < 1e-12
    assert abs(so3["x"]-1) < 1e-12 and abs(so3["y"]-1) < 1e-12
    assert abs(abelian["y"]-1.5) < 1e-12

    fit, st, sy = cfg["ALEPH_x"]
    tests = {}
    for err_transfer in (*cfg["examples_of_additive_y_transfer"], 0.06):
        limit = y_max + err_transfer
        tests[str(err_transfer)] = {"true_y_upper": limit, "excludes_all_proper": limit < sharp,
                    "excludes_all_proper_if_uniform_total_radiation": limit < 3/4}
    return {"W_d_interval": d_interval, "W_integer_dimension": allowed,
            "per_colour_source_limits_for_d_eq4_max_or_d_eq2_min":
               {"d4_max": ratio(b_hi)/4, "d2_min": ratio(b_lo)/2},
            "OPAL_y_upper_two_width": y_max,
            "proper_group_bounds":{"uniform_total_radiation": .75, "all_invariant_metrics": sharp},
            "max_additional_y_error_strict":{"uniform_total_radiation": .75-y_max,
                        "all_invariant_metrics": sharp-y_max}, "transfer_examples": tests,
            "compact_subalgebra_kernel_lower_for_dim_5_to_7":
                [m - (8-m)*(7-m)//2 for m in (5,6,7)],
            "SU3": full, "rotation_SO3": so3, "two_diagonals": abelian, "block_U2": blocks,
            "ALEPH_independent_x_window": [fit-multiple*(st+sy), fit+multiple*(st+sy)]}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    opts = p.add_mutually_exclusive_group()
    opts.add_argument("--write-new", action="store_true")
    opts.add_argument("--check", action="store_true")
    flags = p.parse_args()
    result = solve()
    filename = HERE / "results.json"
    if flags.write_new:
        with filename.open("x", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2)
            handle.write("\n")
    else:
        def compare(a, b):
            if isinstance(a, dict):
                return a.keys() == b.keys() and all(compare(v, b[k]) for k, v in a.items())
            if isinstance(a, list):
                return len(a) == len(b) and all(compare(x, y) for x, y in zip(a,b))
            if isinstance(a, float):
                return math.isclose(a, b, abs_tol=1e-12, rel_tol=1e-12)
            return a == b
        assert compare(result, json.loads(filename.read_text(encoding="utf-8")))
        print("1507: W multiplicity, closed-subgroup threshold and independent jet check agree.")
