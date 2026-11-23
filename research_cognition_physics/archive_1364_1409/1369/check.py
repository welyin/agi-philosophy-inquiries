"""Conditional test of shared matter weights. All running bounds are one-loop."""

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
LOOP = 16 * math.pi ** 2
SLOPES = np.array([-41 / 6, 19 / 6, 7]) / (8 * math.pi ** 2)


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def endpoints(point):
    return [point[0] - point[1], point[0] + point[1]]


def gauge_match(gauge, mz):
    inverse = 1 / (np.array(gauge) ** 2)
    inverse[0] *= 5 / 3  # x_Y = 5/(3 g_1^2)
    functional = np.array([1, -3, 4 / 3])
    t = -np.dot(functional, inverse) / np.dot(functional, SLOPES)
    matched = inverse + t * SLOPES
    assert np.all(matched > 0)
    assert abs(np.dot(functional, matched)) < 1e-13
    kinetic_ratio = 4 * matched[1] / matched[2] - 3
    assert kinetic_ratio > 0
    return dict(t=t, mu_GeV=mz * math.exp(t), wl_wq=kinetic_ratio,
                predicted_flavor_weight=3 / kinetic_ratio)


def p(r):
    return r + 1 / r


def inverse_p(value):
    return (value + math.sqrt(value * value - 4)) / 2


def minimum_weight(rd, re):
    return (math.sqrt(re) - 1 / math.sqrt(re)) / (math.sqrt(rd) - 1 / math.sqrt(rd))


def run():
    data = read(HERE / "inputs.json")
    inherited = read(HERE / data["light_down_lepton"])
    light = inherited["points"]["MZ"]
    centers = np.array(light["central"]) * inherited["units"]
    widths = np.array(light["half_width"]) * inherited["units"]
    low, high = centers - widths, centers + widths
    up_data = read(HERE / data["up_charm"])["MZ"]
    gauge = [data[key] for key in ["g1_GUT", "g2", "g3"]]
    central = gauge_match([v[0] for v in gauge], data["mu_ref_GeV"])
    corners = [gauge_match(v, data["mu_ref_GeV"])
               for v in itertools.product(*(endpoints(g) for g in gauge))]
    box = {k: [min(c[k] for c in corners), max(c[k] for c in corners)] for k in central}
    tmax = box["t"][1]
    assert 10 < box["t"][0] < tmax < 22
    uq0 = (sum(endpoints(x)[1] ** 2 for x in
               [up_data["up"], up_data["charm"], data["yt"], data["yb"]])
           + high[0] ** 2 + high[1] ** 2)
    el0 = high[2] ** 2 + high[3] ** 2 + endpoints(data["ytau"])[1] ** 2
    g20 = endpoints(data["g2"])[0]
    g30 = endpoints(data["g3"])[0]

    def g_sq(g0, b, t):
        return 1 / (g0 ** -2 + b * t / (8 * math.pi ** 2))

    def gauge_damping(t):
        return 16 * g_sq(g30, 7, t) + 4.5 * g_sq(g20, 19 / 6, t)

    l_cap = 0.001
    l_on_barrier = el0 * math.exp((6 + 5 * l_cap) * 22 / LOOP)
    quark_at_10 = uq0 * math.exp(
        ((9 + 2 * l_cap) * 10
         - 16 * (8 * math.pi ** 2 / 7) * math.log1p(7 * g30 ** 2 * 10 / (8 * math.pi ** 2)))
        / LOOP
    )
    assert uq0 < 1 and l_on_barrier < l_cap
    assert gauge_damping(10) > 9 + 2 * l_cap
    assert quark_at_10 < 0.8 and gauge_damping(22) > 9 * 0.8 + 2 * l_cap
    rd_initial_ceiling = high[1] / low[0]
    re_floor = low[3] / high[2]
    integrated_quark_cap = 10 + 0.8 * (tmax - 10)
    rd_ceiling = rd_initial_ceiling * math.exp(1.5 * integrated_quark_cap / LOOP)
    w_floor = minimum_weight(rd_ceiling, re_floor)
    w_high = box["predicted_flavor_weight"][1]
    assert 1 < box["predicted_flavor_weight"][0] and w_high < w_floor
    re_ceiling = inverse_p(2 + w_high ** 2 * (p(rd_ceiling) - 2))
    assert re_ceiling < re_floor
    return dict(
        one_loop_central=central, one_loop_corner_bounds=box,
        barrier=dict(U0_max=uq0, E0_max=el0, E_until_22_upper=l_on_barrier,
                     gauge_damping_lower_at_10=gauge_damping(10), U_at_10_upper=quark_at_10,
                     gauge_damping_lower_at_22=gauge_damping(22)),
        same_scale=dict(rd_initial_upper=rd_initial_ceiling, rd_upper=rd_ceiling,
                        re_lower=re_floor, flavor_weight_needed_at_least=w_floor,
                        predicted_re_at_most=re_ceiling,
                        min_bridge_multiplier=w_floor / w_high,
                        max_preparation_to_kinetic_ratio_of_wl_over_wq=w_high / w_floor,
                        minimum_Q_to_L_hierarchy_deficit=re_floor / re_ceiling)
    )


def same(actual, saved):
    if isinstance(actual, dict):
        return actual.keys() == saved.keys() and all(same(v, saved[k]) for k, v in actual.items())
    if isinstance(actual, list):
        return len(actual) == len(saved) and all(same(a, b) for a, b in zip(actual, saved))
    return math.isclose(actual, saved, rel_tol=1e-10, abs_tol=1e-12)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    output = run()
    target = HERE / "results.json"
    if args.write_new:
        with target.open("x", encoding="utf-8") as out:
            json.dump(output, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        assert same(output, read(target))
        print("PASS: gauge-only crossing, all-corner bounds and common-scale Yukawa test.")

