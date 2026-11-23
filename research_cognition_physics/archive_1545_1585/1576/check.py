"""Check the chiral characters on the whole weight lattice and the matching window."""
import argparse
from collections import Counter
from itertools import combinations, combinations_with_replacement, product
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ZERO = (0, 0, 0, 0)
COLOR = [(1, 0, 0, 0), (0, 1, 0, 0), (-1, -1, 0, 0)]
WEAK = [(0, 0, -1, 0), (0, 0, 1, 0)]


def charged(weights, charge):
    return [(*weight[:3], charge) for weight in weights]


def opposite(weight):
    return tuple(-x for x in weight)


def net_spectrum(weights):
    net = Counter(weights)
    net.subtract(opposite(v) for v in weights)
    return net


def pairs(weights, symmetric):
    choose = combinations_with_replacement if symmetric else combinations
    return [tuple(a + b for a, b in zip(v, u)) for v, u in choose(weights, 2)]


def full_family():
    table = json.loads((HERE / "../1573/inputs.json").read_text())["net_chiral_family"]
    states = []
    for dimension, weak_dimension, y in table.values():
        col = {3: COLOR, -3: [opposite(v) for v in COLOR], 1: [ZERO]}[dimension]
        doublet = WEAK if weak_dimension == 2 else [ZERO]
        states.extend((v[0], v[1], u[2], y) for v, u in product(col, doublet))
    return states


def check_representation(p, r, w, s, hypercharge):
    c, d = p + r, p - r
    e = 1 if c % 2 else 2
    assert not hypercharge or w == c
    col = charged(COLOR, -2 if hypercharge else 0)
    weak = charged(WEAK, 3 if hypercharge else 0)
    f = col * p + [opposite(v) for v in col] * r + [ZERO] * s
    f += weak * (p if hypercharge else w)
    if hypercharge:
        f += [opposite(v) for v in weak] * r
    n = len(f)
    coeff = [e + e * (c - 1) // 2, -e * (c - 1) // 2, e * (4 * c - n)]
    assert coeff[0] * (n - 4) + coeff[1] * (n + 4) + coeff[2] == 0
    actual = Counter()
    for k, rep in zip(coeff, (pairs(f, False), pairs(f, True), f)):
        for weight, mult in net_spectrum(rep).items():
            actual[weight] += k * mult
    expected = full_family()
    if not hypercharge:
        expected = charged(expected, 0)
    target = net_spectrum(expected)
    count = e * d * w
    assert all(actual[v] == count * target[v] for v in set(actual) | set(target))
    return count


def gauge_bounds(config):
    gauge = json.loads((HERE / config["gauge"]).read_text())
    width = config["width_count"]
    delta = config["abs_threshold_each_inverse_coupling"]
    vals = []
    for mu in config["matching_window_GeV"]:
        log = math.log(mu / gauge["mu_ref_GeV"]) / (2 * math.pi)
        for signs in product((-1, 1), repeat=4):
            inv = []
            for i, label in enumerate(("g2", "g3")):
                g, err = gauge[label]
                inv.append(4 * math.pi / (g + signs[i] * width * err) ** 2
                           - config["sm_one_loop_b2_b3"][i] * log + signs[i+2] * delta)
            vals.append(inv[0] / inv[1])
    return min(vals), max(vals)


def calculate():
    config = json.loads((HERE / "inputs.json").read_text())
    samples = 0
    for p, r, w, s in product(range(5), range(5), range(1, 5), (0, 2)):
        if not p + r:
            continue
        check_representation(p, r, w, s, False)
        samples += 1
    diagonal = 0
    for p, r, s in product(range(5), range(5), (0, 2)):
        if p + r:
            check_representation(p, r, p+r, s, True)
            diagonal += 1
    low, high = gauge_bounds(config)
    assert 3/5 < low < 1 < high < 3
    choices = []
    observed = config["three_families"]
    for w in range(1, observed + 1):
        for c in range(1, math.floor(w / low) + 1):
            for d in range(1, c+1):
                e = 1 if c % 2 else 2
                if ((c-d) % 2 == 0 and observed % (e*d*w) == 0
                        and low <= w / c <= high):
                    choices.append([c, w, d, e*d*w])
    assert choices == [[1, 1, 1, 1], [3, 3, 1, 3]]
    return {"bare_weight_lattice_checks": samples, "full_hypercharge_checks": diagonal,
            "one_loop_A2_over_A3_window": [low, high],
            "c_w_abs_d_min_families": choices,
            "SU15_example_net_full_families": check_representation(2, 1, 3, 0, True)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    flags = parser.add_mutually_exclusive_group(required=True)
    flags.add_argument("--check", action="store_true")
    flags.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = calculate()
    dest = HERE / "results.json"
    if args.write_new:
        with dest.open("x", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2)
            handle.write("\n")
    else:
        assert result == json.loads(dest.read_text())
    print(json.dumps(result))
