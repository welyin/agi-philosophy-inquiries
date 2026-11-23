"""Integer four-dimensional checks for anomaly cancellation with repartitioned index parents."""
import argparse
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPECIES = ("Q", "Uc", "Dc", "L", "Ec", "Nc")
# All fields are left-handed. y6=6Y, and T(fund)=1/2.
FIELDS = ((3, 2, 1), (3, 1, -4), (3, 1, 2), (1, 2, -3), (1, 1, 6), (1, 1, 0))


def read(rel):
    return json.loads((HERE / rel).read_text(encoding="utf-8-sig"))


def anomalies(generations):
    totals = [0] * 6
    for (colour, weak, y6), charges in zip(FIELDS, generations):
        for x in charges:
            dim = colour * weak
            coefficients = (weak * int(colour == 3) * x, colour * int(weak == 2) * x,
                            dim * y6 ** 2 * x, dim * y6 * x ** 2, dim * x,
                            dim * x ** 3)
            totals = [a + b for a, b in zip(totals, coefficients)]
    return totals


def signed_squares(cost):
    # Each unit of integral Dynkin-index weight can be split for the bound.
    if cost == 0:
        return {0}
    result = set()
    for n in range(1, cost + 1):
        for value in signed_squares(cost - n):
            result.update([value + n * n, value - n * n])
    return result


def connected_weak_edges():
    cfg = read("../1476/inputs.json")
    working = cfg["edge_work_box"]
    edges = {}
    adjacent = {}
    for edge, value, errors in cfg["edges"]:
        lower = (value - working["widths"] * sum(errors)) / (1 + working["amplitude_transfer"])
        assert lower > 0
        edges[edge] = lower
        a, b = edge
        adjacent.setdefault(a, set()).add(b)
        adjacent.setdefault(b, set()).add(a)
    reached = {"u"}
    for _ in range(5):
        reached |= set().union(*(adjacent[k] for k in reached))
    assert reached == set("uctdsb")
    return edges


def calculate():
    allowed = (-3, -1, 1, 3)
    universal = []
    for charges in itertools.product(allowed, repeat=6):
        if not any(anomalies([[x] * 3 for x in charges])):
            universal.append(charges)
    assert universal == [(-1, 1, 1, 3, -3, -3), (1, -1, -1, -3, 3, 3)]
    partitions = sorted(signed_squares(3))
    assert partitions == [-9, -5, -3, -1, 1, 3, 5, 9]
    minimum_weak = min(abs(3 * q + l) for q in (-9, 9) for l in partitions)
    minimum_colour = min(abs(3 * q + d) for q in (-9, 9) for d in partitions)
    costs = {str(a): min(i for i in range(7) if a in signed_squares(i)) for a in (18,36)}
    assert minimum_weak == minimum_colour == 18 and set(costs.values()) == {6}
    # Distinguish local anomaly freedom from observed quark weak connectivity.
    nonuniversal = [[1, 1, -1], [-1, -1, 1], [-1, -1, 1], [-1] * 3, [1] * 3, [1] * 3]
    assert anomalies(nonuniversal) == [0] * 6
    universal_physical = [[x] * 3 for x in universal[-1]]
    parents = {name: 3 // abs(x) for name, x in zip(SPECIES, universal[-1])}
    assert sum(3 * int(weak == 2) * colour for colour, weak, _ in FIELDS) % 2 == 0
    return {
        "signed_charge_traces_of_arbitrary_three_chiral_modes": partitions,
        "weak_abs_min_Q_still_one_parent": minimum_weak,
        "colour_abs_min_same_Q_U_still_one_parent": minimum_colour,
        "minimum_compensating_gauge_index_cost": costs,
        "family_universal_solutions_in_effective_degree_units": universal,
        "positive_Q_branch_parent_counts": parents,
        "per_three_generation_total_X_anomalies_SU3_SU2_YY_YXX_gravity_cubic": anomalies(universal_physical),
        "nonuniversal_charge_example": nonuniversal,
        "nonuniversal_example_all_anomalies": anomalies(nonuniversal),
        "weak_edges_excluding_nonuniversal_example_unbroken": connected_weak_edges(),
        "example_common_delta_inverse_alpha_1TeV_to_1e16GeV_if_delta_b2_equal_delta_b3_equal_2":
            -math.log(1e13) / math.pi,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = json.loads(json.dumps(calculate()))
    if args.write_new:
        with (HERE / "independent_results.json").open("x", encoding="utf8") as output:
            json.dump(result, output, indent=2)
            output.write("\n")
    else:
        assert json.dumps(result, sort_keys=True) == json.dumps(read("independent_results.json"), sort_keys=True)
    print("Checked index partitions, anomaly-free B-L, and independent connected charged-current data.")
