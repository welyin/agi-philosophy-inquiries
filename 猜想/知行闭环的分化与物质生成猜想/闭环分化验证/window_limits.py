"""Independent checks of action order and distinguishable unfinished commitments."""

import argparse
import json
from itertools import combinations, permutations, product
from pathlib import Path

from collective_maintenance import World, maintain, prepared
from interaction_order import evolve
from joint_action import COUPLER, joint_window, prepare


def apply_in_arrival_order(start, actors):
    lanes = start
    for actor in actors:
        lanes = tuple(evolve(lane, "B" if owner == actor else "A", (0, 1))
                      for owner, lane in zip(COUPLER, lanes))
    return lanes


def check_orders():
    configurations = 0
    groups = {"flip": [], "keep": []}
    for actors in sorted(set(permutations("LLRR"))):
        pairs = sum(a == "L" and b == "R" for a, b in combinations(actors, 2))
        flip = pairs % 2
        for settings, old in product(product((0, 1), repeat=3), (0, 1)):
            start = prepare(settings, old)
            end = apply_in_arrival_order(start, actors)
            assert end == prepare(settings, old ^ flip)
            if actors == tuple("LRLR"):
                assert end == joint_window(start, (1, 1, 1, 1))
            configurations += 1
        groups["flip" if flip else "keep"].append("".join(actors))
    assert groups["flip"] == ["LRLR", "RLRL"]
    return configurations, groups


def check_receipt_capacity():
    # Four distinct single-use slots; an arrival may be duplicated or occur in
    # any order. Evaluate equivalence by testing every possible continuation.
    slots = frozenset(("L1", "L2", "R1", "R2"))
    all_sets = [frozenset(combo) for size in range(5)
                for combo in combinations(sorted(slots), size)]
    signatures = {tuple((seen | suffix) == slots for suffix in all_sets) for seen in all_sets}
    assert len(signatures) == 16
    pairs_checked = 0
    for first, second in combinations(all_sets, 2):
        if not first - second:
            first, second = second, first
        suffix = slots - first
        assert first | suffix == slots and second | suffix != slots
        pairs_checked += 1
    return pairs_checked


def check_private_votes():
    # Identical physical states follow (1) two deliberate private edits and
    # (2) a left edit with a right memory fault. Only history 1 requests change.
    def private_flip(world, who):
        ports = list(world.ports)
        tag, vote = ports[who]
        ports[who] = (tag, 1 - vote)
        return World(world.marks, tuple(ports))

    for old in (0, 1):
        initial = prepared(0, old)
        pending_left = private_flip(initial, 0)
        consented = private_flip(pending_left, 1)
        faulted = World(pending_left.marks, (pending_left.ports[0], (1, 1 - old)))
        assert consented == faulted
        for wire in (0, 1):
            repaired, _ = maintain(faulted, (0, 0), wire)
            assert repaired.marks == (0, 1 - old)


def run():
    count, orders = check_orders()
    check_private_votes()
    return {"order_and_preparation_cases": count, "direct_action_orders": orders,
            "distinguishable_pending_sets": 16,
            "separated_pairs_of_pending_sets": check_receipt_capacity(),
            "one_fault_plus_one_private_edit_can_imitate_joint_vote": True}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    verify = parser.parse_args().check
    result = run()
    path = Path(__file__).with_name("window_limits_results.json")
    if verify:
        assert json.loads(path.read_text(encoding="utf-8")) == result
    else:
        with path.open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
