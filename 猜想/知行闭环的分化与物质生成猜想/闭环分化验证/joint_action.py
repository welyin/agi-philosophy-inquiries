"""Paired local actions change a maintained relation when both take part."""

import argparse
import json
from itertools import product
from pathlib import Path

from collective_maintenance import World, maintain, read_relation, tags
from feedback_routes import BITS, PAIRS, circulate


# A lane contains (setting, target). Only the target's owner can couple it
# to its setting. The other participant can flip the setting at its own turn.
# Lane order: common object, left private claim, right private claim.
COUPLER = ("R", "L", "R")
SLOTS = ("L", "R", "L", "R")


def joint_window(lanes, pulses, fault=None):
    lanes = [list(lane) for lane in lanes]
    for time in range(len(SLOTS) + 1):
        if fault is not None and time == fault[0]:
            _, lane, field = fault
            lanes[lane][field] ^= 1
        if time == len(SLOTS) or not pulses[time]:
            continue
        actor = SLOTS[time]
        for lane, owner in zip(lanes, COUPLER):
            if actor == owner:
                lane[1] ^= lane[0]
            else:
                lane[0] ^= 1
    return tuple(tuple(lane) for lane in lanes)


def prepare(settings, old):
    x, left, right = settings
    return ((x, x ^ old), (left, old), (right, old))


def as_world(lanes, first):
    # Claims are stored by the two participants and travel with their tags.
    by_owner = tuple((owner, lanes[owner + 1][1]) for owner in BITS)
    return World(lanes[0], (by_owner[first], by_owner[1 - first]))


def verify_transactions():
    faults = tuple(product(range(5), range(3), BITS))
    clean = disturbed = 0
    for settings, old, (left, right) in product(product(BITS, repeat=3), BITS, PAIRS):
        start = prepare(settings, old)
        pulses = (left, right, left, right)
        expected = old ^ (left & right)  # Used only in assertions.
        clean_end = joint_window(start, pulses)
        assert clean_end == prepare(settings, expected)
        for fault, wire, first in product((None,) + faults, BITS, BITS):
            actual = as_world(joint_window(start, pulses, fault), first)
            votes = (read_relation(actual), *(claim for _, claim in actual.ports))
            assert sum(v != expected for v in votes) <= 1
            # No action log, desired target or fault metadata enters maintenance.
            restored, _ = maintain(actual, (0, 0), wire)
            assert read_relation(restored) == expected
            assert all(claim == expected for _, claim in restored.ports)
            port_order = (first, 1 - first)
            assert tags(restored) == (port_order[::-1] if expected ^ wire else port_order)
            if fault is None:
                clean += 1
            else:
                disturbed += 1
    return clean, disturbed


def check_limits():
    setup = prepare((0, 1, 0), 0)
    # Both actions are truly asserted in the first run. In the second, R's
    # source is stuck on before it fans out. Downstream evidence is identical.
    intended = (1, 0, 1, 0)
    stuck_right = tuple(1 if side == "R" else on for side, on in zip(SLOTS, intended))
    real_joint = joint_window(setup, (1, 1, 1, 1))
    assert real_joint == joint_window(setup, stuck_right)
    assert real_joint != joint_window(setup, intended)
    # A single unmatched L pulse changes the object's setting and its own
    # private claim. The protected single-lane fault contract no longer holds.
    unmatched = as_world(joint_window(setup, (1, 0, 0, 0)), 0)
    observed = [read_relation(unmatched), *(b for _, b in unmatched.ports)]
    assert observed == [1, 1, 0]
    repaired, _ = maintain(unmatched, (0, 0), 0)
    assert read_relation(repaired) == 1
    # Following a joint change, the old straight return crosses at fixed wire.
    # A real change to the wire can preserve the return while readout changes.
    probes = (1, 2)
    assert circulate(setup[0], probes, 0)[1] == probes
    assert circulate(real_joint[0], probes, 0)[1] == probes[::-1]
    assert circulate(real_joint[0], probes, 1)[1] == probes
    assert read_relation(as_world(setup, 0)) != read_relation(as_world(real_joint, 0))
    return {"shared_source_fault_matches_joint_action": True,
            "unmatched_single_party_votes": observed,
            "physical_rewiring_can_mask_change_in_return_ports": True}


def run():
    clean, disturbed = verify_transactions()
    return {"clean_transactions": clean, "single_lane_fault_transactions": disturbed,
            "injection_boundaries": 5, "limits": check_limits()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    verify = parser.parse_args().check
    results = run()
    path = Path(__file__).with_name("joint_action_results.json")
    if verify:
        assert json.loads(path.read_text(encoding="utf-8")) == results
    else:
        with path.open("x", encoding="utf-8") as output:
            output.write(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results))
