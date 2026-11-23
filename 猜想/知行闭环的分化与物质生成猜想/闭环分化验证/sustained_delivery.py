"""Finite, directly wired two-stage requests, cleared only after joint execution."""

import argparse
import json
from dataclasses import dataclass
from itertools import product
from pathlib import Path

from joint_action import joint_window, prepare


@dataclass(frozen=True)
class Port:
    prep: int
    seal: int
    receipt: int

    def act(self):
        if self.prep == self.receipt:
            return Port(1 - self.prep, self.seal, self.receipt)
        if self.seal == self.receipt:
            return Port(self.prep, 1 - self.seal, self.receipt)
        return self  # Continue driving the outstanding requests.

    def acknowledge(self):
        assert self.prep == self.seal != self.receipt
        return Port(self.prep, self.seal, self.seal)

    def pending(self):
        return (self.prep != self.receipt) + (self.seal != self.receipt)


def visit(ports, lanes, side, offer=True):
    # One local read/act visit. Fixed physical inputs do not queue, drop, or
    # copy old pulses. The last new seal triggers one serialized receiver.
    updated = list(ports)
    if offer:
        updated[side] = updated[side].act()
    left, right = updated
    pulses = (left.prep != left.receipt, right.prep != right.receipt,
              left.seal != left.receipt, right.seal != right.receipt)
    if not all(pulses):
        return tuple(updated), lanes, False
    lanes = joint_window(lanes, pulses)
    return tuple(port.acknowledge() for port in updated), lanes, True


def port_at(phase, depth):
    return Port(phase ^ int(depth >= 1), phase ^ int(depth == 2), phase)


def check_reachability():
    start = (port_at(0, 0), port_at(0, 0))
    lanes = prepare((1, 0, 1), 0)
    seen = {start}
    todo = [start]
    while todo:
        ports = todo.pop()
        for side in range(2):
            reached, _, _ = visit(ports, lanes, side)
            if reached not in seen:
                seen.add(reached)
                todo.append(reached)
    expected = {(port_at(phase, ld), port_at(phase, rd))
                for phase, ld, rd in product(range(2), range(3), range(3))
                if (ld, rd) != (2, 2)}
    assert seen == expected
    return len(seen)


def check_transitions():
    # All stable control states: two receipt phases and eight partial profiles.
    # The full/full profile is transient inside visit's receiver transaction.
    count = 0
    for phase, ld, rd in product(range(2), range(3), range(3)):
        if (ld, rd) == (2, 2):
            continue
        ports = (port_at(phase, ld), port_at(phase, rd))
        for settings, old, side in product(product(range(2), repeat=3), range(2), range(2)):
            lanes = prepare(settings, old)
            ended, result, changed = visit(ports, lanes, side)
            progress = int(ports[side].pending() < 2)
            assert changed == (ld + rd + progress == 4)
            assert result == prepare(settings, old ^ changed)
            for i, (before, after) in enumerate(zip(ports, ended)):
                added = progress if i == side else 0
                # Each consumed batch costs exactly two fresh contacts per side.
                assert after.pending() == before.pending() + added - 2 * changed
                assert after.receipt == before.receipt ^ changed
            count += 1
    return count


def check_absence():
    lanes = prepare((1, 0, 1), 0)
    checks = 0
    for active, other_visits in product(range(2), (False, True)):
        ports = (port_at(0, 0), port_at(0, 0))
        for _ in range(2):
            ports, same, changed = visit(ports, lanes, active)
            assert not changed and same == lanes
        # A participant absent or actively declining every new offer leaves
        # the same fixed point under any repetition of the following visits.
        stopped = (ports, lanes, False)
        assert visit(ports, lanes, active) == stopped
        if other_visits:
            assert visit(ports, lanes, 1 - active, offer=False) == stopped
        checks += 1
    return checks


def compute():
    return {"reachable_handshake_control_states": check_reachability(),
            "control_and_lane_transitions": check_transitions(),
            "one_sided_fixed_points_checked": check_absence()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = compute()
    target = Path(__file__).with_name("sustained_delivery_results.json")
    if args.check:
        assert json.loads(target.read_text(encoding="utf-8")) == result
    else:
        with target.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
