"""Conditional branching from an occupied return port and finite docking resources."""

import argparse
import json
from dataclasses import dataclass, field
from pathlib import Path

from feedback_routes import BITS, circulate, recurrent_classes
from joint_action import prepare
from window_action import Window


@dataclass
class Module:
    # Prepared but unused: no bound participants and no occupied routes.
    serial: int
    lanes: tuple = field(default_factory=lambda: prepare((1, 0, 0), 1))
    window: Window = field(default_factory=Window)
    entrants: tuple = ()


class Growth:
    """Persistent interfaces and their frontier; a path is an address, not a probe."""

    def __init__(self, cells, spares):
        if cells < 0 or spares < 0:
            raise ValueError("nonnegative resource counts required")
        self.unused = [Module(i) for i in range(cells, 0, -1)]
        self.reserves = list(range(spares + 1, 1, -1))
        self.frontier = {(): 1}  # Existing root circuit and its sole moving probe.
        self.junctions = {}
        self.spent = []

    def accessible(self, path):
        # At successive docked interfaces, use the actual switch to test whether
        # a signal returns to the port that admitted it. Visits are serial.
        tag = self.frontier[path]
        for depth, port in enumerate(path):
            parent = self.junctions[path[:depth]]
            signal = (tag, 0) if port == 0 else (0, tag)
            _, returned = circulate(parent.lanes[0], signal)
            if returned != signal:
                return False
        return True

    def split(self, path, consent=(True, True)):
        if path not in self.frontier:
            raise ValueError("can split only an occupied frontier port")
        if not self.unused or not self.reserves or not self.accessible(path):
            return False

        module = self.unused.pop()
        spare = self.reserves.pop()
        original = (self.frontier[path], spare)
        ports = original
        assert len(recurrent_classes(module.lanes[0], ports, False)) == 1
        for stage in BITS:
            _, ports = circulate(module.lanes[0], ports)
            for side, tag in enumerate(original):
                expected_port = side ^ (1 - stage)
                if consent[side] and ports[expected_port] == tag:
                    item = module.window.handles[side].issue(stage)
                    if not module.window.deliver(side, item):
                        raise RuntimeError("a live contact failed to reach its own socket")
        assert ports == original
        module.lanes = module.window.close(module.lanes)
        if len(recurrent_classes(module.lanes[0], ports, False)) != 2:
            # An expired solicitation is not a fresh resource.
            self.spent.append(module)
            self.reserves.append(spare)
            return False
        for side, tag in enumerate(original):
            alone = (tag, 0) if side == 0 else (0, tag)
            assert circulate(module.lanes[0], alone)[1] == alone
        module.entrants = original  # A record of this event, not extra carriers.
        self.junctions[path] = module
        del self.frontier[path]
        for side, tag in enumerate(original):
            self.frontier[path + (side,)] = tag
        return True

    def lineage(self):
        def visit(path):
            if path in self.frontier:
                return self.frontier[path]
            return (visit(path + (0,)), visit(path + (1,)))

        return visit(())

    def relation_lanes(self):
        return {path: cell.lanes for path, cell in self.junctions.items()}


def grow(cells, spares, order="breadth"):
    selectors = {
        "breadth": lambda path: (len(path), path),
        "depth": lambda path: (-len(path), path),
        "interleave": lambda path: (len(path), path[::-1]),
    }
    if order not in selectors:
        raise ValueError("unknown scheduling rule")
    growth = Growth(cells, spares)
    while growth.unused and growth.reserves:
        selected = min(growth.frontier, key=selectors[order])
        if not growth.split(selected):
            raise RuntimeError("prepared split failed")
    return growth


def check_resources():
    count = 0
    for c in range(9):
        for s in range(9):
            for order in ("breadth", "depth", "interleave"):
                growth = grow(c, s, order)
                used = min(c, s)
                assert len(growth.junctions) == used and len(growth.frontier) == used + 1
                assert len(growth.unused) == c - used and len(growth.reserves) == s - used
                assert not growth.spent
                probes = [*growth.frontier.values(), *growth.reserves]
                assert sorted(probes) == list(range(1, s + 2))
                assert not set(growth.junctions).intersection(growth.frontier)
                for path, module in growth.junctions.items():
                    assert module.lanes == prepare((1, 0, 0), 0)
                    for side, entrant in enumerate(module.entrants):
                        below = path + (side,)
                        if below in growth.frontier:
                            assert growth.frontier[below] == entrant
                        else:
                            assert growth.junctions[below].entrants[0] == entrant
                assert all(growth.accessible(path) for path in growth.frontier)
                count += 1
    return count


def check_controls():
    for consent in ((False, False), (False, True), (True, False)):
        growth = Growth(1, 1)
        assert not growth.split((), consent)
        assert growth.lineage() == 1 and growth.reserves == [2] and not growth.unused
        assert len(growth.spent) == 1 and not growth.junctions

    growth = Growth(3, 3)
    assert growth.split(()) and growth.split((0,))
    saved = growth.lineage()
    root = growth.junctions[()]
    original = root.lanes
    x, y = original[0]
    root.lanes = ((x, 1 - y), *original[1:])
    # The genealogy remains, but its inherited access has been disrupted.
    assert growth.lineage() == saved and not growth.split((0, 1))
    assert len(growth.unused) == len(growth.reserves) == 1
    root.lanes = original
    assert growth.split((0, 1))


def run():
    checks = check_resources()
    check_controls()
    schedules = {}
    for order in ("breadth", "depth", "interleave"):
        growth = grow(5, 5, order)
        sides = [sum(path[0] == side for path in growth.frontier) for side in BITS]
        schedules[order] = {"lineage": growth.lineage(), "leaves_per_root_port": sides,
                            "maximum_depth": max(map(len, growth.frontier))}
    assert schedules["breadth"]["leaves_per_root_port"] == [4, 2]
    assert schedules["interleave"]["leaves_per_root_port"] == [3, 3]
    assert schedules["depth"]["maximum_depth"] == 5
    return {"resource_and_schedule_cases": checks, "nonjoint_controls": 3,
            "inherited_route_intervention": True, "five_splits": schedules}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    verify = parser.parse_args().check
    result = run()
    encoded = json.dumps(result, indent=2) + "\n"
    target = Path(__file__).with_name("recursive_growth_results.json")
    if verify:
        assert json.loads(target.read_text(encoding="utf-8")) == json.loads(encoded)
    else:
        with target.open("x", encoding="utf-8") as output:
            output.write(encoded)
    print(json.dumps(result))
