"""Measure repeatable responses of nested instances of the four-contact protocol."""

import argparse
import json
from functools import lru_cache
from itertools import product
from math import comb
from pathlib import Path

from feedback_routes import BITS, circulate
from interaction_order import evolve
from joint_action import SLOTS, joint_window, prepare
from recursive_growth import grow


def read_contact(lanes):
    marks = lanes[0]
    revealed = evolve(marks, "B", BITS)
    pulse = revealed[1]
    assert evolve(revealed, "B", BITS) == marks
    return pulse


class Hierarchy:
    def __init__(self, tree, leaves, formed=None):
        self.lanes = {}
        self.children = {}
        self.leaves = {}
        self.wires = {}
        self.spent = set()

        def connect(fragment, path):
            if isinstance(fragment, int):
                self.leaves[path] = fragment
                self.lanes[path] = prepare((1, 0, 0), leaves[fragment])
                return
            left, right = fragment
            children = (path + (0,), path + (1,))
            self.children[path] = children
            self.lanes[path] = formed[path] if formed is not None else prepare((1, 0, 0), 0)
            self.wires[path] = (0, read_contact(self.lanes[path]))
            connect(left, children[0])
            connect(right, children[1])

        connect(tree, ())
        if formed is not None:
            assert set(formed) == set(self.children)
        self.outward = tuple(sorted(self.children, key=len))

    def actuate(self, path, *, once=False, coadjust=False, disconnected=()):
        if once and path in self.spent:
            return
        children = self.children[path]
        # Each of the four contacts samples the corresponding child's public
        # readout. Values remain stable throughout this one local window.
        pulses = []
        for role in SLOTS:
            child = children[int(role == "R")]
            pulse = 0 if child in disconnected else read_contact(self.lanes[child])
            pulses.append(pulse)
        before = read_contact(self.lanes[path])
        self.lanes[path] = joint_window(self.lanes[path], pulses)
        after = read_contact(self.lanes[path])
        assert all(vote == after for _, vote in self.lanes[path][1:])
        if once and after != before:
            self.spent.add(path)
        if coadjust:
            # Additional wiring actuator, driven by the actual same pulses.
            wire = self.wires[path]
            for role, pulse in zip(SLOTS, pulses):
                if pulse:
                    wire = evolve(wire, "A" if role == "L" else "B", BITS)
            self.wires[path] = wire

    def sweep(self, *, inward=False, once=False, coadjust=False):
        order = reversed(self.outward) if inward else self.outward
        for path in order:
            self.actuate(path, once=once, coadjust=coadjust)

    def value(self, path=()):
        return read_contact(self.lanes[path])

    def still_separate(self, path):
        _, receivers = circulate(self.lanes[path][0], (1, 2), self.wires[path][1])
        return receivers == (1, 2)


@lru_cache(maxsize=None)
def all_trees(first, count):
    if count == 1:
        return (first,)
    trees = []
    for n_left in range(1, count):
        for left in all_trees(first, n_left):
            for right in all_trees(first + n_left, count - n_left):
                trees.append((left, right))
    return tuple(trees)


def indicators(model):
    # Analytical comparator, used only after the initial preparations.
    stats = {p: (0, model.value(p), 1, 0) for p in model.leaves}
    for p in reversed(model.outward):
        left, right = (stats[c] for c in model.children[p])
        k = 1 + (left[0] | right[0])
        leaf_count = left[2] + right[2]
        depth = 1 + max(left[3], right[3])
        assert depth <= k <= leaf_count - 1
        stats[p] = (k, left[1] * right[1], leaf_count, depth)
    return stats


def check_one_use():
    cases = 0
    for n in range(1, 7):
        for tree in all_trees(0, n):
            for bits in product(BITS, repeat=n):
                model = Hierarchy(tree, dict(enumerate(bits)))
                model.sweep(inward=True)
                assert model.value() == int(all(bits))
                cases += 1
    return cases


def least_period(trace):
    for period in range(1, len(trace) // 2 + 1):
        if all(value == trace[i + period] for i, value in enumerate(trace[:-period])):
            return period
    raise AssertionError("too few samples for period check")


def check_time_law():
    cases = 0
    for n in range(1, 7):
        for tree in all_trees(0, n):
            for zero in (None, n // 2):
                leaves = {i: int(i != zero) for i in range(n)}
                model = Hierarchy(tree, leaves)
                expected = indicators(model)
                root_trace = []
                for t in range(17):
                    for path, (k, active, _, _) in expected.items():
                        assert model.value(path) == (comb(t, k) % 2) * active
                    root_trace.append(model.value())
                    if t < 16:
                        model.sweep()
                k, active, _, _ = expected[()]
                period = 1 << k.bit_length() if active else 1
                assert least_period(root_trace) == period
                cases += 1
    return cases


def grown_example(order):
    grown = grow(5, 5, order=order)
    tree = grown.lineage()
    formed = grown.relation_lanes()

    def names(node):
        if isinstance(node, int):
            return [node]
        return names(node[0]) + names(node[1])

    labels = names(tree)
    leaves = dict.fromkeys(labels, 1)
    repeat = Hierarchy(tree, leaves, formed)
    single = Hierarchy(tree, leaves, formed)
    matched = Hierarchy(tree, leaves, formed)
    assert all(repeat.value(p) == 0 for p in repeat.children)
    trajectory, one_event = [], []
    for t in range(17):
        trajectory.append(repeat.value())
        one_event.append(single.value())
        for path in repeat.children:
            assert repeat.value(path) == matched.value(path)
            assert repeat.still_separate(path) == (repeat.value(path) == 0)
            assert matched.still_separate(path)
        if t < 16:
            repeat.sweep()
            single.sweep(once=True)
            matched.sweep(coadjust=True)
    k, _, count, depth = indicators(Hierarchy(tree, leaves, formed))[()]
    assert trajectory.index(1) == k and least_period(trajectory) == 1 << k.bit_length()
    assert one_event == [int(t >= depth) for t in range(len(one_event))]
    return {"leaves_per_side": [len(names(child)) for child in tree],
            "leaf_count": count, "depth": depth, "first_response": k,
            "root_period": least_period(trajectory), "first_nine_readings": trajectory[:9],
            "single_event_stable_from": depth}


def check_physical_link():
    model = Hierarchy(((0, 1), (2, 3)), dict.fromkeys(range(4), 1))
    model.sweep()
    assert [model.value(path) for path in ((), (0,), (1,))] == [0, 1, 1]
    model.actuate((), disconnected=((0,),))
    assert model.value() == 0
    model.actuate(())
    assert model.value() == 1


def run():
    check_physical_link()
    results = {order: grown_example(order) for order in ("breadth", "interleave")}
    assert [results[o]["root_period"] for o in ("breadth", "interleave")] == [8, 4]
    return {"one_use_tree_and_input_cases": check_one_use(),
            "time_law_tree_and_input_cases": check_time_law(),
            "grown_six_leaf_examples": results}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    verify = parser.parse_args().check
    result = run()
    target = Path(__file__).with_name("recursive_response_results.json")
    if verify:
        assert json.loads(target.read_text(encoding="utf-8")) == result
    else:
        with target.open("x", encoding="utf-8") as output:
            output.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
