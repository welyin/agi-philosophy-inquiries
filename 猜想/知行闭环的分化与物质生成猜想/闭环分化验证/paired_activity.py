"""Local probe feedback and the observable effect of asynchronous phase."""

import argparse
import json
from itertools import permutations, product
from pathlib import Path

from feedback_routes import circulate
from recursive_response import Hierarchy


NODES = {"A": (0,), "B": (1,), "P": ()}


class PairedActivity:
    def __init__(self):
        self.device = Hierarchy(((0, 1), (2, 3)), dict.fromkeys(range(4), 1))
        # A local participant sends its next probe from its last return port.
        self.port = dict.fromkeys(self.device.lanes, 0)

    def probe(self, path):
        marks = self.device.lanes[path][0]
        send = (1, 0) if self.port[path] == 0 else (0, 1)
        restored, received = circulate(marks, send)
        assert restored == marks and sum(received) == 1
        self.port[path] = received.index(1)

    def visit(self, node):
        path = NODES[node]
        self.device.actuate(path)
        self.probe(path)
        return self.device.value(path)


def run_blocks(orders):
    activity = PairedActivity()
    root_values = []
    for order in orders:
        assert sorted(order) == ["A", "B", "P"]
        for leaf in activity.device.leaves:
            # The boundary output stays high as its own probe alternates ports.
            previous = activity.port[leaf]
            activity.probe(leaf)
            assert activity.port[leaf] == 1 - previous
        for node in order:
            value = activity.visit(node)
            if node == "P":
                root_values.append(value)
    return root_values


def all_target_changes(length):
    # Construct interleavings for the mathematical scope test only.
    for target in product(range(2), repeat=length):
        orders = []
        expected = []
        root = 0
        for index, change in enumerate(target):
            order = "APB"
            if change:
                order = "ABP" if index % 2 == 0 else "PAB"
            orders.append(order)
            root ^= change
            expected.append(root)
        assert run_blocks(orders) == expected
    return 2 ** length


def run():
    traces = {}
    for perm in permutations("ABP"):
        order = "".join(perm)
        measured = run_blocks([order] * 8)
        root_slot = order.index("P")
        if root_slot == 1:
            expected = [0] * 8
        elif root_slot == 2:
            expected = [1, 1, 0, 0] * 2
        else:
            expected = [0, 1, 1, 0] * 2
        assert measured == expected
        traces[order] = measured
    # The two lower clocks use phases 0, 2 of the same four-position cycle.
    # Moving just the root detector from phase 3 to 1 swaps the two histories.
    phase_orders = {}
    for root_phase in (1, 3):
        phases = {"A": 0, "B": 2, "P": root_phase}
        order = "".join(sorted(phases, key=phases.get))
        phase_orders[root_phase] = order
    assert phase_orders == {1: "APB", 3: "ABP"}
    return {
        "root_values_after_own_actions": traces,
        "visits_per_participant_in_each_trace": 8,
        "active_boundaries_and_internal_loops": [4, 3],
        "distinct_eight_step_root_change_words_verified": all_target_changes(8),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    options = parser.parse_args()
    results = run()
    path = Path(__file__).with_name("paired_activity_results.json")
    if options.check:
        assert json.loads(path.read_text(encoding="utf-8")) == results
    else:
        with path.open("x", encoding="utf-8") as output:
            output.write(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results))
