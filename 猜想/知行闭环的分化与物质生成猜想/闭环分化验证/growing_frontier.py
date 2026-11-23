"""Local read budgets and delayed joint reports on a growing classical tree."""

import argparse
import heapq
import json
from collections import Counter
from itertools import product
from pathlib import Path

from joint_action import joint_window, prepare
from recursive_growth import Growth
from recursive_response import read_contact


def build_frontiers(depth):
    resources = 2 ** depth - 1
    tree = Growth(resources, resources)
    generations = []
    for d in range(depth + 1):
        leaves = tuple(sorted(tree.frontier))
        assert len(leaves) == 2 ** d and all(len(p) == d for p in leaves)
        generations.append(leaves)
        if d < depth:
            for p in leaves:
                assert tree.split(p)
    assert len(tree.junctions) == resources
    return generations


def all_ready(left, right):
    # A new report accumulator starts at zero. The persistent growth relation
    # is on a different lane. The operands are the actual incoming reports.
    lanes = joint_window(prepare((1, 0, 0), 0), (left, right, left, right))
    return read_contact(lanes)


def reports(preparations, birth_interval, hop_interval):
    """One fixed-capacity lane per edge. Times here are model clock units."""
    events = []
    for generation, bits in preparations.items():
        for path, value in bits.items():
            assert len(path) == generation
            heapq.heappush(events, (generation * birth_interval, generation, path, -1, value))
    waiting = {}
    last_send = {}
    results = {}
    messages = Counter()
    while events:
        time, generation, node, side, value = heapq.heappop(events)
        if side >= 0:
            # Pair by physical receiver only. Generation labels assert that the
            # fixed-delay links preserve order; receivers never select by label.
            if node not in waiting:
                waiting[node] = (side, value, generation)
                continue
            other_side, other_value, other_generation = waiting.pop(node)
            assert other_side != side and other_generation == generation
            left, right = (value, other_value) if side == 0 else (other_value, value)
            value = all_ready(left, right)
        if not node:
            results[generation] = (time, value)
            continue
        if node in last_send:
            assert time - last_send[node] >= hop_interval
        last_send[node] = time
        messages[generation] += 1
        heapq.heappush(events, (time + hop_interval, generation, node[:-1], node[-1], value))
    assert not waiting
    return results, messages


def check_snapshots(generations):
    count = 0
    for d, leaves in enumerate(generations[:4]):
        for values in product(range(2), repeat=len(leaves)):
            preparation = dict(zip(leaves, values))
            output, messages = reports({d: preparation}, 1, 1)
            assert output == {d: (2 * d, int(all(values)))}
            assert messages[d] == 2 * (len(leaves) - 1)
            count += 1
    return count


def check_causal_delay(generations):
    depth = len(generations) - 1
    leaves = generations[depth]
    baseline = {d: dict.fromkeys(paths, 1) for d, paths in enumerate(generations)}
    controls = 0
    examples = []
    for birth, hop in ((1, 1), (1, 3), (3, 1)):
        expected = {d: (d * (birth + hop), 1) for d in baseline}
        original, counts = reports(baseline, birth, hop)
        assert original == expected
        assert all(counts[d] == 2 * (2 ** d - 1) for d in baseline)
        for p in leaves:
            changed = dict(baseline)
            changed[depth] = {**baseline[depth], p: 0}
            new, _ = reports(changed, birth, hop)
            assert all(new[d] == original[d] for d in range(depth))
            assert new[depth] == (expected[depth][0], 0)
            controls += 1
        examples.append({"birth_interval": birth, "hop_interval": hop,
                         "frontier_depth": depth, "birth_time": depth * birth,
                         "root_report_time": original[depth][0]})
    return controls, examples


def check_information(leaves):
    n = len(leaves)
    words = list(product(range(2), repeat=n))
    # With no later access, all coordinate questions distinguish all words.
    counts = {
        "all_coordinate_questions": len(set(words)),
        "all_ready_question": len({int(all(word)) for word in words}),
        "population_count_question": len({sum(word) for word in words}),
        "one_seed_and_known_address_rule": len({
            tuple(seed ^ (sum(p) % 2) for p in leaves) for seed in range(2)}),
    }
    assert list(counts.values()) == [2 ** n, 2, n + 1, 2]
    # On all-ones, dropping any inspected coordinate leaves an adversarial
    # preparation with the same other values and a different readiness report.
    for i in range(n):
        almost = [1] * n
        almost[i] = 0
        assert all(almost[:i] + almost[i + 1:]) and not all(almost)
    return counts


def run():
    generations = build_frontiers(6)
    cases = check_snapshots(generations)
    controls, examples = check_causal_delay(generations)
    return {"exhaustive_snapshot_cases": cases, "one_leaf_interventions": controls,
            "distinguishable_answers_for_eight_leaves": check_information(generations[3]),
            "six_generation_clocks": examples}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    target = Path(__file__).with_name("growing_frontier_results.json")
    if args.check:
        assert json.loads(target.read_text(encoding="utf-8")) == result
    else:
        with target.open("x", encoding="utf-8") as output:
            output.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
