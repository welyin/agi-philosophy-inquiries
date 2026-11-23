"""Check what a newly retained interaction can contribute to a feedback loop."""

import argparse
import json
from itertools import combinations, permutations, product
from pathlib import Path


def reachable(edges, start):
    seen = {start}
    pending = [start]
    while pending:
        current = pending.pop()
        for u, v in edges:
            if u == current and v not in seen:
                seen.add(v)
                pending.append(v)
    return seen


def simple_cycles(n):
    cycles = []
    for length in range(2, n + 1):
        for vertices in combinations(range(n), length):
            # Fix the least vertex at the start to count each directed cycle once.
            for tail in permutations(vertices[1:]):
                path = (vertices[0],) + tail
                cycles.append(frozenset(zip(path, path[1:] + path[:1])))
    return cycles


def check_new_edge(n):
    ports = tuple(permutations(range(n), 2))
    cycles = simple_cycles(n)
    containing = {edge: [c for c in cycles if edge in c] for edge in ports}
    closes = 0
    total = 0
    for configuration in product((False, True), repeat=len(ports)):
        edges = {edge for edge, active in zip(ports, configuration) if active}
        returns = [reachable(edges, v) for v in range(n)]
        for u, v in set(ports) - edges:
            enlarged = edges | {(u, v)}
            new_cycle = any(c <= enlarged for c in containing[u, v])
            assert new_cycle == (u in returns[v])
            closes += int(new_cycle)
            total += 1
    return {"additions": total, "with_new_cycle": closes, "without_new_cycle": total - closes}


def encounter(action, observed, retain):
    # A local writable return port latches the relative response when enabled.
    return action ^ observed if retain else None


def use_return_port(relation, claim, feed_back, steps):
    # None means no message. A delivered zero is an actual response.
    messages = []
    for _ in range(steps):
        if relation is None:
            break
        reply = claim ^ relation
        messages.append((claim, reply))
        if not feed_back:
            break
        claim ^= reply
    return messages


def check_encounters():
    cases = 0
    for action, observed, original_claim in product((0, 1), repeat=3):
        relation = encounter(action, observed, True)
        assert relation is not None
        for guess in (original_claim, 1 - original_claim):
            messages = use_return_port(relation, guess, True, 4)
            assert messages == [(guess, guess ^ relation)] + [(relation, 0)] * 3
            # Same preparation without retention, or actual erasure afterwards.
            assert not use_return_port(encounter(action, observed, False), guess, True, 4)
            assert not use_return_port(None, guess, True, 4)
            # The trace still answers once when no answer can trigger a next probe.
            assert use_return_port(relation, guess, False, 4) == messages[:1]
        cases += 1
    old = {(0, 1), (1, 2), (2, 0)}
    new = old | {(2, 1)}
    possible = simple_cycles(3)
    assert sum(c <= old for c in possible) == 1
    assert sum(c <= new for c in possible) == 2
    without_initiator = {edge for edge in new if 0 not in edge}
    assert sum(c <= without_initiator for c in possible) == 1
    return cases


def run():
    return {"encounter_cases": check_encounters(),
            "three_vertices": check_new_edge(3), "four_vertices": check_new_edge(4)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    check = parser.parse_args().check
    result = run()
    path = Path(__file__).with_name("retained_relation_results.json")
    if check:
        assert result == json.loads(path.read_text(encoding="utf-8"))
    else:
        with path.open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
