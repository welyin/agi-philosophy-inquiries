"""One active checker recruits carriers on an anonymous degree-two medium."""

import argparse
import json
from itertools import product
from pathlib import Path


def make_graph(size, closed=True):
    graph = {i: set() for i in range(size)}
    edges = [(i, i + 1) for i in range(size - 1)]
    if closed:
        edges.append((size - 1, 0))
    for u, v in edges:
        graph[u].add(v)
        graph[v].add(u)
    return graph


def probe(traces, belief):
    # Two clocked, read-only accesses. The caller's belief sets the query order.
    first = traces[belief]
    second = traces[1 - belief]
    return first ^ second


def launch(graph, source, traces, belief, free_slots):
    ports = graph[source]
    observation = probe(traces, belief)
    if belief == observation or len(ports) > free_slots:
        return {}, free_slots
    # Every free carrier starts with the old claim; each must check it for itself.
    return {(source, v): belief for v in ports}, free_slots - len(ports)


def tick(graph, source, traces, carriers, glitches):
    moved = {}
    receipts = {}
    for (u, v), claim in carriers.items():
        if (u, v) in glitches:
            claim ^= 1
            glitches.remove((u, v))
        if v == source:
            observation = probe(traces, claim)
            receipts[u] = (claim, observation)
            claim = observation
        for w in graph[v] - {u}:
            assert (v, w) not in moved, "model supports one carrier per directed port"
            moved[v, w] = claim
    return moved, receipts


def travel(graph, source, traces, carriers, ticks, glitches=()):
    reports = []
    pending = set(glitches)
    for _ in range(ticks):
        carriers, arrivals = tick(graph, source, traces, carriers, pending)
        if arrivals:
            reports.append(arrivals)
    assert not pending
    return carriers, reports


def cycles_of_edges(graph):
    remaining = {(u, v) for u in graph for v in graph[u]}
    cycles = []
    while remaining:
        start = next(iter(remaining))
        u, v = start
        orbit = []
        while (u, v) in remaining:
            remaining.remove((u, v))
            orbit.append((u, v))
            w, = graph[v] - {u}
            u, v = v, w
        assert (u, v) == start
        cycles.append(orbit)
    return cycles


def check_ring(graph, source, traces):
    size = len(graph)
    q = traces[0] ^ traces[1]
    assert launch(graph, source, traces, q, 2) == ({}, 2)
    assert launch(graph, source, traces, 1 - q, 1) == ({}, 1)
    launched, reserve = launch(graph, source, traces, 1 - q, 2)
    assert len(launched) == 2 and reserve == 0
    correct, first = travel(graph, source, traces, launched, size)
    assert first == [{p: (1 - q, q) for p in graph[source]}]
    assert correct == {edge: q for edge in launched}

    # Faults are one-use, directed-port interactions, on the first hop of a tour.
    outgoing = list(correct)
    for flags in product((False, True), repeat=2):
        errors = {edge for edge, on in zip(outgoing, flags) if on}
        repaired, reports = travel(graph, source, traces, correct, size, errors)
        assert repaired == correct
        assert len(reports) == 1
        for edge in outgoing:
            return_port, = graph[source] - {edge[1]}
            assert reports[0][return_port] == (q ^ int(edge in errors), q)
        reversed_order = dict(reversed(list(correct.items())))
        assert travel(graph, source, traces, reversed_order, size, errors) == (repaired, reports)

    # Removing a carrier leaves the surviving carrier able to correct a new error.
    for edge in outgoing:
        survivor, reports = travel(graph, source, traces, {edge: q}, size, [edge])
        assert survivor == {edge: q}
        assert list(reports[0].values()) == [(1 - q, q)]

    # An operation on the shared pair changes each carrier's next evidence.
    changed = (traces[0], 1 - traces[1])
    learned, changed_reports = travel(graph, source, changed, correct, size)
    assert learned == {edge: 1 - q for edge in correct}
    assert changed_reports == [{p: (q, 1 - q) for p in graph[source]}]

    # Cut the link between the first two sites; relaunch from the checkpoint on
    # surviving outgoing edges. Dropped carriers cannot traverse the missing link.
    cut = {u: ports.copy() for u, ports in graph.items()}
    cut[0].remove(1)
    cut[1].remove(0)
    surviving = {(u, v): c for (u, v), c in correct.items() if v in cut[u]}
    final, cut_reports = travel(cut, source, traces, surviving, 2 * size)
    assert not final and not cut_reports


def run():
    ring_cases = 0
    topology_cases = 0
    for size in range(3, 13):
        ring = make_graph(size)
        first, second = cycles_of_edges(ring)
        assert len(first) == len(second) == size
        assert set(first) == {(v, u) for u, v in second}
        topology_cases += 1
        for source in ring:
            for traces in product((0, 1), repeat=2):
                check_ring(ring, source, traces)
                ring_cases += 1
        chain = make_graph(size, closed=False)
        source = size // 2
        emitted, _ = launch(chain, source, (0, 1), 0, 2)
        final, reports = travel(chain, source, (0, 1), emitted, 2 * size)
        assert len(emitted) == 2 and not final and not reports
    return {
        "ring_sizes": [3, 12], "distinct_topologies_checked": topology_cases,
        "ring_seed_and_preparation_cases": ring_cases,
        "private_error_patterns_per_case": 4,
        "single_carrier_removals_per_case": 2,
        "shared_trace_changes_per_case": 1,
        "single_link_cuts_per_case": 1,
        "open_chains_with_emission_but_no_return": topology_cases,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    check = parser.parse_args().check
    results = run()
    path = Path(__file__).with_name("return_paths_results.json")
    if check:
        assert results == json.loads(path.read_text(encoding="utf-8"))
    else:
        with path.open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results))
