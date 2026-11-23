"""An internal comparison writes a relation that changes occupied return routes."""

import argparse
import json
from dataclasses import dataclass, replace
from itertools import permutations, product
from pathlib import Path

from interaction_order import evolve


BITS = (0, 1)
PAIRS = tuple(product(BITS, repeat=2))
WRITER = ("compare", "first", "second", "A", "compare")
EMPTY = (0, 0)  # Zero is an empty carrier slot; 1 and 2 are identifiable probes.


@dataclass(frozen=True)
class State:
    marks: tuple
    scratch: int = 0
    ports: tuple = EMPTY
    reserve: tuple = EMPTY


def circulate(marks, ports, crossed=False):
    # Reuse B as a temporary readout; its current y controls a reversible switch.
    revealed = evolve(marks, "B", BITS)
    outgoing = ports[::-1] if revealed[1] else ports
    restored = evolve(revealed, "B", BITS)
    return restored, outgoing[::-1] if crossed else outgoing


def transition(state, operation, traces, prediction, crossed):
    if operation == "compare":
        first, second = traces[prediction], traces[1 - prediction]
        return replace(state, scratch=state.scratch ^ prediction ^ first ^ second)
    if operation == "first":
        operation = "A" if state.scratch else "B"
    elif operation == "second":
        operation = "B" if state.scratch else "A"
    if operation in ("A", "B"):
        return replace(state, marks=evolve(state.marks, operation, BITS))
    if operation == "transfer":
        return replace(state, ports=state.reserve, reserve=state.ports)
    if operation == "route":
        marks, ports = circulate(state.marks, state.ports, crossed)
        return replace(state, marks=marks, ports=ports)
    raise ValueError(operation)


def write_relation(state, traces, prediction):
    for operation in WRITER:
        state = transition(state, operation, traces, prediction, False)
    return state


def clock_program(laps):
    if laps < 2 or laps % 2:
        raise ValueError("the common return window requires a positive even lap count")
    return WRITER + ("transfer",) + ("route",) * laps + ("transfer",) + WRITER[::-1]


def check_mediation():
    cases = 0
    for x, traces, prediction, crossed, start in product(BITS, PAIRS, BITS, BITS, BITS):
        originals = (1, 0) if start == 0 else (0, 1)
        state = write_relation(State((x, 0), reserve=originals), traces, prediction)
        error = prediction ^ traces[0] ^ traces[1]
        assert state.marks == (x, x ^ error) and state.scratch == 0
        state = transition(state, "transfer", traces, prediction, crossed)
        for lap in range(1, 7):
            # The post-handoff transition must not consult the original records.
            alternatives = {
                transition(state, "route", other_traces, other_claim, crossed)
                for other_traces, other_claim in product(PAIRS, BITS)
            }
            assert len(alternatives) == 1
            state = alternatives.pop()
            port = start ^ ((lap % 2) * (error ^ crossed))
            assert state.ports.index(1) == port and state.reserve == EMPTY
        marks = state.marks
        regular = circulate(marks, state.ports, crossed)[1]
        assert circulate((marks[0], 1 - marks[1]), state.ports, crossed)[1] == regular[::-1]
        assert circulate((1 - marks[0], 1 - marks[1]), state.ports, crossed)[1] == regular
        # Disconnect both outgoing ports after one passage. The probe is parked,
        # marks are untouched, and three further route ticks receive no probe.
        once = transition(state, "route", traces, prediction, crossed)
        parked = transition(once, "transfer", traces, prediction, crossed)
        for _ in range(3):
            parked = transition(parked, "route", traces, prediction, crossed)
        assert parked.ports == EMPTY and sum(parked.reserve) == 1 and parked.marks == marks
        initial = State((x, 0), reserve=originals)
        clocked = initial
        for operation in clock_program(6):
            clocked = transition(clocked, operation, traces, prediction, crossed)
        assert clocked == initial
        cases += 1
    dormant = State((0, 0), ports=(1, 0))
    assert circulate(dormant.marks, dormant.ports)[1] == (1, 0)
    awakened = write_relation(dormant, (0, 1), 0)
    assert circulate(awakened.marks, awakened.ports)[1] == (0, 1)
    return cases


def recurrent_classes(marks, occupied, crossed):
    # Follow a single tag to infer the actual orbit of each occupied port.
    classes = set()
    for start, marker in enumerate(occupied):
        if not marker:
            continue
        ports = (1, 0) if start == 0 else (0, 1)
        seen = set()
        while ports.index(1) not in seen:
            seen.add(ports.index(1))
            _, ports = circulate(marks, ports, crossed)
        classes.add(frozenset(seen))
    return classes


def stop_port_one(marks, ports):
    parked = []
    visits_at_zero = []
    for _ in range(4):
        if ports[1]:
            parked.append(ports[1])
            ports = (ports[0], 0)
        visits_at_zero.append(ports[0])
        _, ports = circulate(marks, ports)
    return visits_at_zero, sorted(parked)


def check_occupation():
    cases = 0
    for marks, occupied, crossed in product(PAIRS, (EMPTY, (1, 0), (0, 1), (1, 2)), BITS):
        count = sum(bool(tag) for tag in occupied)
        classes = recurrent_classes(marks, occupied, crossed)
        changes_port = marks[0] ^ marks[1] ^ crossed
        assert len(classes) == (int(count > 0) if changes_port else count)
        cases += 1
    # Same two existing probes. An initially crossed return pattern separates.
    circulating = State((1, 0), ports=(1, 2))
    for expected in ((2, 1), (1, 2)):
        circulating = transition(circulating, "route", (1, 1), 0, False)
        assert circulating.ports == expected
    before = circulating.marks
    separated = write_relation(circulating, (1, 1), 0)
    assert separated.ports == circulating.ports
    after = separated.marks
    assert recurrent_classes(before, (1, 2), False) == {frozenset({0, 1})}
    assert recurrent_classes(after, (1, 2), False) == {frozenset({0}), frozenset({1})}
    assert stop_port_one(before, (1, 2)) == ([1, 0, 0, 0], [1, 2])
    assert stop_port_one(after, (1, 2)) == ([1, 1, 1, 1], [2])
    return cases


def check_reversible_clock():
    # Every instruction acts on the actual state and is its own inverse.
    # A cyclic phase counter selects it; the inverse uses the preceding phase.
    placements = set(permutations((0, 0, 1, 2)))
    cases = 0
    for traces, prediction, crossed in product(PAIRS, BITS, BITS):
        for marks, scratch, slots in product(PAIRS, BITS, placements):
            initial = State(marks, scratch, slots[:2], slots[2:])
            state = initial
            program = clock_program(6)
            phase = 0
            while True:
                operation = program[phase]
                nxt = transition(state, operation, traces, prediction, crossed)
                assert transition(nxt, operation, traces, prediction, crossed) == state
                assert sorted(nxt.ports + nxt.reserve) == [0, 0, 1, 2]
                state = nxt
                phase = (phase + 1) % len(program)
                if phase == 0:
                    break
            assert state == initial
            cases += 1
    return cases


def run():
    return {"internal_feedback_and_interventions": check_mediation(),
            "orbit_occupation_cases": check_occupation(),
            "whole_clock_recurrences": check_reversible_clock(),
            "one_probe_maximum_occupied_orbits": 1,
            "two_probe_example_occupied_orbits_before_after": [1, 2]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    verify = parser.parse_args().check
    result = run()
    target = Path(__file__).with_name("feedback_routes_results.json")
    if verify:
        assert json.loads(target.read_text(encoding="utf-8")) == result
    else:
        with target.open("x", encoding="utf-8") as output:
            output.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
