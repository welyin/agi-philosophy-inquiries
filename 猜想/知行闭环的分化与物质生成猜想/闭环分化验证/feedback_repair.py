"""Correct independently perturbed predictions on routes formed by feedback."""

import argparse
import json
from itertools import product
from pathlib import Path

from feedback_routes import BITS, PAIRS, State, circulate, write_relation
from interaction_order import evolve


def exchange(marks, ports, crossed, reverse=False):
    # Each carrier poses its belief as a binary probe. The return is its
    # disagreement with the actual relation. Retain it in a local work bit.
    readable = evolve(marks, "B", BITS)
    value = readable[1]
    updated = []
    for carrier in ports:
        if carrier is None:
            updated.append(None)
            continue
        tag, belief, record = carrier
        if reverse:
            belief ^= record
            record ^= belief ^ value
        else:
            record ^= belief ^ value
            belief ^= record
        updated.append((tag, belief, record))
    # Port permutation commutes with identical local payload updates and is
    # self-inverse. The marks are restored after each classical read.
    assert evolve(readable, "B", BITS) == marks
    return circulate(marks, tuple(updated), crossed)[1]


def run():
    inverse_cases = 0
    for marks, crossed, left, right in product(PAIRS, BITS, PAIRS, PAIRS):
        ports = ((1, *left), (2, *right))
        after = exchange(marks, ports, crossed)
        assert exchange(marks, after, crossed, reverse=True) == ports
        inverse_cases += 1

    correction_cases = 0
    for crossed in BITS:
        # Start with one coupled orbit and two live carriers, then separate
        # the paths by letting the original comparator choose the operation order.
        original = State((1 - crossed, 0), ports=(1, 2))
        assert circulate(original.marks, original.ports, crossed)[1] == (2, 1)
        separated = write_relation(original, (0, crossed), 0)
        marks = separated.marks
        assert circulate(marks, separated.ports, crossed)[1] == (1, 2)
        distinct_records = set()
        for left, right in PAIRS:
            ports = ((1, left, 0), (2, right, 0))
            after = exchange(marks, ports, crossed)
            assert after == ((1, crossed, left ^ crossed), (2, crossed, right ^ crossed))
            distinct_records.add(tuple(carrier[2] for carrier in after))
            # Remove the other carrier, preserving the actual shared device.
            for port in BITS:
                alone = list(ports)
                alone[1 - port] = None
                recovered = exchange(marks, alone, crossed)
                assert recovered[port] == after[port] and recovered[1 - port] is None
            # Alter the actual relation with claims, switch and wiring fixed:
            # both predictions follow it; the two carriers now cross ports.
            changed = (marks[0], 1 - marks[1])
            shifted = exchange(changed, ports, crossed)
            assert shifted[0][:2] == (2, 1 - crossed)
            assert shifted[1][:2] == (1, 1 - crossed)
            correction_cases += 1
        assert len(distinct_records) == 4
    return {"full_two_carrier_inverse_cases": inverse_cases,
            "split_route_independent_correction_cases": correction_cases,
            "required_distinct_residual_states_per_fixed_relation": 4}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    verify = parser.parse_args().check
    result = run()
    path = Path(__file__).with_name("feedback_repair_results.json")
    if verify:
        assert json.loads(path.read_text(encoding="utf-8")) == result
    else:
        with path.open("x", encoding="utf-8") as output:
            output.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
