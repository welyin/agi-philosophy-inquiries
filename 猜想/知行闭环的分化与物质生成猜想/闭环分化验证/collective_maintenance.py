"""Check reversible maintenance of the same marks that route the two carriers."""

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from pathlib import Path

from feedback_routes import BITS, PAIRS, circulate
from interaction_order import evolve


FAULTS = ("none", "relation", "port0", "port1")


@dataclass(frozen=True)
class World:
    marks: tuple
    ports: tuple  # Each carrier contains a tag and a private belief.


def perturb(world, fault):
    x, y = world.marks
    if fault == "relation":
        return World((x, 1 - y), world.ports)
    ports = list(world.ports)
    if fault in ("port0", "port1"):
        i = 0 if fault == "port0" else 1
        tag, belief = ports[i]
        ports[i] = (tag, 1 - belief)
    elif fault != "none":
        raise ValueError(fault)
    return World(world.marks, tuple(ports))


def xor_evidence(record, ports, observed):
    return tuple(bit ^ belief ^ observed for bit, (_, belief) in zip(record, ports))


def update_beliefs(ports, record, correction):
    return tuple((tag, belief ^ bit ^ correction)
                 for (tag, belief), bit in zip(ports, record))


def change_relation(marks, ports, wire, correction, redeliver):
    if not correction:
        return marks, ports
    changed = (marks[0], 1 - marks[1])
    if redeliver:
        # Exactly one of the two settings physically crosses. If it is the old
        # setting, make the extra passage before changing y; otherwise after.
        observed = evolve(marks, "B", BITS)[1]
        crossing_setting = marks if observed ^ wire else changed
        _, ports = circulate(crossing_setting, ports, wire)
    return changed, ports


def maintain(world, record, wire, *, writable=True, redeliver=True, inverse=False):
    # The pair of comparisons is indexed by arrival port. The coupler is allowed
    # to act only after both real carriers have arrived in this stable window.
    if inverse:
        correction = int(writable) & record[0] & record[1]
        marks_before, ports = change_relation(world.marks, world.ports, wire, correction, redeliver)
        arrivals = update_beliefs(ports, record, correction)
        observed = evolve(marks_before, "B", BITS)[1]
        record_before = xor_evidence(record, arrivals, observed)
        _, departed = circulate(marks_before, arrivals, wire)
        return World(marks_before, departed), record_before

    marks, arrivals = circulate(world.marks, world.ports, wire)
    observed = evolve(marks, "B", BITS)[1]
    assert evolve((marks[0], observed), "B", BITS) == marks
    evidence = xor_evidence(record, arrivals, observed)
    correction = int(writable) & evidence[0] & evidence[1]
    ports = update_beliefs(arrivals, evidence, correction)
    marks, ports = change_relation(marks, ports, wire, correction, redeliver)
    return World(marks, ports), evidence


def tags(world):
    return tuple(tag for tag, _ in world.ports)


def prepared(x, target):
    return World((x, x ^ target), ((0, target), (1, target)))


def read_relation(world):
    return evolve(world.marks, "B", BITS)[1]


def check_inverse():
    cases = 0
    for marks, beliefs, record, wire, writable, redeliver, order in product(
            PAIRS, PAIRS, PAIRS, BITS, BITS, BITS, BITS):
        world = World(marks, ((order, beliefs[0]), (1 - order, beliefs[1])))
        options = dict(writable=writable, redeliver=redeliver)
        output, evidence = maintain(world, record, wire, **options)
        assert maintain(output, evidence, wire, inverse=True, **options) == (world, record)
        cases += 1
    return cases


def check_local_controls():
    required = tuple(a & b for a, b in PAIRS)
    attainable = {tuple(left[a] ^ right[b] for a, b in PAIRS)
                  for left, right in product(PAIRS, repeat=2)}
    assert len(attainable) == 8 and required not in attainable
    return len(attainable)


def check_histories(length):
    cases = 0
    for x, target, wire in product(BITS, repeat=3):
        start = prepared(x, target)
        for faults in product(FAULTS, repeat=length):
            active = only_marks = passive = ideal = start
            object_faults = 0
            for fault in faults:
                object_faults += int(fault == "relation")
                actual = perturb(active, fault)
                _, intermediate = circulate(actual.marks, actual.ports, wire)
                _, expected = circulate(ideal.marks, ideal.ports, wire)
                wrong_transit = tuple(t for t, _ in intermediate) != tuple(t for t, _ in expected)
                assert wrong_transit == (fault == "relation")
                ideal = World(ideal.marks, expected)
                active, _ = maintain(actual, (0, 0), wire)
                only_marks, _ = maintain(perturb(only_marks, fault), (0, 0), wire, redeliver=False)
                passive, _ = maintain(perturb(passive, fault), (0, 0), wire, writable=False)
                assert active == ideal and read_relation(only_marks) == target
                nominal = tags(ideal)
                assert tags(only_marks) == (nominal[::-1] if object_faults % 2 else nominal)
                assert read_relation(passive) == target ^ (object_faults % 2)
                assert all(b == read_relation(passive) for _, b in passive.ports)
            cases += 1
    return cases


def check_lifetime(length):
    # Uniform IID four-symbol stream. Examine every prefix at the actual
    # checkpoint; a recovery at a later checkpoint cannot undo a first failure.
    counts = {"private_only": 0, "joint_write": 0}
    start = prepared(0, 1)
    for faults in product(FAULTS, repeat=length):
        for policy in counts:
            world = start
            for fault in faults:
                world, _ = maintain(perturb(world, fault), (0, 0), 0,
                                    writable=(policy == "joint_write"))
                if read_relation(world) != 1:
                    break
            else:
                counts[policy] += 1
    total = 4 ** length
    assert Fraction(counts["private_only"], total) == Fraction(3, 4) ** length
    assert counts["joint_write"] == total
    return {policy: str(Fraction(count, total)) for policy, count in counts.items()}


def check_limits():
    start = prepared(0, 0)
    both_wrong = perturb(perturb(start, "port0"), "port1")
    other = perturb(prepared(0, 1), "relation")
    assert both_wrong == other  # Same switch, tag order, wiring and available evidence.
    wrong, _ = maintain(both_wrong, (0, 0), 0)
    assert read_relation(wrong) == 1
    # A stable but damaged relation still allows a perfect two-lap echo.
    damaged = perturb(start, "relation")
    _, first = circulate(damaged.marks, damaged.ports, 0)
    _, second = circulate(damaged.marks, first, 0)
    assert first != damaged.ports and second == damaged.ports and read_relation(damaged) == 1


def run():
    check_limits()
    return {"reversible_full_state_cases": check_inverse(),
            "distinct_separate_flip_rules": check_local_controls(),
            "four_round_paired_histories": check_histories(4),
            "six_checkpoint_survival_equal_fault_probabilities": check_lifetime(6),
            "double_private_fault_indistinguishable_from_relation_fault": True,
            "uncorrected_relation_can_pass_two_lap_return_test": True}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    verify = parser.parse_args().check
    results = run()
    path = Path(__file__).with_name("collective_maintenance_results.json")
    if verify:
        assert json.loads(path.read_text(encoding="utf-8")) == results
    else:
        with path.open("x", encoding="utf-8") as output:
            output.write(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results))
