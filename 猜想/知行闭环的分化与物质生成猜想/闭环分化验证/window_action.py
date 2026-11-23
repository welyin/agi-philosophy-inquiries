"""Asynchronous deposits for one isolated object and one finite four-phase window."""

import argparse
import json
from dataclasses import dataclass
from itertools import product
from pathlib import Path

from collective_maintenance import FAULTS, World, maintain, perturb, read_relation, tags
from feedback_routes import BITS, PAIRS, circulate
from joint_action import as_world, joint_window, prepare


@dataclass(frozen=True, eq=False)
class Contact:
    # The exact instance is preassigned to one physical socket. A new object
    # with the same labels is not the assigned contact.
    side: int
    stage: int


class LocalHandle:
    def __init__(self, contacts, in_transit):
        self.remaining = list(contacts)
        self.in_transit = in_transit

    def issue(self, stage):
        if not self.remaining or self.remaining[0].stage != stage:
            return None
        contact = self.remaining.pop(0)
        self.in_transit.add(contact)
        return contact


class Window:
    def __init__(self):
        self.sockets = tuple(tuple(Contact(side, stage) for stage in BITS) for side in BITS)
        self.in_transit = set()
        self.handles = tuple(LocalHandle(row, self.in_transit) for row in self.sockets)
        self.deposits = set()
        self.closed = False

    def deliver(self, ingress, contact):
        if self.closed or contact not in self.in_transit or contact not in self.sockets[ingress]:
            return False
        self.in_transit.remove(contact)
        self.deposits.add(contact)
        return True

    def close(self, lanes, fault=None):
        if self.closed:
            raise ValueError("window already closed")
        self.closed = True
        # A side can release only its own complete pair. The clock feeds the
        # actual L,R,L,R executor; it is not given a common approval or target.
        left, right = (all(c in self.deposits for c in row) for row in self.sockets)
        return joint_window(lanes, (left, right, left, right), fault)


def issue_four(window):
    return tuple(handle.issue(stage) for handle in window.handles for stage in BITS)


def close_history(events, start, fault=None):
    window = Window()
    contacts = issue_four(window)
    for index in events:
        window.deliver(index // 2, contacts[index])
    return window.close(start, fault)


def check_collection():
    start = prepare((0, 1, 0), 0)
    cases = 0
    for length in range(7):
        for events in product(range(4), repeat=length):
            desired = int(set(events) == set(range(4)))
            assert close_history(events, start) == prepare((0, 1, 0), desired)
            cases += 1
    # All transitions from the 16 reachable staging subsets, including repeated
    # submissions, are checked separately, covering arbitrary finite duplicates.
    transitions = 0
    for bits in product(BITS, repeat=4):
        for new in range(4):
            window = Window()
            contacts = issue_four(window)
            present = {c for c, bit in zip(contacts, bits) if bit}
            for c in present:
                assert window.deliver(c.side, c)
            assert window.deliver(new // 2, contacts[new]) == (contacts[new] not in present)
            assert window.deposits == present | {contacts[new]}
            transitions += 1
    return cases, transitions


def check_boundaries():
    settings = (0, 1, 0)
    start = prepare(settings, 0)
    for side in BITS:
        sole = Window()
        handle = sole.handles[side]
        assert sole.deliver(side, handle.issue(0))
        assert sole.deliver(side, handle.issue(1))
        assert sole.close(start) == start  # The other participant never issued.
    current = Window()
    left, right = current.handles
    # Naming even an original socket cannot transfer a contact still held at home.
    assert not current.deliver(0, current.sockets[0][1])
    assert right.issue(1) is None  # Seal cannot be fabricated before preparation.
    left_prep = left.issue(0)
    assert left.issue(0) is None
    assert current.deliver(0, left_prep)
    assert not current.deliver(0, left_prep)
    assert not current.deliver(0, Contact(0, 1))
    assert not current.deliver(1, left_prep)
    # The right side finishes; the left side never delivers its second contact.
    assert current.deliver(1, right.issue(0))
    assert current.deliver(1, right.issue(1))
    assert current.close(start) == start
    left_seal = left.issue(1)
    assert not current.deliver(0, left_seal)
    try:
        current.close(start)
    except ValueError:
        pass
    else:
        raise AssertionError("closed window executed twice")

    for source in ("previous window", "different object"):
        # Distinct live socket families in both cases. Finite labels alone
        # would not distinguish a delayed contact when their names wrap.
        previous = Window()
        new = Window()
        stale = issue_four(previous)[:2]
        fresh = issue_four(new)[2:]
        if source == "previous window":
            previous.close(start)
        for c in stale:
            assert not new.deliver(0, c)
        for c in fresh:
            assert new.deliver(1, c)
        assert new.close(start) == start
        careless = {(c.side, c.stage) for c in stale + fresh}
        assert careless == set(product(BITS, repeat=2))
        assert joint_window(start, (1, 1, 1, 1)) == prepare(settings, 1)


def check_data_faults():
    # Representative of every sender-pair readiness class; collection is
    # independently checked for all prefixes and transition subsets above.
    faults = (None,) + tuple(product(range(5), range(3), BITS))
    clean = damaged = 0
    for settings, old, readiness in product(product(BITS, repeat=3), BITS, PAIRS):
        start = prepare(settings, old)
        events = [3, 1]  # both second contacts arrive before the first ones
        for side in BITS:
            if readiness[side]:
                events.append(2 * side)
        desired = old ^ int(all(readiness))
        for fault in faults:
            lanes = close_history(events, start, fault)
            for wire, first in product(BITS, repeat=2):
                actual = as_world(lanes, first)
                corrected, _ = maintain(actual, (0, 0), wire)
                assert read_relation(corrected) == desired
                assert all(belief == desired for _, belief in corrected.ports)
                order = (first, 1 - first)
                assert tags(corrected) == (order[::-1] if desired ^ wire else order)
                if fault is None:
                    clean += 1
                else:
                    damaged += 1
    return clean, damaged


def check_followup():
    cases = 0
    for old, wire, first, history in product(BITS, BITS, BITS, product(FAULTS, repeat=2)):
        start = prepare((0, 1, 0), old)
        # One real joint change. Neither sender contributes to later maintenance.
        live = as_world(close_history((2, 3, 1, 0), start), first)
        nominal = live
        for fault in history:
            actual = perturb(live, fault)
            live, _ = maintain(actual, (0, 0), wire)
            _, moved = circulate(nominal.marks, nominal.ports, wire)
            nominal = World(nominal.marks, moved)
            assert live == nominal
            assert read_relation(live) == 1 - old
        cases += 1
    return cases


def run():
    arrivals, transitions = check_collection()
    check_boundaries()
    clean, damaged = check_data_faults()
    return {"arrival_words_through_length_six": arrivals,
            "all_staging_subset_transitions": transitions,
            "clean_executions": clean, "single_lane_fault_executions": damaged,
            "post_joint_change_two_check_histories": check_followup(),
            "late_or_other_object_contacts_rejected": True,
            "unscoped_staging_can_merge_two_unilateral_submissions": True}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    verify = parser.parse_args().check
    result = run()
    path = Path(__file__).with_name("window_action_results.json")
    if verify:
        assert json.loads(path.read_text(encoding="utf-8")) == result
    else:
        with path.open("x", encoding="utf-8") as output:
            output.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
