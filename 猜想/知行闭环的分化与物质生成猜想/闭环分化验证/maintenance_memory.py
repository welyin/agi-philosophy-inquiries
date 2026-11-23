"""Independent reversible memory check for maintenance with at most one local fault."""

import argparse
import json
from itertools import product
from pathlib import Path


BITS = (0, 1)
WORDS = range(4)
FAULTS = ("none", "relation", "port0", "port1")


def inject(data, fault):
    # Test fixture only. The downstream repair does not receive the label.
    q, port0, port1 = data
    return q ^ (fault == "relation"), port0 ^ (fault == "port0"), port1 ^ (fault == "port1")


def syndrome(q, port0, port1):
    # High bit refers to actual arrival at port 0, low bit to arrival at port 1.
    return 2 * (q ^ port0) + (q ^ port1)


def maintain(data, work, crossed, template):
    q, port0, port1 = data
    if q ^ crossed:
        port0, port1 = port1, port0
    work ^= syndrome(q, port0, port1)
    if work == 1:
        port1 ^= 1
    elif work == 2:
        port0 ^= 1
    elif work == 3:
        q ^= 1
        # Correct the relation-controlled routing error as well.
        port0, port1 = port1, port0
    return (q, port0, port1), work ^ template


def infer_fault(r, crossed, word):
    if word == 0:
        return "none"
    if word == 3:
        return "relation"
    arrival = 0 if word == 2 else 1
    departure = arrival ^ r ^ crossed
    return "port0" if departure == 0 else "port1"


def check_permutations():
    all_states = set(product(product(BITS, repeat=3), WORDS))
    for crossed, template in product(BITS, WORDS):
        images = {maintain(data, work, crossed, template) for data, work in all_states}
        assert images == all_states
    return len(all_states) * len(BITS) * len(WORDS)


def check_words(horizon):
    counts = {}
    for r, crossed in product(BITS, repeat=2):
        clean = (r, r, r)
        for t in range(horizon + 1):
            memories = set()
            for history in product(FAULTS, repeat=t):
                data = clean
                memory = []
                for fault in history:
                    data, word = maintain(inject(data, fault), 0, crossed, 0)
                    assert data == clean
                    memory.append(word)
                assert tuple(infer_fault(r, crossed, w) for w in memory) == history
                memories.add(tuple(memory))
            assert len(memories) == 4 ** t
            counts[str(t)] = len(memories)
    return counts


def check_learned_templates(repetitions):
    changed = 0
    checks = 0
    for r, crossed, initial_fault in product(BITS, BITS, FAULTS):
        clean = (r, r, r)
        data, template = maintain(inject(clean, initial_fault), 0, crossed, 0)
        assert data == clean
        # First word remains in the template bank. The work bank starts blank.
        work = 0
        for _ in range(repetitions):
            data, work = maintain(inject(data, initial_fault), work, crossed, template)
            assert data == clean and work == 0
            checks += 1
        for fault in FAULTS:
            restored, residual = maintain(inject(data, fault), work, crossed, template)
            assert restored == clean
            assert infer_fault(r, crossed, residual ^ template) == fault
            assert (residual != 0) == (fault != initial_fault)
            changed += int(residual != 0)
    return {"learned_template_repetitions_checked": checks,
            "changed_faults_corrected_with_nonzero_residual": changed}


def run():
    return {"full_reversible_transition_cases": check_permutations(),
            "distinct_memory_states_per_fixed_r_and_w": check_words(5),
            **check_learned_templates(32)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    verify = parser.parse_args().check
    result = run()
    target = Path(__file__).with_name("maintenance_memory_results.json")
    if verify:
        assert json.loads(target.read_text(encoding="utf-8")) == result
    else:
        with target.open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
