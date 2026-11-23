"""Replay two-epoch attribution with individually fallible intervention ports on S3."""

import argparse
import json
from itertools import combinations, combinations_with_replacement, permutations, product
from pathlib import Path


GROUP = tuple(permutations(range(3)))
IDENTITY = (0, 1, 2)
A, B, C = (1, 0, 2), (0, 2, 1), (2, 1, 0)
TURN = (1, 2, 0)
RESULT = Path(__file__).with_name("reflective_probe_drift_results.json")


def compose(left, right):
    return tuple(left[i] for i in right)


def inverse(element):
    return tuple(element.index(i) for i in IDENTITY)


def conjugate(element, transform):
    return compose(compose(inverse(transform), element), transform)


def raw_epoch(world, rules, ports):
    return tuple((compose(rule, world), compose(compose(rule, port), world))
                 for rule, port in zip(rules, ports, strict=True))


def relative_responses(raw):
    return tuple(compose(inverse(baseline), probed) for baseline, probed in raw)


def distance(left, right):
    return sum(a != b for a, b in zip(left, right, strict=True))


def codebook(initial):
    return {change: tuple(conjugate(port, change) for port in initial)
            for change in GROUP}


def minimum_distance(book):
    return min(distance(first, second) for first, second in combinations(book.values(), 2))


def one_fault_inputs(ports):
    yield ports
    for i, original in enumerate(ports):
        for replacement in GROUP:
            if replacement != original:
                yield ports[:i] + (replacement,) + ports[i + 1:]


def check_all_changes(ports):
    epochs = 0
    holdout_checks = 0
    initial_world = TURN
    initial_rules = GROUP[:len(ports)]
    initial = raw_epoch(initial_world, initial_rules, ports)
    book = codebook(relative_responses(initial))
    assert minimum_distance(book) == 3
    for change in GROUP:
        final_world = compose(initial_world, change)
        for changed_ports in one_fault_inputs(ports):
            for rule_shift in GROUP:
                # Different participants may start with and acquire different rules.
                final_rules = tuple(compose(rule_shift, inverse(rule)) for rule in initial_rules)
                final = raw_epoch(final_world, final_rules, changed_ports)
                measured = relative_responses(final)
                compatible = [h for h, word in book.items() if distance(measured, word) <= 1]
                assert compatible == [change]
                decoded = compatible[0]
                for (old_b, _), (new_b, _), old_rule, new_rule in zip(
                    initial, final, initial_rules, final_rules, strict=True
                ):
                    correction = compose(compose(old_b, decoded), inverse(new_b))
                    for fresh_world in GROUP:
                        assert compose(correction, compose(new_rule, fresh_world)) == compose(
                            old_rule, fresh_world
                        )
                        holdout_checks += 1
                epochs += 1
    return {"two_epoch_cases": epochs, "fresh_object_readings": holdout_checks}


def collision_example():
    ports = (A, B, C)
    old_rules = (IDENTITY,) * 3
    old_raw = raw_epoch(IDENTITY, old_rules, ports)
    first_final_ports = (A, B, B)
    second_final_ports = (A, C, C)
    final_raw = raw_epoch(IDENTITY, old_rules, first_final_ports)
    second_final_raw = raw_epoch(A, (A,) * 3, second_final_ports)
    assert final_raw == second_final_raw
    assert distance(ports, first_final_ports) == distance(ports, second_final_ports) == 1
    book = codebook(relative_responses(old_raw))
    observed = relative_responses(final_raw)
    candidates = [h for h, word in book.items() if distance(observed, word) <= 1]
    assert set(candidates) == {IDENTITY, A}
    return {"raw_old": old_raw, "raw_new_in_both_worlds": final_raw,
            "compatible_world_changes": candidates}


def run():
    # For three ports, include identities, repeats and all order types.
    best_three = max(minimum_distance(codebook(ports))
                     for ports in combinations_with_replacement(GROUP, 3))
    assert best_three == 2
    designs = {"three_different_swaps": (A, B, C),
               "three_swaps_and_one_turn": (A, B, C, TURN),
               "five_swaps": (A, A, B, B, C)}
    distances = {name: minimum_distance(codebook(ports)) for name, ports in designs.items()}
    assert tuple(distances.values()) == (2, 3, 3)
    # Check the closed centralizer formula against every pair of codewords.
    for ports in product(GROUP, repeat=3):
        formula = min(sum(compose(p, h) != compose(h, p) for p in ports)
                      for h in GROUP if h != IDENTITY)
        assert formula == minimum_distance(codebook(ports))
    return {"distance_by_design": distances, "best_distance_of_any_three_ports": best_three,
            "all_ordered_three_port_designs_checked": len(GROUP) ** 3,
            "four_port_recovery": check_all_changes(designs["three_swaps_and_one_turn"]),
            "three_port_collision": collision_example()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = json.dumps(run(), ensure_ascii=False, indent=2) + "\n"
    if args.check:
        if json.loads(RESULT.read_text(encoding="utf-8")) != json.loads(result):
            raise SystemExit("Saved drift results differ from recomputation.")
        print("Verified", RESULT.name)
    else:
        RESULT.write_text(result, encoding="utf-8")
