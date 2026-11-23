"""Check which delayed questions protect common corrections after one encoding."""

import argparse
import json
from fractions import Fraction
from itertools import combinations
from pathlib import Path


RESULT = Path(__file__).with_name("recombined_correction_tasks_results.json")


def parity(mask, history):
    return (mask & history).bit_count() % 2


def erasure_check(bits):
    """Enumerate the posterior support for every record under a uniform source."""
    size = 1 << bits
    supports = {}
    for history in range(size):
        for erased in range(bits):
            record = (erased, history & ~(1 << erased))
            supports.setdefault(record, []).append(history)

    fixed_errors = [0] * bits
    random_parity_errors = 0
    adaptive_errors = 0
    whole_history_errors = 0
    for (erased, _), candidates in supports.items():
        assert len(candidates) == 2
        for coordinate in range(bits):
            ones = sum(parity(1 << coordinate, x) for x in candidates)
            fixed_errors[coordinate] += min(ones, 2 - ones)
        for mask in range(size):
            ones = sum(parity(mask, x) for x in candidates)
            random_parity_errors += min(ones, 2 - ones)
        adaptive_errors += min(sum(parity(1 << erased, x) == bit for x in candidates)
                               for bit in (0, 1))
        whole_history_errors += len(candidates) - 1

    events = size * bits
    assert all(Fraction(e, events) == Fraction(1, 2 * bits) for e in fixed_errors)
    assert Fraction(random_parity_errors, events * size) == Fraction(1, 4)
    assert Fraction(adaptive_errors, events) == Fraction(1, 2)
    assert Fraction(whole_history_errors, events) == Fraction(1, 2)
    return {
        "history_bits": bits,
        "record_symbols_including_erasure_index": len(supports),
        "fixed_coordinate_error": str(Fraction(fixed_errors[0], events)),
        "uniform_delayed_parity_error": str(Fraction(random_parity_errors, events * size)),
        "record_dependent_query_error": str(Fraction(adaptive_errors, events)),
        "whole_history_error": str(Fraction(whole_history_errors, events)),
    }


def distance_check(bits):
    """Check the minimum code distance and every possible returned answer word."""
    length = 1 << bits
    words = [sum(parity(q, x) << q for q in range(length)) for x in range(length)]
    min_distance = min((a ^ b).bit_count() for a, b in combinations(words, 2))
    min_errors_on_failure = length
    for answers in range(1 << length):
        distances = [(answers ^ word).bit_count() for word in words]
        decoded = distances.index(min(distances))
        for history, errors in enumerate(distances):
            if history != decoded:
                min_errors_on_failure = min(min_errors_on_failure, errors)
                assert 4 * errors >= length
    assert 2 * min_distance == length
    return {
        "history_bits": bits,
        "parity_questions_including_zero": length,
        "answer_words_exhausted": 1 << length,
        "minimum_code_distance": min_distance,
        "minimum_bit_errors_on_wrong_decoding": min_errors_on_failure,
    }


def quotient_check():
    """A record may merge histories precisely when every permitted answer agrees."""
    # Two origins can lead to the same written record in a noisy classical encoding.
    supports = [{0, 1}, {1, 2}, {3}, {4, 5}]
    truth = [0, 0, 1, 1]

    def conflicts(signatures):
        return [[i, j] for i, j in combinations(range(len(supports)), 2)
                if supports[i] & supports[j] and signatures[i] != signatures[j]]

    original = conflicts(truth)
    refined = conflicts(list(range(len(truth))))
    assert original == [] and refined == [[0, 1]]
    return {"one_question_conflicts": original, "identify_origin_conflicts": refined}


def calculate():
    return {
        "exact_task_quotient": quotient_check(),
        "single_random_erasure": [erasure_check(bits) for bits in (1, 2, 3, 4, 6, 8)],
        "whole_history_via_parity": [distance_check(bits) for bits in (1, 2, 3)],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="compare with saved JSON")
    args = parser.parse_args()
    result = calculate()
    if args.check:
        if result != json.loads(RESULT.read_text(encoding="utf-8")):
            raise SystemExit("saved result differs")
        print("PASS: recombined correction tasks")
    else:
        with RESULT.open("x", encoding="utf-8") as out:
            out.write(json.dumps(result, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
