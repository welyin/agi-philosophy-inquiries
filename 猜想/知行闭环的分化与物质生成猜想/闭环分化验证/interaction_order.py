"""Compare two histories of the same reversible actions with exact probabilities."""

import argparse
import json
from collections import defaultdict
from fractions import Fraction
from itertools import product
from pathlib import Path


STATES = tuple(product((0, 1), repeat=2))


def evolve(state, actions, response):
    x, y = state
    for action in actions:
        if action == "A":
            x ^= 1
        elif action == "B":
            y ^= response[x]
        else:
            raise ValueError(action)
    return x, y


def push(distribution, operation):
    output = defaultdict(Fraction)
    for state, mass in distribution.items():
        output[operation(state)] += mass
    return dict(output)


def distance(left, right):
    return sum(abs(left.get(s, 0) - right.get(s, 0)) for s in left.keys() | right.keys()) / 2


def single_end_noise(distribution, error):
    return {(x, y): (1 - error) * distribution[x, y] + error * distribution[x, 1 - y]
            for x, y in STATES}


def run():
    biases = [Fraction(i, 4) for i in range(5)]
    examples = []
    cases = 0
    for response in STATES:  # Truth table (f(0), f(1)).
        depends_on_setting = response[0] != response[1]
        for state in STATES:
            first = evolve(state, "ABA", response)
            second = evolve(state, "BAA", response)
            assert first[0] == second[0] == state[0]
            assert first[1] ^ second[1] == int(depends_on_setting)
            assert evolve(state, "BB", response) == state
        for p in biases:
            initial = {(x, y): (p if y == 0 else 1 - p) / 2 for x, y in STATES}
            formed = [push(initial, lambda s: evolve(s, word, response)) for word in ("ABA", "BAA")]
            for error in biases:
                left, right = formed
                for steps in range(7):
                    expected = int(depends_on_setting) * abs(2 * p - 1) * abs(1 - 2 * error) ** steps
                    assert distance(left, right) == expected
                    for index in (0, 1):
                        assert push(left, lambda s: s[index]) == push(right, lambda s: s[index])
                    left_parity = push(left, lambda s: s[0] ^ s[1])
                    right_parity = push(right, lambda s: s[0] ^ s[1])
                    assert distance(left_parity, right_parity) == expected
                    read_y = lambda s: evolve(s, "B", response)[1]
                    assert distance(push(left, read_y), push(right, read_y)) == expected
                    if response == (0, 1) and p == 1 and error == Fraction(1, 4) and steps < 4:
                        examples.append(str(expected))
                    left = single_end_noise(left, error)
                    right = single_end_noise(right, error)
                    cases += 1
    return {"truth_tables": 4, "exact_probability_comparisons": cases,
            "example": {"f": [0, 1], "initial_probability_y_zero": 1,
                        "single_end_flip_probability": "1/4", "TV_after_0_to_3_steps": examples}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    verify = parser.parse_args().check
    result = run()
    path = Path(__file__).with_name("interaction_order_results.json")
    if verify:
        assert json.loads(path.read_text(encoding="utf-8")) == result
    else:
        with path.open("x", encoding="utf-8") as output:
            output.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
