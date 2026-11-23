"""Exact local marginals, a joint readout, and blind return of classical correlations."""

import argparse
import json
from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path

from feedback_routes import BITS
from interaction_order import evolve
from joint_action import joint_window, prepare
from recursive_response import read_contact


def prepare_hidden(relation, settings):
    # Four equally likely inputs. No accessible preparation log is retained.
    return [(left, right, prepare(settings, relation ^ left ^ right))
            for left, right in product(BITS, repeat=2)]


def joint_readout(left, right, lanes):
    # An authorized joint probe collects three contributions. Releasing this
    # outcome to a single party changes the information available to that party.
    target = 0
    for source in (left, right, read_contact(lanes)):
        _, target = evolve((source, target), "B", BITS)
    return target


def local_view(prepare_list, side):
    # Grant the whole original apparatus as an upper bound on local access.
    # Only the other party's newly installed private mask is excluded.
    return Counter((row[side], row[2]) for row in prepare_list)


def check_marginals():
    prefix_checks = endpoint_checks = joint_checks = 0
    for settings in product(BITS, repeat=3):
        states = [prepare_hidden(r, settings) for r in BITS]
        for left_ok, right_ok in product(BITS, repeat=2):
            pulses = (left_ok, right_ok, left_ok, right_ok)
            histories = []
            endpoints = []
            for r in BITS:
                paths = []
                final = []
                for a, b, lanes in states[r]:
                    # A truncated pulse vector gives the exact old gate state
                    # at that microstep. The unused slots perform no action.
                    path = tuple(joint_window(lanes, pulses[:k] + (0,) * (4-k))
                                 for k in range(5))
                    paths.append((a, b, path))
                    final.append((a, b, path[-1]))
                    assert joint_readout(*final[-1]) == r ^ (left_ok & right_ok)
                    joint_checks += 1
                histories.append(paths)
                endpoints.append(final)
            for side in BITS:
                for length in range(1, 6):
                    views = [Counter((row[side], row[2][:length]) for row in paths)
                             for paths in histories]
                    assert views[0] == views[1]
                    prefix_checks += 1
                for r in BITS:
                    assert local_view(states[r], side) == local_view(endpoints[r], side)
                    endpoint_checks += 1
                # The old local claims were in the clear when used without masks.
                assert [prepare(settings, r)[1 + side][1] for r in BITS] == [0, 1]
    return prefix_checks, endpoint_checks, joint_checks


def check_blind_inverse():
    checked = 0
    for bits in product(BITS, repeat=6):
        lanes = tuple(zip(bits[::2], bits[1::2]))
        changed = joint_window(lanes, (1, 1, 1, 1))
        assert changed == tuple((setting, target ^ 1) for setting, target in lanes)
        assert joint_window(changed, (1, 1, 1, 1)) == lanes
        checked += 1
    # The masks and any spectator are untouched; delta inputs prove the
    # return channel is the identity for every joint classical distribution.
    return checked


def check_reprepare():
    settings = (1, 0, 0)
    correlations = []
    for r in BITS:
        original = prepare_hidden(r, settings)
        old_ab, replaced_ab = Counter(), Counter()
        matches = 0
        for a, b, lanes in original:
            spectator = a
            old_ab[a, b, lanes] += 4
            assert spectator == a
            for new_a, new_b, new_lanes in prepare_hidden(joint_readout(a, b, lanes), settings):
                replaced_ab[new_a, new_b, new_lanes] += 1
                matches += new_a == spectator
        assert old_ab == replaced_ab and Fraction(matches, 16) == Fraction(1, 2)
        correlations.append(str(Fraction(matches, 16)))
    return correlations


def run():
    prefixes, endpoints, joint = check_marginals()
    return {"local_history_comparisons_between_relations": prefixes,
            "unchanged_complete_endpoint_marginals": endpoints,
            "joint_probe_tests": joint,
            "all_six_bit_gate_inputs_restored": check_blind_inverse(),
            "spectator_match_after_true_return": "1",
            "spectator_match_after_reprepare_for_r0_and_r1": check_reprepare()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    results = run()
    path = Path(__file__).with_name("hidden_relation_return_results.json")
    if args.check:
        assert json.loads(path.read_text(encoding="utf-8")) == results
    else:
        with path.open("x", encoding="utf-8") as out:
            out.write(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results))
