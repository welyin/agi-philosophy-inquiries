"""Exhaustive two-observer reflection test on the six ordinary permutations."""

import argparse
import json
from itertools import permutations, product
from pathlib import Path


PERMS = list(permutations(range(3)))
ELEMENTS = range(len(PERMS))
IDENTITY = PERMS.index((0, 1, 2))
MUL = [[PERMS.index(tuple(left[i] for i in right)) for right in PERMS]
       for left in PERMS]
INV = [next(j for j in ELEMENTS if MUL[i][j] == IDENTITY) for i in ELEMENTS]
PROBES = (PERMS.index((1, 0, 2)), PERMS.index((0, 2, 1)))
RESULT = Path(__file__).with_name("reflective_attribution_results.json")


def conjugate(element, by):
    return MUL[MUL[INV[by]][element]][by]


def observe(state, probes):
    world, *readers = state
    return tuple((MUL[reader][world], MUL[MUL[reader][probe]][world])
                 for reader, probe in zip(readers, probes, strict=True))


def world_signature(record):
    return tuple(MUL[INV[base]][intervened] for base, intervened in record)


def possible_changes(before, after):
    return [h for h in ELEMENTS if tuple(conjugate(x, h) for x in before) == after]


def run():
    states = list(product(ELEMENTS, repeat=3))
    records = {s: observe(s, PROBES) for s in states}
    signatures = {s: world_signature(record) for s, record in records.items()}
    # This lookup sees only the observed signatures. Probe values are hidden.
    seen = set(signatures.values())
    decoders = {(old, new): possible_changes(old, new) for old, new in product(seen, repeat=2)}
    for (old, new), candidates in decoders.items():
        assert len(candidates) == 1
        for i in range(2):
            assert len(possible_changes((old[i],), (new[i],))) == 2

    corrected_responses = 0
    for old, new in product(states, repeat=2):
        change, = decoders[signatures[old], signatures[new]]
        assert change == MUL[INV[old[0]]][new[0]]
        for i in range(2):
            old_base, new_base = records[old][i][0], records[new][i][0]
            correction = MUL[MUL[old_base][change]][INV[new_base]]
            assert (correction == IDENTITY) == (old[i + 1] == new[i + 1])
            for fresh in ELEMENTS:
                assert MUL[correction][MUL[new[i + 1]][fresh]] == MUL[old[i + 1]][fresh]
                corrected_responses += 1

    # The same local transcript can require correcting one's reader or leaving it alone.
    for i, probe in enumerate(PROBES):
        reader_only = (IDENTITY, probe, probe)
        world_only = (probe, IDENTITY, IDENTITY)
        assert records[reader_only][i] == records[world_only][i]
        assert records[reader_only][1 - i] != records[world_only][1 - i]

    # A passive prediction can remain right after the world and both readers change.
    initial = records[(IDENTITY,) * 3]
    cancelled = records[(PROBES[0],) * 3]
    assert all(old[0] == new[0] for old, new in zip(initial, cancelled, strict=True))
    assert initial[0] == cancelled[0] and initial[1][1] != cancelled[1][1]

    drift_checks = 0
    for (world, a, b), gauge in product(states, ELEMENTS):
        changed = (MUL[gauge][world], MUL[a][INV[gauge]], MUL[b][INV[gauge]])
        drifting_probes = tuple(conjugate(g, INV[gauge]) for g in PROBES)
        assert observe(changed, drifting_probes) == records[(world, a, b)]
        drift_checks += 1
    return {
        "group_size": len(ELEMENTS),
        "epoch_states": len(states),
        "two_epoch_histories": len(states) ** 2,
        "candidates_for_world_change_single_observer": 2,
        "candidates_for_world_change_together": 1,
        "corrected_responses_on_all_next_inputs": corrected_responses,
        "simultaneous_changes_hidden_from_all_passive_readings": True,
        "joint_probe_drift_raw_transcript_checks": drift_checks,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.check:
        if result != json.loads(RESULT.read_text(encoding="utf-8")):
            raise SystemExit("Saved answers differ from recomputation.")
        print("Verified", RESULT.name)
    else:
        RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
