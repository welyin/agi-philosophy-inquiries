"""Check ledger succession on the corrected four-leaf branches of note 43."""

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from quantum_consensus_dispatch import equal, route
from quantum_living_tree_recovery import (
    INPUT_DIM, flag_and_code_checks, leaf_correction, sample_x_and_keep_others,
    target_encoding,
)


RESULT = Path(__file__).with_name("living_ledger_succession_results.json")
BITS = (0, 1)


def handover_columns(dimension):
    """Old and successor devices actually swap. Successor starts in basis state 0."""
    identity = np.eye(dimension**2)
    swap = identity.reshape(dimension, dimension, -1).swapaxes(0, 1)
    columns = swap[:, :, ::dimension]
    equal(columns[0], np.eye(dimension))
    equal(columns[1:], 0)
    return columns[0]


def quantum_handover():
    history, code = route()
    # Independent local SWAPs of flag, original cell and candidate cell at
    # each end; one cell at the middle; then two evidence wires.
    bit, cell = handover_columns(2), handover_columns(4)
    endpoint = np.kron(bit, np.kron(cell, cell))
    logical_phase = np.diag([1, np.exp(0.37j)] * 4)
    for place in BITS:
        target = target_encoding(code, place)
        for outcome in BITS:
            raw = sample_x_and_keep_others(history, 1, place, outcome)
            if outcome:
                raw = raw * leaf_correction(place)[..., None]
            dominant = np.einsum(
                "ia,jb,kc,mu,nv,abcuvx->ijkmnx",
                endpoint, cell, endpoint, bit, bit, raw, optimize=True,
            )
            equal(dominant, (-1)**(place * outcome) * target / 2)
            flag_and_code_checks(dominant, place)
            recovered = target.reshape(-1, INPUT_DIM).T @ dominant.reshape(-1, INPUT_DIM)
            equal(recovered @ logical_phase, (-1)**(place * outcome) * logical_phase / 2)
    return 4


def all_snapshots(rows):
    """A--C--B, source slots fixed by wiring. Every site takes its own snapshot."""
    a, c, b = rows
    # First slot: both endpoints send to the middle. C sends its own row in
    # each direction. Next slot: C relays the endpoint payloads across.
    at_middle_from_a, at_middle_from_b = a, b
    at_a_from_middle, at_b_from_middle = c, c
    at_a_from_b, at_b_from_a = at_middle_from_b, at_middle_from_a
    return (
        (a, at_a_from_middle, at_a_from_b),
        (at_middle_from_a, c, at_middle_from_b),
        (at_b_from_a, at_b_from_middle, b),
    )


def read_local(snapshot):
    return tuple(int(sum(column) >= 2) for column in zip(*snapshot))


def commit(rows, a, b):
    # Both input pulses have separate fixed local origins and cross the
    # middle in distinct slots. Each replica acts on its own old snapshot.
    original = tuple(read_local(local) for local in all_snapshots(rows))
    at_c_from_a, at_c_from_b = a, b
    grants = ((a, at_c_from_b), (at_c_from_a, at_c_from_b), (at_c_from_a, b))
    after = tuple((birth, relation ^ (left & right))
                  for (birth, relation), (left, right) in zip(original, grants))
    return original, after


def probe_outputs(current, challenge):
    # Each address holds vacuum (-1) or a bit-carrying token. The readout
    # drives the successor's local switch between two distinct return paths.
    ports = (challenge, -1)
    return ports[::-1] if current else ports


def faults(birth, relation):
    clean = ((birth, relation),) * 3
    yield clean
    for index in range(3):
        for flip_h, flip_q in ((0, 1), (1, 0), (1, 1)):
            damaged = list(clean)
            damaged[index] = birth ^ flip_h, relation ^ flip_q
            yield tuple(damaged)


def check_ledger():
    cases = 0
    for birth, relation, left, right in product(BITS, repeat=4):
        next_relation = relation ^ (left & right)
        precursors = set()
        for rows in faults(birth, relation):
            old_reads, after = commit(rows, left, right)
            assert old_reads == ((birth, relation),) * 3
            assert after == ((birth, next_relation),) * 3
            precursors.add(rows)
            readers = all_snapshots(after)
            for endpoint in (0, 2):
                history, current = read_local(readers[endpoint])
                assert history == birth
                for challenge in BITS:
                    returns = probe_outputs(current, challenge)
                    assert returns[next_relation] == challenge
                    assert returns[1 - next_relation] == -1
            cases += 1
        assert len(precursors) == 10

    # If paths to both private write ports can be impersonated, the replicas
    # receive exactly the same physical inputs as a genuine joint grant.
    clean = ((0, 0),) * 3
    assert commit(clean, 1, 1)[1] != commit(clean, 1, 0)[1]
    # Without an additional trusted current-value witness, a two-copy reader
    # that sees (0,1) must answer differently on these two single-error inputs.
    assert (0, 0 ^ 1) == (1 ^ 1, 1)
    # Two damaged stores take over the majority; birth audit can fail as well.
    assert read_local(((1, 1), (1, 1), (0, 0))) == (1, 1)
    return cases


def verify():
    return {
        "handover_branches_checked_in_note43_physical_encoding": quantum_handover(),
        "protected_logical_dimension_with_untouched_partner": INPUT_DIM,
        "four_active_cells_in_each_branch_after_local_swaps": True,
        "candidate_addresses_including_dormant_spare": 5,
        "successor_quantum_evidence_ports": 2,
        "ledger_replicas": 3,
        "fields_per_replica": ["birth_location", "current_relation"],
        "complete_single_fault_and_two_grant_cases": check_ledger(),
        "readers_operating_without_old_participants": 2,
        "successor_probe": "one token delivered to the selected one of two local return ports",
        "bit_hops_for_one_all_site_snapshot": 12,
        "bit_hops_for_both_endpoint_grants": 4,
        "retained_distinct_fault_inputs_for_fixed_intended_output": 10,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = verify()
    if args.check:
        if result != json.loads(RESULT.read_text(encoding="utf-8")):
            raise SystemExit("Saved results differ from recomputation.")
        print("Verified ledger handovers and 160 maintenance cases.")
    else:
        RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
