"""Can the returns distinguish a requested change from the same unrequested edit?"""

import argparse
import json
from itertools import product
from pathlib import Path

from return_paths import launch, make_graph, travel


def flip_recent(traces):
    return (traces[0], 1 - traces[1])


def shared_transition(traces, requested, keep_receipt):
    # Both internal events have the same action on the traces. Only a retained
    # source record can affect future queries in this finite example.
    return flip_recent(traces), int(requested) if keep_receipt else None


def act_on_disagreement(traces, visits, receipt):
    if len(visits) != 2 or len(set(visits.values())) != 1:
        raise ValueError("requires two concordant completed probes")
    prior, observed = next(iter(visits.values()))
    if prior == observed or receipt not in (0, 1):
        raise ValueError("requires a changed reading and a known source record")
    return traces if receipt else flip_recent(traces)


def run():
    graph = make_graph(5)
    actions = set()
    blinded = []
    accepted = repaired = rechecked = 0
    for traces in product((0, 1), repeat=2):
        old = traces[0] ^ traces[1]
        launched, _ = launch(graph, 0, traces, 1 - old, 2)
        calibrated, _ = travel(graph, 0, traces, launched, len(graph))
        paired_futures = []
        for requested in (False, True):
            changed, absent = shared_transition(traces, requested, False)
            learned, reports = travel(graph, 0, changed, calibrated, len(graph))
            assert absent is None and len(reports) == 1
            visits = reports[0]
            paired_futures.append((changed, learned, visits))
            blinded.append((old, requested))
            actions.add(tuple(sorted(visits.values())))

            same_changed, receipt = shared_transition(traces, requested, True)
            assert same_changed == changed
            acted = act_on_disagreement(changed, visits, receipt)
            desired = changed if requested else traces
            assert acted == desired
            # Both carriers revisit the traces after a real keep/repair operation.
            again, _ = travel(graph, 0, acted, learned, len(graph))
            assert all(value == (desired[0] ^ desired[1]) for value in again.values())
            rechecked += 1
            accepted += int(requested)
            repaired += int(not requested)
        assert paired_futures[0] == paired_futures[1]
    assert len(actions) == 2
    # A binary choice for each of the two possible disagreement transcripts.
    scores = [
        sum(bool(accept[old]) == requested for old, requested in blinded)
        for accept in product((0, 1), repeat=2)
    ]
    assert scores == [4] * 4
    return {"paired_indistinguishable_histories": 4, "observable_transcripts": len(actions),
            "receipt_free_policies_correct_out_of_eight": scores,
            "with_receipt_changes_accepted": accepted, "with_receipt_edits_repaired": repaired,
            "post_action_dual_rechecks": rechecked}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    check = parser.parse_args().check
    results = run()
    path = Path(__file__).with_name("record_stability_results.json")
    if check:
        assert results == json.loads(path.read_text(encoding="utf-8"))
    else:
        with path.open("x", encoding="utf-8") as output:
            output.write(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results))
