"""A closed finite controller testing whether two internally retained traces agree."""
import argparse
import json
from itertools import product
from pathlib import Path


HERE = Path(__file__).resolve().parent
ACTIONS = ("wait", "exchange", "contrast")
HISTORIES = tuple(product((0, 1), repeat=2))


def step(traces, action):
    older, recent = traces
    if action == "wait":
        return traces
    if action == "exchange":
        return recent, older
    if action == "contrast":
        return older, recent ^ 1
    raise ValueError(action)


def read(traces, convention):
    return traces[0] ^ convention


def predict_step(reading, relation, action):
    if action == "exchange":
        return reading ^ relation, relation
    if action == "contrast":
        return reading, relation ^ 1
    if action == "wait":
        return reading, relation
    raise ValueError(action)


def probe(traces, convention, exchanges):
    before = read(traces, convention)
    for _ in range(exchanges):
        traces = step(traces, "exchange")
    after = read(traces, convention)
    return traces, before ^ after


def self_explore(seed, convention):
    # One controller. The comparison and the decision to change a trace are internal.
    traces = (seed, seed)
    before = read(traces, convention)
    after_wait = read(step(traces, "wait"), convention)
    traces, ledger = probe(traces, convention, 1)
    initial_ledger = ledger
    if before != after_wait:
        raise RuntimeError("wait disagreed with the stored reading")
    # Return to the original preparation; test a different word before altering it.
    traces = step(traces, "exchange")
    expected = ledger
    traces, observed = probe(traces, convention ^ 1, 3)
    if observed == expected:
        traces = step(traces, "contrast")
    # The stale claim is committed before the next probe.
    stale = ledger
    traces, observed = probe(traces, convention, 1)
    correction = (observed != stale)
    ledger = observed
    # Use the corrected record to choose which probes to perform.
    for count in ((3, 2) if ledger else (2, 3)):
        expected = ledger if count % 2 else 0
        traces, observed = probe(traces, convention ^ 1, count)
        if expected != observed:
            raise RuntimeError("held-out prediction failed")
    # Overwrite the claim alone. The same two physical traces are read again.
    ledger = 0
    traces, recovered = probe(traces, convention, 1)
    overwritten_claim_failed = ledger != recovered
    ledger = recovered
    return [initial_ledger, correction, overwritten_claim_failed, ledger]


def check_future_records(max_length):
    words = [w for n in range(max_length+1) for w in product(ACTIONS, repeat=n)]
    cases = 0
    for traces, convention, word in product(HISTORIES, (0, 1), words):
        current = traces
        raw = read(traces, convention)
        relation = traces[0] ^ traces[1]
        for action in word:
            current = step(current, action)
            raw, relation = predict_step(raw, relation, action)
            assert raw == read(current, convention)
            assert relation == current[0] ^ current[1]
        cases += 1
    return cases


def check_attribution(max_length):
    # eta is the response of a readout fault to an exchange. A contrast flips the
    # trace parity k, and neither action supplies a privileged reference for eta.
    # Registers here are in response coordinates, not a second observer.
    def transcript(k, eta, word):
        response, displayed = k, 0
        output = [displayed]
        for action in word:
            if action == "exchange":
                displayed ^= response ^ eta
            elif action == "contrast":
                response ^= 1
            output.append(displayed)
        return output

    cases = 0
    for n in range(max_length+1):
        for word in product(ACTIONS, repeat=n):
            assert transcript(1, 0, word) == transcript(0, 1, word)
            cases += 1
    return cases


def check_reference_preparation():
    # A persistent readout offset flips at every exchange if eta is set.
    def compare_pair(traces, offset, eta):
        first = read(traces, offset)
        traces = step(traces, "exchange")
        offset ^= eta
        return traces, offset, first ^ read(traces, offset)

    cases = 0
    for (x,y), offset, eta in product(HISTORIES, (0,1), (0,1)):
        traces, offset, original = compare_pair((x,y), offset, eta)
        matched = (traces[0], traces[0])
        _, _, reference = compare_pair(matched, offset, eta)
        assert reference == eta
        assert original ^ reference == x ^ y
        cases += 1
    return cases


def reachable_ledgers(n, independent):
    start = (0,)*n
    found, todo = {start}, [start]
    masks = [tuple(int(i == j) for i in range(n)) for j in range(n)]
    if not independent:
        masks = [(1,)*n]
    while todo:
        current = todo.pop()
        for mask in masks:
            nxt = tuple(bit ^ flip for bit, flip in zip(current, mask))
            if nxt not in found:
                found.add(nxt)
                todo.append(nxt)
    return len(found)


def run():
    for seed, convention in product((0, 1), repeat=2):
        assert self_explore(seed, convention) == [0, True, True, 1]
    for convention, x in product((0, 1), repeat=2):
        assert read((x, x), convention) == read((x, 1-x), convention)
        assert probe((x, x), convention, 1)[1] != probe((x, 1-x), convention, 1)[1]
    for state, convention in product(HISTORIES, (0, 1)):
        # Wait and exchange generate this full orbit since two exchanges give identity.
        orbit = {state, step(state, "exchange")}
        can_reverse = any(read(t,convention) != read(state,convention) for t in orbit)
        assert can_reverse == bool(state[0] ^ state[1])
    equivariant = leaves_diagonal = 0
    for outputs in product(HISTORIES, repeat=4):
        transform = dict(zip(HISTORIES, outputs))
        if all(transform[(y, x)] == transform[(x, y)][::-1] for x,y in HISTORIES):
            equivariant += 1
            leaves_diagonal += any(transform[(i, i)][0] != transform[(i, i)][1] for i in (0, 1))
    assert equivariant == 16 and leaves_diagonal == 0
    counts = [[n, reachable_ledgers(n, True), reachable_ledgers(n, False)] for n in range(1, 6)]
    assert all((a,b) == (2**n,2) for n,a,b in counts)
    return {"symmetric_seed_runs": 4, "example_initial_final_ledger": [0,1],
            "equivariant_updates_of_256": equivariant, "equivariant_updates_leaving_diagonal": leaves_diagonal,
            "future_record_cases_through_length_6": check_future_records(6),
            "readout_ambiguity_words_through_length_6": check_attribution(6),
            "calibration_preparations_verified": check_reference_preparation(),
            "copies_independently_addressed_or_collective": counts}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    data = run()
    path = HERE/"results.json"
    if parser.parse_args().check:
        assert data == json.loads(path.read_text(encoding="utf-8"))
        print("PASS: closed controller, necessity, calibration ambiguity, compositional reachability")
    else:
        with path.open("x", encoding="utf-8") as output:
            output.write(json.dumps(data, indent=2)+"\n")
