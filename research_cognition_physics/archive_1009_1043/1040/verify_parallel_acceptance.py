"""Read-only receipt verification; --write creates the first receipt exclusively.

Checks provenance and dependency composition, not a substitute for scientific
proof review or each round's separate numerical verifier. Old receipts remain
immutable; the prior batch is referenced rather than duplicated.
"""
from pathlib import Path
from urllib.parse import unquote
import argparse
import hashlib
import json
import re
import runpy

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
BASE = STAGE.parent
ROOT = BASE.parent
OUT = HERE / 'parallel_acceptance_1039_1040.json'
PRIOR = STAGE / '1038/parallel_acceptance_1036_1038.json'
DELTA = HERE / 'dependency_delta_1036_1040.json'
REVIEWS = {
    n: STAGE / str(n) / 'independent_review_by_dynamics_agent.md'
    for n in (1032, 1033, 1034, 1039)
}
REVIEWS[1040] = HERE / 'independent_review_by_species_agent.md'
NAV = [BASE / 'README.md', BASE / 'research_direction.md',
       BASE / 'RESEARCH_STATE.md', STAGE / 'README.md', STAGE / '文件索引.md']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text('utf8'))


def relative(path):
    return path.resolve().relative_to(BASE).as_posix()


def check_map(parent, mapping):
    result = {}
    for rel, digest in mapping.items():
        p = (parent / rel).resolve()
        assert p.is_relative_to(BASE), str(p)
        assert sha(p) == digest, str(p)
        result[relative(p)] = digest
    return result


def check_links(path, prospective=False):
    body = re.sub(r'```.*?```|\$\$.*?\$\$', '', path.read_text('utf8'), flags=re.S)
    total = 0
    for raw in re.findall(r'(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)', body):
        target = unquote(raw.strip().strip('<>').split('#', 1)[0])
        if not target or re.match(r'^[a-zA-Z][a-zA-Z0-9+.-]*:', target):
            continue
        p = (path.parent / target).resolve()
        if prospective and p == OUT:
            continue
        assert p.exists(), (str(path), target)
        total += 1
    return total


def evidence():
    previous = read(PRIOR)
    frozen_previous = dict(previous)
    frozen_previous.pop('navigation_links_at_acceptance')
    previous_module = runpy.run_path(str(PRIOR.parent / 'verify_parallel_acceptance.py'))
    assert frozen_previous == previous_module['evidence'](), 'Prior acceptance changed'
    assert previous['cumulative_after'] == 3815

    science, history, receipts = {}, {}, {}
    for n in (1032, 1039, 1040):
        receipt = STAGE / str(n) / f'research_round_{n}_checks.json'
        d = read(receipt)
        receipts[relative(receipt)] = sha(receipt)
        if n == 1032:
            assert d['local_delivery_checks_passed']
            science.update(check_map(ROOT, d['source_sha256']))
            history.update(check_map(BASE, d['historical_source_sha256']))
        elif n == 1039:
            assert d['passed'] and d['scientific_groups_proposed'] == 1
            science.update(check_map(STAGE, d['artifact_sha256']))
            history.update(check_map(BASE, d['input_sha256']))
        else:
            assert d['all_checks_passed'] and d['same_round_science_calibration_groups'] == 1
            science.update(check_map(HERE, d['owned_sha256']))
            history.update(check_map(BASE, d['historical_sha256']))

    delta = read(DELTA)
    base_path = BASE / delta['base_ledger']
    assert sha(base_path) == delta['base_sha256']
    base = read(base_path)
    assert len(base['rows']) == delta['base_row_count'] == 27
    appended = delta['rows_append']
    rows = base['rows'] + appended
    assert len(rows) == delta['resolved_row_count'] == 32
    assert len({row['id'] for row in rows}) == len(rows)
    assert [r['source_rounds'] for r in appended] == [[n] for n in range(1036, 1041)]
    for row in rows:
        for rel in row['source_paths']:
            assert (BASE / rel).exists(), rel
    for row in appended:
        assert row['requires'] and row['not_proved'] and row['object']
        for flag in ('cognitive_origin_closed', 'full_parent_process_certified',
                     'physical_implementation_certified'):
            assert row[flag] is False
    assert not delta['goal_completed'] and not delta['full_common_parent_model_established']
    assert delta['new_adopted_cognitive_axioms'] == 0

    supplements = [Path(__file__), DELTA, HERE / 'parallel_integration_1039_1040.md',
                   STAGE / '1038/rf2_selection_admission.md',
                   HERE / 'review_1040_by_species.py',
                   HERE / 'review_1040_by_species_results.json']
    for n, path in REVIEWS.items():
        body = path.read_text('utf8')
        assert str(n) in body and '通过' in body, str(path)
        supplements.append(path)
    all_science = set(previous['frozen_science_sha256']) | set(science)
    all_history = set(previous['historical_input_sha256']) | set(history)
    return dict(
        schema='append_only_parallel_acceptance_v2', date='2026-10-08',
        rounds_accepted=[1039, 1040], cumulative_before=3815,
        cumulative_by_round={'1039': 3816, '1040': 3817}, cumulative_after=3817,
        new_scientific_calibration_groups=2, integration_calibration_groups=0,
        new_adopted_cognitive_axioms=0,
        prior_acceptance=relative(PRIOR), prior_acceptance_sha256=sha(PRIOR),
        independent_agent_reviews={str(n): relative(p) for n, p in REVIEWS.items()},
        review_backlog_closed=[1032, 1033, 1034],
        independent_review_pending_in_range_1032_1040=[],
        historical_author_receipts_preserved=True, human_peer_review=False,
        dependency_base=relative(base_path), dependency_base_sha256=sha(base_path),
        dependency_delta=relative(DELTA), dependency_base_rows=27,
        dependency_added_rows=5, dependency_resolved_rows=32,
        full_common_parent_model_established=False,
        physical_implementation_certified=False, goal_completed=False,
        science_file_count_including_prior=len(all_science),
        historical_input_count_including_prior=len(all_history),
        additional_frozen_science_sha256=science, additional_historical_sha256=history,
        additional_author_receipt_sha256=receipts,
        supplemental_sha256={relative(p): sha(p) for p in supplements})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = evidence()
    links = sum(check_links(p, args.write) for p in NAV)
    for p in [HERE / 'parallel_integration_1039_1040.md', *REVIEWS.values()]:
        check_links(p, args.write)
    if args.write:
        with OUT.open('x', encoding='utf8') as stream:
            json.dump(dict(result, navigation_links_at_acceptance=links), stream,
                      ensure_ascii=False, indent=2, allow_nan=False)
            stream.write('\n')
    else:
        prior = read(OUT)
        prior.pop('navigation_links_at_acceptance')
        assert prior == result, 'Frozen integration evidence changed'
    print(json.dumps(dict(passed=True, accepted=result['rounds_accepted'], cumulative=3817,
                          dependency_rows=32, live_navigation_links=links,
                          science_files=result['science_file_count_including_prior'],
                          historical_inputs=result['historical_input_count_including_prior']),
                     ensure_ascii=False))
