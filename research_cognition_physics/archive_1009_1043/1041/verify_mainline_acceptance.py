"""Verify appended provenance and dependency composition, without rewriting history.

Scientific replay is provided by verify_round1041.py. This entry checks the
frozen author files, separate reviews, prior acceptance and the dependency delta.
Only --write creates the first receipt, exclusively; the default is read-only.
"""
from pathlib import Path
import argparse
import json
import runpy

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
BASE = STAGE.parent
OUT = HERE / 'mainline_acceptance.json'
PRIOR = STAGE / '1040/parallel_acceptance_1039_1040.json'
DELTA = HERE / 'dependency_delta_1041.json'
AUTHOR = HERE / 'research_round_1041_checks.json'
REVIEWS = [HERE / 'independent_review_by_dynamics_agent.md',
           HERE / 'operational_scope_review_by_species.md']
ADMISSIONS = [STAGE / '1040' / name for name in (
    'next_generation_admission.md', 'next_matter_admission.md',
    'next_dynamics_admission.md')]
NAV = [BASE / 'README.md', BASE / 'research_direction.md',
       BASE / 'RESEARCH_STATE.md', STAGE / 'README.md', STAGE / '文件索引.md']
OLD = runpy.run_path(str(PRIOR.parent / 'verify_parallel_acceptance.py'))
sha, read, relative, check_map = [OLD[n] for n in ('sha', 'read', 'relative', 'check_map')]


def evidence():
    assert sha(PRIOR) == '75c72c3142ac05c63393bcfd1558213fdf31fee3d4860a0ea9eb90fc4048c336'
    previous = read(PRIOR)
    comparison = dict(previous)
    comparison.pop('navigation_links_at_acceptance')
    assert comparison == OLD['evidence'](), 'Prior frozen evidence changed'
    assert previous['cumulative_after'] == 3817

    author = read(AUTHOR)
    assert author['all_checks_passed'] and author['science_read_only_compare_passed']
    assert author['same_round_science_groups'] == 1
    assert author['additional_science_groups'] == 0
    science = check_map(HERE, author['owned_sha256'])
    history = check_map(BASE, author['historical_sha256'])
    assert len(science) == 8 and len(history) == 17

    delta = read(DELTA)
    base_delta_path = BASE / delta['base_delta']
    assert sha(base_delta_path) == delta['base_delta_sha256']
    base_delta = read(base_delta_path)
    base_ledger_path = BASE / base_delta['base_ledger']
    assert sha(base_ledger_path) == base_delta['base_sha256']
    inherited = read(base_ledger_path)['rows'] + base_delta['rows_append']
    assert len(inherited) == delta['base_resolved_rows'] == 32
    appended = delta['rows_append']
    assert len(appended) == 1 and appended[0]['source_rounds'] == [233, 1041]
    rows = inherited + appended
    assert len(rows) == delta['resolved_row_count'] == 33
    assert len({r['id'] for r in rows}) == len(rows)
    for row in rows:
        for source in row['source_paths']:
            assert (BASE / source).exists(), source
    for row in appended:
        assert row['requires'] and row['not_proved'] and row['object']
        for flag in ('cognitive_origin_closed', 'full_parent_process_certified',
                     'physical_implementation_certified'):
            assert row[flag] is False
    assert delta['new_adopted_cognitive_axioms'] == 0 and not delta['goal_completed']

    for review in REVIEWS:
        body = review.read_text('utf8')
        assert '1041' in body and '通过' in body
        assert author['owned_sha256']['../research_note_1041.md'] in body
    supplements = [Path(__file__), DELTA, HERE / 'mainline_integration.md',
                   *REVIEWS, *ADMISSIONS]
    return dict(
        schema='append_only_mainline_acceptance_v1', date='2026-10-08',
        rounds_accepted=[1041], cumulative_before=3817, cumulative_after=3818,
        new_scientific_calibration_groups=1, integration_calibration_groups=0,
        new_adopted_cognitive_axioms=0,
        prior_acceptance=relative(PRIOR), prior_acceptance_sha256=sha(PRIOR),
        author_receipt=relative(AUTHOR), author_receipt_sha256=sha(AUTHOR),
        independent_agent_reviews=[relative(p) for p in REVIEWS],
        historical_author_receipts_preserved=True, human_peer_review=False,
        dependency_delta=relative(DELTA), dependency_base_rows=32,
        dependency_added_rows=1, dependency_resolved_rows=33,
        full_common_parent_model_established=False,
        physical_implementation_certified=False, goal_completed=False,
        additional_frozen_science_sha256=science,
        additional_historical_sha256=history,
        supplemental_sha256={relative(p): sha(p) for p in supplements})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = evidence()
    # Reuse the established link checker, with this receipt as the sole allowed
    # prospective target on its first write. Navigation remains live, not frozen.
    checker = OLD['check_links']
    checker.__globals__['OUT'] = OUT
    links = sum(checker(p, args.write) for p in NAV)
    for p in [HERE / 'mainline_integration.md', *REVIEWS, *ADMISSIONS]:
        checker(p, args.write)
    if args.write:
        with OUT.open('x', encoding='utf8') as stream:
            json.dump(dict(result, navigation_links_at_acceptance=links), stream,
                      ensure_ascii=False, indent=2, allow_nan=False)
            stream.write('\n')
    else:
        stored = read(OUT)
        stored.pop('navigation_links_at_acceptance')
        assert stored == result, 'Frozen mainline evidence changed'
    print(json.dumps(dict(passed=True, accepted=[1041], cumulative=3818,
                          dependency_rows=33, live_navigation_links=links,
                          new_author_assets=len(result['additional_frozen_science_sha256']),
                          separate_agent_reviews=len(REVIEWS)), ensure_ascii=False))
