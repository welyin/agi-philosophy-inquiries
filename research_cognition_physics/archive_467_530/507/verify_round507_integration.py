"""Read-only integration audit for round 507, preserving the complete baseline."""
import argparse
import ast
import json
from pathlib import Path
import re
import verify_round506_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE/'round507_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE/'research_round_507_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 506
    assert checked['fresh_tests'] == dict(run=6, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert len(checked['preserved_draft_hashes']) == 2
    for group in ['new_file_hashes', 'preserved_draft_hashes']:
        for name, sha in checked[group].items():
            assert core.digest(HERE/name) == sha, name
    saved = core.read(HERE/'distributed_role_reader_results.json')
    assert saved['scope'] == checked['scope']
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    note = HERE/'research_note_507.md'
    assert core.text_checks(note) == checked['text_checks']
    assert note.read_bytes() == (HERE/'round507_drafts/research_note_507.txt').read_bytes()
    for name in ['distributed_role_reader.py', 'verify_distributed_role_reader_round.py', Path(__file__).name]:
        ast.parse((HERE/name).read_text(encoding='utf-8'))
    links = 0
    for doc in [note, HERE/'three_dimensional_four_conditions_review.md']:
        for link in core.link_parser()(doc.read_text(encoding='utf-8')):
            dest = (doc.parent/link).resolve()
            assert dest.exists() or (pending and dest == TARGET), link
            links += 1
    assert result.pop('science_hashes_verified_231_506') == 829
    assert result.pop('unchanged_prior_science_hashes_231_505') == 826
    assert result['stage_saved_tests'] == 2468
    assert result['total_protected_evidence_hashes'] == 912
    direction = (HERE.parent/'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 507 and int(latest[2]) >= 2474
    assert '## 143. 第507轮' in (HERE/'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    assert '## 50. 507' in (HERE/'three_dimensional_four_conditions_review.md').read_text(encoding='utf-8')
    result.update(date='2026-09-27', rounds=[507],
        execution_mode='distributed finite-time cat reader with scale-uniform complete instrument error; baseline 506',
        scientific_base_through_round_by_round={507:506}, fresh_tests_by_round={507:6},
        fresh_tests=6, stage_saved_tests=2474,
        science_hashes_verified_231_507=832, unchanged_prior_science_hashes_231_506=829,
        total_protected_evidence_hashes=917, unchanged_prior_evidence_hashes=912,
        newly_preserved_draft_files=2, newly_protected_unnumbered_files=0,
        new_scientific_display_formulas_checked=checked['text_checks']['display_formulas'],
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        latest_round_scope=saved['scope'], cognitive_axioms_silently_strengthened=False,
        physical_positions_or_dimension_generated=False, full_cognition_to_gr_refuted=False,
        phase_closure_triggered=False, old_scientific_experiments_rerun=False,
        independent_final_agent_review_completed=True, independent_analytic_scope_review_completed=True,
        parent_analytic_and_reproduction_review_completed=True, all_reported_checks_passed=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ['rounds', 'stage_saved_tests',
        'total_protected_evidence_hashes', 'all_reported_checks_passed']}, ensure_ascii=False))
