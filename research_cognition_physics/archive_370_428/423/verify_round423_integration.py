"""Integrate round 423 without rerunning prior science or changing frozen files."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round421_422_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round423_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_423_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 422
    assert checked['fresh_tests'] == dict(run=9, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert len(checked['preserved_draft_hashes']) == 2
    for key in ('new_file_hashes', 'preserved_draft_hashes'):
        for name, sha in checked[key].items():
            assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'finite_sampling_source_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (9, 0, 0)
    for key in ('three_dimensional_space_derived', 'full_GR_goal_completed',
                'full_cognitive_countermodel_completed', 'phase_closure_triggered',
                'restricted_bias_same_bound_claimed', 'original_round413_algorithm_invalidated'):
        assert not saved['scope'][key], key
    note = HERE / 'research_note_423.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 14
    for name in ('finite_sampling_source_audit.py', 'verify_finite_sampling_source_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_422') == 577
    assert result.pop('unchanged_prior_science_hashes_231_420') == 571
    assert result['stage_saved_tests'] == 1951
    assert result['total_protected_evidence_hashes'] == 610
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)',
                       (HERE.parent / 'research_direction.md').read_text(encoding='utf-8'))
    assert latest and int(latest[1]) >= 423 and int(latest[2]) >= 1960
    result.update(date='2026-09-24', rounds=[423],
        execution_mode='one complete source-capacity round with independent final review',
        scientific_base_through_round_by_round={423: 422},
        fresh_tests_by_round={423: 9}, fresh_tests=9,
        stage_saved_tests=1960, science_hashes_verified_231_423=580,
        unchanged_prior_science_hashes_231_422=577, total_protected_evidence_hashes=615,
        unchanged_prior_evidence_hashes=610, newly_preserved_draft_files=2,
        new_scientific_display_formulas_checked=14,
        local_links_checked=result['local_links_checked'] + links,
        navigation_and_notes_checked=result['navigation_and_notes_checked'] + 1,
        full_bias_and_local_bias_sharp_source_dimensions_proved=True,
        physical_positions_or_dimension_generated=False, full_cognition_to_gr_refuted=False,
        phase_closure_triggered=False, old_scientific_experiments_rerun=False,
        independent_final_agent_review_completed=True, all_reported_checks_passed=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
