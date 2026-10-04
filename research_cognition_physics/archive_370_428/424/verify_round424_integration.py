"""Integrate round 424, retaining all historical science and checks."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round423_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round424_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_424_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 423
    assert checked['fresh_tests'] == dict(run=9, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'single_contraction_coordinate_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (9, 0, 0)
    for key in ('three_dimensional_space_unconditionally_derived', 'full_GR_goal_completed',
                'full_cognitive_countermodel_completed', 'phase_closure_triggered',
                'old_budget_power_conclusion_retained', 'all_contractions_claimed_physically_implemented'):
        assert not saved['scope'][key], key
    assert saved['scope']['single_contracting_endpoint_automorphism_suffices']
    note = HERE / 'research_note_424.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 14
    for name in ('single_contraction_coordinate_audit.py', 'verify_single_contraction_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_423') == 580
    assert result.pop('unchanged_prior_science_hashes_231_422') == 577
    assert result['stage_saved_tests'] == 1960
    assert result['total_protected_evidence_hashes'] == 615
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)',
                       (HERE.parent / 'research_direction.md').read_text(encoding='utf-8'))
    assert latest and int(latest[1]) >= 424 and int(latest[2]) >= 1969
    result.update(date='2026-09-24', rounds=[424],
        execution_mode='one complete premise-reduction round with independent final review',
        scientific_base_through_round_by_round={424: 423},
        fresh_tests_by_round={424: 9}, fresh_tests=9,
        stage_saved_tests=1969, science_hashes_verified_231_424=583,
        unchanged_prior_science_hashes_231_423=580, total_protected_evidence_hashes=618,
        unchanged_prior_evidence_hashes=615, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=14,
        local_links_checked=result['local_links_checked'] + links,
        navigation_and_notes_checked=result['navigation_and_notes_checked'] + 1,
        all_scale_action_removed_from_conditional_three_dimensional_route=True,
        exact_budget_homogeneity_removed_from_conditional_coordinate_route=True,
        actual_single_contraction_source_still_open=True,
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
