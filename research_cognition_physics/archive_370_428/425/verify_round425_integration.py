"""Integrate round 425 without changing historical scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round424_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round425_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_425_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 424
    assert checked['fresh_tests'] == dict(run=8, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'local_halving_coordinate_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (8, 0, 0)
    for key in ('three_dimensional_space_unconditionally_derived', 'full_GR_goal_completed',
                'full_cognitive_countermodel_completed', 'phase_closure_triggered',
                'additive_displacement_coordinates_inherited',
                'full_round424_contract_claimed_logically_weakened',
                'witness_claimed_complete_cognitive_realization'):
        assert not saved['scope'][key], key
    assert saved['scope']['local_halving_implies_NSS']
    note = HERE / 'research_note_425.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 16
    for name in ('local_halving_coordinate_audit.py', 'verify_local_halving_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_424') == 583
    assert result.pop('unchanged_prior_science_hashes_231_423') == 580
    assert result['stage_saved_tests'] == 1969
    assert result['total_protected_evidence_hashes'] == 618
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)',
                       (HERE.parent / 'research_direction.md').read_text(encoding='utf-8'))
    assert latest and int(latest[1]) >= 425 and int(latest[2]) >= 1977
    result.update(date='2026-09-24', rounds=[425],
        execution_mode='one complete local-coordinate bridge round with independent final review',
        scientific_base_through_round_by_round={425: 424},
        fresh_tests_by_round={425: 8}, fresh_tests=8,
        stage_saved_tests=1977, science_hashes_verified_231_425=586,
        unchanged_prior_science_hashes_231_424=583, total_protected_evidence_hashes=621,
        unchanged_prior_evidence_hashes=618, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=16,
        local_links_checked=result['local_links_checked'] + links,
        navigation_and_notes_checked=result['navigation_and_notes_checked'] + 1,
        consistent_local_halving_to_NSS_bridge_verified=True,
        alternative_local_Lie_route_allows_nonabelian_groups=True,
        full_previous_contract_claimed_strictly_weaker=False,
        actual_endpoint_and_local_halving_source_still_open=True,
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
