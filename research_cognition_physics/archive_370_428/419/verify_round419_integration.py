"""Integrate round 419 without rerunning old scientific experiments."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round418_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round419_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE/'research_round_419_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 418
    assert checked['fresh_tests'] == dict(run=8,failures=0,errors=0)
    assert len(checked['new_file_hashes']) == 3
    for name,sha in checked['new_file_hashes'].items():
        assert core.digest(HERE/name) == sha,name
    saved = core.read(HERE/'internal_rate_window_audit_results.json')
    assert saved['checks'] == checked['fresh_tests']
    assert saved['deterministic_time_window_proved']
    assert saved['rate_preparation_and_distribution_assumed_as_counted_initial_resources']
    for key in ('random_read_time_choice_required','same_12_state_local_rule_unchanged',
                'fixed_local_dimension_for_all_windows','arbitrary_active_readout_schedule_proved',
                'online_new_subject_admission_proved','all_late_time_guarantee',
                'full_cognitive_countermodel_completed','physical_dimension_derived','phase_closure_triggered'):
        assert not saved[key],key
    assert core.text_checks(HERE/'research_note_419.md') == checked['text_checks']
    links = 0
    for name in ('research_note_419.md','spatial_premise_closure_audit.md'):
        for link in core.link_parser()((HERE/name).read_text(encoding='utf-8')):
            dest = (HERE/link).resolve()
            assert dest.exists() or (pending and dest == TARGET),(name,link)
            links += 1
    for name in ('internal_rate_window_audit.py','verify_internal_rate_round.py',Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding='utf-8'))
    assert result.pop('science_hashes_verified_231_418') == 565
    assert result.pop('unchanged_prior_science_hashes_231_417') == 562
    assert result['stage_saved_tests'] == 1915
    assert result['total_protected_evidence_hashes'] == 590
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)',
                       (HERE.parent/'research_direction.md').read_text(encoding='utf-8'))
    assert latest and int(latest[1]) >= 419 and int(latest[2]) >= 1923
    result.update(date='2026-09-24',rounds=[419],
        execution_mode='one complete internal-rate-window round with independent final review',
        scientific_base_through_round_by_round={419:418},fresh_tests_by_round={419:8},
        fresh_tests=8,stage_saved_tests=1923,science_hashes_verified_231_419=568,
        unchanged_prior_science_hashes_231_418=565,total_protected_evidence_hashes=593,
        unchanged_prior_evidence_hashes=590,new_scientific_display_formulas_checked=9,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+2,
        finite_internal_rate_window_proved=True,random_read_time_required=False,
        autonomous_stopping_or_online_admission_proved=False,
        physical_positions_or_dimension_generated=False,full_cognition_to_gr_refuted=False,
        phase_closure_triggered=False,old_scientific_experiments_rerun=False,
        independent_final_agent_review_completed=True,all_reported_checks_passed=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args = parser.parse_args()
    result = verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
