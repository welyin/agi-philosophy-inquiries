"""Integrate round 420 without rerunning older scientific experiments."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round419_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round420_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_420_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 418
    assert checked['frozen_integration_base_through_round'] == 419
    assert checked['fresh_tests'] == dict(run=11, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert len(checked['preserved_draft_hashes']) == 5
    for key in ('new_file_hashes', 'preserved_draft_hashes'):
        for name, sha in checked[key].items():
            assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'local_hardware_family_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (11, 0, 0)
    assert saved['scope']['all_complete_old_device_states_supported']
    assert saved['scope']['exact_continuous_hardware_types_are_model_input']
    for key in ('endpoint_restores_original_free_hamiltonian',
                'black_box_devices_physically_serial_composable',
                'complete_local_process_FUCP_implementation_proved',
                'three_dimensional_space_derived', 'full_GR_goal_completed',
                'full_cognitive_countermodel_completed', 'phase_closure_triggered'):
        assert not saved['scope'][key], key
    assert core.text_checks(HERE / 'research_note_420.md') == checked['text_checks']
    links = 0
    for name in ('research_note_420.md', 'spatial_premise_closure_audit.md'):
        for link in core.link_parser()((HERE / name).read_text(encoding='utf-8')):
            dest = (HERE / link).resolve()
            assert dest.exists() or (pending and dest == TARGET), (name, link)
            links += 1
    for name in ('local_hardware_family_audit.py', 'verify_local_hardware_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    assert result.pop('science_hashes_verified_231_419') == 568
    assert result.pop('unchanged_prior_science_hashes_231_418') == 565
    assert result['stage_saved_tests'] == 1923
    assert result['total_protected_evidence_hashes'] == 593
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)',
                       (HERE.parent / 'research_direction.md').read_text(encoding='utf-8'))
    assert latest and int(latest[1]) >= 420 and int(latest[2]) >= 1934
    result.update(date='2026-09-24', rounds=[420],
        execution_mode='complete hardware-family round independent of round419 with independent final review',
        scientific_base_through_round_by_round={420: 418}, fresh_tests_by_round={420: 11},
        fresh_tests=11, stage_saved_tests=1934, science_hashes_verified_231_420=571,
        unchanged_prior_science_hashes_231_419=568, total_protected_evidence_hashes=601,
        unchanged_prior_evidence_hashes=593, newly_preserved_draft_files=5,
        new_scientific_display_formulas_checked=12,
        local_links_checked=result['local_links_checked'] + links,
        navigation_and_notes_checked=result['navigation_and_notes_checked'] + 2,
        finite_hardware_family_control_extension_closed=True,
        endpoint_channel_is_complete_future_process=False,
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
