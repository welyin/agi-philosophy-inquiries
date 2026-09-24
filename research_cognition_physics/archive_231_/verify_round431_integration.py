"""Integrate round 431 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round430_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round431_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_431_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 430
    assert checked['fresh_tests'] == dict(run=6, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'exchange_reader_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    for key in ('exact_programmable_processor_required', 'permanent_record_proved',
                'complete_readout_amplifier_constructed', 'physical_clock_derived',
                'full_cognitive_implementation_completed', 'three_dimensional_space_unconditionally_derived',
                'full_GR_goal_completed', 'phase_closure_triggered'):
        assert not saved['scope'][key], key
    for key in ('internal_pair_record_transfer_constructed', 'exact_integer_and_rational_signal_certificate',
                'finite_read_window_proved', 'same_continuous_exchange_rule_only',
                'unknown_full_joint_reference_information_preserved', 'singlet_source_and_coupling_pattern_are_inputs'):
        assert saved['scope'][key], key
    note = HERE / 'research_note_431.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 10
    for name in ('exchange_reader_audit.py', 'verify_exchange_reader_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_430') == 601
    assert result.pop('unchanged_prior_science_hashes_231_429') == 598
    assert result['stage_saved_tests'] == 2008
    assert result['total_protected_evidence_hashes'] == 636
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 431 and int(latest[2]) >= 2014
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 59. 第431轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[431],
        execution_mode='one complete internal-reader round with exact rational certificate and independent spectral verification',
        scientific_base_through_round_by_round={431: 430}, fresh_tests_by_round={431: 6}, fresh_tests=6,
        stage_saved_tests=2014, science_hashes_verified_231_431=604,
        unchanged_prior_science_hashes_231_430=601, total_protected_evidence_hashes=639,
        unchanged_prior_evidence_hashes=636, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=10,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        three_body_exchange_reader_block_verified=True,
        finite_read_window_certified_by_exact_series_and_analytic_bound=True,
        reader_pair_distinction_without_single_qubit_signal_verified=True,
        continuous_natural_evolution_distinguished_from_exact_processor=True,
        original_229_theorem_not_rewritten=True,
        physical_positions_or_dimension_generated=False, full_cognition_to_gr_refuted=False,
        phase_closure_triggered=False, old_scientific_experiments_rerun=False,
        independent_final_agent_review_completed=False, all_reported_checks_passed=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
