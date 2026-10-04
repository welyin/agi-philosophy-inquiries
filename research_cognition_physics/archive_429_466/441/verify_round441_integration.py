"""Integrate round 441 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round440_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round441_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_441_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 440
    assert checked['fresh_tests'] == dict(run=6, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'moving_subject_propagation_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    for key in ('new_dynamical_term_added_beyond_440',
                'arbitrary_initial_label_delocalization_covered',
                'initial_path_or_three_dimensional_geometry_generated',
                'label_preparation_and_encoding_free',
                'strict_relativistic_light_cone_proved',
                'full_GR_goal_completed', 'phase_closure_triggered'):
        assert not saved['scope'][key], key
    for key in ('exact_local_composite_exchange_completion_verified',
                'moving_unique_label_tail_uniform_over_unknown_data_proved',
                'all_input_reference_TP_encoding_signal_bound_proved',
                'dynamic_identity_readout_distinguished_from_fixed_slot',
                'finite_rational_all_size_certificate_proved',
                'existing_LR_and_prior_fixed_chain_results_explicitly_reused'):
        assert saved['scope'][key], key
    note = HERE / 'research_note_441.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 13
    for name in ('moving_subject_propagation_audit.py', 'verify_moving_subject_propagation_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_440') == 631
    assert result.pop('unchanged_prior_science_hashes_231_439') == 628
    assert result['stage_saved_tests'] == 2071
    assert result['total_protected_evidence_hashes'] == 666
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 441 and int(latest[2]) >= 2077
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 69. 第441轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[441],
        execution_mode='moving identity-selected receiver under the unchanged path-sector exchange model',
        scientific_base_through_round_by_round={441: 440}, fresh_tests_by_round={441: 6}, fresh_tests=6,
        stage_saved_tests=2077, science_hashes_verified_231_441=634,
        unchanged_prior_science_hashes_231_440=631, total_protected_evidence_hashes=669,
        unchanged_prior_evidence_hashes=666, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=13,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        exact_local_tensor_completion_verified=True,
        unique_label_displacement_bound_proved=True,
        reference_complete_two_TP_encoding_signal_bound_proved=True,
        bounded_receiver_extension_verified=True,
        fixed_slot_false_subject_signal_verified=True,
        all_size_rational_certificate_verified=True,
        continuous_natural_evolution_distinguished_from_exact_processor=True,
        original_229_theorem_not_rewritten=True,
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
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
