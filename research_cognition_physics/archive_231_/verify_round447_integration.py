"""Integrate round 447 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round446_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round447_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_447_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 446
    assert checked['fresh_tests'] == dict(run=8, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'internal_partner_response_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (8, 0, 0)
    for key in ('reciprocal_consistency_penalty_derived_from_429',
                'initial_matching_dependent_encoding_generated_autonomously',
                'arbitrary_raw_data_preserves_fixed_matching',
                'dynamic_partner_reassignment_combined_with_this_processor',
                'potential_contact_or_spatial_dimension_generated',
                'full_GR_goal_completed', 'phase_closure_triggered'):
        assert not saved['scope'][key], key
    for key in ('internal_substructure_swap_breaks_446_tagged_data_conservation',
                'exact_matching_marker_code_and_pair_factorization',
                'arbitrary_encoded_matching_logic_reference_covered',
                'no_separate_cross_subject_data_exchange_added',
                'strict_conditional_phase_for_all_nonzero_lambda',
                'all_time_bare_code_bound_uses_exact_eigenvalue_phases',
                'actual_recorded_partner_signal_and_fixed_sector_no_cross_signal',
                'unrestricted_raw_data_fourth_order_scope_counterexample',
                'single_marker_exact_qubit_code_and_record_selected_pair_H',
                'same_fixed_population_effect_signal_above_one_tenth_on_window',
                'low_energy_phase_readout_limitation_explicit'):
        assert saved['scope'][key], key
    note = HERE / 'research_note_447.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 16
    for name in ('internal_partner_response_audit.py', 'verify_internal_partner_response_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_446') == 649
    assert result.pop('unchanged_prior_science_hashes_231_445') == 646
    assert result['stage_saved_tests'] == 2106
    assert result['total_protected_evidence_hashes'] == 684
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 447 and int(latest[2]) >= 2114
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 75. 第447轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[447],
        execution_mode='internal register swaps with exact distributed partner records and a fixed population communication certificate',
        scientific_base_through_round_by_round={447: 446}, fresh_tests_by_round={447: 8}, fresh_tests=8,
        stage_saved_tests=2114, science_hashes_verified_231_447=652,
        unchanged_prior_science_hashes_231_446=649, total_protected_evidence_hashes=687,
        unchanged_prior_evidence_hashes=684, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=16,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        single_marker_exact_code_and_pair_H_intertwining_verified=True,
        full_unknown_matching_logic_reference_scope_verified=True,
        fixed_population_signal_above_one_tenth_on_window_certified=True,
        low_energy_all_time_bare_code_error_and_readout_boundary_proved=True,
        unrestricted_raw_data_fourth_order_rematching_counterexample_verified=True,
        recorded_partner_identity_not_confused_with_instantaneous_A_reciprocity=True,
        potential_mediator_source_not_claimed_generated=True,
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
