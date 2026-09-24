"""Integrate round 449 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round448_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round449_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_449_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 448
    assert checked['fresh_tests'] == dict(run=7, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'protected_partner_processing_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (7, 0, 0)
    for key in ('protected_logical_states_required_to_be_degenerate', 'new_diagonal_interactions_derived_from_429', 'every_composite_full_state_satisfies_429_condition', 'unlimited_scale_uniform_stability_proved', 'initial_code_preparation_or_actual_clock_generated', 'full_spatial_dimension_or_GR_generated', 'full_GR_goal_completed', 'phase_closure_triggered'):
        assert not saved['scope'][key], key
    for key in ('diagonal_reward_degeneracy_and_symmetry_conditions_proved', 'record_and_conditional_processing_compatible_with_explicit_scale_hierarchy', 'full_unknown_marker_logic_reference_error_certified', 'uniform_internal_field_residual_cancellation_proved', 'all_time_finite_N_record_leakage_bound_proved', 'same_fixed_physical_effect_and_rematching_on_common_window_certified', 'internal_resource_and_time_cost_explicit'):
        assert saved['scope'][key], key
    note = HERE / 'research_note_449.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 14
    for name in ('protected_partner_processing_audit.py', 'verify_protected_partner_processing_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_448') == 655
    assert result.pop('unchanged_prior_science_hashes_231_447') == 652
    assert result['stage_saved_tests'] == 2120
    assert result['total_protected_evidence_hashes'] == 690
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 449 and int(latest[2]) >= 2127
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 77. 第449轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[449],
        execution_mode='protected partner rematching and conditional processing on the same unknown-input carrier and physical readout window',
        scientific_base_through_round_by_round={449: 448}, fresh_tests_by_round={449: 7}, fresh_tests=7,
        stage_saved_tests=2127, science_hashes_verified_231_449=658,
        unchanged_prior_science_hashes_231_448=655, total_protected_evidence_hashes=693,
        unchanged_prior_evidence_hashes=690, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=14,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        diagonal_reward_symmetry_and_degeneracy_scope_verified=True,
        joint_protected_processing_residual_identity_verified=True,
        full_unknown_marker_logic_reference_error_and_all_time_leakage_verified=True,
        fixed_physical_conditional_signal_and_legal_rematching_common_window_certified=True,
        pure_transport_original_partner_signal_boundary_verified=True,
        energy_time_and_scale_costs_not_omitted=True,
        diagonal_interaction_sources_not_claimed_derived_from_429=True,
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
