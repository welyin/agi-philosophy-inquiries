"""Integrate round 445 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round444_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round445_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_445_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 444
    assert checked['fresh_tests'] == dict(run=6, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'reciprocal_exchange_dynamics_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    for key in ('consistency_energy_and_address_access_derived_from_429',
                'existing_path_locality_selected_by_this_consistency_model',
                'arbitrary_unknown_matching_state_has_uniform_change_probability',
                'data_conditional_transport_implemented_in_this_model',
                'arbitrary_inconsistent_input_autonomously_repaired',
                'bounded_capacity_implies_uniform_interaction_budget',
                'three_dimensional_space_generated',
                'full_GR_goal_completed', 'phase_closure_triggered'):
        assert not saved['scope'][key], key
    for key in ('explicit_two_subject_continuous_address_model',
                'all_even_N_second_order_matching_flip_formula_proved',
                'finite_window_bare_code_error_and_probability_certificate',
                'arbitrary_unknown_matching_coherence_reference_error_covered',
                'all_time_initial_code_leakage_bound_proved',
                'known_perfect_matching_flip_graph_identified',
                'direct_two_body_virtual_process_evades_exact_code_support_obstruction'):
        assert saved['scope'][key], key
    note = HERE / 'research_note_445.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 12
    for name in ('reciprocal_exchange_dynamics_audit.py', 'verify_reciprocal_exchange_dynamics_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_444') == 643
    assert result.pop('unchanged_prior_science_hashes_231_443') == 640
    assert result['stage_saved_tests'] == 2095
    assert result['total_protected_evidence_hashes'] == 678
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 445 and int(latest[2]) >= 2101
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 73. 第445轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[445],
        execution_mode='autonomous two-subject address exchange with reciprocal consistency energy, effective rematching and finite-window certificate',
        scientific_base_through_round_by_round={445: 444}, fresh_tests_by_round={445: 6}, fresh_tests=6,
        stage_saved_tests=2101, science_hashes_verified_231_445=646,
        unchanged_prior_science_hashes_231_444=643, total_protected_evidence_hashes=681,
        unchanged_prior_evidence_hashes=678, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=12,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        all_even_N_second_order_matching_coefficients_proved=True,
        exact_fraction_duhamel_and_window_certificate_verified=True,
        bare_unknown_code_and_reference_scope_verified=True,
        all_time_initial_code_leakage_bound_proved=True,
        finite_N4_actual_rematching_probability_above_39_over_50=True,
        old_path_selection_failure_of_this_candidate_verified=True,
        extra_consistency_energy_and_potential_access_explicit=True,
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
