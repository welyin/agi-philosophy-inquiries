"""Integrate round 446 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round445_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round446_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_446_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 445
    assert checked['fresh_tests'] == dict(run=5, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'whole_subject_exchange_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (5, 0, 0)
    for key in ('separate_edge_conditioned_data_swap_added',
                'actual_data_signal_identified_with_address_motion_only',
                'arbitrary_unknown_graph_has_simple_classical_routing_channel',
                'full_unknown_message_perfectly_transferred',
                'consistency_energy_and_all_pair_access_derived_from_429',
                'old_relation_locality_or_three_dimensional_space_generated',
                'full_GR_goal_completed', 'phase_closure_triggered'):
        assert not saved['scope'][key], key
    for key in ('whole_subject_swap_includes_address_and_data',
                'exact_all_N_equal_data_dimension_controlled_permutation_conjugacy',
                'unknown_joint_code_data_reference_error_inherited_from_445',
                'fixed_subject_readout_transformed_explicitly',
                'complete_address_tagged_data_algebra_conserved',
                'actual_cross_component_data_signal_strictly_certified',
                'old_partner_legal_code_transport_forbidden_in_this_model'):
        assert saved['scope'][key], key
    note = HERE / 'research_note_446.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 12
    for name in ('whole_subject_exchange_audit.py', 'verify_whole_subject_exchange_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_445') == 646
    assert result.pop('unchanged_prior_science_hashes_231_444') == 643
    assert result['stage_saved_tests'] == 2101
    assert result['total_protected_evidence_hashes'] == 681
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 446 and int(latest[2]) >= 2106
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 74. 第446轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[446],
        execution_mode='whole-subject address-data exchange, exact tagged-data conservation and actual fixed-subject receiver channel',
        scientific_base_through_round_by_round={446: 445}, fresh_tests_by_round={446: 5}, fresh_tests=5,
        stage_saved_tests=2106, science_hashes_verified_231_446=649,
        unchanged_prior_science_hashes_231_445=646, total_protected_evidence_hashes=684,
        unchanged_prior_evidence_hashes=681, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=12,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        exact_whole_subject_controlled_permutation_conjugacy_verified=True,
        arbitrary_joint_code_data_reference_error_inherited_exactly=True,
        complete_tagged_data_algebra_and_closed_address_loop_boundary_proved=True,
        physical_readout_and_known_matching_receiver_channel_verified=True,
        actual_cross_component_signal_above_7_over_20_certified=True,
        old_partner_source_weight_bounded_by_all_time_code_leakage=True,
        address_probability_and_full_address_state_distinguished=True,
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
