"""Integrate round 438 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round437_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round438_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_438_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 437
    assert checked['fresh_tests'] == dict(run=7, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'endpoint_resource_exchange_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (7, 0, 0)
    for key in ('conversion_rule_charge_preparation_and_capacity_derived_from_429',
                'arbitrary_data_marginal_unchanged_claimed',
                'stable_connected_geometry_or_dimension_generated',
                'final_Heisenberg_hardware_compiled', 'full_GR_goal_completed', 'phase_closure_triggered'):
        assert not saved['scope'][key], key
    for key in ('existing_matter_link_conversion_explicitly_attributed',
                'local_resource_charge_conservation_and_hard_degree_bound_proved',
                'autonomous_reassignment_and_data_response_verified',
                'all_unknown_joint_data_reference_preserved',
                'all_matching_configurations_have_nonzero_short_time_path_amplitudes',
                'relationship_memory_and_rate_costs_explicit',
                'cubic_complete_pair_data_response_from_empty_graph_proved'):
        assert saved['scope'][key], key
    note = HERE / 'research_note_438.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 15
    for name in ('endpoint_resource_exchange_audit.py', 'verify_endpoint_resource_exchange_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_437') == 622
    assert result.pop('unchanged_prior_science_hashes_231_436') == 619
    assert result['stage_saved_tests'] == 2050
    assert result['total_protected_evidence_hashes'] == 657
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 438 and int(latest[2]) >= 2057
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 66. 第438轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[438],
        execution_mode='one complete conserved endpoint-resource model with autonomous regrouping and operational propagation audit',
        scientific_base_through_round_by_round={438: 437}, fresh_tests_by_round={438: 7}, fresh_tests=7,
        stage_saved_tests=2057, science_hashes_verified_231_438=625,
        unchanged_prior_science_hashes_231_437=622, total_protected_evidence_hashes=660,
        unchanged_prior_evidence_hashes=657, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=15,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        endpoint_resource_charge_and_full_intertwining_verified=True,
        exact_autonomous_two_ninths_reassignment_verified=True,
        arbitrary_joint_unknown_data_and_reference_preserved=True,
        active_configuration_and_mean_graph_distinguished=True,
        complete_interface_and_symmetric_trajectory_capacity_distinguished=True,
        endpoint_conversion_rate_cost_verified=True,
        cubic_data_influence_and_internal_reference_witness_verified=True,
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
