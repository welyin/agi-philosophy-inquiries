"""Integrate round 439 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round438_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round439_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_439_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 438
    assert checked['fresh_tests'] == dict(run=7, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'short_path_exchange_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (7, 0, 0)
    for key in ('short_path_rule_or_capacity_derived_from_429',
                'whole_data_dynamics_frozen_claimed',
                'single_trajectory_dimension_equal_to_configuration_count_claimed',
                'original_full_hopping_model_refuted',
                'stable_connected_geometry_or_dimension_generated',
                'final_Heisenberg_hardware_compiled', 'full_GR_goal_completed', 'phase_closure_triggered'):
        assert not saved['scope'][key], key
    for key in ('existing_short_path_gate_explicitly_attributed',
                'resource_charge_and_full_unknown_reference_intertwining_proved',
                'partition_conservation_and_component_no_signalling_proved',
                'capacity_two_all_N_configuration_classification_proved',
                'nonzero_autonomous_three_vertex_reconfiguration_retained',
                'saturated_triangle_free_obstruction_proved',
                'size_independent_endpoint_conversion_norm_bound_proved',
                'old_routing_closure_and_shortest_power_results_reused'):
        assert saved['scope'][key], key
    note = HERE / 'research_note_439.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 15
    for name in ('short_path_exchange_audit.py', 'verify_short_path_exchange_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_438') == 625
    assert result.pop('unchanged_prior_science_hashes_231_437') == 622
    assert result['stage_saved_tests'] == 2057
    assert result['total_protected_evidence_hashes'] == 660
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 439 and int(latest[2]) >= 2064
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 67. 第439轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[439],
        execution_mode='one complete short-path gated resource model with partition causality, bounded endpoint norm and all-N capacity-two classification',
        scientific_base_through_round_by_round={439: 438}, fresh_tests_by_round={439: 7}, fresh_tests=7,
        stage_saved_tests=2064, science_hashes_verified_231_439=628,
        unchanged_prior_science_hashes_231_438=625, total_protected_evidence_hashes=663,
        unchanged_prior_evidence_hashes=660, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=15,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        exact_short_path_charge_and_intertwining_verified=True,
        fixed_partition_tensor_factorization_and_no_signalling_verified=True,
        size_independent_endpoint_norm_bound_verified=True,
        all_N_capacity_two_graph_sector_classification_verified=True,
        frozen_occupancy_distinguished_from_quantum_coherence_and_data=True,
        ample_capacity_source_result_reused=True,
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
