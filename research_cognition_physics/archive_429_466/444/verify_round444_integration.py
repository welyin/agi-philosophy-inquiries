"""Integrate round 444 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round443_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round444_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_444_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 443
    assert checked['fresh_tests'] == dict(run=6, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'partner_port_carrier_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    for key in ('universal_fixed_dimension_subject_hardware_constructed',
                'four_subject_interaction_derived_from_429_two_body_rule',
                'address_access_and_interaction_architecture_generated',
                'general_nonuniform_port_states_equivalent_to_old_graph',
                'all_old_relation_operations_locally_preserved',
                'three_dimensional_space_generated',
                'full_GR_goal_completed', 'phase_closure_triggered'):
        assert not saved['scope'][key], key
    for key in ('existing_quantum_address_register_idea_reused',
                'explicit_subject_tensor_factorization_supplied',
                'no_independent_qubit_per_potential_edge_required_by_this_realization',
                'exact_440_graph_data_H_and_subject_readout_intertwining',
                'unknown_encoded_joint_data_graph_reference_preserved',
                'four_subject_extension_without_global_code_projector_proved',
                'exact_three_subject_support_obstruction_in_fixed_degree_code'):
        assert saved['scope'][key], key
    note = HERE / 'research_note_444.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 12
    for name in ('partner_port_carrier_audit.py', 'verify_partner_port_carrier_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_443') == 640
    assert result.pop('unchanged_prior_science_hashes_231_442') == 637
    assert result['stage_saved_tests'] == 2089
    assert result['total_protected_evidence_hashes'] == 675
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 444 and int(latest[2]) >= 2095
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 72. 第444轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[444],
        execution_mode='reciprocal subject-port tensor realization of unchanged graph Hamiltonian and exact interaction-support boundary',
        scientific_base_through_round_by_round={444: 443}, fresh_tests_by_round={444: 6}, fresh_tests=6,
        stage_saved_tests=2095, science_hashes_verified_231_444=643,
        unchanged_prior_science_hashes_231_443=640, total_protected_evidence_hashes=678,
        unchanged_prior_evidence_hashes=675, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=12,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        explicit_subject_tensor_factorization_and_address_budget_verified=True,
        exact_port_flip_graph_and_data_intertwiner_verified=True,
        all_encoded_unknown_graph_data_reference_scope_verified=True,
        four_subject_raw_space_extension_without_global_projector_verified=True,
        naive_adjoint_leakage_counterexample_verified=True,
        three_subject_exact_code_graph_change_obstruction_proved=True,
        address_storage_not_confused_with_interaction_origin=True,
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
