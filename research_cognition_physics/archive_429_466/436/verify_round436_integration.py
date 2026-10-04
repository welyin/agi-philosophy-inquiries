"""Integrate round 436 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round435_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round436_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_436_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 435
    assert checked['fresh_tests'] == dict(run=6, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'static_exchange_network_bridge_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    for key in ('final_Heisenberg_coupling_list_or_hardware_run_completed',
                'exact_finite_resource_universal_implementation_claimed',
                'natural_architecture_or_actual_positions_generated',
                'three_dimensional_space_unconditionally_derived', 'full_GR_goal_completed', 'phase_closure_triggered'):
        assert not saved['scope'][key], key
    for key in ('full_finite_network_static_exchange_simulation_exists_by_cited_theorems',
                'explicit_first_no_Y_layer_checked_on_all_64_input_columns',
                'complete_unknown_reference_dynamical_contract_proved',
                'target_signal_and_simulation_budget_certified',
                'known_universal_Hamiltonian_theorems_not_claimed_as_new',
                'target_dependent_resources_and_nonuniform_weights_still_inputs'):
        assert saved['scope'][key], key
    note = HERE / 'research_note_436.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 12
    for name in ('static_exchange_network_bridge_audit.py', 'verify_static_exchange_network_bridge_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_435') == 616
    assert result.pop('unchanged_prior_science_hashes_231_434') == 613
    assert result['stage_saved_tests'] == 2038
    assert result['total_protected_evidence_hashes'] == 651
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 436 and int(latest[2]) >= 2044
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 64. 第436轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[436],
        execution_mode='one complete full-network static-simulation bridge with source theorem audit and explicit first-layer verification',
        scientific_base_through_round_by_round={436: 435}, fresh_tests_by_round={436: 6}, fresh_tests=6,
        stage_saved_tests=2044, science_hashes_verified_231_436=619,
        unchanged_prior_science_hashes_231_435=616, total_protected_evidence_hashes=654,
        unchanged_prior_evidence_hashes=651, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=12,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        full_finite_network_static_exchange_existence_theorem_applied=True,
        full_first_layer_64_columns_and_signal_certificate_verified=True,
        block_partition_gap_and_coupling_strengths_still_inputs=True,
        local_static_simulation_unknown_reference_contract_verified=True,
        final_Heisenberg_coupling_list_generated=False,
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
