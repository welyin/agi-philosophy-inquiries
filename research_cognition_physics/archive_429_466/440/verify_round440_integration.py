"""Integrate round 440 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round439_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round440_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_440_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 439
    assert checked['fresh_tests'] == dict(run=7, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'simultaneous_edge_exchange_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (7, 0, 0)
    for key in ('graph_flip_generator_selected_by_429',
                'individual_degrees_conserved_by_combined_model',
                'classical_mixing_claimed_for_closed_unitary',
                'physical_dimension_or_new_subject_growth_derived',
                'final_Heisenberg_hardware_compiled', 'full_GR_goal_completed', 'phase_closure_triggered'):
        assert not saved['scope'][key], key
    for key in ('existing_quantum_graphity_flip_explicitly_attributed',
                'simultaneous_rewiring_without_extra_slots_proved',
                'endpoint_degree_charge_and_partition_invariants_proved',
                'unknown_data_reference_and_relational_frame_verified',
                'combined_C2_all_N_configuration_classification_proved',
                'size_independent_endpoint_flip_norm_bound_proved',
                'minimum_four_edge_transition_support_proved'):
        assert saved['scope'][key], key
    note = HERE / 'research_note_440.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 14
    for name in ('simultaneous_edge_exchange_audit.py', 'verify_simultaneous_edge_exchange_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_439') == 628
    assert result.pop('unchanged_prior_science_hashes_231_438') == 625
    assert result['stage_saved_tests'] == 2064
    assert result['total_protected_evidence_hashes'] == 663
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 440 and int(latest[2]) >= 2071
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 68. 第440轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[440],
        execution_mode='one complete simultaneous flip model with unknown data, degree constraints and exact relational frame',
        scientific_base_through_round_by_round={440: 439}, fresh_tests_by_round={440: 7}, fresh_tests=7,
        stage_saved_tests=2071, science_hashes_verified_231_440=631,
        unchanged_prior_science_hashes_231_439=628, total_protected_evidence_hashes=666,
        unchanged_prior_evidence_hashes=663, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=14,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        known_quantum_graphity_flip_normalization_verified=True,
        no_temporary_slots_and_endpoint_flip_bound_verified=True,
        combined_capacity_two_configuration_count_proved=True,
        all_unknown_data_reference_transfer_certificate_verified=True,
        complete_controlled_frame_and_label_readouts_verified=True,
        minimum_four_edge_support_obstruction_proved=True,
        markov_mixing_not_transferred_to_closed_quantum_model=True,
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
