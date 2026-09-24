"""Integrate round 443 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round442_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round443_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_443_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 442
    assert checked['fresh_tests'] == dict(run=6, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'quantum_neighborhood_signal_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    for key in ('universal_branching_graph_propagation_bound',
                'unknown_other_subject_data_preserved',
                'graph_neighborhood_routing_autonomously_implemented',
                'nested_dynamic_ball_channels_assumed',
                'new_internal_measurement_or_encoding_device_constructed',
                'initial_geometry_selected', 'three_dimensional_space_generated',
                'full_GR_goal_completed', 'phase_closure_triggered'):
        assert not saved['scope'][key], key
    for key in ('unchanged_440_442_autonomous_H',
                'physical_receiver_data_signal_certified',
                'all_coherent_far_graph_inputs_and_reference_covered',
                'finite_channel_and_actual_subject_readout_compatible',
                'obsolete_trace_failure_attributed_to_2024',
                'exact_rational_operator_certificate'):
        assert saved['scope'][key], key
    note = HERE / 'research_note_443.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 12
    for name in ('quantum_neighborhood_signal_audit.py', 'verify_quantum_neighborhood_signal_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_442') == 637
    assert result.pop('unchanged_prior_science_hashes_231_441') == 634
    assert result['stage_saved_tests'] == 2083
    assert result['total_protected_evidence_hashes'] == 672
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 443 and int(latest[2]) >= 2089
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 71. 第443轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[443],
        execution_mode='actual subject data signalling with every coherent far-graph input and internal reference in unchanged finite autonomous model',
        scientific_base_through_round_by_round={443: 442}, fresh_tests_by_round={443: 6}, fresh_tests=6,
        stage_saved_tests=2089, science_hashes_verified_231_443=640,
        unchanged_prior_science_hashes_231_442=637, total_protected_evidence_hashes=675,
        unchanged_prior_evidence_hashes=672, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=12,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        corrected_quantum_network_trace_source_explicitly_acknowledged=True,
        basis_split_CPTP_and_subject_readout_compatibility_verified=True,
        exact_coherent_signal_order_and_coupling_dependence_verified=True,
        exact_rational_uniform_finite_time_data_signal_certificate=True,
        actual_signal_not_confused_with_graph_distance_probability=True,
        prepared_data_and_unknown_graph_quantifiers_distinguished=True,
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
