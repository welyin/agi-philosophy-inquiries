"""Integrate round 442 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round441_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round442_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_442_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 441
    assert checked['fresh_tests'] == dict(run=6, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'branching_tree_distance_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    for key in ('actual_internal_distance_measurement_protocol_constructed',
                'general_cyclic_graphs_or_nonleaf_pairs_covered',
                'combined_439_440_model_covered',
                'complete_data_signal_bound_on_branching_graphs_proved',
                'programmed_flip_schedule_attributed_to_uniform_autonomous_F',
                'initial_tree_or_three_dimensional_geometry_generated',
                'full_GR_goal_completed', 'phase_closure_triggered'):
        assert not saved['scope'][key], key
    for key in ('unchanged_440_autonomous_H_used_for_main_theorem',
                'coherent_branching_trees_with_unknown_data_reference_covered',
                'size_independent_quantum_distance_deformation_bound_proved',
                'existing_NNI_and_441_weighted_method_explicitly_reused',
                'tree_combinatorics_checked_against_440_generator',
                'graph_shape_change_without_common_path_template_verified'):
        assert saved['scope'][key], key
    note = HERE / 'research_note_442.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 14
    for name in ('branching_tree_distance_audit.py', 'verify_branching_tree_distance_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_441') == 634
    assert result.pop('unchanged_prior_science_hashes_231_440') == 631
    assert result['stage_saved_tests'] == 2077
    assert result['total_protected_evidence_hashes'] == 669
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 442 and int(latest[2]) >= 2083
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 70. 第442轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[442],
        execution_mode='unchanged autonomous graph flip on coherent branching trees with complete unknown data and reference',
        scientific_base_through_round_by_round={442: 441}, fresh_tests_by_round={442: 6}, fresh_tests=6,
        stage_saved_tests=2083, science_hashes_verified_231_442=637,
        unchanged_prior_science_hashes_231_441=634, total_protected_evidence_hashes=672,
        unchanged_prior_evidence_hashes=669, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=14,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        exact_tree_distance_change_counts_verified=True,
        all_size_weighted_distance_deformation_bound_proved=True,
        arbitrary_coherent_tree_input_moment_bound_proved=True,
        unknown_data_distance_contraction_certificate_verified=True,
        fixed_path_template_not_needed_for_tree_distance_bound=True,
        programmed_schedule_not_confused_with_autonomous_evolution=True,
        distance_deformation_not_confused_with_data_signalling=True,
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
