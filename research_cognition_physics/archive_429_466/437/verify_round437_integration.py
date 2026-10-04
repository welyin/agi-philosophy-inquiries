"""Integrate round 437 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round436_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round437_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_437_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 436
    assert checked['fresh_tests'] == dict(run=6, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'exchange_participation_capacity_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    for key in ('all_exchange_based_sparse_mechanisms_refuted',
                'degree_three_or_physical_space_generated',
                'final_Heisenberg_coupling_list_or_hardware_run_completed',
                'full_GR_goal_completed', 'phase_closure_triggered'):
        assert not saved['scope'][key], key
    for key in ('known_singlet_monogamy_reused_not_claimed_as_new',
                'exchange_correlation_budget_distinguished_from_participation_capacity',
                'all_N_all_unknown_data_activation_lower_bound_proved',
                'exact_joint_dense_graph_counterexample_with_reference_proved',
                'additional_relation_hardware_and_initial_energy_accounted',
                'naive_433_automatic_sparse_limit_refuted'):
        assert saved['scope'][key], key
    note = HERE / 'research_note_437.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 12
    for name in ('exchange_participation_capacity_audit.py', 'verify_exchange_participation_capacity_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_436') == 619
    assert result.pop('unchanged_prior_science_hashes_231_435') == 616
    assert result['stage_saved_tests'] == 2044
    assert result['total_protected_evidence_hashes'] == 654
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 437 and int(latest[2]) >= 2050
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 65. 第437轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[437],
        execution_mode='one complete all-scale participation-capacity audit, exact dense-sector counterexample and internal resource account',
        scientific_base_through_round_by_round={437: 436}, fresh_tests_by_round={437: 6}, fresh_tests=6,
        stage_saved_tests=2050, science_hashes_verified_231_437=622,
        unchanged_prior_science_hashes_231_436=619, total_protected_evidence_hashes=657,
        unchanged_prior_evidence_hashes=654, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=12,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        all_input_scale_uniform_activation_bound_verified=True,
        exact_joint_dense_graph_distribution_verified=True,
        complete_symmetric_data_and_reference_retained=True,
        internal_relation_carriers_and_initial_energy_accounted=True,
        all_exchange_sparsity_routes_refuted=False,
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
