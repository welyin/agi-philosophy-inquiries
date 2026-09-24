"""Integrate round 434 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round433_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round434_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_434_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 433
    assert checked['fresh_tests'] == dict(run=6, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'encoded_exchange_response_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    for key in ('derived_from_429_exchange_alone', 'full_cognitive_implementation_completed',
                'new_fundamental_flip_or_controlled_exchange_required_in_this_example',
                'exact_raw_433_hamiltonian_implemented', 'complete_433_three_edge_model_implemented',
                'three_dimensional_space_unconditionally_derived', 'full_GR_goal_completed', 'phase_closure_triggered'):
        assert not saved['scope'][key], key
    for key in ('continuous_exchange_only_conditional_response_constructed',
                'uniform_finite_time_error_and_reference_contract_proved',
                'known_singlet_code_and_effective_interaction_tools_attributed',
                'partition_couplings_gap_and_code_preparation_are_inputs'):
        assert saved['scope'][key], key
    note = HERE / 'research_note_434.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 12
    for name in ('encoded_exchange_response_audit.py', 'verify_encoded_exchange_response_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_433') == 610
    assert result.pop('unchanged_prior_science_hashes_231_432') == 607
    assert result['stage_saved_tests'] == 2026
    assert result['total_protected_evidence_hashes'] == 645
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 434 and int(latest[2]) >= 2032
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 62. 第434轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[434],
        execution_mode='one complete continuous encoded-exchange response round with uniform error proof',
        scientific_base_through_round_by_round={434: 433}, fresh_tests_by_round={434: 6}, fresh_tests=6,
        stage_saved_tests=2032, science_hashes_verified_231_434=613,
        unchanged_prior_science_hashes_231_433=610, total_protected_evidence_hashes=648,
        unchanged_prior_evidence_hashes=645, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=12,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        encoded_exchange_only_conditional_response_verified=True,
        block_partition_gap_and_coupling_strengths_still_inputs=True,
        uniform_finite_time_error_for_unknown_encoded_states_verified=True,
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
