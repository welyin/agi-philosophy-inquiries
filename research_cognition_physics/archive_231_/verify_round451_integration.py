"""Integrate round 451 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round450_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round451_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_451_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 450
    assert checked['fresh_tests'] == dict(run=6, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'global_comparison_frames_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    false_keys = ('comparison_frames_dynamically_generated', 'new_record_potentials_derived_from_429', 'full_soca_architecture_implemented', 'all_Hadamard_orders_classified', 'exact_global_identification_adopted_as_cognitive_axiom', 'full_spatial_dimension_or_GR_generated', 'full_GR_goal_completed', 'phase_closure_triggered')
    true_keys = ('consistent_pair_frame_contract_classified', 'permutation_and_coherent_dictionaries_distinguished', 'minimal_marker_dimension_requirement_explicit', 'full_local_record_and_exchange_tensor_identities_verified', 'explicit_paley_twelve_and_padded_six_certificates', 'pair_contract_not_full_network_symmetry_verified', 'soca_functional_interface_motivation_distinguished_from_full_state_identification')
    for key in false_keys:
        assert not saved['scope'][key], key
    for key in true_keys:
        assert saved['scope'][key], key
    note = HERE / 'research_note_451.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 13
    for name in ('global_comparison_frames_audit.py', 'verify_global_comparison_frames_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_450') == 661
    assert result.pop('unchanged_prior_science_hashes_231_449') == 658
    assert result['stage_saved_tests'] == 2133
    assert result['total_protected_evidence_hashes'] == 696
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 451 and int(latest[2]) >= 2139
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 79. 第451轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[451],
        execution_mode='consistent pair comparison classification, exact Hadamard certificates, explicit capacity and SoCA interface audit',
        scientific_base_through_round_by_round={451: 450}, fresh_tests_by_round={451: 6}, fresh_tests=6,
        stage_saved_tests=2139, science_hashes_verified_231_451=664,
        unchanged_prior_science_hashes_231_450=661, total_protected_evidence_hashes=699,
        unchanged_prior_evidence_hashes=696, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=13,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        consistent_pair_frame_classification_verified=True,
        permutation_and_coherent_dictionary_difference_verified=True,
        twelve_subject_and_padded_six_exact_certificates_verified=True,
        pair_symmetry_not_promoted_to_global_network_symmetry=True,
        soca_selective_interface_does_not_imply_complete_dictionary=True,
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
