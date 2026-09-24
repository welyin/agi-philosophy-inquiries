"""Integrate round 450 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round449_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round450_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_450_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 449
    assert checked['fresh_tests'] == dict(run=6, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'composite_agreement_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    for key in ('weak_agreement_adopted_as_new_cognitive_axiom', 'strong_429_theorem_invalidated', 'new_record_potentials_derived_from_429', 'all_encodings_or_effective_low_body_implementations_excluded', 'multi_subject_reference_composition_selected', 'full_spatial_dimension_or_GR_generated', 'full_GR_goal_completed', 'phase_closure_triggered'):
        assert not saved['scope'][key], key
    for key in ('weak_agreement_Hamiltonian_iff_classification_proved', 'microscopic_exchange_weak_contract_closed_under_composition', 'strong_contract_failure_for_independent_composite_Bell_inputs_proved', 'exact_all_time_raw_support_boundary_explicit', 'raw_record_comparison_dictionary_obstruction_and_two_subject_alternative_verified', 'internal_and_cross_subject_correlations_distinguished', 'fixed_identification_source_not_erased'):
        assert saved['scope'][key], key
    note = HERE / 'research_note_450.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 13
    for name in ('composite_agreement_audit.py', 'verify_composite_agreement_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_449') == 658
    assert result.pop('unchanged_prior_science_hashes_231_448') == 655
    assert result['stage_saved_tests'] == 2127
    assert result['total_protected_evidence_hashes'] == 693
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 450 and int(latest[2]) >= 2133
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 78. 第450轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[450],
        execution_mode='composite agreement contract classification, independent internally correlated input counterexample and comparison dictionary audit',
        scientific_base_through_round_by_round={450: 449}, fresh_tests_by_round={450: 6}, fresh_tests=6,
        stage_saved_tests=2133, science_hashes_verified_231_450=661,
        unchanged_prior_science_hashes_231_449=658, total_protected_evidence_hashes=696,
        unchanged_prior_evidence_hashes=693, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=13,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        weak_agreement_iff_swap_symmetry_and_derivative_necessity_verified=True,
        independent_composite_Bell_counterexample_verified=True,
        weak_contract_closure_under_microscopic_composition_verified=True,
        exact_raw_support_boundary_not_extended_to_effective_models=True,
        raw_record_dictionary_failure_and_twisted_identification_verified=True,
        strong_429_theorem_preserved_and_weak_condition_not_adopted=True,
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
