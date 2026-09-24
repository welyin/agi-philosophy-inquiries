"""Integrate round 452 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round451_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round452_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_452_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 451
    assert checked['fresh_tests'] == dict(run=6, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'macro_interface_closure_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    false_keys = ('all_arbitrary_quantum_compressions_excluded', 'complete_macroscopic_agent_or_SoCA_implemented', 'macro_full_tomography_or_control_derived', 'new_rule_derived_from_cognitive_principles', 'full_spatial_dimension_or_GR_generated', 'full_GR_goal_completed', 'phase_closure_triggered')
    true_keys = ('existing_physical_readout_interface_used', 'exact_minimal_invariant_Cstar_algebra_classified', 'matching_coherence_future_physical_signal_certified', 'actual_raw_finite_window_gap_transferred', 'linear_statistics_and_quantum_algebra_distinguished')
    for key in false_keys:
        assert not saved['scope'][key], key
    for key in true_keys:
        assert saved['scope'][key], key
    note = HERE / 'research_note_452.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 13
    for name in ('macro_interface_closure_audit.py', 'verify_macro_interface_closure_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_451') == 664
    assert result.pop('unchanged_prior_science_hashes_231_450') == 661
    assert result['stage_saved_tests'] == 2139
    assert result['total_protected_evidence_hashes'] == 699
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 452 and int(latest[2]) >= 2145
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 80. 第452轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[452],
        execution_mode='physical output algebra classification and finite-window relation-coherence certificate, independent science baseline 451',
        scientific_base_through_round_by_round={452: 451}, fresh_tests_by_round={452: 6}, fresh_tests=6,
        stage_saved_tests=2145, science_hashes_verified_231_452=667,
        unchanged_prior_science_hashes_231_451=664, total_protected_evidence_hashes=702,
        unchanged_prior_evidence_hashes=699, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=13,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        physical_output_algebra_and_parameter_boundary_verified=True,
        same_classical_relation_summary_future_signal_verified=True,
        raw_window_gap_exceeds_one_tenth_verified=True,
        stationary_indistinguishability_limits_tomography=True,
        effective_algebra_not_promoted_to_exact_raw_compression=True,
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
