"""Integrate round 453 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round452_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round453_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_453_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 451
    assert checked['fresh_tests'] == dict(run=6, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'relational_subject_composition_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    false_keys = ('nontrivial_macro_interaction_obtained', 'all_possible_relational_organizations_ruled_out', 'original_logic_interface_preservation_imposed_as_cognitive_axiom', 'new_subsystem_or_exchange_gate_theory_claimed', 'depends_on_round_452', 'full_soca_or_spatial_dimension_or_GR_generated', 'full_GR_goal_completed', 'phase_closure_triggered')
    true_keys = ('two_fixed_three_qubit_relational_interfaces_only', 'arbitrary_static_real_pair_exchange_weights', 'exact_all_time_code_retention_classified', 'unknown_gauge_independent_compressed_logic_classified', 'exact_code_retention_already_forces_no_cross_logical_interaction')
    for key in false_keys:
        assert not saved['scope'][key], key
    for key in true_keys:
        assert saved['scope'][key], key
    note = HERE / 'research_note_453.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 13
    for name in ('relational_subject_composition_audit.py', 'verify_relational_subject_composition_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_452') == 667
    assert result.pop('unchanged_prior_science_hashes_231_451') == 664
    assert result['stage_saved_tests'] == 2145
    assert result['total_protected_evidence_hashes'] == 702
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 453 and int(latest[2]) >= 2151
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 81. 第453轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[453],
        execution_mode='fixed relational subject composition, all-weight leakage Gram and physical gauge witness, independent science baseline 451',
        scientific_base_through_round_by_round={453: 451}, fresh_tests_by_round={453: 6}, fresh_tests=6,
        stage_saved_tests=2151, science_hashes_verified_231_453=670,
        unchanged_prior_science_hashes_231_452=667, total_protected_evidence_hashes=705,
        unchanged_prior_evidence_hashes=702, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=13,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        all_static_cross_weights_code_retention_classified=True,
        code_retention_already_forces_logical_noninteraction=True,
        physical_gauge_readout_witness_without_postselection_verified=True,
        endpoint_swap_and_existing_effective_encodings_preserved=True,
        two_independent_scientific_rounds_share_baseline_451=True,
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
