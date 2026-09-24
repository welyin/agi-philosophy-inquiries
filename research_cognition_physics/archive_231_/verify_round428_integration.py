"""Integrate round 428 without modifying any historical scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round427_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round428_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_428_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 427
    assert checked['fresh_tests'] == dict(run=6, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'internal_termination_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    for key in ('three_dimensional_space_unconditionally_derived', 'full_GR_goal_completed',
                'full_cognitive_implementation_completed', 'phase_closure_triggered',
                'common_deadline_or_fixed_hamiltonian_required',
                'exact_universal_almost_sure_processor_under_contract',
                'candidate_principle_itself_refuted', 'known_phase_processor_claimed_as_new_discovery'):
        assert not saved['scope'][key], key
    for key in ('internal_generation_is_user_authorized_candidate',
                'countable_first_termination_reduction_proved',
                'separable_normal_program_contract_explicit',
                'all_target_information_accounted_at_initial_interface',
                'growing_memory_and_random_retries_allowed'):
        assert saved['scope'][key], key
    note = HERE / 'research_note_428.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 11
    for name in ('internal_termination_audit.py', 'verify_internal_termination_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_427') == 592
    assert result.pop('unchanged_prior_science_hashes_231_426') == 589
    assert result['stage_saved_tests'] == 1990
    assert result['total_protected_evidence_hashes'] == 627
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 428 and int(latest[2]) >= 1996
    assert '统一内生演化' in direction
    assert '## 56. 第428轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[428],
        execution_mode='one complete internal-termination research round with independent proof and reproduction review',
        scientific_base_through_round_by_round={428: 427}, fresh_tests_by_round={428: 6}, fresh_tests=6,
        stage_saved_tests=1996, science_hashes_verified_231_428=595,
        unchanged_prior_science_hashes_231_427=592, total_protected_evidence_hashes=630,
        unchanged_prior_evidence_hashes=627, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=11,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        candidate_principle_distinguished_from_specific_law=True,
        internal_random_termination_reduced_to_normal_processor=True,
        known_stochastic_phase_program_explicitly_attributed=True,
        no_new_cognitive_axiom_inferred_from_implementation_conditions=True,
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
