"""Integrate round 427 while retaining historical scientific hashes."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round426_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round427_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_427_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 426
    assert checked['fresh_tests'] == dict(run=5, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'finite_probe_topology_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (5, 0, 0)
    for key in ('three_dimensional_space_unconditionally_derived', 'full_GR_goal_completed',
                'full_cognitive_countermodel_completed', 'phase_closure_triggered',
                'Baire_derived_from_cognitive_principles', 'preparation_and_readout_autonomously_generated',
                'finite_test_list_is_exact_global_tomography', 'sampled_inputs_used_to_prove_diamond_bound'):
        assert not saved['scope'][key], key
    assert saved['scope']['test_Baire_topology_recovers_uniform_topology']
    assert saved['scope']['compact_budget_finite_test_certificate_proved']
    note = HERE / 'research_note_427.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 11
    for name in ('finite_probe_topology_audit.py', 'verify_finite_probe_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_426') == 589
    assert result.pop('unchanged_prior_science_hashes_231_425') == 586
    assert result['stage_saved_tests'] == 1985
    assert result['total_protected_evidence_hashes'] == 624
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)',
                       (HERE.parent / 'research_direction.md').read_text(encoding='utf-8'))
    assert latest and int(latest[1]) >= 427 and int(latest[2]) >= 1990
    result.update(date='2026-09-24', rounds=[427],
        execution_mode='one complete finite-test research round with independent proof and reproduction audit',
        scientific_base_through_round_by_round={427: 426}, fresh_tests_by_round={427: 5}, fresh_tests=5,
        stage_saved_tests=1990, science_hashes_verified_231_427=592,
        unchanged_prior_science_hashes_231_426=589, total_protected_evidence_hashes=627,
        unchanged_prior_evidence_hashes=624, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=11,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        finite_probability_certificate_includes_infinite_tail=True,
        test_topology_Baire_to_uniform_topology_bridge_verified=True,
        additional_coherent_testing_permissions_explicit=True,
        cognitive_conjecture_retested_as_new_round=False,
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
