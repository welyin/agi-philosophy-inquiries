"""Integrate round 401; preserve earlier science and check snapshots."""
import argparse
import ast
import json
from pathlib import Path
import verify_round400_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round401_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_401_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 400
    assert checked['batch_scientific_dependencies'] == []
    assert len(checked['new_file_hashes']) == 3
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'persistent_prefix_history_audit_results.json')
    assert saved['checks'] == {'run': 8, 'failures': 0, 'errors': 0}
    assert saved['same_fixed_hamiltonian_for_all_finite_prefixes']
    assert saved['infinite_internal_workspace_assumed']
    assert saved['frozen_output_contract_required']
    for key in ('postselection_used', 'finite_time_exact_completion_claimed',
                'arbitrary_repeated_readout_guaranteed', 'whole_state_convergence_claimed',
                'full_sustained_cognition_model_constructed', 'spatial_dimension_generated'):
        assert not saved[key], key
    links = 0
    for path in (HERE / 'research_note_401.md', HERE / 'spatial_premise_closure_audit.md'):
        content = path.read_text(encoding='utf-8')
        for link in core.link_parser()(content):
            destination = (path.parent / link).resolve()
            assert destination.exists() or (pending and destination == TARGET), (path, link)
            links += 1
    formulas = core.text_checks(HERE / 'research_note_401.md')['display_formulas']
    assert formulas == 11
    for name in ('persistent_prefix_history_audit.py', 'verify_persistent_prefix_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    assert result.pop('science_hashes_verified_231_400') == 511
    assert result.pop('unchanged_prior_science_hashes_231_399') == 508
    assert result['stage_saved_tests'] == 1781
    assert result['total_protected_evidence_hashes'] == 532
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    assert '完成231—401轮' in direction
    assert '401／1789' in direction
    result.update(rounds=[401], execution_mode='persistent finite archives under one fixed internal Hamiltonian',
                  scientific_base_through_round_by_round={401: 400},
                  fresh_tests_by_round={401: 8}, fresh_tests=8,
                  stage_saved_tests=1789, science_hashes_verified_231_401=514,
                  unchanged_prior_science_hashes_231_400=511,
                  total_protected_evidence_hashes=535, unchanged_prior_evidence_hashes=532,
                  new_scientific_display_formulas_checked=formulas,
                  local_links_checked=result['local_links_checked']+links,
                  navigation_and_notes_checked=result['navigation_and_notes_checked']+2,
                  infinite_internal_workspace_assumed=True,
                  full_sustained_cognition_model_constructed=False,
                  spatial_dimension_generated=False,
                  independent_final_agent_review_completed=False,
                  all_reported_checks_passed=True)
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
