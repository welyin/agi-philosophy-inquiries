"""Integrate round 396; preserve earlier science and check snapshots."""
import argparse
import ast
import json
from pathlib import Path
import verify_round395_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round396_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_396_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 395
    assert checked['batch_scientific_dependencies'] == []
    assert len(checked['new_file_hashes']) == 3
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'timed_position_response_audit_results.json')
    assert saved['checks'] == {'run': 6, 'failures': 0, 'errors': 0}
    assert saved['timing']['exact_deadline_case']['both_on_time']
    assert saved['timing']['positive_processing_case']['both_on_time']
    assert not saved['timing']['zero_slack_rejects_positive_delay']['both_on_time']
    assert not saved['trusted_unique_local_marker_disproved']
    assert not saved['quantum_identity_authentication_disproved']
    assert not saved['full_position_generation_completed']
    links = 0
    for path in (HERE / 'research_note_396.md', HERE / 'spatial_premise_closure_audit.md'):
        content = path.read_text(encoding='utf-8')
        for link in core.link_parser()(content):
            destination = (path.parent / link).resolve()
            assert destination.exists() or (pending and destination == TARGET), (path, link)
            links += 1
    formulas = core.text_checks(HERE / 'research_note_396.md')['display_formulas']
    assert formulas == 11
    for name in ('timed_position_response_audit.py', 'verify_timed_position_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    assert result.pop('science_hashes_verified_231_395') == 496
    assert result.pop('unchanged_prior_science_hashes_231_394') == 493
    assert result['stage_saved_tests'] == 1746
    assert result['total_protected_evidence_hashes'] == 517
    result.update(rounds=[396], execution_mode='timed position-response counterexample',
                  scientific_base_through_round_by_round={396: 395},
                  fresh_tests_by_round={396: 6}, fresh_tests=6,
                  stage_saved_tests=1752, science_hashes_verified_231_396=499,
                  unchanged_prior_science_hashes_231_395=496,
                  total_protected_evidence_hashes=520, unchanged_prior_evidence_hashes=517,
                  new_scientific_display_formulas_checked=formulas,
                  local_links_checked=result['local_links_checked']+links,
                  navigation_and_notes_checked=result['navigation_and_notes_checked']+2,
                  trusted_unique_local_marker_disproved=False,
                  quantum_identity_authentication_disproved=False,
                  full_position_generation_completed=False,
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
