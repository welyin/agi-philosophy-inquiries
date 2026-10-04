"""Integrate round 399; preserve earlier science and check snapshots."""
import argparse
import ast
import json
from pathlib import Path
import verify_round398_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round399_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_399_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 398
    assert checked['batch_scientific_dependencies'] == []
    assert len(checked['new_file_hashes']) == 3
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'operational_region_certificate_audit_results.json')
    assert saved['checks'] == {'run': 8, 'failures': 0, 'errors': 0}
    assert not saved['assumed_spatial_distance']
    assert not saved['exact_support_gap_required']
    assert saved['all_background_states_covered_by_proof']
    assert saved['arbitrary_outside_control_covered_by_proof']
    assert not saved['postselected_conditional_state_guarantee']
    assert not saved['finite_geometry_or_dimension_generated']
    links = 0
    for path in (HERE / 'research_note_399.md', HERE / 'spatial_premise_closure_audit.md'):
        content = path.read_text(encoding='utf-8')
        for link in core.link_parser()(content):
            destination = (path.parent / link).resolve()
            assert destination.exists() or (pending and destination == TARGET), (path, link)
            links += 1
    formulas = core.text_checks(HERE / 'research_note_399.md')['display_formulas']
    assert formulas == 13
    for name in ('operational_region_certificate_audit.py', 'verify_operational_region_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    assert result.pop('science_hashes_verified_231_398') == 505
    assert result.pop('unchanged_prior_science_hashes_231_397') == 502
    assert result['stage_saved_tests'] == 1765
    assert result['total_protected_evidence_hashes'] == 526
    result.update(rounds=[399], execution_mode='finite-time operational region certificates',
                  scientific_base_through_round_by_round={399: 398},
                  fresh_tests_by_round={399: 8}, fresh_tests=8,
                  stage_saved_tests=1773, science_hashes_verified_231_399=508,
                  unchanged_prior_science_hashes_231_398=505,
                  total_protected_evidence_hashes=529, unchanged_prior_evidence_hashes=526,
                  new_scientific_display_formulas_checked=formulas,
                  local_links_checked=result['local_links_checked']+links,
                  navigation_and_notes_checked=result['navigation_and_notes_checked']+2,
                  finite_geometry_or_dimension_generated=False,
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
