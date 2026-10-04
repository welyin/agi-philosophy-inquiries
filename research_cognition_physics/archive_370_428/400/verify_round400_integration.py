"""Integrate round 400; preserve earlier science and check snapshots."""
import argparse
import ast
import json
from pathlib import Path
import verify_round399_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round400_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_400_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 399
    assert checked['batch_scientific_dependencies'] == []
    assert len(checked['new_file_hashes']) == 3
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'joint_region_composition_audit_results.json')
    assert saved['checks'] == {'run': 8, 'failures': 0, 'errors': 0}
    assert saved['near_unitary_reference_required']
    assert saved['actual_common_channel_required']
    assert not saved['merely_matching_marginals_sufficient']
    assert not saved['independent_receiver_noise_assumed']
    assert not saved['ideal_local_references_assumed_to_commute']
    assert not saved['uniform_infinite_receiver_bound_proved']
    assert not saved['spatial_dimension_generated']
    links = 0
    for path in (HERE / 'research_note_400.md', HERE / 'spatial_premise_closure_audit.md'):
        content = path.read_text(encoding='utf-8')
        for link in core.link_parser()(content):
            destination = (path.parent / link).resolve()
            assert destination.exists() or (pending and destination == TARGET), (path, link)
            links += 1
    formulas = core.text_checks(HERE / 'research_note_400.md')['display_formulas']
    assert formulas == 13
    for name in ('joint_region_composition_audit.py', 'verify_joint_region_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    assert result.pop('science_hashes_verified_231_399') == 508
    assert result.pop('unchanged_prior_science_hashes_231_398') == 505
    assert result['stage_saved_tests'] == 1773
    assert result['total_protected_evidence_hashes'] == 529
    result.update(rounds=[400], execution_mode='joint receiver composition of operational regions',
                  scientific_base_through_round_by_round={400: 399},
                  fresh_tests_by_round={400: 8}, fresh_tests=8,
                  stage_saved_tests=1781, science_hashes_verified_231_400=511,
                  unchanged_prior_science_hashes_231_399=508,
                  total_protected_evidence_hashes=532, unchanged_prior_evidence_hashes=529,
                  new_scientific_display_formulas_checked=formulas,
                  local_links_checked=result['local_links_checked']+links,
                  navigation_and_notes_checked=result['navigation_and_notes_checked']+2,
                  spatial_dimension_generated=False,
                  uniform_infinite_receiver_bound_proved=False,
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
