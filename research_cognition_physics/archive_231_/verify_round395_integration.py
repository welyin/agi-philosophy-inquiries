"""Integrate round 395; preserve earlier science and check snapshots."""
import argparse
import ast
import json
from pathlib import Path
import verify_round394_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round395_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_395_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 394
    assert checked['batch_scientific_dependencies'] == []
    assert len(checked['new_file_hashes']) == 3
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'multiscale_readout_carrier_audit_results.json')
    assert saved['checks'] == {'run': 5, 'failures': 0, 'errors': 0}
    assert saved['noisy_case']['minimum_integer_dimension'] == 11
    assert saved['preparation_and_measurement_independence_required']
    assert not saved['memory_dimension_identified_with_spatial_dimension']
    assert not saved['full_position_generation_completed']
    links = 0
    for path in (HERE / 'research_note_395.md', HERE / 'spatial_premise_closure_audit.md'):
        content = path.read_text(encoding='utf-8')
        for link in core.link_parser()(content):
            destination = (path.parent / link).resolve()
            assert destination.exists() or (pending and destination == TARGET), (path, link)
            links += 1
    formulas = core.text_checks(HERE / 'research_note_395.md')['display_formulas']
    assert formulas == 12
    for name in ('multiscale_readout_carrier_audit.py', 'verify_multiscale_readout_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    assert result.pop('science_hashes_verified_231_394') == 493
    assert result.pop('unchanged_prior_science_hashes_231_393') == 490
    assert result['stage_saved_tests'] == 1741
    assert result['total_protected_evidence_hashes'] == 514
    result.update(rounds=[395], execution_mode='finite multiscale readout carrier witness',
                  scientific_base_through_round_by_round={395: 394},
                  fresh_tests_by_round={395: 5}, fresh_tests=5,
                  stage_saved_tests=1746, science_hashes_verified_231_395=496,
                  unchanged_prior_science_hashes_231_394=493,
                  total_protected_evidence_hashes=517, unchanged_prior_evidence_hashes=514,
                  new_scientific_display_formulas_checked=formulas,
                  local_links_checked=result['local_links_checked']+links,
                  navigation_and_notes_checked=result['navigation_and_notes_checked']+2,
                  memory_dimension_identified_with_spatial_dimension=False,
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
