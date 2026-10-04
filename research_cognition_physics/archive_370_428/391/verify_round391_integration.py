"""Integrate round 391 without rerunning or rewriting old scientific results."""
import argparse
import ast
import json
from pathlib import Path
import verify_round390_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round391_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_391_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 390
    assert checked['batch_scientific_dependencies'] == []
    assert len(checked['new_file_hashes']) == 3
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'shared_sharp_record_audit_results.json')
    assert saved['checks'] == {'run': 9, 'failures': 0, 'errors': 0}
    paths = [HERE / 'research_note_391.md', HERE / 'spatial_premise_closure_audit.md']
    links = 0
    for path in paths:
        text = path.read_text(encoding='utf-8')
        assert not any(ord(c) < 32 and c not in '\r\n\t' for c in text)
        for link in core.link_parser()(text):
            destination = (path.parent / link).resolve()
            assert destination.exists() or (pending and destination == TARGET), (path, link)
            links += 1
    formulas = core.text_checks(paths[0])['display_formulas']
    assert formulas == 12
    for name in ('shared_sharp_record_audit.py', 'verify_shared_sharp_record_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    assert result.pop('science_hashes_verified_231_390') == 481
    assert result.pop('unchanged_prior_science_hashes_231_389') == 478
    assert result['stage_saved_tests'] == 1707
    assert result['total_protected_evidence_hashes'] == 502
    result.update(rounds=[391], execution_mode='complete shared sharp record round on frozen round-390 baseline',
                  scientific_base_through_round_by_round={391: 390},
                  fresh_tests_by_round={391: 9}, fresh_tests=9,
                  stage_saved_tests=1716, science_hashes_verified_231_391=484,
                  unchanged_prior_science_hashes_231_390=481,
                  total_protected_evidence_hashes=505, unchanged_prior_evidence_hashes=502,
                  new_scientific_display_formulas_checked=formulas,
                  local_links_checked=result['local_links_checked']+links,
                  navigation_and_notes_checked=result['navigation_and_notes_checked']+2,
                  common_record_algebra_identified_with_physical_space=False,
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
