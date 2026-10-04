"""Integrate round 390 without rerunning or rewriting old scientific results."""
import argparse
import ast
import json
from pathlib import Path
import verify_round389_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round390_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_390_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 389
    assert checked['batch_scientific_dependencies'] == []
    assert len(checked['new_file_hashes']) == 3
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'finite_time_trace_audit_results.json')
    assert saved['checks'] == {'run': 9, 'failures': 0, 'errors': 0}
    paths = [HERE / 'research_note_390.md', HERE / 'spatial_premise_closure_audit.md']
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
    for name in ('finite_time_trace_audit.py', 'verify_finite_time_trace_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    assert result.pop('science_hashes_verified_231_389') == 478
    assert result.pop('unchanged_prior_science_hashes_231_388') == 475
    assert result['stage_saved_tests'] == 1698
    assert result['total_protected_evidence_hashes'] == 499
    result.update(rounds=[390], execution_mode='complete finite-time trace round on frozen round-389 baseline',
                  scientific_base_through_round_by_round={390: 389},
                  fresh_tests_by_round={390: 9}, fresh_tests=9,
                  stage_saved_tests=1707, science_hashes_verified_231_390=481,
                  unchanged_prior_science_hashes_231_389=478,
                  total_protected_evidence_hashes=502, unchanged_prior_evidence_hashes=499,
                  new_scientific_display_formulas_checked=formulas,
                  local_links_checked=result['local_links_checked']+links,
                  navigation_and_notes_checked=result['navigation_and_notes_checked']+2,
                  ontology_and_overwrite_candidate_proved_as_physical_law=False,
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
