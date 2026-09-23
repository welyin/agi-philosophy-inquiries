"""Integrate round 389; reuse old hash checks without rerunning old science."""
import argparse
import ast
import json
from pathlib import Path
import verify_round388_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round389_integration_checks.json'


def verify(pending=False):
    # Only change the loaded module's pending-report destination, never its file.
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_389_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 388
    assert checked['batch_scientific_dependencies'] == []
    assert len(checked['new_file_hashes']) == 3
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'additive_history_cone_audit_results.json')
    assert saved['checks'] == {'run': 10, 'failures': 0, 'errors': 0}
    paths = [HERE / 'research_note_389.md', HERE / 'spatial_premise_closure_audit.md']
    links = 0
    for path in paths:
        text = path.read_text(encoding='utf-8')
        assert not any(ord(c) < 32 and c not in '\r\n\t' for c in text)
        for link in core.link_parser()(text):
            target = (path.parent / link).resolve()
            assert target.exists() or (pending and target == TARGET), (path, link)
            links += 1
    formulas = core.text_checks(paths[0])['display_formulas']
    assert formulas == 12
    for name in ('additive_history_cone_audit.py', 'verify_additive_history_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    old_science_count = result.pop('science_hashes_verified_231_388')
    old_science_unchanged = result.pop('unchanged_prior_science_hashes_231_387')
    assert old_science_count == 475 and old_science_unchanged == 472
    assert result['stage_saved_tests'] == 1688
    assert result['total_protected_evidence_hashes'] == 496
    result.update(rounds=[389], execution_mode='complete additive history order round on frozen round-388 baseline',
                  scientific_base_through_round_by_round={389: 388},
                  fresh_tests_by_round={389: 10}, fresh_tests=10,
                  stage_saved_tests=1698, science_hashes_verified_231_389=478,
                  unchanged_prior_science_hashes_231_388=475,
                  total_protected_evidence_hashes=499, unchanged_prior_evidence_hashes=496,
                  new_scientific_display_formulas_checked=formulas,
                  local_links_checked=result['local_links_checked']+links,
                  navigation_and_notes_checked=result['navigation_and_notes_checked']+2,
                  unnumbered_four_function_review_counted_as_science=False,
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
