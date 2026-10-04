"""Integrate round 394; preserve earlier science and check snapshots."""
import argparse
import ast
import json
from pathlib import Path
import verify_round393_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round394_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_394_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 393
    assert checked['batch_scientific_dependencies'] == []
    assert len(checked['new_file_hashes']) == 3
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'static_counterflow_audit_results.json')
    assert saved['checks'] == {'run': 5, 'failures': 0, 'errors': 0}
    assert saved['static_finite_range_transport_generator_excluded_analytically']
    assert not saved['prepared_returning_clock_excluded']
    assert not saved['full_computing_rule_static_generator_excluded']
    assert not saved['spatial_dimension_selected']
    links = 0
    for path in (HERE / 'research_note_394.md', HERE / 'spatial_premise_closure_audit.md'):
        content = path.read_text(encoding='utf-8')
        for link in core.link_parser()(content):
            destination = (path.parent / link).resolve()
            assert destination.exists() or (pending and destination == TARGET), (path, link)
            links += 1
    formulas = core.text_checks(HERE / 'research_note_394.md')['display_formulas']
    assert formulas == 10
    for name in ('static_counterflow_audit.py', 'verify_static_counterflow_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    assert result.pop('science_hashes_verified_231_393') == 490
    assert result.pop('unchanged_prior_science_hashes_231_392') == 487
    assert result['stage_saved_tests'] == 1736
    assert result['total_protected_evidence_hashes'] == 511
    result.update(rounds=[394], execution_mode='static counterflow generator obstruction',
                  scientific_base_through_round_by_round={394: 393},
                  fresh_tests_by_round={394: 5}, fresh_tests=5,
                  stage_saved_tests=1741, science_hashes_verified_231_394=493,
                  unchanged_prior_science_hashes_231_393=490,
                  total_protected_evidence_hashes=514, unchanged_prior_evidence_hashes=511,
                  new_scientific_display_formulas_checked=formulas,
                  local_links_checked=result['local_links_checked']+links,
                  navigation_and_notes_checked=result['navigation_and_notes_checked']+2,
                  pure_transport_static_finite_range_generator_excluded=True,
                  prepared_returning_clock_excluded=False,
                  full_computing_rule_static_generator_excluded=False,
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
