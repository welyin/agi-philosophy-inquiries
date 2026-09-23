"""Integrate round 393; preserve earlier science and check snapshots."""
import argparse
import ast
import json
from pathlib import Path
import verify_round392_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round393_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_393_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 392
    assert checked['batch_scientific_dependencies'] == []
    assert len(checked['new_file_hashes']) == 3
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'program_transport_index_audit_results.json')
    assert saved['checks'] == {'run': 10, 'failures': 0, 'errors': 0}
    assert saved['counterflow']['combined_index'] == '1'
    assert not saved['constant_autonomous_hamiltonian_derived']
    links = 0
    for path in (HERE / 'research_note_393.md', HERE / 'spatial_premise_closure_audit.md'):
        content = path.read_text(encoding='utf-8')
        for link in core.link_parser()(content):
            destination = (path.parent / link).resolve()
            assert destination.exists() or (pending and destination == TARGET), (path, link)
            links += 1
    formulas = core.text_checks(HERE / 'research_note_393.md')['display_formulas']
    assert formulas == 12
    for name in ('program_transport_index_audit.py', 'verify_program_transport_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    assert result.pop('science_hashes_verified_231_392') == 487
    assert result.pop('unchanged_prior_science_hashes_231_391') == 484
    assert result['stage_saved_tests'] == 1726
    assert result['total_protected_evidence_hashes'] == 508
    result.update(rounds=[393], execution_mode='program transport index and auxiliary completion',
                  scientific_base_through_round_by_round={393: 392},
                  fresh_tests_by_round={393: 10}, fresh_tests=10,
                  stage_saved_tests=1736, science_hashes_verified_231_393=490,
                  unchanged_prior_science_hashes_231_392=487,
                  total_protected_evidence_hashes=511, unchanged_prior_evidence_hashes=508,
                  new_scientific_display_formulas_checked=formulas,
                  local_links_checked=result['local_links_checked']+links,
                  navigation_and_notes_checked=result['navigation_and_notes_checked']+2,
                  constant_autonomous_hamiltonian_derived=False,
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
