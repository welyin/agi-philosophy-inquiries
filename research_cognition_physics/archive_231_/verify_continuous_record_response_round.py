"""Read-only round 458 science and historical-evidence verification."""
import argparse
import ast
import importlib.util
import json
from pathlib import Path

import verify_round457_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'research_round_458_checks.json'
NAMES = ['continuous_record_response_audit.py', 'continuous_record_response_audit_results.json', 'research_note_458.md']


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        frozen = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert frozen['science_hashes_verified_231_457'] == 682
    assert frozen['total_protected_evidence_hashes'] == 717
    if TARGET.exists():
        for name, sha in core.read(TARGET)['new_file_hashes'].items():
            assert core.digest(HERE / name) == sha, name
    revised = HERE/'research_note_456_text_v2.md'
    old_note = (HERE/'research_note_456.md').read_text(encoding='utf-8')
    revised_note = revised.read_text(encoding='utf-8').split('\n\n',1)[1]
    assert revised_note == old_note.replace(r'\frac1{24}quad\Longrightarrow',r'\frac1{24}\quad\Longrightarrow',1)
    assert core.text_checks(revised)['display_formulas'] == 11
    if TARGET.exists():
        for name, sha in core.read(TARGET)['preserved_draft_hashes'].items():
            assert core.digest(HERE/name) == sha, name
    source = HERE / NAMES[0]
    ast.parse(source.read_text(encoding='utf-8'))
    spec = importlib.util.spec_from_file_location('continuous_record_response_audit_checked', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    reproduced = module.run()
    saved = core.read(HERE / NAMES[1])
    assert saved == json.loads(json.dumps(reproduced))
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    assert (saved['round'], saved['baseline_round']) == (458, 457)
    checked_text = core.text_checks(HERE / NAMES[2])
    assert checked_text['display_formulas'] == 13
    links = 0
    for link in core.link_parser()((HERE / NAMES[2]).read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    return dict(date='2026-09-24', round=458,
        scientific_base_through_round=457, frozen_integration_base_through_round=457,
        fresh_tests=dict(run=6, failures=0, errors=0), saved_results_reproduced=True,
        scientific_results_rewritten=False, previous_scientific_file_hashes_verified=682,
        previous_protected_evidence_hashes_verified=717, text_checks=checked_text,
        local_links_checked=links, broken_links=0,
        new_file_hashes={name: core.digest(HERE / name) for name in NAMES},
        preserved_draft_hashes={'research_note_456_text_v2.md': core.digest(HERE/'research_note_456_text_v2.md')}, visual_rendering_performed=False, legacy_science_tests_rerun=False,
        independent_review='parent analytic review plus independent agent review of the physical memory sector, reference-complete response channel, all-state activation bound and rational simultaneous response certificate',
        scope=saved['scope'], all_reported_checks_passed=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
