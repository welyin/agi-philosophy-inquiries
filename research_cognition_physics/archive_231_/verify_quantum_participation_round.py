"""Read-only round 433 science and historical-evidence verification."""
import argparse
import ast
import importlib.util
import json
from pathlib import Path

import verify_round432_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'research_round_433_checks.json'
NAMES = ['quantum_participation_audit.py', 'quantum_participation_audit_results.json', 'research_note_433.md']


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        frozen = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert frozen['science_hashes_verified_231_432'] == 607
    assert frozen['total_protected_evidence_hashes'] == 642
    if TARGET.exists():
        for name, sha in core.read(TARGET)['new_file_hashes'].items():
            assert core.digest(HERE / name) == sha, name
    source = HERE / NAMES[0]
    ast.parse(source.read_text(encoding='utf-8'))
    spec = importlib.util.spec_from_file_location('quantum_participation_checked', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    reproduced = module.run()
    saved = core.read(HERE / NAMES[1])
    assert saved == json.loads(json.dumps(reproduced))
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    assert (saved['round'], saved['baseline_round']) == (433, 432)
    checked_text = core.text_checks(HERE / NAMES[2])
    assert checked_text['display_formulas'] == 10
    links = 0
    for link in core.link_parser()((HERE / NAMES[2]).read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    return dict(date='2026-09-24', round=433,
        scientific_base_through_round=432, frozen_integration_base_through_round=432,
        fresh_tests=dict(run=6, failures=0, errors=0), saved_results_reproduced=True,
        scientific_results_rewritten=False, previous_scientific_file_hashes_verified=607,
        previous_protected_evidence_hashes_verified=642, text_checks=checked_text,
        local_links_checked=links, broken_links=0,
        new_file_hashes={name: core.digest(HERE / name) for name in NAMES},
        preserved_draft_hashes={}, visual_rendering_performed=False, legacy_science_tests_rerun=False,
        independent_review='parent exact short-time operator identities, independent spectral response and new-primitive scope audit; no separate agent review claimed',
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
