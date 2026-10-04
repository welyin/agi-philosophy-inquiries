"""Read-only science/evidence verification for round 487."""
import argparse
import ast
import importlib.util
import json
import sys
from pathlib import Path

import verify_round485_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'research_round_487_checks.json'
NAMES = ['zero_velocity_acceleration_audit.py',
         'zero_velocity_acceleration_audit_results.json', 'research_note_487.md']


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        frozen = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert frozen['science_hashes_verified_231_485'] == 766
    assert frozen['total_protected_evidence_hashes'] == 811
    if TARGET.exists():
        for group in ['new_file_hashes', 'preserved_draft_hashes']:
            for name, sha in core.read(TARGET)[group].items():
                assert core.digest(HERE / name) == sha, name
    assert (HERE / NAMES[2]).read_bytes() == (HERE / 'round487_drafts/research_note_487.txt').read_bytes()
    source = HERE / NAMES[0]
    ast.parse(source.read_text(encoding='utf-8'))
    spec = importlib.util.spec_from_file_location('zero_velocity_acceleration_audit_checked', source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    reproduced = module.run()
    saved = core.read(HERE / NAMES[1])
    assert saved == json.loads(json.dumps(reproduced))
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    assert (saved['round'], saved['baseline_round']) == (487, 485)
    assert not saved['scope']['full_GR_goal_completed']
    assert not saved['scope']['phase_closure_triggered']
    checked_text = core.text_checks(HERE / NAMES[2])
    assert checked_text['display_formulas'] == 16
    links = 0
    for link in core.link_parser()((HERE / NAMES[2]).read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    return dict(date='2026-09-27', round=487,
        scientific_base_through_round=485, frozen_integration_base_through_round=485,
        fresh_tests=dict(run=6, failures=0, errors=0), saved_results_reproduced=True,
        scientific_results_rewritten=False, previous_scientific_file_hashes_verified=766,
        previous_protected_evidence_hashes_verified=811, text_checks=checked_text,
        local_links_checked=links, broken_links=0,
        new_file_hashes={name: core.digest(HERE / name) for name in NAMES},
        preserved_draft_hashes={name: core.digest(HERE / name) for name in ['round487_drafts/research_note_487.txt']}, visual_rendering_performed=False,
        legacy_science_tests_rerun=False,
        independent_review='Parent and independent reviewer completed the 16-formula whole-family dichotomy, exact acceleration moments, all-control bound, frozen 479 degenerate-plane application, initial-noise scope, actual readout, reproduction and final hashes before first save',
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
