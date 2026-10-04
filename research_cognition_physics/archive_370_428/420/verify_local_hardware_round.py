"""Verify round 420; preserve all earlier science and the five draft files."""
import argparse
import ast
import importlib.util
import io
import json
from pathlib import Path
import unittest

import verify_round419_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'research_round_420_checks.json'
NAMES = ['local_hardware_family_audit.py', 'local_hardware_family_audit_results.json',
         'research_note_420.md']
DRAFTS = ['round420_drafts/' + name for name in (
    'local_hardware_family_audit.py', 'local_hardware_family_audit_results.json',
    'local_hardware_family_audit_results_v2.json',
    'local_hardware_family_audit_results_v3.json', 'research_note_420.md')]


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        frozen = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert frozen['science_hashes_verified_231_419'] == 568
    assert frozen['total_protected_evidence_hashes'] == 593
    if TARGET.exists():
        checked = core.read(TARGET)
        for key in ('new_file_hashes', 'preserved_draft_hashes'):
            for name, sha in checked[key].items():
                assert core.digest(HERE / name) == sha, name
    source = HERE / NAMES[0]
    ast.parse(source.read_text(encoding='utf-8'))
    spec = importlib.util.spec_from_file_location('local_hardware_checked', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    output = io.StringIO()
    tests = unittest.TextTestRunner(stream=output).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(module.Audit))
    assert tests.wasSuccessful() and tests.testsRun == 11, output.getvalue()
    saved = core.read(HERE / NAMES[1])
    assert saved['observations'] == json.loads(json.dumps(module.OBS))
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (11, 0, 0)
    assert (saved['round'], saved['baseline_round']) == (420, 418)
    checked_text = core.text_checks(HERE / NAMES[2])
    assert checked_text['display_formulas'] == 12
    links = 0
    for link in core.link_parser()((HERE / NAMES[2]).read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    return dict(date='2026-09-24', round=420, scientific_base_through_round=418,
        frozen_integration_base_through_round=419,
        additional_frozen_dependency_rounds=[222, 229, 230, 346, 348, 406, 409, 418],
        fresh_tests=dict(run=11, failures=0, errors=0), saved_results_reproduced=True,
        scientific_results_rewritten=False, previous_scientific_file_hashes_verified=568,
        previous_protected_evidence_hashes_verified=593, text_checks=checked_text,
        local_links_checked=links, broken_links=0,
        new_file_hashes={name: core.digest(HERE / name) for name in NAMES},
        preserved_draft_hashes={name: core.digest(HERE / name) for name in DRAFTS},
        visual_rendering_performed=False, legacy_science_tests_rerun=False,
        independent_review='round419_completion_candidate reviewed the final eleven-check implementation and scope',
        complete_future_process_distinguished_from_endpoint=True,
        scope=saved['scope'], all_reported_checks_passed=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
