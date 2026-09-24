"""Read-only science verification of round 423, retaining all frozen evidence."""
import argparse
import ast
import importlib.util
import io
import json
from pathlib import Path
import unittest

import verify_round421_422_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'research_round_423_checks.json'
NAMES = ['finite_sampling_source_audit.py', 'finite_sampling_source_audit_results.json',
         'research_note_423.md']
DRAFTS = ['round423_drafts/finite_sampling_source_audit_before_local_window.py',
          'round423_drafts/finite_sampling_source_audit_results_before_local_window.json']


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        frozen = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert frozen['science_hashes_verified_231_422'] == 577
    assert frozen['total_protected_evidence_hashes'] == 610
    if TARGET.exists():
        checked = core.read(TARGET)
        for key in ('new_file_hashes', 'preserved_draft_hashes'):
            for name, sha in checked[key].items():
                assert core.digest(HERE / name) == sha, name
    source = HERE / NAMES[0]
    ast.parse(source.read_text(encoding='utf-8'))
    spec = importlib.util.spec_from_file_location('sampling_checked', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    output = io.StringIO()
    tests = unittest.TextTestRunner(stream=output).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(module.Audit))
    assert tests.wasSuccessful() and tests.testsRun == 9, output.getvalue()
    saved = core.read(HERE / NAMES[1])
    assert saved['observations'] == json.loads(json.dumps(module.OBS))
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (9, 0, 0)
    assert (saved['round'], saved['baseline_round']) == (423, 422)
    checked_text = core.text_checks(HERE / NAMES[2])
    assert checked_text['display_formulas'] == 14
    links = 0
    for link in core.link_parser()((HERE / NAMES[2]).read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    return dict(date='2026-09-24', round=423,
        scientific_base_through_round=422, frozen_integration_base_through_round=422,
        fresh_tests=dict(run=9, failures=0, errors=0), saved_results_reproduced=True,
        scientific_results_rewritten=False, previous_scientific_file_hashes_verified=577,
        previous_protected_evidence_hashes_verified=610, text_checks=checked_text,
        local_links_checked=links, broken_links=0,
        new_file_hashes={name: core.digest(HERE / name) for name in NAMES},
        preserved_draft_hashes={name: core.digest(HERE / name) for name in DRAFTS},
        visual_rendering_performed=False, legacy_science_tests_rerun=False,
        independent_review='round419_completion_candidate final read-only review; formal_completion_scope literature and dedup audit',
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
