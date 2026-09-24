"""Verify round 419 and the frozen evidence through round 418."""
import argparse
import ast
import importlib.util
import io
import json
from pathlib import Path
import unittest

import verify_round418_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'research_round_419_checks.json'
NAMES = ['internal_rate_window_audit.py', 'internal_rate_window_audit_results.json',
         'research_note_419.md']


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        frozen = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert frozen['science_hashes_verified_231_418'] == 565
    assert frozen['total_protected_evidence_hashes'] == 590
    if TARGET.exists():
        for name, sha in core.read(TARGET)['new_file_hashes'].items():
            assert core.digest(HERE / name) == sha, name
    source = HERE / NAMES[0]
    ast.parse(source.read_text(encoding='utf-8'))
    spec = importlib.util.spec_from_file_location('internal_rate_checked', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    saved = core.read(HERE / NAMES[1])
    assert {k:v for k,v in saved.items() if k not in ('checks','runtime')} == json.loads(json.dumps(module.report()))
    output = io.StringIO()
    tests = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromModule(module))
    assert tests.wasSuccessful() and tests.testsRun == 8, output.getvalue()
    assert saved['checks'] == dict(run=8,failures=0,errors=0)
    checked = core.text_checks(HERE / NAMES[2])
    assert checked['display_formulas'] == 9
    links = 0
    for link in core.link_parser()((HERE/NAMES[2]).read_text(encoding='utf-8')):
        dest = (HERE/link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    return dict(date='2026-09-24',round=419,scientific_base_through_round=418,
        additional_frozen_dependency_rounds=[230,346,399,401,409,418],
        fresh_tests=dict(run=8,failures=0,errors=0),saved_results_reproduced=True,
        scientific_results_rewritten=False,previous_scientific_file_hashes_verified=565,
        previous_protected_evidence_hashes_verified=590,text_checks=checked,
        local_links_checked=links,broken_links=0,
        new_file_hashes={name:core.digest(HERE/name) for name in NAMES},
        visual_rendering_performed=False,legacy_science_tests_rerun=False,
        independent_review='round419_completion_candidate reviewed actual note and code and independently ran eight checks',
        internal_random_rate_with_uniform_finite_window=True,
        autonomous_stopping_or_online_admission_proved=False,
        full_cognitive_countermodel_completed=False,phase_closure_triggered=False,
        scope=saved['scope'],all_reported_checks_passed=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args = parser.parse_args()
    result = verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
