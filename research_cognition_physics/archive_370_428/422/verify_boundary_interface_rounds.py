"""Read-only verification of independent rounds 421 and 422, based on frozen 420."""
import argparse
import ast
import importlib.util
import io
import json
from pathlib import Path
import unittest

import verify_round420_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
CONFIG = {421: ('boundary_interaction_energy_audit', 'Audit', 8, 10),
          422: ('translation_clock_interface_audit', 'TranslationClockTests', 9, 18)}


def verify(number, pending=False):
    target = HERE / f'research_round_{number}_checks.json'
    old_target = previous.TARGET
    previous.TARGET = target
    try:
        frozen = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert frozen['science_hashes_verified_231_420'] == 571
    assert frozen['total_protected_evidence_hashes'] == 601
    stem, test_class, count, formulas = CONFIG[number]
    names = [stem + '.py', stem + '_results.json', f'research_note_{number}.md']
    drafts = ['round422_drafts/' + n for n in names] if number == 422 else []
    if target.exists():
        checked = core.read(target)
        for key in ('new_file_hashes', 'preserved_draft_hashes'):
            for name, sha in checked[key].items():
                assert core.digest(HERE / name) == sha, name
    source = HERE / names[0]
    ast.parse(source.read_text(encoding='utf-8'))
    spec = importlib.util.spec_from_file_location(stem + '_checked', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    output = io.StringIO()
    tests = unittest.TextTestRunner(stream=output).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(getattr(module, test_class)))
    assert tests.wasSuccessful() and tests.testsRun == count, output.getvalue()
    saved = core.read(HERE / names[1])
    assert saved['observations'] == json.loads(json.dumps(module.OBS))
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (count, 0, 0)
    assert (saved['round'], saved['baseline_round']) == (number, 420)
    checked_text = core.text_checks(HERE / names[2])
    assert checked_text['display_formulas'] == formulas
    links = 0
    for link in core.link_parser()((HERE / names[2]).read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == target), link
        links += 1
    return dict(date='2026-09-24', round=number,
        scientific_base_through_round=420, frozen_integration_base_through_round=420,
        fresh_tests=dict(run=count, failures=0, errors=0), saved_results_reproduced=True,
        scientific_results_rewritten=False, previous_scientific_file_hashes_verified=571,
        previous_protected_evidence_hashes_verified=601, text_checks=checked_text,
        local_links_checked=links, broken_links=0,
        new_file_hashes={name: core.digest(HERE / name) for name in names},
        preserved_draft_hashes={name: core.digest(HERE / name) for name in drafts},
        visual_rendering_performed=False, legacy_science_tests_rerun=False,
        independent_review='round419_completion_candidate reviewed 421; parent independently reviewed and ran 422',
        scope=saved['scope'], all_reported_checks_passed=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round', type=int, choices=CONFIG)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.round, args.write_checks)
    if args.write_checks:
        with (HERE / f'research_round_{args.round}_checks.json').open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
