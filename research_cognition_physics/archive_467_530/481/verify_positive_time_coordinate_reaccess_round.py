"""Read-only science/evidence verification for round 481."""
import argparse
import ast
import importlib.util
import json
import sys
from pathlib import Path

import verify_round480_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'research_round_481_checks.json'
NAMES = ['positive_time_coordinate_reaccess.py',
         'positive_time_coordinate_reaccess_results.json', 'research_note_481.md']


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        frozen = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert frozen['science_hashes_verified_231_480'] == 751
    assert frozen['total_protected_evidence_hashes'] == 790
    if TARGET.exists():
        for group in ['new_file_hashes', 'preserved_draft_hashes']:
            for name, sha in core.read(TARGET)[group].items():
                assert core.digest(HERE / name) == sha, name
    assert (HERE / NAMES[2]).read_bytes() == (HERE / 'round481_drafts/research_note_481.txt').read_bytes()
    source = HERE / NAMES[0]
    ast.parse(source.read_text(encoding='utf-8'))
    spec = importlib.util.spec_from_file_location('positive_time_coordinate_reaccess_checked', source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    reproduced = module.run()
    saved = core.read(HERE / NAMES[1])
    assert saved == json.loads(json.dumps(reproduced))
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    assert (saved['round'], saved['baseline_round']) == (481, 480)
    assert not saved['scope']['full_GR_goal_completed']
    assert not saved['scope']['phase_closure_triggered']
    checked_text = core.text_checks(HERE / NAMES[2])
    assert checked_text['display_formulas'] >= 1
    links = 0
    for link in core.link_parser()((HERE / NAMES[2]).read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    return dict(date='2026-09-27', round=481,
        scientific_base_through_round=480, frozen_integration_base_through_round=480,
        fresh_tests=dict(run=6, failures=0, errors=0), saved_results_reproduced=True,
        scientific_results_rewritten=False, previous_scientific_file_hashes_verified=751,
        previous_protected_evidence_hashes_verified=790, text_checks=checked_text,
        local_links_checked=links, broken_links=0,
        new_file_hashes={name: core.digest(HERE / name) for name in NAMES},
        preserved_draft_hashes={'round481_drafts/research_note_481.txt': core.digest(HERE / 'round481_drafts/research_note_481.txt')}, visual_rendering_performed=False,
        legacy_science_tests_rerun=False,
        independent_review='parent and independent full-note/code review and actual reproduction; finite recurrence certificate, current-center Brouwer cover, complete-process accumulation and actual data-readout budgets',
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
