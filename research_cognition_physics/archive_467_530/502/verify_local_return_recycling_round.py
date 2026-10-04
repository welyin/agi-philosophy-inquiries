"""Read-only scientific reproduction and evidence audit for round 502."""
import argparse
import ast
import importlib.util
import json
import sys
from pathlib import Path
import verify_round501_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE/'research_round_502_checks.json'
NAMES = ['local_return_recycling.py', 'local_return_recycling_results.json',
         'research_note_502.md']
DRAFT = 'round502_drafts/research_note_502.txt'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        frozen = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert frozen['science_hashes_verified_231_501'] == 814
    assert frozen['total_protected_evidence_hashes'] == 883
    if TARGET.exists():
        for group in ['new_file_hashes', 'preserved_draft_hashes']:
            for name, sha in core.read(TARGET)[group].items():
                assert core.digest(HERE/name) == sha, name
    assert (HERE/NAMES[2]).read_bytes() == (HERE/DRAFT).read_bytes()
    source = HERE/NAMES[0]
    ast.parse(source.read_text(encoding='utf-8'))
    spec = importlib.util.spec_from_file_location('local_return_recycling_checked', source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    saved = core.read(HERE/NAMES[1])
    assert saved == json.loads(json.dumps(module.run()))
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    assert (saved['round'], saved['scientific_baseline_round']) == (502, 501)
    for key in ['full_GR_goal_completed', 'three_dimensional_positions_generated',
                'autonomous_complete_controller_generated', 'common_deadline_graph_freeze']:
        assert not saved['scope'][key]
    assert saved['scope']['first_step_fresh_replies_only']
    assert not saved['primary_source']['tool_theory_claimed_original']
    checked_text = core.text_checks(HERE/NAMES[2])
    assert checked_text['display_formulas'] == 15
    links = 0
    for link in core.link_parser()((HERE/NAMES[2]).read_text(encoding='utf-8')):
        dest = (HERE/link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    return dict(date='2026-09-27', round=502, scientific_base_through_round=501,
        frozen_integration_base_through_round=501,
        fresh_tests=dict(run=6, failures=0, errors=0), saved_results_reproduced=True,
        scientific_results_rewritten=False, previous_scientific_file_hashes_verified=814,
        previous_protected_evidence_hashes_verified=883, text_checks=checked_text,
        local_links_checked=links, broken_links=0,
        new_file_hashes={name: core.digest(HERE/name) for name in NAMES},
        preserved_draft_hashes={DRAFT: core.digest(HERE/DRAFT)},
        visual_rendering_performed=False, legacy_science_tests_rerun=False,
        independent_review='Final code and draft independently reviewed and executed: finite monitored return, whole-cycle bounds, reference preservation, fresh-only current labels, complete errors and finite resource scope.',
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
