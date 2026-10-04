"""Read-only scientific reproduction and evidence audit for round 501."""
import argparse
import ast
import importlib.util
import json
import sys
from pathlib import Path
import verify_round500_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE/'research_round_501_checks.json'
NAMES = ['local_neighbor_receipt.py', 'local_neighbor_receipt_results.json',
         'research_note_501.md']
DRAFT = 'round501_drafts/research_note_501.txt'
CANDIDATE = ['neighbor_return_candidate_drafts/'+name for name in
             ['README.md', 'local_return_probe_candidate.py', 'local_return_probe_candidate_results.json']]


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        frozen = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert frozen['science_hashes_verified_231_500'] == 811
    assert frozen['total_protected_evidence_hashes'] == 876
    if TARGET.exists():
        for group in ['new_file_hashes', 'preserved_draft_hashes', 'preserved_candidate_hashes']:
            for name, sha in core.read(TARGET)[group].items():
                assert core.digest(HERE/name) == sha, name
    assert (HERE/NAMES[2]).read_bytes() == (HERE/DRAFT).read_bytes()
    source = HERE/NAMES[0]
    ast.parse(source.read_text(encoding='utf-8'))
    spec = importlib.util.spec_from_file_location('local_neighbor_receipt_checked', source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    saved = core.read(HERE/NAMES[1])
    assert saved == json.loads(json.dumps(module.run()))
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    assert (saved['round'], saved['scientific_baseline_round']) == (501, 500)
    assert not saved['scope']['full_GR_goal_completed']
    assert not saved['scope']['per_label_posterior_bound']
    assert not saved['scope']['original_qubit_H_implemented']
    checked_text = core.text_checks(HERE/NAMES[2])
    assert checked_text['display_formulas'] == 18
    links = 0
    for link in core.link_parser()((HERE/NAMES[2]).read_text(encoding='utf-8')):
        dest = (HERE/link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    return dict(date='2026-09-27', round=501, scientific_base_through_round=500,
        frozen_integration_base_through_round=500,
        fresh_tests=dict(run=6, failures=0, errors=0), saved_results_reproduced=True,
        scientific_results_rewritten=False, previous_scientific_file_hashes_verified=811,
        previous_protected_evidence_hashes_verified=876, text_checks=checked_text,
        local_links_checked=links, broken_links=0,
        new_file_hashes={name: core.digest(HERE/name) for name in NAMES},
        preserved_draft_hashes={DRAFT: core.digest(HERE/DRAFT)},
        preserved_candidate_hashes={name: core.digest(HERE/name) for name in CANDIDATE},
        visual_rendering_performed=False, legacy_science_tests_rerun=False,
        independent_review='Actual final code and note independently reviewed and executed: current-relation bound, root-only complete instrument, rare-success error, no-reset counterexample, full carrier and resource scope.',
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
