"""Read-only scientific reproduction and evidence audit for round 500."""
import argparse
import ast
import importlib.util
import json
import sys
from pathlib import Path
import verify_round499_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE/'research_round_500_checks.json'
NAMES = ['current_neighbor_detection.py', 'current_neighbor_detection_results.json',
         'research_note_500.md']
DRAFT = 'round500_drafts/research_note_500.txt'
DRAFTS = [DRAFT, 'round500_drafts/reviewed_before_link_fix.txt']
LITERATURE = ['correlation_geometry_bridge_checks.py',
              'correlation_geometry_bridge_results.json',
              'correlation_geometry_bridge_review.md',
              'correlation_geometry_bridge_audit.json']


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        frozen = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert frozen['science_hashes_verified_231_499'] == 808
    assert frozen['total_protected_evidence_hashes'] == 867
    if TARGET.exists():
        for group in ['new_file_hashes', 'preserved_draft_hashes',
                      'new_unnumbered_literature_file_hashes']:
            for name, sha in core.read(TARGET)[group].items():
                assert core.digest(HERE/name) == sha, name
    for name, sha in core.read(HERE/LITERATURE[-1])['file_sha256'].items():
        assert core.digest(HERE/name) == sha, name
    assert (HERE/NAMES[2]).read_bytes() == (HERE/DRAFT).read_bytes()
    source = HERE/NAMES[0]
    ast.parse(source.read_text(encoding='utf-8'))
    spec = importlib.util.spec_from_file_location('current_neighbor_detection_checked', source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    reproduced = module.run()
    saved = core.read(HERE/NAMES[1])
    assert saved == json.loads(json.dumps(reproduced))
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    assert (saved['round'], saved['scientific_baseline_round']) == (500, 499)
    assert not saved['scope']['full_GR_goal_completed']
    assert not saved['scope']['individual_label_posterior_guarantee']
    checked_text = core.text_checks(HERE/NAMES[2])
    assert checked_text['display_formulas'] == 19
    links = 0
    allowed_pending = {TARGET, HERE/'round500_integration_checks.json'}
    for link in core.link_parser()((HERE/NAMES[2]).read_text(encoding='utf-8')):
        dest = (HERE/link).resolve()
        assert dest.exists() or (pending and dest in allowed_pending), link
        links += 1
    return dict(date='2026-09-27', round=500, scientific_base_through_round=499,
        frozen_integration_base_through_round=499,
        fresh_tests=dict(run=6, failures=0, errors=0), saved_results_reproduced=True,
        scientific_results_rewritten=False, previous_scientific_file_hashes_verified=808,
        previous_protected_evidence_hashes_verified=867, text_checks=checked_text,
        local_links_checked=links, broken_links=0,
        new_file_hashes={name: core.digest(HERE/name) for name in NAMES},
        preserved_draft_hashes={name: core.digest(HERE/name) for name in DRAFTS},
        final_note_link_only_change_reviewed_by_parent=True,
        new_unnumbered_literature_file_hashes={name: core.digest(HERE/name) for name in LITERATURE},
        unnumbered_literature_checks_not_added_to_science_count=True,
        visual_rendering_performed=False, legacy_science_tests_rerun=False,
        independent_review='Parent and independent reviewer checked the actual final code and note: uniform current-relation cq bound, unknown root degree, stopped instrument, delayed fixed label, complete process errors, resources and scope.',
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
