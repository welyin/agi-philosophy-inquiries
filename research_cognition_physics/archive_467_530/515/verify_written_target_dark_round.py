"""Read-only evidence, reproduction and scope checks for round 515."""
import argparse
import ast
import importlib.util
import json
from pathlib import Path
import sys
sys.setrecursionlimit(max(sys.getrecursionlimit(), 10000))
import verify_round514_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE/'research_round_515_checks.json'
NAMES = ['written_target_dark_obstruction.py', 'written_target_dark_obstruction_results.json',
         'research_note_515.md']
DRAFTS = ['round515_drafts/research_note_515.txt',
          'round515_drafts/written_target_dark_obstruction_before_source_phase.txt']


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        base = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert base['stage_saved_tests'] == 2519
    assert base['science_hashes_verified_231_514'] == 853
    assert base['total_protected_evidence_hashes'] == 960
    additional = core.read(HERE/'record_monitoring_scope_checks.json')
    hashes = additional['new_unnumbered_evidence_hashes']
    assert len(hashes) == 5
    assert additional['total_protected_including_this_review'] == 965
    for name, sha in hashes.items():
        assert core.digest(HERE/name) == sha, name
    if TARGET.exists():
        for group in ['new_file_hashes','preserved_draft_hashes','preserved_additional_hashes']:
            for name, sha in core.read(TARGET)[group].items():
                assert core.digest(HERE/name) == sha, name
    assert (HERE/NAMES[2]).read_bytes() == (HERE/DRAFTS[0]).read_bytes()
    for name in [NAMES[0], Path(__file__).name]:
        ast.parse((HERE/name).read_text(encoding='utf8'))
    spec = importlib.util.spec_from_file_location('written_target_dark_checked', HERE/NAMES[0])
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    saved = core.read(HERE/NAMES[1])
    assert saved == json.loads(json.dumps(module.run()))
    assert (saved['round'], saved['scientific_base_through_round']) == (515, 514)
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    for name, sha in saved['dependency_sha256'].items():
        assert core.digest(HERE/name) == sha, name
    true = ['same_H_as_actual_record_formation', 'actual_514_source_with_strictly_positive_probability',
            'target_only_arbitrary_schedule_obstruction_proved']
    false = ['all_sources_impossible', 'arbitrary_active_routing_excluded',
             'one_failed_spatial_model_refutes_cognitive_axioms', 'autonomous_readout_or_timing_generated',
             'physical_position_or_dimension_generated', 'full_GR_goal_completed']
    assert set(saved['scope']) == set(true+false)
    assert all(saved['scope'][k] is True for k in true)
    assert all(saved['scope'][k] is False for k in false)
    text_checked = core.text_checks(HERE/NAMES[2])
    assert text_checked['display_formulas'] == 13
    links = 0
    for link in core.link_parser()((HERE/NAMES[2]).read_text(encoding='utf8')):
        assert (HERE/link).resolve().exists(), link
        links += 1
    return dict(date='2026-09-28', round=515, scientific_base_through_round=514,
        frozen_integration_base_through_round=514, fresh_tests=dict(run=6, failures=0, errors=0),
        saved_results_reproduced=True, scientific_results_rewritten=False,
        previous_scientific_file_hashes_verified=853,
        previous_protected_evidence_hashes_verified=960+len(hashes),
        text_checks=text_checked, local_links_checked=links, broken_links=0,
        new_file_hashes={name: core.digest(HERE/name) for name in NAMES},
        preserved_draft_hashes={name: core.digest(HERE/name) for name in DRAFTS},
        preserved_additional_hashes=hashes,
        visual_rendering_performed=False, legacy_science_tests_rerun=False,
        independent_review='Independent integer dark projector, actual 514 source, rational finite-time bounds, instrument scope, code and final draft review.',
        scope=saved['scope'], all_reported_checks_passed=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x', encoding='utf8', newline='\n') as f:
            f.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ('round','fresh_tests','all_reported_checks_passed')}, ensure_ascii=False))
