"""Read-only science, certificate and scope audit for round 514."""
import argparse
import ast
import importlib.util
import json
import sys
from pathlib import Path
sys.setrecursionlimit(max(sys.getrecursionlimit(), 10000))
import verify_round513_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE/'research_round_514_checks.json'
NAMES = ['written_record_persistence.py', 'written_record_persistence_results.json',
         'research_note_514.md']
DRAFTS = ['round514_drafts/research_note_514.txt',
          'round514_drafts/written_record_persistence_before_scope_clarification.txt']


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        frozen = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    assert frozen['science_hashes_verified_231_513'] == 850
    assert frozen['total_protected_evidence_hashes'] == 953
    if TARGET.exists():
        for group in ['new_file_hashes', 'preserved_draft_hashes']:
            for name, sha in core.read(TARGET)[group].items():
                assert core.digest(HERE/name) == sha, name
    assert (HERE/NAMES[2]).read_bytes() == (HERE/DRAFTS[0]).read_bytes()
    assert core.digest(HERE/NAMES[0]) == '41d3db3b6a6ffb0aff85f6a29ef024da98b8a7814720449344d5d5da5c051351'
    assert core.digest(HERE/NAMES[2]) == '6a81b2e793ac33572828cdd0acc19f1d4a285dfa62cb4ee1a20b5faf26a11aec'
    source = HERE/NAMES[0]
    ast.parse(source.read_text(encoding='utf-8'))
    spec = importlib.util.spec_from_file_location('written_record_persistence_checked', source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    saved = core.read(HERE/NAMES[1])
    assert saved == json.loads(json.dumps(module.run()))
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    assert (saved['round'], saved['scientific_baseline_round']) == (514, 513)
    for name, sha in saved['dependency_sha256'].items():
        assert core.digest(HERE/name) == sha, name
    true_keys = ['same_H_for_formation_and_retention', 'true_remote_Z_record_source', 'size_uniform_source_and_retention_bounds', 'arbitrary_unknown_graph_reference_covered', 'all_failure_branches_and_conditioning_error_accounted', 'retention_clock_starts_at_source_instrument_end', 'unchanged_real_memory_coupling']
    false_keys = ['deterministic_or_high_probability_source', 'source_success_for_each_specified_remote_site', 'permanent_record_or_entire_member_table_protected', 'unmeasured_natural_history_proved', 'autonomous_readout_controller_generated', 'source_record_collection_deadline_proved', 'dimension_three_generated', 'full_GR_goal_completed', 'phase_closure_triggered']
    assert set(saved['scope']) == set(true_keys+false_keys)
    assert all(saved['scope'][k] is True for k in true_keys)
    assert all(saved['scope'][k] is False for k in false_keys)
    checked_text = core.text_checks(HERE/NAMES[2])
    assert checked_text['display_formulas'] == 16
    links = 0
    for link in core.link_parser()((HERE/NAMES[2]).read_text(encoding='utf-8')):
        assert (HERE/link).resolve().exists(), link
        links += 1
    return dict(date='2026-09-28', round=514, scientific_base_through_round=513,
        frozen_integration_base_through_round=513,
        fresh_tests=dict(run=6, failures=0, errors=0), saved_results_reproduced=True,
        scientific_results_rewritten=False, previous_scientific_file_hashes_verified=850,
        previous_protected_evidence_hashes_verified=953, text_checks=checked_text,
        local_links_checked=links, broken_links=0,
        new_file_hashes={name: core.digest(HERE/name) for name in NAMES},
        preserved_draft_hashes={name: core.digest(HERE/name) for name in DRAFTS},
        visual_rendering_performed=False, legacy_science_tests_rerun=False,
        independent_review='Independent size-uniform source, same-H retention, complete-instrument scope, code and final 16-equation draft review; read-only reproduction passed.',
        scope=saved['scope'], all_reported_checks_passed=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(dict(round=result['round'], fresh_tests=result['fresh_tests'],
        all_reported_checks_passed=result['all_reported_checks_passed']), ensure_ascii=False))
