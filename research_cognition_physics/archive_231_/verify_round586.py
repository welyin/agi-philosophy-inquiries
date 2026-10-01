"""Verify round 586 with frozen historical evidence and fresh diagnostics."""
import argparse
import ast
import json
from pathlib import Path
import joint_record_source_compression as model
import verify_interaction_rounds as core
import verify_round585 as previous
import cognition_forward_bridge_review_585_audit as historical_review

HERE = Path(__file__).resolve().parent
TARGET = HERE / "research_round_586_checks.json"


def verify():
    base = previous.verify()
    assert (base['round'], base['cumulative_numbered_tests'],
            base['cumulative_numbered_scientific_files'],
            base['cumulative_unique_protected_evidence_files']) == (585, 2984, 1067, 1630)
    assert historical_review.verify()['all_documentation_checks_passed']
    result = core.read(model.TARGET)
    assert result == model.run()
    assert (result['tests_run'], result['failures'], result['errors']) == (5, 0, 0)
    for name, digest in result['dependency_hashes'].items():
        assert core.digest(HERE/name) == digest, name
    new = {name: core.digest(HERE/name) for name in (
        'research_note_586.md', 'joint_record_source_compression.py',
        'joint_record_source_compression_results.json')}
    preserved = {name: core.digest(HERE/name) for name in (
        'unified_physics_condition_ledger_586.md', 'round586_drafts/STATUS.md',
        'round586_drafts/entry_scope_records.md', 'round586_drafts/research_note_586_draft.md',
        'round586_drafts/final_review.txt', 'cognition_forward_bridge_review_585.md',
        'cognition_forward_bridge_review_585_audit.py', 'cognition_forward_bridge_review_585_inventory.json',
        'cognition_forward_bridge_review_585_checks.json')}
    assert len(new) == 3 and len(preserved) == 9
    assert (HERE/'research_note_586.md').read_bytes() == (HERE/'round586_drafts/research_note_586_draft.md').read_bytes()
    review = (HERE/'round586_drafts/final_review.txt').read_text('utf8')
    for name in ('research_note_586.md', 'joint_record_source_compression.py',
                 'joint_record_source_compression_results.json', 'unified_physics_condition_ledger_586.md'):
        assert core.digest(HERE/name) in review, name
    text = core.text_checks(HERE/'research_note_586.md')
    assert text['display_formulas'] == 14
    links = 0
    for name in ('research_note_586.md', 'unified_physics_condition_ledger_586.md'):
        for link in core.link_parser()((HERE/name).read_text('utf8')):
            p = (HERE/link).resolve()
            assert p.exists() or p == TARGET.resolve(), (name, link)
            links += 1
    for name in ('joint_record_source_compression.py', Path(__file__).name):
        ast.parse((HERE/name).read_text('utf8'))
    if TARGET.exists():
        old = core.read(TARGET)
        assert old['new_file_hashes'] == new and old['preserved_draft_hashes'] == preserved
    return dict(date='2026-10-01', round=586, scientific_base_through_round=585,
                fresh_tests=dict(run=5, failures=0, errors=0), cumulative_numbered_tests=2989,
                cumulative_numbered_scientific_files=1070, unchanged_prior_evidence_files=1630,
                cumulative_unique_protected_evidence_files=1642,
                new_file_hashes=new, preserved_draft_hashes=preserved, text_checks=text,
                local_links_checked=links, broken_links=0, saved_results_reproduced=True,
                previous_results_unchanged=True, full_historical_science_rerun=False,
                direct_round584_585_dependencies_reproduced=True,
                historical_review_and_manifest_verified=True,
                inherited_theorems_not_counted_as_new=True,
                primary_code_and_note_review_completed=True,
                independent_final_code_and_draft_review_completed=False,
                visual_checks_performed=False, active_goal_unchanged=True,
                scope=result['scope'], all_reported_checks_passed=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args(); result = verify()
    if args.write_checks:
        with TARGET.open('x', encoding='utf8', newline='\n') as f:
            f.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    else:
        assert core.read(TARGET) == result
    print(json.dumps({k: result[k] for k in ('round', 'cumulative_numbered_tests',
          'cumulative_unique_protected_evidence_files', 'all_reported_checks_passed')}))
