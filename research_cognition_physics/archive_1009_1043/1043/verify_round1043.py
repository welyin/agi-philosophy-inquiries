"""Verify the stage synthesis, old provenance and bounded regression evidence.

Does not infer mathematical truth from labels. Scientific composition and goal
scope are separately reviewed. Default is read-only; --write creates only the
first final receipt. --check permits checking before that receipt is created.
"""
from pathlib import Path
import argparse
import json
import os
import re
import runpy
import subprocess
import sys

import build_chain_manifest as manifest

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
BASE = STAGE.parent
ROOT = BASE.parent
OUT = HERE / 'research_round_1043_checks.json'
PRIOR = STAGE / '1042/mainline_acceptance.json'
REVIEWS = ['independent_goal_scope_review.md', 'independent_scientific_review.md',
           'independent_replay_review.md']
OWN = ['../research_note_1043.md', 'README.md', 'NEXT.md', 'completion_audit.md',
       'build_chain_manifest.py', 'chain_manifest.json', 'verify_round1043.py',
       'replay_stage.py', 'stage_replay_results.json',
       'legacy_1009_scope_check.py', 'legacy_1009_scope_results.json',
       'review_science.py', 'review_science_results.json', *REVIEWS]
NAV = [BASE/'README.md', BASE/'research_direction.md', BASE/'RESEARCH_STATE.md',
       STAGE/'README.md', STAGE/'文件索引.md']
sha, read, rel = manifest.sha, manifest.read, manifest.rel


def replay_local(script):
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1', PYTHONUTF8='1')
    run = subprocess.run([sys.executable, '-B', '-X', 'utf8', str(HERE/script)],
                         cwd=HERE, env=env, capture_output=True, text=True,
                         encoding='utf8', timeout=120)
    assert run.returncode == 0, (script, run.stdout[-3000:], run.stderr[-3000:])


def evidence():
    index = read(HERE/'chain_manifest.json')
    assert index == manifest.build(), 'Cited synthesis or historical source changed'
    old_module = runpy.run_path(str(PRIOR.parent/'verify_mainline_acceptance.py'))
    old_receipt = read(PRIOR)
    old_receipt.pop('navigation_links_at_acceptance')
    assert old_receipt == old_module['evidence'](), 'Previous frozen mainline changed'
    raw = read(HERE/'stage_replay_results.json')
    assert raw['rounds'] == list(range(1009, 1043))
    assert raw['passed_rounds'] == list(range(1010, 1043))
    assert raw['failed_rounds'] == [1009] and not raw['all_regressions_passed']
    assert raw['historical_files_changed'] == []
    assert raw['historical_snapshot_sha256_before'] == raw['historical_snapshot_sha256_after']
    replay = runpy.run_path(str(HERE/'replay_stage.py'))
    assert replay['manifest_digest'](replay['historical_snapshot']()) == raw['historical_snapshot_sha256_after']
    assert sha(HERE/'replay_stage.py') == raw['replay_script_sha256']
    for row in raw['results']:
        for name in ('verifier', 'receipt'):
            assert sha(ROOT/row[name]) == row[name+'_sha256_before'] == row[name+'_sha256_after']
        assert row['verifier_and_receipt_unchanged'] and not row['timeout']
        assert row['exit_code'] == (1 if row['round'] == 1009 else 0)
    replay_local('legacy_1009_scope_check.py')
    replay_local('review_science.py')
    diagnosis = read(HERE/'legacy_1009_scope_results.json')
    assert diagnosis['original_1009_default_verifier_passed'] is False
    assert diagnosis['run_false_full_scientific_and_source_checks_passed']
    assert diagnosis['all_seven_frozen_original_assets_unchanged']
    assert diagnosis['preserved_raw_stage_replay_sha256'] == sha(HERE/'stage_replay_results.json')
    audit = (HERE/'completion_audit.md').read_text('utf8')
    assert re.findall(r'^\|(R\d\d) ', audit, re.M) == [f'R{i:02d}' for i in range(1, 17)]
    # Reviews are frozen as evidence, not converted into truth by string tests.
    for filename in OWN:
        assert (HERE/filename).is_file(), filename
    return {
        'schema': 'stage_synthesis_acceptance_v1', 'round': 1043, 'date': '2026-10-08',
        'stage_deliverables': 'conditional_generation_chain_and_stratified_remaining_freedom',
        'cumulative_before': 3819, 'cumulative_after': 3819,
        'new_scientific_calibration_groups': 0, 'new_adopted_cognitive_axioms': 0,
        'full_physics_or_cognitive_generation_claimed': False,
        'prior_acceptance': rel(PRIOR), 'prior_acceptance_sha256': sha(PRIOR),
        'resolved_dependency_rows': 34, 'input_categories': 14, 'original_phenomena': 16,
        'goal_requirements_reviewed': 16,
        'old_default_verifiers_passed': list(range(1010, 1043)),
        'legacy_default_failure_preserved': [1009],
        'legacy_science_and_frozen_assets_separately_reverified': [1009],
        'historical_snapshot_sha256': raw['historical_snapshot_sha256_after'],
        'old_science_replay_does_not_add_experiments': True,
        'independent_agent_reviews': [rel(HERE/p) for p in REVIEWS],
        'human_peer_review': False, 'machine_check_is_not_semantic_completion_proof': True,
        'application_goal_status_not_set_by_this_script': True,
        'source_index_sha256': sha(HERE/'chain_manifest.json'),
        'current_artifacts_sha256': {rel(HERE/p): sha(HERE/p) for p in OWN},
        'new_app_task_created': False, 'new_schedule_created': False, 'image_checks': False,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = evidence()
    link_module = runpy.run_path(str(STAGE/'1040/verify_parallel_acceptance.py'))
    check_links = link_module['check_links']
    check_links.__globals__['OUT'] = OUT
    documents = [HERE/p for p in OWN if p.endswith('.md')]
    document_links = sum(check_links(p, args.write or args.check) for p in documents)
    navigation_links = sum(check_links(p, args.write or args.check) for p in NAV)
    if args.write:
        with OUT.open('x', encoding='utf8') as stream:
            json.dump(dict(result, document_links_at_acceptance=document_links,
                           navigation_links_at_acceptance=navigation_links), stream,
                      ensure_ascii=False, indent=2, allow_nan=False)
            stream.write('\n')
    elif not args.check:
        saved = read(OUT)
        saved.pop('document_links_at_acceptance')
        saved.pop('navigation_links_at_acceptance')
        assert saved == result, 'Frozen synthesis evidence changed'
    print(json.dumps({'passed': True, 'round': 1043, 'cumulative': 3819,
                      'original_default_passes': 33, 'legacy_scope_reverified': 1,
                      'dependency_rows': 34, 'frozen_artifacts': len(OWN),
                      'document_links': document_links, 'navigation_links': navigation_links,
                      'full_physics_generated': False}, ensure_ascii=False))
