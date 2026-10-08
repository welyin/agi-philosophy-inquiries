"""Reproduce round 1009; check frozen inputs and current round links."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import re
import sys

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
ROOT=BASE.parent
OUT=HERE/'research_round_1009_checks.json'
sys.path.insert(0,str(ROOT/'scripts'))
from organize_research_231_775 import links
import remaining_freedom


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(writing):
    result=remaining_freedom.run()
    saved=json.loads((HERE/'remaining_freedom_results.json').read_text('utf8'))
    remaining_freedom.compare(result,saved)
    for name,value in saved['source_sha256'].items():
        assert sha(ROOT/name)==value,name
    docs=list(HERE.parent.rglob('*.md'))
    checked=0
    for path in docs:
        for _,_,target,local in links(path.read_text('utf-8-sig')):
            resolved=(path.parent/local.replace('\\','/')).resolve()
            assert resolved.exists() or (writing and resolved==OUT),(path,target)
            checked+=1
    ledger=(HERE/'input_dependency_ledger_v0_1.md').read_text('utf8')
    categories=re.findall(r'^\|(I\d\d) ',ledger,re.M)
    assert categories==[f'I{i:02d}' for i in range(1,15)],categories
    goal=json.loads((HERE/'goal_start.json').read_text('utf8'))
    assert goal['goal']['status']=='active'
    assert '生成性理论' in goal['goal']['objective']
    for name in ('README.md','research_direction.md','RESEARCH_STATE.md'):
        prose=(BASE/name).read_text('utf-8-sig')
        assert 'archive_1009_' in prose[:4000] and '1009' in prose[:4000]
    # Historical migration receipts describe their publication-time navigation;
    # they are not rewritten to pretend that live navigation is frozen forever.
    receipt=BASE/'_migration/layout_764_1008_20261008/replay_checks.json'
    prior=json.loads(receipt.read_text('utf8'))
    for name,value in prior['verified_runtime_sha256'].items():
        assert sha(ROOT/name)==value,name
    evidence=[p for p in HERE.parent.rglob('*') if p.is_file() and p!=OUT]
    return dict(round=1009,date='2026-10-08',all_delivery_checks_passed=True,
        scientific_result_reproduced=True,input_categories=14,new_stage_local_links_checked=checked,
        fresh_model_calibration_groups=1,cumulative_test_groups=3787,
        archived_scientific_inputs_unchanged=True,previous_migration_manifest_unchanged=True,
        live_navigation_updated=True,old_publication_navigation_reverified=False,
        source_sha256={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in evidence},
        goal_status='active',goal_completed=False,new_app_task_created=False,
        visual_checks_performed=False,proof_scope='analytic note; numerical checks calibrate the stated finite material family')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();result=run(args.write)
    if args.write:
        with OUT.open('x',encoding='utf8') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2);stream.write('\n')
    else:
        assert result==json.loads(OUT.read_text('utf8'))
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},ensure_ascii=False,indent=2))
