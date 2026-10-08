"""Verify the fixed 1010 derivation artifacts without freezing live navigation."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import sys

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
ROOT=BASE.parent
OUT=HERE/'research_round_1010_checks.json'
sys.path.insert(0,str(ROOT/'scripts'))
from organize_research_231_775 import links
import hypercharge_neutrino_selection as current


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(writing=False):
    fresh=current.main_result()
    saved=json.loads((HERE/'hypercharge_neutrino_selection_results.json').read_text('utf8'))
    assert fresh==saved
    for name,digest in saved['source_sha256'].items():
        assert sha(BASE/name)==digest,name
    previous=HERE.parent/'1009/remaining_freedom.py'
    spec=importlib.util.spec_from_file_location('r1009_frozen',previous)
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    old.compare(old.run(),json.loads((previous.parent/'remaining_freedom_results.json').read_text('utf8')))
    note=HERE.parent/'research_note_1010.md'
    docs=[note]+sorted(HERE.glob('*.md'))
    count=0
    for p in docs:
        for _,_,target,local in links(p.read_text('utf-8-sig')):
            path=(p.parent/local.replace('\\','/')).resolve()
            assert path.exists() or (writing and path==OUT),(p,target)
            count+=1
    assert not any(saved['claim_boundaries'].values())
    evidence=[note]+sorted(p for p in HERE.iterdir() if p.is_file() and p!=OUT)
    result=dict(round=1010,date='2026-10-08',all_delivery_checks_passed=True,
        saved_result_reproduced=True,previous_1009_core_reproduced=True,
        frozen_inputs_verified=len(saved['source_sha256']),round_local_links_checked=count,
        new_calibration_groups=1,cumulative_research_groups=3788,
        numerical_scope='exact rational anomaly checks and finite complex matrices; analytic completeness remains in note',
        source_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in evidence},
        live_navigation_part_of_frozen_receipt=False,goal_completed=False,
        new_cognitive_axiom=False,visual_checks_performed=False)
    if not writing:
        assert result==json.loads(OUT.read_text('utf8'))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();result=run(args.write)
    if args.write:
        with OUT.open('x',encoding='utf8') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2);stream.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},ensure_ascii=False,indent=2))
