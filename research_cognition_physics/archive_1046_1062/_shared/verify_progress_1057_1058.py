"""Verify accepted process rounds, pending drafts, and unchanged roadmap scope.

--freeze saves an exclusive snapshot. Default verifies frozen assets; later
navigation additions are allowed but original requirements remain protected.
"""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import re
import runpy
import subprocess
import sys

HERE=Path(__file__).resolve().parent
PHASE=HERE.parent
ROOT=PHASE.parent
RECEIPT=HERE/'progress_1057_1058_checks.json'
BASE=runpy.run_path(str(HERE/'verify_adoptions_after1055.py'))
ROADMAP_SHA='0572fd4469c75ffa86d4e08162ca11784dfdbaa983762c520076a5e152f10fa4'


def digest(path): return sha256(path.read_bytes()).hexdigest()


def local_links(path):
    body=path.read_text(encoding='utf-8')
    if path.parent==ROOT and path.name in BASE['TAILS']:
        body=body.split(BASE['TAILS'][path.name][0],1)[0]
    body=re.sub(r'\\\[[\s\S]*?\\\]','',body)
    body=re.sub(r'\$\$[\s\S]*?\$\$','',body)
    count=0
    for target in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)',body):
        target=target.strip().strip('<>')
        if target.startswith(('http://','https://','#')): continue
        assert (path.parent/target.split('#',1)[0]).is_file(),(path,target)
        count+=1
    return count


def run(freeze):
    accepted=[]; author_only=[]
    owned=['integration_1057_1058.md','verify_progress_1057_1058.py',
           '../_admission/after1057_selection/next_direction_scope.md',
           '../_admission/after1057_selection/independent_math_gate.md']
    for n in (1057,1058):
        folder=PHASE/str(n)
        proc=subprocess.run([sys.executable,'-B','-X','utf8',str(folder/f'verify_round{n}.py')],
                            capture_output=True,text=True,encoding='utf-8',check=True)
        assert json.loads(proc.stdout)['author_checks_passed']
        receipt=folder/f'research_round_{n}_checks.json'
        owned.append(f'../{n}/{receipt.name}')
        acceptance=folder/'mainline_acceptance.json'
        if not acceptance.exists():
            author_only.append(n); continue
        a=json.loads(acceptance.read_text(encoding='utf-8'))
        assert a['accepted'] and a['independent_review_status']=='passed'
        assert a['author_receipt_sha256']==digest(receipt)
        assert a['independent_review_sha256']==digest(folder/'independent_review.md')
        assert a['independent_checks_sha256']==digest(folder/'review_checks.json')
        assert not a['full_GR_or_roadmap_completed']
        accepted.append(n)
        owned += [f'../{n}/{x}' for x in ('mainline_acceptance.json','independent_review.md','review_checks.json')]
    assert accepted and accepted[0]==1057
    for name,(marker,expected) in BASE['TAILS'].items():
        raw=(ROOT/name).read_bytes()
        assert sha256(raw[raw.index(marker.encode('utf-8')):]).hexdigest()==expected
    raw=(ROOT/'ROADMAP.md').read_bytes(); marker='## M1 检验认知要求的独立选择力'.encode('utf-8')
    assert sha256(raw[raw.index(marker):]).hexdigest()==ROADMAP_SHA
    assets={name:digest(HERE/name) for name in owned}
    links=sum(local_links(HERE/name) for name in owned if name.endswith('.md'))
    navlinks=sum(local_links(p) for p in BASE['NAV'])
    record={'date':'2026-10-09','author_rounds_checked':[1057,1058],
            'accepted_rounds_in_batch':accepted,'pending_independent_review':author_only,
            'latest_accepted_round':max(accepted),'cumulative_scientific_calibrations':3832+len(accepted),
            'stage_groups':11+len(accepted),'new_cognitive_axioms':0,
            'assets_sha256':assets,'asset_links_checked':links,'all_checks_passed':True,
            'original_roadmap_requirements_sha256':ROADMAP_SHA,'whole_roadmap_complete':False,
            'user_continuation_authorized':True,'goal_tool_status_observed':'blocked',
            'assistant_goal_status_mutation':False}
    if freeze:
        record['navigation_snapshot']={'sha256':{str(p.relative_to(ROOT)):digest(p) for p in BASE['NAV']},
                                       'local_links_checked':navlinks}
        with RECEIPT.open('x',encoding='utf-8') as out:
            json.dump(record,out,ensure_ascii=False,indent=2); out.write('\n')
    elif RECEIPT.exists():
        saved=json.loads(RECEIPT.read_text(encoding='utf-8')); saved.pop('navigation_snapshot')
        assert saved==record,'Frozen batch changed; retain history and write a new acceptance batch.'
    else:
        raise RuntimeError('Run --freeze once after the review state is finalized.')
    return record,navlinks


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--freeze',action='store_true')
    args=parser.parse_args(); r,links=run(args.freeze)
    print(json.dumps({'all_checks_passed':True,'accepted_rounds':r['accepted_rounds_in_batch'],
                      'pending_review':r['pending_independent_review'],'assets':len(r['assets_sha256']),
                      'current_navigation_links':links,'whole_roadmap_complete':False},ensure_ascii=False))
