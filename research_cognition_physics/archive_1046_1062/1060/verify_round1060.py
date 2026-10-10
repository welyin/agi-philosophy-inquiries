"""Freeze/recheck round1060 evidence. No edits to frozen scientific files."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import runpy
import subprocess
import sys
from urllib.parse import unquote

HERE=Path(__file__).resolve().parent
PHASE=HERE.parent
ROOT=PHASE.parent
SHARED=PHASE/'_shared'
RECEIPT=HERE/'research_round_1060_checks.json'
ACCEPT=HERE/'mainline_acceptance.json'
BASE=runpy.run_path(str(SHARED/'verify_adoptions_after1058.py'))
NAV=BASE['BASE']['NAV']
AUTHOR=['../research_note_1060.md','proof.md','selection.md','sources.md',
        'dependency_update.md','global_anomaly_pullback.py','results.json']
REVIEW=['independent_review.md','independent_check.py','independent_checks.json']
HISTORY=['archive_585_628/research_note_614.md',
         'archive_585_628/614/joint_spinor_subgroup_mass.py',
         'archive_585_628/research_note_627.md','archive_585_628/research_note_628.md',
         'archive_1009_1043/research_note_1032.md',
         'archive_1009_1043/research_note_1034.md',
         'archive_1009_1043/research_note_1035.md',
         'archive_1009_1043/1035/drafts/joint_selection_derivation.md',
         'archive_1046_/1059/mainline_acceptance.json',
         'archive_1046_/_shared/finite_source_readout_checks.json']


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def load(p):
    return json.loads(p.read_text(encoding='utf-8'))


def links(p, freeze=False):
    text=p.read_text(encoding='utf-8')
    if p.parent==ROOT and p.name in BASE['BASE']['TAILS']:
        text=text.split(BASE['BASE']['TAILS'][p.name][0],1)[0]
    text=re.sub(r'\\\[[\s\S]*?\\\]|\$\$[\s\S]*?\$\$','',text)
    n=0
    for target in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)',text):
        target=target.strip().strip('<>')
        if re.match(r'[a-zA-Z][\w+.-]*://',target) or target.startswith('#'):
            continue
        dest=(p.parent/unquote(target.split('#',1)[0])).resolve()
        assert dest.is_file() or (freeze and dest in {RECEIPT.resolve(),ACCEPT.resolve()}), (p,target)
        n+=1
    return n


def run(freeze):
    old=load(SHARED/'finite_source_readout_checks.json')
    assert old['all_checks_passed'] and old['latest_accepted_round_at_batch']==1059
    for rel,h in old['assets_sha256'].items():
        assert digest(SHARED/rel)==h, rel
    for rel,h in old['historical_sha256'].items():
        assert digest(ROOT/rel)==h, rel
    for name,(marker,h) in BASE['BASE']['TAILS'].items():
        raw=(ROOT/name).read_bytes()
        assert hashlib.sha256(raw[raw.index(marker.encode('utf-8')):]).hexdigest()==h
    raw=(ROOT/'ROADMAP.md').read_bytes()
    marker='## M1 检验认知要求的独立选择力'.encode('utf-8')
    assert hashlib.sha256(raw[raw.index(marker):]).hexdigest()==BASE['ROADMAP_SHA']
    review=(HERE/'independent_review.md').read_text(encoding='utf-8')
    assert '通过' in review
    for rel in AUTHOR:
        assert digest(HERE/rel) in review, rel
    for script in ['global_anomaly_pullback.py','independent_check.py']:
        cp=subprocess.run([sys.executable,'-B','-X','utf8',str(HERE/script)],
                          capture_output=True,text=True,encoding='utf-8',check=True)
        data=json.loads(cp.stdout)
        assert data.get('passed',data.get('all_checks_passed',False)), script
    files=AUTHOR+REVIEW+['verify_round1060.py']
    local=sum(links(HERE/f,freeze) for f in files if f.endswith('.md'))
    navlinks=sum(links(p,freeze) for p in NAV)
    stable={'date':'2026-10-09','round':1060,'all_checks_passed':True,
            'author_assets_sha256':{f:digest(HERE/f) for f in AUTHOR},
            'review_assets_sha256':{f:digest(HERE/f) for f in REVIEW},
            'verifier_sha256':digest(Path(__file__)),
            'historical_sha256':{f:digest(ROOT/f) for f in HISTORY},
            'original_roadmap_requirements_sha256':BASE['ROADMAP_SHA'],
            'asset_local_links_checked':local,
            'new_project_conditional_groups':1,'new_cognitive_axioms':0,
            'new_empirical_groups':0,'full_declared_fermion_anomaly_pullback_proved':True,
            'whole_cognitive_axiom_countermodel_proved':False,
            'physical_Spin10_gauge_group_added':False,
            'whole_roadmap_completed':False,'goal_objective_or_status_mutated':False}
    if freeze:
        record=dict(stable)
        record['navigation_snapshot']={'sha256':{str(p.relative_to(ROOT)):digest(p) for p in NAV},
                                       'local_links_checked':navlinks}
        with RECEIPT.open('x',encoding='utf-8') as f:
            json.dump(record,f,ensure_ascii=False,indent=2);f.write('\n')
    else:
        record=load(RECEIPT);record.pop('navigation_snapshot')
        assert record==stable
        acceptance=load(ACCEPT)
        assert acceptance['accepted'] and acceptance['round']==1060
        assert acceptance['receipt_sha256']==digest(RECEIPT)
        assert acceptance['cumulative_scientific_calibrations']==3836
        assert acceptance['stage_groups']==15
        assert not acceptance['whole_roadmap_completed']
    print(json.dumps({'all_checks_passed':True,'round':1060,
                      'author_assets':len(AUTHOR),'review_assets':len(REVIEW),
                      'history':len(HISTORY),'navigation_links':navlinks,
                      'whole_roadmap_completed':False}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--freeze',action='store_true')
    run(p.parse_args().freeze)
