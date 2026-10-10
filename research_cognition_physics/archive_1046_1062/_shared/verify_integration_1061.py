"""Verify accepted 1061, frozen history and current navigation; no default writes."""
from pathlib import Path
from hashlib import sha256
from urllib.parse import unquote
import argparse
import json
import re
import runpy
import subprocess
import sys

HERE = Path(__file__).resolve().parent
PHASE = HERE.parent
ROOT = PHASE.parent
ROUND = PHASE/'1061'
RECEIPT = HERE/'integration_1061_checks.json'
PREVIOUS = runpy.run_path(str(HERE/'verify_feedback_dense_matter_after1060.py'))
NAV = PREVIOUS['NAV']
TAILS = PREVIOUS['BASE']['BASE']['BASE']['TAILS']


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text(encoding='utf-8'))


def links(path, freeze):
    text = path.read_text(encoding='utf-8')
    if path.parent == ROOT and path.name in TAILS:
        text = text.split(TAILS[path.name][0],1)[0]
    text = re.sub(r'\\\[[\s\S]*?\\\]|\$\$[\s\S]*?\$\$', '', text)
    count=0
    for ref in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)',text):
        ref=ref.strip().strip('<>')
        if re.match(r'[a-zA-Z][\w+.-]*://',ref) or ref.startswith('#'):
            continue
        p=(path.parent/unquote(ref.split('#',1)[0])).resolve()
        assert p.is_file() or (freeze and p == RECEIPT.resolve()), (path,ref)
        count+=1
    return count


def run(freeze=False):
    # The previous verifier reads all current nav links, so do not point nav
    # at a pending new receipt until this receipt has actually been written.
    PREVIOUS['run'](False)
    acceptance=load(ROUND/'mainline_acceptance.json')
    assert acceptance['accepted'] and acceptance['round']==1061
    assert acceptance['cumulative_scientific_calibrations']==3837
    assert acceptance['stage_groups']==16
    assert not acceptance['whole_roadmap_completed']
    author=load(ROUND/'research_round_1061_checks.json')
    assert digest(ROUND/'research_round_1061_checks.json')==acceptance['receipt_sha256']
    for field in ('author_assets_sha256','review_assets_sha256'):
        for name,expected in acceptance[field].items():
            assert digest(ROUND/name)==expected, name
    assert author['owned_sha256']==acceptance['author_assets_sha256']
    for review in ('statistical_review.md','physical_review.md'):
        body=(ROUND/review).read_text(encoding='utf-8')
        assert '通过' in body and acceptance['receipt_sha256'] in body
        assert all(value in body for value in author['owned_sha256'].values())
    for code in ('verify_round1061.py','independent_check.py'):
        cp=subprocess.run([sys.executable,'-B','-X','utf8',str(ROUND/code)],
                          capture_output=True,text=True,encoding='utf-8',check=True)
        assert cp.stdout.strip(), code
    paths=[ROUND/name for name in acceptance['author_assets_sha256']]
    paths += [ROUND/name for name in acceptance['review_assets_sha256']]
    paths += [ROUND/'research_round_1061_checks.json',ROUND/'mainline_acceptance.json',
              HERE/'integration_1061.md',Path(__file__)]
    own_links=sum(links(p,freeze) for p in paths if p.suffix=='.md')
    nav_count=sum(links(p,freeze) for p in NAV)
    stable={
        'round':1061,'date':'2026-10-09','all_checks_passed':True,
        'assets_sha256':{str(p.relative_to(ROOT)):digest(p) for p in paths},
        'baseline_receipts_sha256':{
            'round1060':digest(PHASE/'1060/mainline_acceptance.json'),
            'feedback_dense_matter':digest(HERE/'feedback_dense_matter_after1060_checks.json')},
        'original_roadmap_requirements_sha256':
            load(HERE/'feedback_dense_matter_after1060_checks.json')['original_roadmap_requirements_sha256'],
        'asset_links_checked':own_links,
        'independent_statistical_and_physical_reviews_passed':True,
        'latest_accepted_round':1061,'cumulative_scientific_calibrations':3837,
        'stage_groups':16,'new_project_empirical_groups':1,'new_cognitive_axioms':0,
        'whole_roadmap_completed':False,'goal_objective_or_status_mutated':False,
    }
    if freeze:
        frozen=dict(stable)
        frozen['navigation_snapshot']={'local_links_checked':nav_count,
            'sha256':{str(p.relative_to(ROOT)):digest(p) for p in NAV}}
        with RECEIPT.open('x',encoding='utf-8') as out:
            json.dump(frozen,out,ensure_ascii=False,indent=2);out.write('\n')
    else:
        record=load(RECEIPT);record.pop('navigation_snapshot')
        assert record==stable
    print(json.dumps({'round':1061,'all_checks_passed':True,'assets':len(paths),
                      'navigation_links':nav_count,'whole_roadmap_completed':False}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--freeze',action='store_true')
    run(p.parse_args().freeze)
