"""Author delivery verification; default read-only, first --write exclusive."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
import subprocess
import sys
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / 'research_round_1042_checks.json'
OWN = ['../research_note_1042.md', 'complementary_entropy.py',
       'complementary_entropy_results.json', 'selection_audit.md',
       'input_dependency_update.md', 'NEXT.md', 'review.md', 'verify_round1042.py']
EXTRA_HISTORY = ['archive_629_652/research_note_636.md',
                 'archive_819_853/research_note_837.md',
                 'archive_819_853/research_note_839.md']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def evidence():
    d = json.loads((HERE/'complementary_entropy_results.json').read_text('utf8'))
    run = subprocess.run([sys.executable,'-B','-X','utf8',str(HERE/'complementary_entropy.py')],
                         check=True,capture_output=True,text=True,encoding='utf8')
    assert json.loads(run.stdout)['passed']
    assert d['new_scientific_calibration_groups']==1
    assert d['new_adopted_cognitive_axioms']==0 and not d['goal_completed']
    assert not d['geometry_certified'] and not d['full_common_parent_model_established']
    history = d['historical_sha256'].copy()
    for p in EXTRA_HISTORY:
        history[p] = sha(BASE/p)
    for p,h in history.items():
        assert sha(BASE/p)==h, p
    # Independent scalar checks of the analytic two-resource comparison.
    s_half=math.log(2)
    s_quarter=2*math.log(2)-3*math.log(3)/4
    assert abs(d['central_entropy_gap']-(s_half-s_quarter))<1e-13
    assert abs(d['physical_edge_trace_distance']-1/4)<1e-13
    cross_d=math.log(4/3)/2
    assert cross_d>0
    links=0
    for name in OWN:
        p=(HERE/name).resolve()
        if p.suffix!='.md': continue
        for raw in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)',p.read_text('utf8')):
            target=unquote(raw.strip().strip('<>').split('#',1)[0])
            if not target or re.match(r'^[a-zA-Z][a-zA-Z0-9+.-]*:',target): continue
            assert (p.parent/target).resolve().exists(), (str(p),target)
            links+=1
    return dict(round=1042,date='2026-10-08',passed=True,author_science_replay_passed=True,
                scientific_groups_proposed=1,additional_author_check_groups=0,
                independent_agent_review_certified_here=False,
                historical_sha256=history,owned_sha256={p:sha(HERE/p) for p in OWN},
                independent_scalar_cross_encoding_relative_entropy=cross_d,
                local_links=links,freeze_scope='author_assets_only')


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--write',action='store_true')
    args=parser.parse_args(); result=evidence()
    if args.write:
        with OUT.open('x',encoding='utf8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False); f.write('\n')
    else:
        assert result==json.loads(OUT.read_text('utf8')), 'Author evidence changed'
    print(json.dumps(dict(passed=True,round=1042,author_assets=len(OWN),
                         historical_inputs=len(result['historical_sha256']),links=result['local_links'])))
