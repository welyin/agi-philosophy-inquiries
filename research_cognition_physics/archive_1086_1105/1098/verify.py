"""Read-only acceptance verification for the finite 1098 task bridge."""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    a=json.loads((HERE/'acceptance.json').read_text(encoding='utf8'))
    assert a['status']=='PASS_CONDITIONAL_INFLUENCE_BRIDGE'
    assert a['new_adopted_axioms']==0 and a['scientific_count_increment']==0
    assert not a['entire_conjecture_decided'] and not a['Lorentz_derived']
    assert not a['all_physical_influences_certified'] and not a['new_OC_or_MC_adopted']
    for name,value in a['assets_sha256'].items():assert sha(STAGE/name)==value,name
    for name,value in a['prior_sources_sha256'].items():assert sha(ROOT/name)==value,name
    evidence=json.loads((HERE/'review_evidence.json').read_text(encoding='utf8'))
    assert len(evidence['reviews'])==2
    for review in evidence['reviews']:
        assert review['status']=='PASS_CONDITIONAL_INFLUENCE_BRIDGE'
        assert review['delivery']=='read_only_collaboration_message'
        assert not review['reviewer_files_created']
        assert review['actual_recomputation']['groups']==6
        for name,value in review['source_sha256'].items():assert sha(STAGE/name)==value,name
    links=0
    for name in a['assets_sha256']:
        p=STAGE/name
        if p.suffix!='.md':continue
        content=p.read_text(encoding='utf8')
        assert not any(ord(c)<32 and c not in '\n\r\t' for c in content)
        assert not re.search(r'\\[\[\]]',content)
        assert len(re.findall(r'(?m)^\$\$\s*$',content))%2==0
        for target in re.findall(r'\]\(([^)]+)\)',content):
            if not target.startswith(('http:','https:','#')):
                assert (p.parent/target.split('#',1)[0]).exists(),(p,target)
                links+=1
    p=subprocess.run([sys.executable,'-B','-X','utf8',str(HERE/'check.py')],
        capture_output=True,text=True,encoding='utf8')
    assert p.returncode==0,p.stdout+p.stderr
    check=json.loads(p.stdout);assert check['status']=='PASS' and check['groups']==6
    print(json.dumps(dict(round=1098,status='PASS',frozen_assets=len(a['assets_sha256']),
        prior_sources=len(a['prior_sources_sha256']),local_links=links,
        independent_read_only_reviews=2,recomputation=check,
        entire_conjecture_decided=False,goal_status='active'),ensure_ascii=False))


if __name__=='__main__':main()
