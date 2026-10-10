"""Read-only verification of round 1094 and its explicitly scoped evidence."""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    evidence=json.loads((HERE/'acceptance.json').read_text(encoding='utf8'))
    assert evidence['status']=='PASS_CONDITIONAL_CLASSIFICATION'
    assert evidence['new_adopted_axioms']==0
    assert not evidence['entire_conjecture_decided'] and not evidence['Lorentz_derived']
    assert not evidence['all_six_protocols_implemented']
    for path,value in evidence['assets_sha256'].items():
        assert sha(STAGE/path)==value,path
    for path,value in evidence['prior_sources_sha256'].items():
        assert sha(ROOT/path)==value,path
    reviews=json.loads((HERE/'review_evidence.json').read_text(encoding='utf8'))
    assert len(reviews['reviews'])==2
    for review in reviews['reviews']:
        assert review['status']=='PASS_CONDITIONAL_CLASSIFICATION'
        assert review['delivery']=='read_only_collaboration_message'
        assert review['reviewer_files_created'] is False
        assert review['actual_recomputation']['groups']==7
        for path,value in review['source_sha256'].items():
            assert sha(STAGE/path)==value,path
    links=0
    for name in evidence['assets_sha256']:
        p=STAGE/name
        if p.suffix!='.md': continue
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
    result=json.loads(p.stdout)
    assert result['status']=='PASS' and result['groups']==7
    print(json.dumps(dict(round=1094,status='PASS',frozen_assets=len(evidence['assets_sha256']),
        prior_sources=len(evidence['prior_sources_sha256']),independent_read_only_reviews=2,
        local_links=links,recomputation=result,entire_conjecture_decided=False,
        goal_status='active'),ensure_ascii=False))


if __name__=='__main__': main()
