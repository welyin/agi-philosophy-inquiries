"""Read-only verification of frozen round1093 evidence and finite recomputation."""
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
    assert evidence['status']=='PASS_CONDITIONAL_BRIDGE'
    assert evidence['new_adopted_axioms']==0
    assert not evidence['entire_conjecture_decided'] and not evidence['Lorentz_derived']
    for path,value in evidence['assets_sha256'].items():
        assert sha(STAGE/path)==value, path
    for path,value in evidence['prior_sources_sha256'].items():
        assert sha(ROOT/path)==value, path
    reviews=json.loads((HERE/'review_evidence.json').read_text(encoding='utf8'))
    assert len(reviews['reviews'])==2
    for review in reviews['reviews']:
        assert review['status']=='PASS_CONDITIONAL_BRIDGE'
        assert review['delivery']=='read_only_collaboration_message'
        assert review['reviewer_files_created'] is False
        assert review['actual_recomputation']['groups']==9
        for path,value in review['source_sha256'].items():
            assert sha(STAGE/path)==value, path
    links=0
    for name in evidence['assets_sha256']:
        p=STAGE/name
        if p.suffix!='.md': continue
        text=p.read_text(encoding='utf8')
        assert not any(ord(c)<32 and c not in '\n\r\t' for c in text)
        assert not re.search(r'\\[\[\]]',text)
        assert len(re.findall(r'(?m)^\$\$\s*$',text))%2==0
        for target in re.findall(r'\]\(([^)]+)\)',text):
            if not target.startswith(('http:','https:','#')):
                assert (p.parent/target.split('#',1)[0]).exists(), (p,target)
                links+=1
    process=subprocess.run([sys.executable,'-B','-X','utf8',str(HERE/'check.py')],
        capture_output=True,text=True,encoding='utf8')
    assert process.returncode==0, process.stdout+process.stderr
    recalculation=json.loads(process.stdout)
    assert recalculation['status']=='PASS'
    print(json.dumps(dict(round=1093,status='PASS',frozen_assets=len(evidence['assets_sha256']),
        prior_sources=len(evidence['prior_sources_sha256']),independent_read_only_reviews=2,
        local_links=links,recomputation=recalculation,entire_conjecture_decided=False,
        goal_status='active'),ensure_ascii=False))


if __name__=='__main__': main()
