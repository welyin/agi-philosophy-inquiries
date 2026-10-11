"""Publish reviewed, unchanged round cores; never treat this as scientific testing."""
from pathlib import Path
import argparse
import hashlib
import json
import re

STAGE=Path(__file__).resolve().parents[1]
sha=lambda b:hashlib.sha256(b).hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('round',type=int);args=ap.parse_args()
    n=args.round; folder=STAGE/str(n); acceptance=folder/'acceptance.json'
    assert not acceptance.exists(),'Already accepted; do not overwrite historical acceptance'
    review=(folder/'independent_review.md').read_bytes()
    result=json.loads((folder/'results.json').read_text(encoding='utf-8'))
    assert result['status'].startswith('PASS')
    groups=result.get('numerical_check_groups',result.get('check_groups'))
    assert isinstance(groups,int) and groups>0,'Missing check group count'
    core={p:sha((folder/p).read_bytes()) for p in ['proof.md','check.py','results.json','NEXT.md',f'../research_note_{n}.md']}
    for p,h in core.items():assert h in review.decode('utf-8-sig'),('Unreviewed core',p)
    assert 'PASS' in review.decode('utf-8-sig')
    lp=STAGE/'_shared/dependency_ledger.json'; old=lp.read_bytes(); ledger=json.loads(old)
    assert n in ledger['in_progress_rounds'] and n not in ledger['completed_rounds']
    ledger['in_progress_rounds'].remove(n);ledger['completed_rounds'].append(n);ledger['completed_rounds'].sort()
    ledger['latest_completed_round']=max(ledger['completed_rounds'])
    ledger.setdefault('round_summary',{})[str(n)]=dict(status=result['status'],new_axioms=0,
        check_groups=groups,scientific_calibration=0,
        independent_review='independent_review.md')
    idx=STAGE/'文件索引.md';before=idx.read_bytes();txt=before.decode('utf-8-sig').replace('\r\n','\n')
    pattern=rf'^\|{n}\|([^|]+)\|[^|]+\|$'
    changed,count=re.subn(pattern,rf'|{n}|[\1](research_note_{n}.md)|已验收；[证明]({n}/proof.md) · [结果]({n}/results.json) · [独审]({n}/independent_review.md)|',txt,flags=re.M)
    assert count==1
    record=dict(round=n,status='PASS',date='2026-10-11',core_sha256=core,
        review_sha256=sha(review),mainline_reviewer='/root',
        new_adopted_axioms=0,scientific_calibration_increment=0,scientific_count_total=3860,
        numerical_check_groups=groups,goal_status='active',goal_complete=False,
        scope=result['status'],external_human_review=False)
    assert lp.read_bytes()==old and idx.read_bytes()==before,'Concurrent navigation change'
    acceptance.write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    lp.write_text(json.dumps(ledger,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    idx.write_text(changed,encoding='utf-8')
    print(json.dumps(dict(status='PASS',accepted_round=n,check_groups=groups,goal_complete=False)))

if __name__=='__main__':main()
