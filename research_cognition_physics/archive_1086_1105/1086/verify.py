"""Read-only verification of round 1086 author assets and review persistence."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import subprocess
import sys
from urllib.parse import unquote

BASE=Path(__file__).resolve().parent
STAGE=BASE.parent
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--recompute',action='store_true')
a=p.parse_args()
record=json.loads((BASE/'acceptance.json').read_text(encoding='utf-8'))
assert record['science_increment']==0
for rel,digest in record['author_assets_sha256'].items():
    assert hashlib.sha256((STAGE/rel).read_bytes()).hexdigest()==digest,rel
if a.recompute:
    subprocess.run([sys.executable,'-B','-X','utf8',str(BASE/'check.py')],check=True)
review_saved=(BASE/'independent_review.md').exists() and (BASE/'independent_review_checks.json').exists()
assert review_saved == record['independent_review_files_saved']
if review_saved:
    review=json.loads((BASE/'independent_review_checks.json').read_text(encoding='utf-8'))
    assert review['review_status']=='PASS'
    for rel,digest in review['author_assets_sha256'].items():
        assert hashlib.sha256((BASE/rel).read_bytes()).hexdigest()==digest,rel
links=0
pending=[]
for f in STAGE.rglob('*.md'):
    for target in re.findall(r'\]\(([^)]+)\)',f.read_text(encoding='utf-8-sig')):
        target=unquote(target.strip().strip('<>').split('#')[0])
        if not target or re.match(r'^(?:https?://|mailto:|codex:)',target):
            continue
        q=Path(target)
        if not q.is_absolute():
            q=f.parent/q
        if not q.exists():
            if q.resolve()==(BASE/'independent_review.md').resolve() and not review_saved:
                pending.append(str(q))
                continue
            raise AssertionError((str(f),target))
        links+=1
print(json.dumps({'round':1086,'author_assets_passed':True,'local_links':links,
    'independent_review_received':'PASS','review_files_saved':review_saved,
    'pending_review_links':pending,'recomputed':a.recompute,
    'lorentz_from_cognitive_axioms_proved':False},ensure_ascii=False))
