"""Verify frozen round assets; --recompute additionally runs both science codes."""
from pathlib import Path
import argparse,hashlib,json,re,subprocess,sys
from urllib.parse import unquote
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--recompute',action='store_true')
args=parser.parse_args()
base=Path(__file__).resolve().parent
root=base.parents[2]
r=json.loads((base/'mainline_acceptance.json').read_text(encoding='utf-8'))
assert r['accepted'] and not r['spatial_generation_goal_completed']
for name,digest in r['frozen_assets'].items():
    assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,name
out=json.loads((base/'results.json').read_text(encoding='utf-8'))
for name,digest in out['source_hashes'].items():
    assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,name
if args.recompute:
    subprocess.run([sys.executable,'-B','-X','utf8',str(base/'check.py'),'--check'],check=True,stdout=subprocess.DEVNULL)
    actual=json.loads(subprocess.check_output([sys.executable,'-B','-X','utf8',str(base/'independent_check.py')],text=True,encoding='utf-8'))
    assert actual==json.loads((base/'independent_results.json').read_text(encoding='utf-8'))
links=0
for name in r['frozen_assets']:
    p=root/name
    if p.suffix!='.md':continue
    for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8-sig')):
        target=unquote(target.strip().strip('<>').split('#')[0])
        if not target or re.match(r'^(?:https?://|mailto:|codex:)',target):continue
        q=Path(target)
        if not q.is_absolute():q=p.parent/q
        assert q.exists(),(str(p),target)
        links+=1
print(json.dumps(dict(round=1084,passed=True,frozen_assets=len(r['frozen_assets']),
     frozen_sources=len(out['source_hashes']),local_links=links,
     science_recomputed_in_this_call=args.recompute,
     spatial_generation_goal_completed=False)))
