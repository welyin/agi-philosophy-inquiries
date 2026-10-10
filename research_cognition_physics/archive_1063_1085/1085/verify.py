"""Verify frozen synthesis assets; optional --recompute repeats finite model checks."""
from pathlib import Path
import argparse,hashlib,json,re,subprocess,sys
from urllib.parse import unquote
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--recompute',action='store_true')
a=p.parse_args()
base=Path(__file__).resolve().parent
root=base.parents[2]
m=json.loads((base/'mainline_acceptance.json').read_text(encoding='utf-8'))
assert m['accepted'] and m['new_scientific_calibrations']==0
for rel,digest in m['frozen_assets'].items():
    assert hashlib.sha256((root/rel).read_bytes()).hexdigest()==digest,rel
r=json.loads((base/'results.json').read_text(encoding='utf-8'))
for rel,digest in r['source_hashes'].items():
    assert hashlib.sha256((root/rel).read_bytes()).hexdigest()==digest,rel
if a.recompute:
    for f in ['check.py','independent_model/check.py']:
        subprocess.run([sys.executable,'-B','-X','utf8',str(base/f),'--check'],check=True,stdout=subprocess.DEVNULL)
links=0
for rel in m['frozen_assets']:
    f=root/rel
    if f.suffix!='.md':continue
    for target in re.findall(r'\]\(([^)]+)\)',f.read_text(encoding='utf-8-sig')):
        target=unquote(target.strip().strip('<>').split('#')[0])
        if not target or re.match(r'^(?:https?://|mailto:|codex:)',target):continue
        q=Path(target)
        if not q.is_absolute():q=f.parent/q
        assert q.exists(),(rel,target)
        links+=1
print(json.dumps({'round':1085,'passed':True,'frozen_assets':len(m['frozen_assets']),
    'frozen_sources':len(r['source_hashes']),'local_links':links,'finite_checks_recomputed':a.recompute,
    'conditional_classification_complete':True,'actual_cognitive_source_proved':False,'goal_status':'active'}))
