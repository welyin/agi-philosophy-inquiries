"""Read-only reproducibility and scope checks for a counted-zero bridge."""
from pathlib import Path
import json,hashlib,re,subprocess,sys
from urllib.parse import unquote
B=Path(__file__).resolve().parent
S=B.parents[1];ROOT=S.parents[1]
def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert abs(a-b)<2e-12,(a,b)
    else:assert a==b,(a,b)
for src,dst in (("check.py","results.json"),("independent_check.py","independent_results.json")):
    current=json.loads(subprocess.check_output([sys.executable,"-B","-X","utf8",str(B/src)],encoding="utf-8"))
    compare(current,json.loads((B/dst).read_text(encoding="utf-8")))
prior={}
for n in (1077,1079,1080):
    a=json.loads((S/str(n)/"mainline_acceptance.json").read_text(encoding="utf-8"))
    for name,sha in a["frozen_assets"].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,name
    prior[str(n)]=len(a["frozen_assets"])
links=0
for f in B.glob("*.md"):
    for target in re.findall(r"\]\(([^)]+)\)",f.read_text(encoding="utf-8")):
        if re.match(r"^(https?://|#)",target):continue
        q=(f.parent/unquote(target.split("#")[0])).resolve();assert q.exists(),(f,target)
        links+=1
print(json.dumps({"passed":True,"formal_rounds_added":0,"author_recomputed":True,
  "independent_recomputed":True,"prior_frozen_assets_verified":prior,"local_links":links,
  "spatial_generation_goal_completed":False}))
