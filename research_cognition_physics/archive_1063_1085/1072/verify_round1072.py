"""Read-only receipt, recomputation and local-link verification for round 1072."""
from pathlib import Path
import hashlib,json,re,subprocess,sys
from urllib.parse import unquote
base=Path(__file__).resolve().parent
root=base.parents[2]
receipt=json.loads((base/"mainline_acceptance.json").read_text(encoding="utf-8"))
assert receipt["accepted"] and not receipt["spatial_generation_goal_completed"]
for name,sha in receipt["frozen_assets"].items():
    assert hashlib.sha256((root/name).read_bytes()).hexdigest()==sha,name

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert abs(a-b)<2e-12,(a,b)
    else:assert a==b,(a,b)
for script,result in [("check.py","results.json"),("independent_check.py","independent_results.json")]:
    output=subprocess.check_output([sys.executable,"-B","-X","utf8",str(base/script)],text=True,encoding="utf-8")
    compare(json.loads(output),json.loads((base/result).read_text(encoding="utf-8")))
links=0
for name in receipt["frozen_assets"]:
    p=root/name
    if p.suffix!=".md" or not (p.parent==base or p.name=="research_note_1072.md"):continue
    for target in re.findall(r"\]\(([^)]+)\)",p.read_text(encoding="utf-8-sig")):
        target=target.strip().strip("<>")
        if re.match(r"^(?:https?://|mailto:|#|codex:)",target):continue
        target=unquote(target.split("#")[0])
        if not target:continue
        q=Path(target)
        if not q.is_absolute():q=p.parent/q
        assert q.exists(),(str(p),target)
        links+=1
print(json.dumps({"round":1072,"passed":True,"frozen_assets":len(receipt["frozen_assets"]),
 "local_links":links,"author_recomputed":True,"independent_recomputed":True,
 "spatial_generation_goal_completed":False,"whole_roadmap_completed":False}))
