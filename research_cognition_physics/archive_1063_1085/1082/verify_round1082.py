"""Read-only round receipt, exact recomputation and local-link checks."""
from pathlib import Path
import hashlib,json,re,subprocess,sys
from urllib.parse import unquote
base=Path(__file__).resolve().parent
root=base.parents[2]
receipt=json.loads((base/"mainline_acceptance.json").read_text(encoding="utf-8"))
assert receipt["accepted"] and not receipt["spatial_generation_goal_completed"]
for name,sha in receipt["frozen_assets"].items():
    assert hashlib.sha256((root/name).read_bytes()).hexdigest()==sha,name
for script,result in [("check.py","results.json"),("independent_check.py","independent_results.json")]:
    output=subprocess.check_output([sys.executable,"-B","-X","utf8",str(base/script)],text=True,encoding="utf-8")
    assert json.loads(output)==json.loads((base/result).read_text(encoding="utf-8")),script
links=0
for name in receipt["frozen_assets"]:
    p=root/name
    if p.suffix!=".md":continue
    for target in re.findall(r"\]\(([^)]+)\)",p.read_text(encoding="utf-8-sig")):
        target=target.strip().strip("<>")
        if re.match(r"^(?:https?://|mailto:|#|codex:)",target):continue
        target=unquote(target.split("#")[0])
        if not target:continue
        q=Path(target)
        if not q.is_absolute():q=p.parent/q
        assert q.exists(),(str(p),target)
        links+=1
print(json.dumps(dict(round=1082,passed=True,frozen_assets=len(receipt["frozen_assets"]),
 local_links=links,author_recomputed_exact=True,independent_recomputed_exact=True,
 spatial_generation_goal_completed=False,whole_roadmap_completed=False)))
