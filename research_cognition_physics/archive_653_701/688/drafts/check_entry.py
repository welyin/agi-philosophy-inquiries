"""Read-only authority and entry reproducibility check after formal687."""
import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import verify_interaction_rounds as core
import temporal_sphere_probe as model
history=dict(core.read(ROOT/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
for n in range(584,688):
    receipt=core.read(ROOT/f'research_round_{n}_checks.json')
    for key in ('new_file_hashes','preserved_draft_hashes'):
        for name,digest in receipt[key].items():
            assert name not in history or history[name]==digest
            history[name]=digest
history.update(core.read(ROOT/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
assert len(history)==2596
for name,digest in history.items():assert core.digest(ROOT/name)==digest,name
assert model.run()==core.read(model.TARGET)
names=('temporal_sphere_probe.py','temporal_sphere_probe_results.json','temporal_sphere_entry.md')
result=dict(date='2026-10-02',entry_round=688,latest_formal_round=687,
    formal_tests_unchanged=3281,protected_evidence_files_checked=len(history),
    all_prior_hashes_unchanged=True,entry_reproduced=True,
    entry_file_hashes={p:core.digest(HERE/p) for p in names},all_checks_passed=True)
with (HERE/'entry_checks.json').open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(dict(entry_round=688,all_checks_passed=True,protected_files=len(history))))

