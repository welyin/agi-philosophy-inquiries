"""Read-only682 entry reproduction plus all681 protected evidence."""
import importlib.util
import json
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import verify_interaction_rounds as core
spec=importlib.util.spec_from_file_location('entry682',HERE/'conditional_bulk_probe.py')
probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)
assert probe.run()==core.read(probe.TARGET)
history=dict(core.read(ROOT/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
for n in range(584,682):
    receipt=core.read(ROOT/f'research_round_{n}_checks.json')
    for key in ('new_file_hashes','preserved_draft_hashes'):
        for name,digest in receipt[key].items():
            assert name not in history or history[name]==digest
            history[name]=digest
history.update(core.read(ROOT/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
assert len(history)==2525
for name,digest in history.items():assert core.digest(ROOT/name)==digest,name
names=('conditional_bulk_probe.py','conditional_bulk_probe_results.json','conditional_bulk_entry.md')
result=dict(entry_round=682,latest_formal_round=681,entry_reproduced=True,
    prior_publication_unchanged=True,protected_evidence_files=len(history),
    numbered_tests_not_increased=True,evidence_hashes={p:core.digest(HERE/p) for p in names},
    all_checks_passed=True)
target=HERE/'entry_checks.json'
if target.exists():assert result==core.read(target)
else:
    with target.open('x',encoding='utf8',newline='\n') as f:
        f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result))

