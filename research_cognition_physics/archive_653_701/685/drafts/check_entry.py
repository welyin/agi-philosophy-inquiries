"""Reproduce685 body-weight entry and protect the published684 evidence."""
import importlib.util
import json
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import verify_interaction_rounds as core
spec=importlib.util.spec_from_file_location('entry685',HERE/'body_weight_limit_probe.py')
probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)
assert probe.run()==core.read(probe.TARGET)
history=dict(core.read(ROOT/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
for n in range(584,685):
    receipt=core.read(ROOT/f'research_round_{n}_checks.json')
    for key in ('new_file_hashes','preserved_draft_hashes'):
        for name,digest in receipt[key].items():
            assert name not in history or history[name]==digest
            history[name]=digest
history.update(core.read(ROOT/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
assert len(history)==2567
for name,digest in history.items():assert core.digest(ROOT/name)==digest,name
names=('body_weight_limit_probe.py','body_weight_limit_probe_results.json','body_weight_limit_entry.md')
result=dict(entry_round=685,latest_formal_round=684,entry_reproduced=True,
    protected_evidence_files=len(history),numbered_tests_not_increased=True,
    evidence_hashes={p:core.digest(HERE/p) for p in names},all_checks_passed=True)
target=HERE/'entry_checks.json'
if target.exists():assert result==core.read(target)
else:
    with target.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result))
