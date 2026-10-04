"""Read-only reproduce680 entry and verify679 immutable publication."""
import hashlib
import json
from pathlib import Path
import mass_congruence_probe as model
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
saved=json.loads(model.TARGET.read_text('utf8'))
assert model.run()==saved
receipt=json.loads((ROOT/'research_round_679_checks.json').read_text('utf8'))
for key in ('new_file_hashes','preserved_draft_hashes'):
    for p,h in receipt[key].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
post=json.loads((ROOT/'round679_postpublication_checks.json').read_text('utf8'))
assert post['all_checks_passed'] and post['protected_evidence_files']==2501
assert receipt['cumulative_numbered_tests']==3265
files=('mass_congruence_probe.py','mass_congruence_probe_results.json','mass_congruence_entry.md')
result=dict(entry_round=680,formal_latest_round=679,entry_reproduced=True,
    parent_publication_unchanged=True,not_counted_as_formal_round=True,
    all_checks_passed=True,
    hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files})
with (HERE/'entry_checks.json').open('x',encoding='utf8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
print(json.dumps(result))
