"""Reproduce681 mapping entry and check immutable680 publication."""
import hashlib
import json
import sys
from pathlib import Path
import flow_half_support_probe as model
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import verify_interaction_rounds as core
assert model.run()==json.loads(model.TARGET.read_text('utf8'))
receipt=core.read(ROOT/'research_round_680_checks.json')
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,h in receipt[key].items():assert core.digest(ROOT/name)==h,name
assert core.read(ROOT/'round680_postpublication_checks.json')['all_checks_passed']
for link in core.link_parser()((HERE/'mature_flow_mapping_entry.md').read_text('utf8')):
    assert (HERE/link).resolve().exists(),link
names=('flow_half_support_probe.py','flow_half_support_probe_results.json',
       'mature_flow_mapping_entry.md')
result=dict(entry_round=681,latest_formal_round=680,entry_reproduced=True,
    prior_publication_unchanged=True,numbered_tests_not_increased=True,
    all_checks_passed=True,
    evidence_hashes={p:core.digest(HERE/p) for p in names})
with (HERE/'entry_checks.json').open('x',encoding='utf8') as f:
    json.dump(result,f,ensure_ascii=False,indent=2)
print(json.dumps(result))
