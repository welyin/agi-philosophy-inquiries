"""Reproduce692 entry without changing formal691 science or frozen evidence."""
import hashlib
import json
from pathlib import Path
import odd_gauss_source_probe as probe

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
TARGET=HERE/'entry_checks.json'
assert not TARGET.exists()
saved=json.loads(probe.TARGET.read_text('utf8'))
assert probe.run()==saved
receipt=json.loads((ROOT/'research_round_691_checks.json').read_text('utf8'))
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in receipt[key].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
for name,digest in saved['dependency_hashes'].items():
    assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
import verify_interaction_rounds as core
links=0
for name in ('odd_gauss_source_entry.md','STATUS.md'):
    text=(HERE/name).read_text('utf8')
    assert text.count('$$')%2==0
    for link in core.link_parser()(text):
        assert (HERE/link).resolve().exists() or (HERE/link).resolve()==TARGET
        links+=1
names=('odd_gauss_source_probe.py','odd_gauss_source_probe_results.json',
       'odd_gauss_source_entry.md','check_entry.py')
out=dict(date='2026-10-02',entry_round=692,latest_formal_round=691,
    saved_results_reproduced=True,prior691_evidence_unchanged=True,
    dependency_hashes_unchanged=True,links_checked=links,broken_links=0,
    complete_averaged_RP_not_evaluated=True,new_formal_round_claimed=False,
    old_space_contracts_and_goal_unchanged=True,all_checks_passed=True,
    artifact_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names})
with TARGET.open('x',encoding='utf8',newline='\n') as f:
    f.write(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(out,ensure_ascii=False))
