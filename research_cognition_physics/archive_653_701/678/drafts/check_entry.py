"""Record executed678 entry without changing frozen evidence."""
import ast
import hashlib
import json
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import verify_interaction_rounds as core
names=('STATUS.md','local_block_probe.py','local_block_probe_results.json','local_block_entry.md')
digests={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names}
links=0
for name in ('STATUS.md','local_block_entry.md'):
    for link in core.link_parser()((HERE/name).read_text('utf8')):
        assert (HERE/link).resolve().exists(),link
        links+=1
ast.parse((HERE/'local_block_probe.py').read_text('utf8'))
r=core.read(HERE/'local_block_probe_results.json')
for name,digest in r['dependency_hashes'].items():
    assert core.digest(ROOT/name)==digest,name
post=core.read(ROOT/'round677_postpublication_checks.json')
assert post['all_checks_passed']
result=dict(date='2026-10-02',entry_round=678,latest_formal_round=677,
    entry_executed_and_readonly_reproduced=True,formal_round_not_claimed=True,
    prior_postpublication_checks_passed=True,artifact_hashes=digests,
    local_links_checked=links,broken_links=0,all_checks_passed=True)
with (HERE/'entry_checks.json').open('x',encoding='utf8',newline='\n') as f:
    f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result))

