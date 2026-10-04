"""Read-only evidence/navigation verification after664 publication."""
import json
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
history=dict(core.read(HERE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
for n in range(584,665):
    receipt=core.read(HERE/f'research_round_{n}_checks.json')
    for key in ('new_file_hashes','preserved_draft_hashes'):
        for name,digest in receipt[key].items():
            assert name not in history or history[name]==digest
            history[name]=digest
history.update(core.read(HERE/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
assert len(history)==2306
for name,digest in history.items():assert core.digest(HERE/name)==digest,name
snap=HERE/'navigation_before_round664_20261002'
manifest=core.read(snap/'manifest.json');links=0
for name,entry in manifest.items():
    assert core.digest(snap/entry['snapshot'])==entry['sha256']
    p=ROOT/name;content=p.read_text('utf-8-sig')
    assert '**第664轮完成：**' in content
    for link in core.link_parser()(content):
        assert (p.parent/link).resolve().exists(),(name,link)
        links+=1
result=dict(date='2026-10-02',round=664,protected_evidence_files=len(history),
    navigation_snapshots=len(manifest),navigation_links=links,all_hashes_unchanged=True,
    broken_links=0,all_checks_passed=True)
with (HERE/'round664_postpublication_checks.json').open('x',encoding='utf8') as f:
    json.dump(result,f,ensure_ascii=False,indent=2)
print(json.dumps(result))
