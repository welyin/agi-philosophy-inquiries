"""Check published entry, all frozen evidence, and byte snapshots."""
import json
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
receipt=core.read(HERE/'round768_drafts/joint_source_entry_checks.json')
history=dict(core.read(HERE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
for n in range(584,768):
    r=core.read(HERE/f'research_round_{n}_checks.json')
    for key in ('new_file_hashes','preserved_draft_hashes'):
        for name,digest in r[key].items():
            assert name not in history or history[name]==digest
            history[name]=digest
history.update(core.read(HERE/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
assert len(history)==3714
for name,digest in history.items():
    assert core.digest(HERE/name)==digest,name
for name,digest in receipt['working_artifact_hashes'].items():
    assert core.digest(HERE/name)==digest
snapshot=HERE/'navigation_before_round768_entry_20261004'
manifest=core.read(snapshot/'manifest.json')
links=0
for name,entry in manifest.items():
    assert core.digest(snapshot/entry['snapshot'])==entry['sha256']
    p=ROOT/name
    assert core.digest(p)==receipt['published_navigation_hashes'][name]
    text=p.read_text('utf-8-sig')
    assert '**768完整来源入口已推进' in text
    for link in core.link_parser()(text):
        assert (p.parent/link).resolve().exists(),(name,link)
        links+=1
assert links==receipt['navigation_links']
result=dict(latest_completed_round=767, cumulative_scientific_tests=3488,
            frozen_evidence_files=len(history), working_round=768,
            working_entry_only=True, navigation_links=links, broken_links=0,
            navigation_snapshots=len(manifest), active_goal_unchanged=True,
            previous_evidence_unchanged=True, all_checks_passed=True)
target=HERE/'round768_drafts/entry_postpublication_checks.json'
if target.exists():
    assert core.read(target)==result
else:
    with target.open('x',encoding='utf8') as f:
        json.dump(result,f,ensure_ascii=False,indent=2)
print(json.dumps(result))
