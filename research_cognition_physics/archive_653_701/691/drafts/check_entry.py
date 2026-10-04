"""Read-only science reproduction plus new691 entry receipt; no round publication."""
import hashlib
import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import verify_interaction_rounds as core
import projected_cut_current_probe as probe
result=probe.run()
assert result==core.read(probe.TARGET)
for name,digest in result['dependency_hashes'].items():assert core.digest(ROOT/name)==digest,name
base=core.read(ROOT/'round690_postpublication_checks.json')
assert base['all_checks_passed'] and base['round']==690
history=dict(core.read(ROOT/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
for n in range(584,691):
    receipt=core.read(ROOT/f'research_round_{n}_checks.json')
    for key in ('new_file_hashes','preserved_draft_hashes'):
        for name,digest in receipt[key].items():
            assert name not in history or history[name]==digest
            history[name]=digest
history.update(core.read(ROOT/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
assert len(history)==2628
for name,digest in history.items():assert core.digest(ROOT/name)==digest,name
note=HERE/'projected_cut_current_entry.md'
assert core.text_checks(note)['display_formulas']==4
target=HERE/'entry_checks.json'
for link in core.link_parser()(note.read_text('utf8')):
    assert (HERE/link).resolve().exists() or (HERE/link).resolve()==target.resolve(),link
names=('projected_cut_current_probe.py','projected_cut_current_probe_results.json','projected_cut_current_entry.md','check_entry.py')
out=dict(date='2026-10-02',entry_round=691,latest_formal_round=690,not_formal_round=True,
    actual_original_matrices_and_source_reproduced=True,protected_prior_hashes_checked=len(history),
    mature_result_not_claimed_new=True,physical_charge_operator_not_yet_identified=True,
    historical_results_unchanged=True,active_goal_unchanged=True,visual_checks=False,
    artifact_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in names},all_checks_passed=True)
if target.exists():assert core.read(target)==out
else:
    with target.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(dict(entry_round=691,all_checks_passed=True,prior_evidence_files=len(history))))
