"""Verify joint pointer premeasurement and terminal-record energy boundary."""
import argparse
import ast
import json
from pathlib import Path
import verify_interaction_rounds as core
import verify_round593 as previous
import joint_pointer_energy_boundary as model
HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_594_checks.json'


def verify():
    base=previous.verify()
    assert (base['round'],base['cumulative_numbered_tests'],base['cumulative_numbered_scientific_files'])==(593,3024,1091)
    historical=dict(core.read(HERE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
    for n in range(584,594):
        receipt=core.read(HERE/f'research_round_{n}_checks.json')
        for key in ('new_file_hashes','preserved_draft_hashes'):
            for name,digest in receipt[key].items():
                assert name not in historical or historical[name]==digest
                historical[name]=digest
    assert len(historical)==1724
    for name,digest in historical.items():assert core.digest(HERE/name)==digest,name
    result=core.read(model.TARGET)
    assert result==model.run()
    assert (result['tests_run'],result['failures'],result['errors'])==(4,0,0)
    for name,digest in result['dependency_hashes'].items():assert core.digest(HERE/name)==digest,name
    entry=core.read(HERE/'round594_drafts/premeasurement_energy_entry_checks.json')
    assert entry['checks_passed'] and entry['saved_results_reproduced']
    for name,digest in entry['hashes'].items():assert core.digest(HERE/name)==digest,name
    new={name:core.digest(HERE/name) for name in ('research_note_594.md','joint_pointer_energy_boundary.py','joint_pointer_energy_boundary_results.json')}
    names=('unified_physics_condition_ledger_594.md','round594_drafts/premeasurement_energy_entry.md',
           'round594_drafts/premeasurement_energy_entry.py','round594_drafts/premeasurement_energy_entry_results.json',
           'round594_drafts/premeasurement_energy_entry_checks.json','round594_drafts/research_note_594_draft.md',
           'round594_drafts/final_review.txt','round595_drafts/STATUS.md')
    preserved={name:core.digest(HERE/name) for name in names}
    assert not set(new)&set(historical) and not set(preserved)&set(historical)
    assert len(historical|new|preserved)==1735
    assert (HERE/'research_note_594.md').read_bytes()==(HERE/'round594_drafts/research_note_594_draft.md').read_bytes()
    review=(HERE/'round594_drafts/final_review.txt').read_text('utf8')
    for name in (*new,'unified_physics_condition_ledger_594.md'):assert core.digest(HERE/name) in review,name
    text=core.text_checks(HERE/'research_note_594.md');assert text['display_formulas']==10
    links=0
    for name in ('research_note_594.md','unified_physics_condition_ledger_594.md','round595_drafts/STATUS.md'):
        for link in core.link_parser()((HERE/name).read_text('utf8')):
            path=(HERE/name).parent.joinpath(link).resolve()
            assert path.exists() or path==TARGET.resolve(),(name,link)
            links+=1
    for name in ('joint_pointer_energy_boundary.py',Path(__file__).name):ast.parse((HERE/name).read_text('utf8'))
    if TARGET.exists():
        old=core.read(TARGET)
        assert old['new_file_hashes']==new and old['preserved_draft_hashes']==preserved
    return dict(date='2026-10-01',round=594,scientific_base_through_round=593,
        fresh_tests=dict(run=4,failures=0,errors=0),cumulative_numbered_tests=3028,
        cumulative_numbered_scientific_files=1094,unchanged_prior_evidence_files=1724,
        cumulative_unique_protected_evidence_files=1735,new_file_hashes=new,preserved_draft_hashes=preserved,
        text_checks=text,local_links_checked=links,broken_links=0,saved_results_reproduced=True,
        previous_results_unchanged=True,full_historical_science_rerun=False,
        relevant_recent_dependencies_reproduced=True,primary_code_and_note_review_completed=True,
        independent_final_code_and_draft_review_completed=False,inherited_results_not_claimed_as_new_theorems=True,
        visual_checks_performed=False,active_goal_unchanged=True,scope=result['scope'],all_reported_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-checks',action='store_true');args=parser.parse_args()
    result=verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert core.read(TARGET)==result
    print(json.dumps({k:result[k] for k in ('round','cumulative_numbered_tests','cumulative_unique_protected_evidence_files','all_reported_checks_passed')}))
