"""Verify positive-energy record dilation, joint clock, and apparatus force."""
import argparse
import ast
import json
from pathlib import Path
import verify_interaction_rounds as core
import verify_round594 as previous
import joint_battery_clock_record as model
HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_595_checks.json'


def verify():
    base=previous.verify()
    assert (base['round'],base['cumulative_numbered_tests'],base['cumulative_numbered_scientific_files'])==(594,3028,1094)
    historical=dict(core.read(HERE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
    for n in range(584,595):
        receipt=core.read(HERE/f'research_round_{n}_checks.json')
        for key in ('new_file_hashes','preserved_draft_hashes'):
            for name,digest in receipt[key].items():
                assert name not in historical or historical[name]==digest
                historical[name]=digest
    assert len(historical)==1735
    for name,digest in historical.items():assert core.digest(HERE/name)==digest,name
    result=core.read(model.TARGET)
    assert result==model.run()
    assert (result['tests_run'],result['failures'],result['errors'])==(4,0,0)
    for name,digest in result['dependency_hashes'].items():assert core.digest(HERE/name)==digest,name
    entry=core.read(HERE/'round595_drafts/positive_battery_entry_checks.json')
    assert entry['checks_passed'] and entry['saved_results_reproduced']
    for name,digest in entry['hashes'].items():assert core.digest(HERE/name)==digest,name
    new={name:core.digest(HERE/name) for name in ('research_note_595.md','joint_battery_clock_record.py','joint_battery_clock_record_results.json')}
    names=('unified_physics_condition_ledger_595.md','round595_drafts/positive_battery_entry.md',
           'round595_drafts/positive_battery_entry.py','round595_drafts/positive_battery_entry_results.json',
           'round595_drafts/positive_battery_entry_checks.json','round595_drafts/research_note_595_draft.md',
           'round595_drafts/final_review.txt','round596_drafts/STATUS.md')
    preserved={name:core.digest(HERE/name) for name in names}
    assert not set(new)&set(historical) and not set(preserved)&set(historical)
    assert len(historical|new|preserved)==1746
    assert (HERE/'research_note_595.md').read_bytes()==(HERE/'round595_drafts/research_note_595_draft.md').read_bytes()
    review=(HERE/'round595_drafts/final_review.txt').read_text('utf8')
    for name in (*new,'unified_physics_condition_ledger_595.md'):assert core.digest(HERE/name) in review,name
    text=core.text_checks(HERE/'research_note_595.md');assert text['display_formulas']==14
    links=0
    for name in ('research_note_595.md','unified_physics_condition_ledger_595.md','round596_drafts/STATUS.md'):
        for link in core.link_parser()((HERE/name).read_text('utf8')):
            path=(HERE/name).parent.joinpath(link).resolve()
            assert path.exists() or path==TARGET.resolve(),(name,link)
            links+=1
    for name in ('joint_battery_clock_record.py',Path(__file__).name):ast.parse((HERE/name).read_text('utf8'))
    if TARGET.exists():
        old=core.read(TARGET)
        assert old['new_file_hashes']==new and old['preserved_draft_hashes']==preserved
    return dict(date='2026-10-01',round=595,scientific_base_through_round=594,
        fresh_tests=dict(run=4,failures=0,errors=0),cumulative_numbered_tests=3032,
        cumulative_numbered_scientific_files=1097,unchanged_prior_evidence_files=1735,
        cumulative_unique_protected_evidence_files=1746,new_file_hashes=new,preserved_draft_hashes=preserved,
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
