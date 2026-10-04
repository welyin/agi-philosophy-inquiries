"""Freeze same-model current, classical matching, and record fluctuations."""
import argparse
import ast
import json
from pathlib import Path
import verify_interaction_rounds as core
import verify_round586 as previous
import joint_matter_energy_current as model
import closed_time_path_bridge_entry as bridge

HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_587_checks.json'


def verify():
    base=previous.verify()
    assert (base['round'],base['cumulative_numbered_tests'],base['cumulative_numbered_scientific_files'])==(586,2989,1070)
    historical=dict(core.read(HERE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
    for n in (584,585,586):
        receipt=core.read(HERE/f'research_round_{n}_checks.json')
        for key in ('new_file_hashes','preserved_draft_hashes'):
            for name,digest in receipt[key].items():
                assert name not in historical or historical[name]==digest
                historical[name]=digest
    assert len(historical)==1642
    for name,digest in historical.items():assert core.digest(HERE/name)==digest,name
    # The earlier entry inventories every then-existing 587 draft. New draft
    # files extend that inventory; they must not change any frozen evidence.
    bridge_fresh=bridge.run();bridge_saved=core.read(bridge.TARGET)
    fresh_sources=bridge_fresh.pop('source_hashes');saved_sources=bridge_saved.pop('source_hashes')
    assert bridge_fresh==bridge_saved
    assert all(fresh_sources[name]==digest for name,digest in saved_sources.items())
    bridge_receipt=core.read(HERE/'closed_time_path_bridge_review_checks.json')
    for name,digest in bridge_receipt['artifact_hashes'].items():assert core.digest(HERE/name)==digest,name
    result=core.read(model.TARGET)
    assert result==model.run()
    for name,digest in result['dependency_hashes'].items():assert core.digest(HERE/name)==digest,name
    assert (result['tests_run'],result['failures'],result['errors'])==(5,0,0)
    new={name:core.digest(HERE/name) for name in ('research_note_587.md','joint_matter_energy_current.py','joint_matter_energy_current_results.json')}
    preserved_names=('unified_physics_condition_ledger_587.md','round587_drafts/research_note_587_draft.md',
        'round587_drafts/final_review.txt','round587_drafts/STATUS.md','round587_drafts/local_current_derivation.md',
        'round587_drafts/local_current_entry.py','round587_drafts/local_current_entry_results.json',
        'closed_time_path_bridge_review_586.md','closed_time_path_bridge_entry.py','closed_time_path_bridge_entry_results.json',
        'publish_closed_time_path_bridge.py','closed_time_path_bridge_review_checks.json','round588_drafts/STATUS.md')
    preserved={name:core.digest(HERE/name) for name in preserved_names}
    assert not set(new)&set(historical) and not set(preserved)&set(historical)
    assert len(historical|new|preserved)==1658
    assert (HERE/'research_note_587.md').read_bytes()==(HERE/'round587_drafts/research_note_587_draft.md').read_bytes()
    review=(HERE/'round587_drafts/final_review.txt').read_text('utf8')
    for name in (*new,'unified_physics_condition_ledger_587.md'):assert core.digest(HERE/name) in review,name
    text=core.text_checks(HERE/'research_note_587.md');assert text['display_formulas']==14
    links=0
    for name in ('research_note_587.md','unified_physics_condition_ledger_587.md','round588_drafts/STATUS.md'):
        for link in core.link_parser()((HERE/name).read_text('utf8')):
            p=(HERE/name).parent.joinpath(link).resolve()
            assert p.exists() or p==TARGET.resolve(),(name,link)
            links+=1
    for name in ('joint_matter_energy_current.py',Path(__file__).name):ast.parse((HERE/name).read_text('utf8'))
    if TARGET.exists():
        old=core.read(TARGET)
        assert old['new_file_hashes']==new and old['preserved_draft_hashes']==preserved
    return dict(date='2026-10-01',round=587,scientific_base_through_round=586,
        fresh_tests=dict(run=5,failures=0,errors=0),cumulative_numbered_tests=2994,
        cumulative_numbered_scientific_files=1073,unchanged_prior_evidence_files=1642,
        cumulative_unique_protected_evidence_files=1658,new_file_hashes=new,preserved_draft_hashes=preserved,
        text_checks=text,local_links_checked=links,broken_links=0,saved_results_reproduced=True,
        previous_results_unchanged=True,full_historical_science_rerun=False,
        relevant_584_586_and_nonnumbered_bridge_reproduced=True,
        primary_code_and_note_review_completed=True,independent_final_code_and_draft_review_completed=False,
        inherited_results_not_claimed_as_new_theorems=True,visual_checks_performed=False,
        active_goal_unchanged=True,scope=result['scope'],all_reported_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-checks',action='store_true');args=parser.parse_args()
    result=verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert core.read(TARGET)==result
    print(json.dumps({k:result[k] for k in ('round','cumulative_numbered_tests','cumulative_unique_protected_evidence_files','all_reported_checks_passed')}))
