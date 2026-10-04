"""Verify full fixed spatial metric, same matter, and actual records."""
import argparse
import ast
import json
from pathlib import Path
import verify_interaction_rounds as core
import verify_round588 as previous
import joint_full_spatial_metric as model
HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_589_checks.json'

def verify():
    base=previous.verify()
    assert (base['round'],base['cumulative_numbered_tests'],base['cumulative_numbered_scientific_files'])==(588,2999,1076)
    historical=dict(core.read(HERE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
    for n in (584,585,586,587,588):
        receipt=core.read(HERE/f'research_round_{n}_checks.json')
        for key in ('new_file_hashes','preserved_draft_hashes'):
            for name,digest in receipt[key].items():
                assert name not in historical or historical[name]==digest
                historical[name]=digest
    assert len(historical)==1668
    for name,digest in historical.items():assert core.digest(HERE/name)==digest,name
    result=core.read(model.TARGET)
    assert result==model.run()
    assert (result['tests_run'],result['failures'],result['errors'])==(6,0,0)
    for name,digest in result['dependency_hashes'].items():assert core.digest(HERE/name)==digest,name
    new={name:core.digest(HERE/name) for name in ('research_note_589.md','joint_full_spatial_metric.py','joint_full_spatial_metric_results.json')}
    names=('unified_physics_condition_ledger_589.md','round589_drafts/full_metric_magnetic_entry.md',
           'round589_drafts/full_metric_magnetic_entry.py','round589_drafts/full_metric_magnetic_entry_results.json',
           'round589_drafts/research_note_589_draft.md','round589_drafts/final_review.txt','round590_drafts/STATUS.md')
    preserved={name:core.digest(HERE/name) for name in names}
    assert not set(new)&set(historical) and not set(preserved)&set(historical)
    assert len(historical|new|preserved)==1678
    assert (HERE/'research_note_589.md').read_bytes()==(HERE/'round589_drafts/research_note_589_draft.md').read_bytes()
    review=(HERE/'round589_drafts/final_review.txt').read_text('utf8')
    for name in (*new,'unified_physics_condition_ledger_589.md'):assert core.digest(HERE/name) in review,name
    text=core.text_checks(HERE/'research_note_589.md');assert text['display_formulas']==14
    links=0
    for name in ('research_note_589.md','unified_physics_condition_ledger_589.md','round590_drafts/STATUS.md'):
        for link in core.link_parser()((HERE/name).read_text('utf8')):
            path=(HERE/name).parent.joinpath(link).resolve()
            assert path.exists() or path==TARGET.resolve(),(name,link)
            links+=1
    for name in ('joint_full_spatial_metric.py',Path(__file__).name):ast.parse((HERE/name).read_text('utf8'))
    if TARGET.exists():
        old=core.read(TARGET)
        assert old['new_file_hashes']==new and old['preserved_draft_hashes']==preserved
    return dict(date='2026-10-01',round=589,scientific_base_through_round=588,
        fresh_tests=dict(run=6,failures=0,errors=0),cumulative_numbered_tests=3005,
        cumulative_numbered_scientific_files=1079,unchanged_prior_evidence_files=1668,
        cumulative_unique_protected_evidence_files=1678,new_file_hashes=new,preserved_draft_hashes=preserved,
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
