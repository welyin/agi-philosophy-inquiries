"""Freeze round 585 with preserved history and fresh scientific reproduction."""
import argparse
import ast
import json
from pathlib import Path
import joint_local_lapse_source as model
import verify_interaction_rounds as core
import verify_round584 as previous

HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_585_checks.json'


def verify():
    base=previous.verify()
    assert (base['round'],base['cumulative_numbered_tests'],
        base['cumulative_numbered_scientific_files'],
        base['cumulative_unique_protected_evidence_files'])==(584,2979,1064,1620)
    result=core.read(model.TARGET)
    assert result==model.run()
    assert (result['tests_run'],result['failures'],result['errors'])==(5,0,0)
    for name,sha in result['dependency_hashes'].items():assert core.digest(HERE/name)==sha,name
    new={name:core.digest(HERE/name) for name in ('research_note_585.md',
        'joint_local_lapse_source.py','joint_local_lapse_source_results.json')}
    preserved={name:core.digest(HERE/name) for name in ('unified_physics_condition_ledger_585.md',
        'round585_drafts/STATUS.md','round585_drafts/entry_scope.md',
        'round585_drafts/joint_local_lapse_source_before_resource_check.py',
        'round585_drafts/joint_local_lapse_source_before_resource_check_results.json',
        'round585_drafts/research_note_585_draft.md','round585_drafts/final_review.txt')}
    assert len(new)==3 and len(preserved)==7
    review=(HERE/'round585_drafts/final_review.txt').read_text('utf8')
    for name in ('research_note_585.md','joint_local_lapse_source.py',
        'joint_local_lapse_source_results.json','unified_physics_condition_ledger_585.md'):
        assert core.digest(HERE/name) in review
    assert (HERE/'research_note_585.md').read_bytes()==(HERE/'round585_drafts/research_note_585_draft.md').read_bytes()
    text=core.text_checks(HERE/'research_note_585.md');assert text['display_formulas']==14
    links=0
    for name in ('research_note_585.md','unified_physics_condition_ledger_585.md'):
        for link in core.link_parser()((HERE/name).read_text('utf8')):
            path=(HERE/link).resolve();assert path.exists() or path==TARGET.resolve(),link
            links+=1
    for name in ('joint_local_lapse_source.py',Path(__file__).name):ast.parse((HERE/name).read_text('utf8'))
    if TARGET.exists():
        old=core.read(TARGET)
        assert old['new_file_hashes']==new and old['preserved_draft_hashes']==preserved
    return dict(date='2026-10-01',round=585,scientific_base_through_round=584,
        fresh_tests=dict(run=5,failures=0,errors=0),cumulative_numbered_tests=2984,
        cumulative_numbered_scientific_files=1067,unchanged_prior_evidence_files=1620,
        cumulative_unique_protected_evidence_files=1630,new_file_hashes=new,preserved_draft_hashes=preserved,
        text_checks=text,local_links_checked=links,broken_links=0,saved_results_reproduced=True,
        previous_results_unchanged=True,full_historical_science_rerun=False,
        direct_round584_dependency_reproduced=True,
        inherited_theorems_not_counted_as_new=True,primary_code_and_note_review_completed=True,
        independent_final_code_and_draft_review_completed=False,visual_checks_performed=False,
        active_goal_unchanged=True,scope=result['scope'],all_reported_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args();result=verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert core.read(TARGET)==result
    print(json.dumps({k:result[k] for k in ('round','cumulative_numbered_tests',
        'cumulative_unique_protected_evidence_files','all_reported_checks_passed')}))
