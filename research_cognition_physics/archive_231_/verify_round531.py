"""Verify round 531, the joint-condition baseline, and all frozen predecessors."""
import argparse
import ast
import json
from pathlib import Path
import joint_gauge_matter_constraints as model
import unified_physics_condition_audit as inventory
import verify_interaction_rounds as core
import verify_round530 as previous

HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_531_checks.json'


def verify():
    if TARGET.exists():
        base=previous.verify()
    else:
        import verify_round388_integration as anchor
        original=anchor.verify
        def pending_once(pending_report=False):
            old_target=anchor.TARGET
            anchor.TARGET=TARGET
            try:
                return original(True)
            finally:
                anchor.TARGET=old_target
        anchor.verify=pending_once
        try:
            base=previous.verify()
        finally:
            anchor.verify=original
    assert (base['round'],base['cumulative_numbered_tests'],base['cumulative_numbered_scientific_files'],
            base['cumulative_unique_protected_evidence_files']) == (530,2622,902,1114)
    result=core.read(model.TARGET)
    assert result==model.run()
    assert (result['tests_run'],result['failures'],result['errors']) == (9,0,0)
    assert inventory.run()==core.read(inventory.TARGET)
    for name, sha in result['dependency_hashes'].items():
        assert core.digest(HERE/name)==sha,name
    publication=core.read(HERE/'joint_conditions_navigation_checks.json')
    for name,sha in publication['artifact_hashes'].items():
        assert core.digest(HERE/name)==sha,name
    for name,sha in core.read(inventory.TARGET)['preserved_unfinished_candidate_hashes'].items():
        assert core.digest(HERE/name)==sha,name
    note=HERE/'research_note_531.md'
    checked=core.text_checks(note)
    assert checked['display_formulas']==12
    new={name:core.digest(HERE/name) for name in ('research_note_531.md',
         'joint_gauge_matter_constraints.py','joint_gauge_matter_constraints_results.json')}
    extras=set(publication['artifact_hashes'])|set(core.read(inventory.TARGET)['preserved_unfinished_candidate_hashes'])
    extras.update(('joint_conditions_navigation_checks.json','research_goal_joint_physics_20260930.txt',
                   'unified_physics_condition_ledger_531.md',
                   'round531_joint_drafts/independent_review.txt','round531_joint_drafts/final_review.txt'))
    preserved={name:core.digest(HERE/name) for name in sorted(extras)}
    assert not set(new)&set(preserved)
    if TARGET.exists():
        old=core.read(TARGET)
        assert old['new_file_hashes']==new and old['preserved_draft_hashes']==preserved
    links=0
    for doc in (note,HERE/'unified_physics_condition_ledger_531.md'):
        for link in core.link_parser()(doc.read_text('utf8')):
            target=(doc.parent/link).resolve()
            assert target.exists() or target==TARGET.resolve(),(doc,link)
            links+=1
    for name in ('joint_gauge_matter_constraints.py',Path(__file__).name):
        ast.parse((HERE/name).read_text('utf8'))
    return dict(date='2026-09-30',round=531,scientific_base_through_round=530,
        fresh_tests=dict(run=9,failures=0,errors=0),cumulative_numbered_tests=2631,
        cumulative_numbered_scientific_files=905,unchanged_prior_evidence_files=1114,
        cumulative_unique_protected_evidence_files=1114+len(new)+len(preserved),
        new_file_hashes=new,preserved_draft_hashes=preserved,text_checks=checked,
        local_links_checked=links,broken_links=0,saved_results_reproduced=True,
        previous_results_unchanged=True,inherited_theorems_not_counted_as_new=True,
        independent_final_code_and_draft_review_completed=True,
        old_anchor_candidate_not_claimed_completed=True,visual_checks_performed=False,
        active_goal_unchanged=True,scope=result['scope'],all_reported_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args()
    result=verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert core.read(TARGET)==result
    print(json.dumps({k:result[k] for k in ('round','cumulative_numbered_tests',
        'cumulative_unique_protected_evidence_files','all_reported_checks_passed')}))
