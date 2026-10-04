"""Reproduce732 and check every protected preceding artifact."""
import argparse
import ast
import json
from pathlib import Path
import verify_interaction_rounds as core
import joint_relative_source_development as model
HERE=Path(__file__).resolve().parent;TARGET=HERE/'research_round_732_checks.json'


def verify():
    base=core.read(HERE/'research_round_731_checks.json');assert base['all_reported_checks_passed']
    history=dict(core.read(HERE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
    for n in range(584,732):
        receipt=core.read(HERE/f'research_round_{n}_checks.json')
        for key in ('new_file_hashes','preserved_draft_hashes'):
            for name,digest in receipt[key].items():
                assert name not in history or history[name]==digest
                history[name]=digest
    history.update(core.read(HERE/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
    assert len(history)==3242
    for name,digest in history.items():assert core.digest(HERE/name)==digest,name
    import joint_source_constraint_response as previous
    assert previous.run()==core.read(previous.TARGET)
    result=model.run();assert result==core.read(model.TARGET)
    assert (result['tests_run'],result['failures'],result['errors'])==(3,0,0)
    for name,digest in result['dependency_hashes'].items():assert core.digest(HERE/name)==digest,name
    for name,digest in core.read(HERE/'round732_drafts/entry_checks.json')['artifact_hashes'].items():
        assert core.digest(HERE/'round732_drafts'/name)==digest,name
    main=('research_note_732.md','joint_relative_source_development.py','joint_relative_source_development_results.json')
    names=('unified_physics_condition_ledger_732.md','round732_drafts/research_note_732_draft.md',
           'round732_drafts/final_review.txt','round732_drafts/literature_scope_audit.json',
           'round732_drafts/scope_and_dedup_review.md','round733_drafts/STATUS.md',
           'round732_drafts/source_feedback_entry.py','round732_drafts/source_feedback_entry_results.json',
           'round732_drafts/source_feedback_entry.md','round732_drafts/check_and_publish_entry.py',
           'round732_drafts/entry_checks.json')
    new={name:core.digest(HERE/name) for name in main};preserved={name:core.digest(HERE/name) for name in names}
    assert not (set(new)|set(preserved))&set(history)
    assert len(history|new|preserved)==3256
    assert (HERE/main[0]).read_bytes()==(HERE/'round732_drafts/research_note_732_draft.md').read_bytes()
    review=(HERE/'round732_drafts/final_review.txt').read_text('utf8')
    for name in (*main,'unified_physics_condition_ledger_732.md'):assert core.digest(HERE/name) in review
    text=core.text_checks(HERE/main[0]);assert text['display_formulas']==18
    links=0
    for name in (main[0],names[0],'round733_drafts/STATUS.md'):
        for link in core.link_parser()((HERE/name).read_text('utf8')):
            p=(HERE/name).parent.joinpath(link).resolve()
            assert p.exists() or p==TARGET.resolve(),(name,link)
            links+=1
    for name in (main[1],Path(__file__).name,'prepare_round732.py'):ast.parse((HERE/name).read_text('utf8'))
    return dict(date='2026-10-04',round=732,scientific_base_through_round=731,
        fresh_tests=dict(run=3,failures=0,errors=0),cumulative_numbered_tests=3401,
        cumulative_numbered_scientific_files=1508,unchanged_prior_evidence_files=3242,
        cumulative_unique_protected_evidence_files=3256,new_file_hashes=new,preserved_draft_hashes=preserved,
        text_checks=text,local_links_checked=links,broken_links=0,saved_results_reproduced=True,
        previous_results_unchanged=True,full_historical_science_rerun=False,
        previous_round_receipt_and_all_evidence_hashes_checked=True,relevant_recent_dependencies_reproduced=True,
        primary_code_and_note_review_completed=True,independent_final_code_and_draft_review_completed=False,
        inherited_results_not_claimed_as_new_theorems=True,visual_checks_performed=False,
        active_goal_unchanged=True,scope=result['scope'],all_reported_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-checks',action='store_true');args=p.parse_args();r=verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==core.read(TARGET)
    print(json.dumps({k:r[k] for k in ('round','cumulative_numbered_tests','all_reported_checks_passed')}))
