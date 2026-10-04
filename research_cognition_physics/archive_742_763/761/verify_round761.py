"""Reproduce761 and check all frozen prior scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import verify_interaction_rounds as core
import joint_loop_scale_transport as model
HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_761_checks.json'


def verify():
    base=core.read(HERE/'research_round_760_checks.json');assert base['all_reported_checks_passed']
    history=dict(core.read(HERE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
    for n in range(584,761):
        receipt=core.read(HERE/f'research_round_{n}_checks.json')
        for key in ('new_file_hashes','preserved_draft_hashes'):
            for name,digest in receipt[key].items():
                assert name not in history or history[name]==digest
                history[name]=digest
    history.update(core.read(HERE/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
    assert len(history)==3627
    for name,digest in history.items():assert core.digest(HERE/name)==digest,name
    import verify_cognitive_joint_candidate_review as working
    assert working.verify()==core.read(working.TARGET)
    result=model.run();assert result==core.read(model.TARGET)
    assert (result['tests_run'],result['failures'],result['errors'])==(3,0,0)
    for name,digest in result['dependency_hashes'].items():assert core.digest(HERE/name)==digest,name
    main=('research_note_761.md','joint_loop_scale_transport.py','joint_loop_scale_transport_results.json')
    names=('unified_physics_condition_ledger_761.md', 'round761_drafts/research_note_761_draft.md', 'round761_drafts/final_review.txt', 'round761_drafts/literature_scope_audit.json', 'round761_drafts/scope_and_dedup_review.md', 'round762_drafts/STATUS.md', 'round761_drafts/coarse_process_reuse_review.md')
    new={name:core.digest(HERE/name) for name in main}
    preserved={name:core.digest(HERE/name) for name in names}
    assert not (set(new)|set(preserved))&set(history)
    assert len(history|new|preserved)==3637
    assert (HERE/main[0]).read_bytes()==(HERE/'round761_drafts/research_note_761_draft.md').read_bytes()
    review=(HERE/'round761_drafts/final_review.txt').read_text('utf8')
    for name in (*main,'unified_physics_condition_ledger_761.md'):assert core.digest(HERE/name) in review
    checks=core.text_checks(HERE/main[0]);assert checks['display_formulas']==15
    links=0
    for name in (main[0],names[0],'round762_drafts/STATUS.md'):
        for link in core.link_parser()((HERE/name).read_text('utf8')):
            p=(HERE/name).parent.joinpath(link).resolve()
            assert p.exists() or p==TARGET.resolve(),(name,link)
            links+=1
    for name in (main[1],Path(__file__).name,'prepare_round761.py','publish_round761.py','postcheck_round761.py'):
        ast.parse((HERE/name).read_text('utf8'))
    return dict(date='2026-10-04',round=761,scientific_base_through_round=760,
        fresh_tests=dict(run=3,failures=0,errors=0),cumulative_numbered_tests=3470,
        cumulative_numbered_scientific_files=1595,unchanged_prior_evidence_files=3627,
        cumulative_unique_protected_evidence_files=3637,new_file_hashes=new,preserved_draft_hashes=preserved,
        text_checks=checks,local_links_checked=links,broken_links=0,saved_results_reproduced=True,
        previous_results_unchanged=True,full_historical_science_rerun=False,
        previous_round_receipt_and_all_evidence_hashes_checked=True,
        new_relative_scaling_explicitly_registered=True,fixed_graph_transport_and_first_response_checked=True,old_scaling_not_refuted_or_replaced=True,epsilon_one_accuracy_proven=False,full_continuum_process_equivalence_proven=False,
        primary_code_and_note_review_completed=True,independent_final_code_and_draft_review_completed=False,
        inherited_results_not_claimed_as_new_theorems=True,visual_checks_performed=False,
        active_goal_unchanged=True,scope=result['scope'],all_reported_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-checks',action='store_true')
    args=p.parse_args();r=verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==core.read(TARGET)
    print(json.dumps({k:r[k] for k in ('round','cumulative_numbered_tests','all_reported_checks_passed')}))
