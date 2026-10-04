"""Reproduce742 and check all frozen prior scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import verify_interaction_rounds as core
import joint_record_gaussian_boundary as model
HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_742_checks.json'


def verify():
    base=core.read(HERE/'research_round_741_checks.json');assert base['all_reported_checks_passed']
    history=dict(core.read(HERE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
    for n in range(584,742):
        receipt=core.read(HERE/f'research_round_{n}_checks.json')
        for key in ('new_file_hashes','preserved_draft_hashes'):
            for name,digest in receipt[key].items():
                assert name not in history or history[name]==digest
                history[name]=digest
    history.update(core.read(HERE/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
    assert len(history)==3373
    for name,digest in history.items():assert core.digest(HERE/name)==digest,name
    entry=core.read(HERE/'round742_drafts/entry_checks.json')
    for name,digest in entry['artifact_hashes'].items():assert core.digest(HERE/'round742_drafts'/name)==digest,name
    result=model.run();assert result==core.read(model.TARGET)
    assert (result['tests_run'],result['failures'],result['errors'])==(1,0,0)
    for name,digest in result['dependency_hashes'].items():assert core.digest(HERE/name)==digest,name
    main=('research_note_742.md','joint_record_gaussian_boundary.py','joint_record_gaussian_boundary_results.json')
    names=('unified_physics_condition_ledger_742.md','round742_drafts/research_note_742_draft.md',
           'round742_drafts/final_review.txt','round742_drafts/literature_scope_audit.json',
           'round742_drafts/scope_and_dedup_review.md','round743_drafts/STATUS.md',
           'round742_drafts/gaussian_record_entry.py','round742_drafts/gaussian_record_entry_results.json',
           'round742_drafts/research_note_742_working.md','round742_drafts/check_and_publish_entry.py',
           'round742_drafts/entry_checks.json')
    new={name:core.digest(HERE/name) for name in main}
    preserved={name:core.digest(HERE/name) for name in names}
    assert not (set(new)|set(preserved))&set(history)
    assert len(history|new|preserved)==3387
    assert (HERE/main[0]).read_bytes()==(HERE/'round742_drafts/research_note_742_draft.md').read_bytes()
    review=(HERE/'round742_drafts/final_review.txt').read_text('utf8')
    for name in (*main,'unified_physics_condition_ledger_742.md'):assert core.digest(HERE/name) in review
    checks=core.text_checks(HERE/main[0]);assert checks['display_formulas']==12
    links=0
    for name in (main[0],names[0],'round743_drafts/STATUS.md'):
        for link in core.link_parser()((HERE/name).read_text('utf8')):
            p=(HERE/name).parent.joinpath(link).resolve()
            assert p.exists() or p==TARGET.resolve(),(name,link)
            links+=1
    for name in (main[1],Path(__file__).name,'prepare_round742.py','publish_round742.py','postcheck_round742.py'):
        ast.parse((HERE/name).read_text('utf8'))
    return dict(date='2026-10-04',round=742,scientific_base_through_round=741,
        fresh_tests=dict(run=1,failures=0,errors=0),cumulative_numbered_tests=3422,
        cumulative_numbered_scientific_files=1538,unchanged_prior_evidence_files=3373,
        cumulative_unique_protected_evidence_files=3387,new_file_hashes=new,preserved_draft_hashes=preserved,
        text_checks=checks,local_links_checked=links,broken_links=0,saved_results_reproduced=True,
        previous_results_unchanged=True,full_historical_science_rerun=False,
        previous_round_receipt_and_all_evidence_hashes_checked=True,
        original_reference_record_partner_and_positive_Gaussian_gap_checked=True,
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
