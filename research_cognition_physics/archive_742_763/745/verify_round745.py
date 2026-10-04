"""Reproduce745 and check all frozen prior scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import verify_interaction_rounds as core
import joint_native_history_certificate as model
HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_745_checks.json'


def verify():
    base=core.read(HERE/'research_round_744_checks.json');assert base['all_reported_checks_passed']
    history=dict(core.read(HERE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
    for n in range(584,745):
        receipt=core.read(HERE/f'research_round_{n}_checks.json')
        for key in ('new_file_hashes','preserved_draft_hashes'):
            for name,digest in receipt[key].items():
                assert name not in history or history[name]==digest
                history[name]=digest
    history.update(core.read(HERE/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
    assert len(history)==3421
    for name,digest in history.items():assert core.digest(HERE/name)==digest,name
    result=model.run();assert result==core.read(model.TARGET)
    assert (result['tests_run'],result['failures'],result['errors'])==(2,0,0)
    for name,digest in result['dependency_hashes'].items():assert core.digest(HERE/name)==digest,name
    main=('research_note_745.md','joint_native_history_certificate.py','joint_native_history_certificate_results.json')
    names=('unified_physics_condition_ledger_745.md','round745_drafts/research_note_745_draft.md',
           'round745_drafts/final_review.txt','round745_drafts/literature_scope_audit.json',
           'round745_drafts/scope_and_dedup_review.md','round746_drafts/STATUS.md',
           'round745_drafts/history_mechanism_probe.py','round745_drafts/history_mechanism_probe_results.json',
           'round745_drafts/exact_history_density.py','round745_drafts/exact_history_density_results.json',
           'round745_drafts/certify_history_integral.py','round745_drafts/certify_history_integral_results.json')
    new={name:core.digest(HERE/name) for name in main}
    preserved={name:core.digest(HERE/name) for name in names}
    assert not (set(new)|set(preserved))&set(history)
    assert len(history|new|preserved)==3436
    assert (HERE/main[0]).read_bytes()==(HERE/'round745_drafts/research_note_745_draft.md').read_bytes()
    review=(HERE/'round745_drafts/final_review.txt').read_text('utf8')
    for name in (*main,'unified_physics_condition_ledger_745.md'):assert core.digest(HERE/name) in review
    checks=core.text_checks(HERE/main[0]);assert checks['display_formulas']==16
    links=0
    for name in (main[0],names[0],'round746_drafts/STATUS.md'):
        for link in core.link_parser()((HERE/name).read_text('utf8')):
            p=(HERE/name).parent.joinpath(link).resolve()
            assert p.exists() or p==TARGET.resolve(),(name,link)
            links+=1
    for name in (main[1],Path(__file__).name,'prepare_round745.py','publish_round745.py','postcheck_round745.py'):
        ast.parse((HERE/name).read_text('utf8'))
    return dict(date='2026-10-04',round=745,scientific_base_through_round=744,
        fresh_tests=dict(run=2,failures=0,errors=0),cumulative_numbered_tests=3428,
        cumulative_numbered_scientific_files=1547,unchanged_prior_evidence_files=3421,
        cumulative_unique_protected_evidence_files=3436,new_file_hashes=new,preserved_draft_hashes=preserved,
        text_checks=checks,local_links_checked=links,broken_links=0,saved_results_reproduced=True,
        previous_results_unchanged=True,full_historical_science_rerun=False,
        previous_round_receipt_and_all_evidence_hashes_checked=True,
        exact_native_coefficients_and_strict_integral_sign_checked=True,
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
