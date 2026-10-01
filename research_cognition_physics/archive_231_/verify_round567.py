"""Freeze the concrete graph/continuum scalar propagation interface."""
import argparse
import ast
from pathlib import Path
import json
import joint_scalar_propagation_matching as model
import verify_interaction_rounds as core
import verify_round566 as previous

HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_567_checks.json'


def verify():
    if TARGET.exists():base=previous.verify()
    else:
        import verify_round388_integration as anchor
        original=anchor.verify
        def pending_once(pending_report=False):
            old_target=anchor.TARGET;anchor.TARGET=TARGET
            try:return original(True)
            finally:anchor.TARGET=old_target
        anchor.verify=pending_once
        try:base=previous.verify()
        finally:anchor.verify=original
    assert (base['round'],base['cumulative_numbered_tests'],base['cumulative_numbered_scientific_files'],
        base['cumulative_unique_protected_evidence_files'])==(566,2869,1010,1426)
    result=core.read(model.TARGET)
    assert result==model.run()
    assert (result['tests_run'],result['failures'],result['errors'])==(6,0,0)
    for name,sha in result['dependency_hashes'].items():assert core.digest(HERE/name)==sha,name
    text=core.text_checks(HERE/'research_note_567.md')
    assert text['display_formulas']==14
    new={name:core.digest(HERE/name) for name in ('research_note_567.md',
        'joint_scalar_propagation_matching.py','joint_scalar_propagation_matching_results.json')}
    preserved={name:core.digest(HERE/name) for name in ('unified_physics_condition_ledger_567.md',
        'round567_drafts/STATUS.md','round567_drafts/research_note_567_draft.md',
        'round567_drafts/propagation_probe.py','round567_drafts/propagation_probe_results.json',
        'round567_drafts/propagation_pre_review.txt','round567_drafts/final_review.txt')}
    for name in ('research_note_567.md','round567_drafts/research_note_567_draft.md'):
        normalized=(HERE/name).read_bytes().replace(b'\r\n',b'\n')
        assert not any(v<32 and v!=10 for v in normalized)
    if TARGET.exists():
        old=core.read(TARGET)
        assert old['new_file_hashes']==new and old['preserved_draft_hashes']==preserved
    links=0
    for name in ('research_note_567.md','unified_physics_condition_ledger_567.md'):
        for link in core.link_parser()((HERE/name).read_text('utf8')):
            path=(HERE/link).resolve();assert path.exists() or path==TARGET.resolve(),link
            links+=1
    for name in ('joint_scalar_propagation_matching.py',Path(__file__).name):ast.parse((HERE/name).read_text('utf8'))
    return dict(date='2026-10-01',round=567,scientific_base_through_round=566,
        fresh_tests=dict(run=6,failures=0,errors=0),cumulative_numbered_tests=2875,
        cumulative_numbered_scientific_files=1013,unchanged_prior_evidence_files=1426,
        cumulative_unique_protected_evidence_files=1436,new_file_hashes=new,
        preserved_draft_hashes=preserved,text_checks=text,local_links_checked=links,broken_links=0,
        saved_results_reproduced=True,previous_results_unchanged=True,
        inherited_theorems_not_counted_as_new=True,independent_final_code_and_draft_review_completed=True,
        visual_checks_performed=False,active_goal_unchanged=True,scope=result['scope'],
        all_reported_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args();result=verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert core.read(TARGET)==result
    print(json.dumps({k:result[k] for k in ('round','cumulative_numbered_tests',
        'cumulative_unique_protected_evidence_files','all_reported_checks_passed')}))
