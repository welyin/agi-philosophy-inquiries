"""Freeze the gauged target connection and common stress reduction."""
import argparse
import ast
from pathlib import Path
import json
import joint_gauged_effective_geometry as model
import verify_interaction_rounds as core
import verify_round581 as previous

HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_582_checks.json'


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
        base['cumulative_unique_protected_evidence_files'])==(581,2964,1055,1580)
    result=core.read(model.TARGET)
    assert result==model.run()
    assert (result['tests_run'],result['failures'],result['errors'])==(6,0,0)
    for name,sha in result['dependency_hashes'].items():assert core.digest(HERE/name)==sha,name
    text=core.text_checks(HERE/'research_note_582.md')
    assert text['display_formulas']==14
    new={name:core.digest(HERE/name) for name in ('research_note_582.md',
        'joint_gauged_effective_geometry.py','joint_gauged_effective_geometry_results.json')}
    preserved={name:core.digest(HERE/name) for name in ('unified_physics_condition_ledger_582.md', 'round582_drafts/STATUS.md', 'round582_drafts/research_note_582_draft.md', 'round582_drafts/entry_progress.md', 'round582_drafts/gauged_target_pre_review.md', 'round582_drafts/gauged_target_entry_probe.py', 'round582_drafts/gauged_target_entry_results.json', 'round582_drafts/gauge_stress_compression_review.md', 'round582_drafts/final_review.txt')}
    assert len(new)==3 and len(preserved)==9
    review=(HERE/'round582_drafts/final_review.txt').read_text('utf8')
    for name in ('joint_gauged_effective_geometry.py','joint_gauged_effective_geometry_results.json',
                 'round582_drafts/research_note_582_draft.md','unified_physics_condition_ledger_582.md'):
        assert core.digest(HERE/name) in review,('reviewed_hash_missing',name)
    assert (HERE/'research_note_582.md').read_bytes()==(HERE/'round582_drafts/research_note_582_draft.md').read_bytes()
    for name in ('research_note_582.md','round582_drafts/research_note_582_draft.md'):
        normalized=(HERE/name).read_bytes().replace(b'\r\n',b'\n')
        assert not any(v<32 and v!=10 for v in normalized)
    if TARGET.exists():
        old=core.read(TARGET)
        assert old['new_file_hashes']==new and old['preserved_draft_hashes']==preserved
    links=0
    for name in ('research_note_582.md','unified_physics_condition_ledger_582.md'):
        for link in core.link_parser()((HERE/name).read_text('utf8')):
            path=(HERE/link).resolve();assert path.exists() or path==TARGET.resolve(),link
            links+=1
    for name in ('joint_gauged_effective_geometry.py',Path(__file__).name):ast.parse((HERE/name).read_text('utf8'))
    return dict(date='2026-10-01',round=582,scientific_base_through_round=581,
        fresh_tests=dict(run=6,failures=0,errors=0),cumulative_numbered_tests=2970,
        cumulative_numbered_scientific_files=1058,unchanged_prior_evidence_files=1580,
        cumulative_unique_protected_evidence_files=1592,new_file_hashes=new,
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
