"""Freeze the joint bare-coefficient result without changing earlier evidence."""
import argparse
import ast
import json
from pathlib import Path
import joint_singlet_common_mass_rg as model
import verify_interaction_rounds as core
import verify_round543 as previous

HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_544_checks.json'


def verify():
    if TARGET.exists():
        base=previous.verify()
    else:
        import verify_round388_integration as anchor
        original=anchor.verify
        def pending_once(pending_report=False):
            old_target=anchor.TARGET; anchor.TARGET=TARGET
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
            base['cumulative_unique_protected_evidence_files'])==(543,2724,941,1215)
    result=core.read(model.TARGET)
    assert result==model.run()
    assert (result['tests_run'],result['failures'],result['errors'])==(8,0,0)
    for name,sha in result['dependency_hashes'].items():
        assert core.digest(HERE/name)==sha,name
    text=core.text_checks(HERE/'research_note_544.md')
    assert text['display_formulas']==10
    new={name:core.digest(HERE/name) for name in ('research_note_544.md',
        'joint_singlet_common_mass_rg.py','joint_singlet_common_mass_rg_results.json')}
    preserved={name:core.digest(HERE/name) for name in ('unified_physics_condition_ledger_544.md',
        'round544_drafts/independent_review.txt','round544_drafts/final_review.txt',
        'round544_drafts/STATUS.md','round544_drafts/explore_common_mass.py')}
    if TARGET.exists():
        old=core.read(TARGET)
        assert old['new_file_hashes']==new and old['preserved_draft_hashes']==preserved
    links=0
    for name in ('research_note_544.md','unified_physics_condition_ledger_544.md'):
        for link in core.link_parser()((HERE/name).read_text('utf8')):
            path=(HERE/link).resolve()
            assert path.exists() or path==TARGET.resolve(),link
            links+=1
    for name in ('joint_singlet_common_mass_rg.py',Path(__file__).name):
        ast.parse((HERE/name).read_text('utf8'))
    return dict(date='2026-09-30',round=544,scientific_base_through_round=543,
        fresh_tests=dict(run=8,failures=0,errors=0),cumulative_numbered_tests=2732,
        cumulative_numbered_scientific_files=944,unchanged_prior_evidence_files=1215,
        cumulative_unique_protected_evidence_files=1223,new_file_hashes=new,
        preserved_draft_hashes=preserved,text_checks=text,local_links_checked=links,broken_links=0,
        saved_results_reproduced=True,previous_results_unchanged=True,
        inherited_theorems_not_counted_as_new=True,independent_final_code_and_draft_review_completed=True,
        visual_checks_performed=False,active_goal_unchanged=True,scope=result['scope'],
        all_reported_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args(); result=verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert core.read(TARGET)==result
    print(json.dumps({k:result[k] for k in ('round','cumulative_numbered_tests',
        'cumulative_unique_protected_evidence_files','all_reported_checks_passed')}))
