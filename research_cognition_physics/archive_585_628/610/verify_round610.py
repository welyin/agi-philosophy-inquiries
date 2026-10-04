"""Verify projected process completion, full histories and common sources."""
import argparse
import ast
import json
import re
from pathlib import Path
import verify_interaction_rounds as core
import verify_round609 as previous
import joint_gapped_link_locality as model
HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_610_checks.json'


def verify():
    base=previous.verify()
    assert (base['round'],base['cumulative_numbered_tests'],base['cumulative_numbered_scientific_files'])==(609,3072,1139)
    historical=dict(core.read(HERE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
    for n in range(584,610):
        receipt=core.read(HERE/f'research_round_{n}_checks.json')
        for key in ('new_file_hashes','preserved_draft_hashes'):
            for name,digest in receipt[key].items():
                assert name not in historical or historical[name]==digest
                historical[name]=digest
    supplement=core.read(HERE/'cognitive_foundation_bridge_605_navigation.json')
    for name,digest in supplement['supplementary_artifact_hashes'].items():
        assert name not in historical or historical[name]==digest
        historical[name]=digest
    assert len(historical)==1856
    for name,digest in historical.items():assert core.digest(HERE/name)==digest,name
    result=core.read(model.TARGET)
    assert result==model.run()
    assert (result['tests_run'],result['failures'],result['errors'])==(3,0,0)
    for name,digest in result['dependency_hashes'].items():assert core.digest(HERE/name)==digest,name
    new={name:core.digest(HERE/name) for name in ('research_note_610.md','joint_gapped_link_locality.py','joint_gapped_link_locality_results.json')}
    names=('unified_physics_condition_ledger_610.md','round610_drafts/research_note_610_draft.md',
           'round610_drafts/final_review.txt','round611_drafts/STATUS.md')
    preserved={name:core.digest(HERE/name) for name in names}
    assert not set(new)&set(historical) and not set(preserved)&set(historical)
    assert len(historical|new|preserved)==1863
    assert (HERE/'research_note_610.md').read_bytes()==(HERE/'round610_drafts/research_note_610_draft.md').read_bytes()
    review=(HERE/'round610_drafts/final_review.txt').read_text('utf8')
    for name in (*new,'unified_physics_condition_ledger_610.md'):assert core.digest(HERE/name) in review,name
    text=core.text_checks(HERE/'research_note_610.md');assert text['display_formulas']==12
    assert not re.search(r'(?<![A-Za-z\\])(qquad|quad)\b', (HERE/'research_note_610.md').read_text('utf8'))
    links=0
    for name in ('research_note_610.md','unified_physics_condition_ledger_610.md','round611_drafts/STATUS.md'):
        for link in core.link_parser()((HERE/name).read_text('utf8')):
            path=(HERE/name).parent.joinpath(link).resolve()
            assert path.exists() or path==TARGET.resolve(),(name,link)
            links+=1
    for name in ('joint_gapped_link_locality.py',Path(__file__).name):ast.parse((HERE/name).read_text('utf8'))
    if TARGET.exists():
        old=core.read(TARGET)
        assert old['new_file_hashes']==new and old['preserved_draft_hashes']==preserved
    return dict(date='2026-10-01',round=610,scientific_base_through_round=609,
        fresh_tests=dict(run=3,failures=0,errors=0),cumulative_numbered_tests=3075,
        cumulative_numbered_scientific_files=1142,unchanged_prior_evidence_files=1856,
        cumulative_unique_protected_evidence_files=1863,new_file_hashes=new,preserved_draft_hashes=preserved,
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
