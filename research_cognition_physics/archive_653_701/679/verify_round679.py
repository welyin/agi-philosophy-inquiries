"""Verify physical-source reflection, singular extension and preserved evidence."""
import argparse
import ast
import json
import re
from pathlib import Path
import verify_interaction_rounds as core
import joint_physical_source_reflection as model
HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_679_checks.json'


def verify():
    base=core.read(HERE/'research_round_678_checks.json')
    assert base['all_reported_checks_passed']
    assert (base['round'],base['cumulative_numbered_tests'],base['cumulative_numbered_scientific_files'])==(678,3263,1346)
    historical=dict(core.read(HERE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
    for n in range(584,679):
        receipt=core.read(HERE/f'research_round_{n}_checks.json')
        for key in ('new_file_hashes','preserved_draft_hashes'):
            for name,digest in receipt[key].items():
                assert name not in historical or historical[name]==digest
                historical[name]=digest
    supplement=core.read(HERE/'cognitive_foundation_bridge_605_navigation.json')
    for name,digest in supplement['supplementary_artifact_hashes'].items():
        assert name not in historical or historical[name]==digest
        historical[name]=digest
    assert len(historical)==2485
    for name,digest in historical.items():assert core.digest(HERE/name)==digest,name
    import joint_local_source_lift as prior_model
    assert prior_model.run()==core.read(prior_model.TARGET)
    result=core.read(model.TARGET)
    assert result==model.run()
    assert (result['tests_run'],result['failures'],result['errors'])==(2,0,0)
    for name,digest in result['dependency_hashes'].items():assert core.digest(HERE/name)==digest,name
    new={name:core.digest(HERE/name) for name in ('research_note_679.md','joint_physical_source_reflection.py','joint_physical_source_reflection_results.json')}
    names=('unified_physics_condition_ledger_679.md', 'round679_drafts/research_note_679_draft.md', 'round679_drafts/final_review.txt', 'round680_drafts/STATUS.md', 'round679_drafts/physical_reflection_probe.py', 'round679_drafts/physical_reflection_probe_results.json', 'round679_drafts/physical_reflection_entry.md', 'round679_drafts/check_entry.py', 'round679_drafts/entry_checks.json', 'round679_drafts/reflection_algebra_probe.py', 'round679_drafts/reflection_algebra_probe_results.json', 'round679_drafts/source_shift_probe.py', 'round679_drafts/source_shift_probe_results.json')
    preserved={name:core.digest(HERE/name) for name in names}
    assert not set(new)&set(historical) and not set(preserved)&set(historical)
    assert len(historical|new|preserved)==2501
    assert (HERE/'research_note_679.md').read_bytes()==(HERE/'round679_drafts/research_note_679_draft.md').read_bytes()
    review=(HERE/'round679_drafts/final_review.txt').read_text('utf8')
    for name in (*new,'unified_physics_condition_ledger_679.md'):assert core.digest(HERE/name) in review,name
    text=core.text_checks(HERE/'research_note_679.md');assert text['display_formulas']==16
    assert not re.search(r'(?<!\\)\b(?:qquad|quad|longrightarrow|operatorname|stackrel)\b', (HERE/'research_note_679.md').read_text('utf8'))
    assert all(c>=32 or c in (10,9) for c in (HERE/'research_note_679.md').read_bytes())
    assert not re.search(r'\\[ \t]+(?:stackrel|frac|sqrt|begin|end|tag|operatorname)\b',(HERE/'research_note_679.md').read_text('utf8'))
    links=0
    for name in ('research_note_679.md','unified_physics_condition_ledger_679.md','round680_drafts/STATUS.md'):
        for link in core.link_parser()((HERE/name).read_text('utf8')):
            path=(HERE/name).parent.joinpath(link).resolve()
            assert path.exists() or path==TARGET.resolve(),(name,link)
            links+=1
    for name in ('joint_physical_source_reflection.py',Path(__file__).name):ast.parse((HERE/name).read_text('utf8'))
    if TARGET.exists():
        old=core.read(TARGET)
        assert old['new_file_hashes']==new and old['preserved_draft_hashes']==preserved
    return dict(date='2026-10-02',round=679,scientific_base_through_round=678,
        fresh_tests=dict(run=2,failures=0,errors=0),cumulative_numbered_tests=3265,
        cumulative_numbered_scientific_files=1349,unchanged_prior_evidence_files=2485,
        cumulative_unique_protected_evidence_files=2501,new_file_hashes=new,preserved_draft_hashes=preserved,
        text_checks=text,local_links_checked=links,broken_links=0,saved_results_reproduced=True,
        previous_results_unchanged=True,full_historical_science_rerun=False,
        previous_round_receipt_and_all_evidence_hashes_checked=True,
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
