"""Read-only science verification, with optional exclusive creation of receipt."""
import argparse
import ast
import json
from pathlib import Path
import verify_interaction_rounds as core
import joint_block_fluctuation_contract as model

HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_708_checks.json'


def verify():
    base=core.read(HERE/'research_round_707_checks.json')
    assert base['all_reported_checks_passed']
    historical=dict(core.read(HERE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
    for n in range(584,708):
        receipt=core.read(HERE/f'research_round_{n}_checks.json')
        for key in ('new_file_hashes','preserved_draft_hashes'):
            for name,digest in receipt[key].items():
                assert name not in historical or historical[name]==digest
                historical[name]=digest
    historical.update(core.read(HERE/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
    assert len(historical)==2900
    for name,digest in historical.items():assert core.digest(HERE/name)==digest,name
    # Only relevant new computation and immediate original full-model dependency.
    import joint_geodesic_spatial_block as prior
    assert prior.run()==core.read(prior.TARGET)
    result=model.run();assert result==core.read(model.TARGET)
    assert (result['tests_run'],result['failures'],result['errors'])==(3,0,0)
    for name,digest in result['dependency_hashes'].items():assert core.digest(HERE/name)==digest,name
    entry=core.read(HERE/'round708_drafts/entry_checks.json')
    for name,digest in entry['artifact_hashes'].items():assert core.digest(HERE/'round708_drafts'/name)==digest,name
    main=('research_note_708.md','joint_block_fluctuation_contract.py','joint_block_fluctuation_contract_results.json')
    names=('unified_physics_condition_ledger_708.md','round708_drafts/research_note_708_draft.md',
           'round708_drafts/final_review.txt','round708_drafts/literature_scope_audit.json',
           'round708_drafts/scope_and_dedup_review.md','round709_drafts/STATUS.md',
           'round708_drafts/nested_mass_entry.py','round708_drafts/nested_mass_entry_results.json',
           'round708_drafts/nested_mass_entry.md','round708_drafts/prepare_entry_publication.py',
           'round708_drafts/check_and_publish_entry.py','round708_drafts/entry_checks.json')
    new={name:core.digest(HERE/name) for name in main}
    preserved={name:core.digest(HERE/name) for name in names}
    assert not (set(new)|set(preserved))&set(historical)
    assert len(historical|new|preserved)==2915
    assert (HERE/main[0]).read_bytes()==(HERE/'round708_drafts/research_note_708_draft.md').read_bytes()
    review=(HERE/'round708_drafts/final_review.txt').read_text('utf8')
    for name in (*main,'unified_physics_condition_ledger_708.md'):assert core.digest(HERE/name) in review
    text=core.text_checks(HERE/main[0]);assert text['display_formulas']==17
    links=0
    for name in (main[0],names[0],'round709_drafts/STATUS.md'):
        for link in core.link_parser()((HERE/name).read_text('utf8')):
            p=(HERE/name).parent.joinpath(link).resolve()
            assert p.exists() or p==TARGET.resolve(),(name,link)
            links+=1
    for name in (main[1],Path(__file__).name,'prepare_round708.py'):
        ast.parse((HERE/name).read_text('utf8'))
    return dict(date='2026-10-03',round=708,scientific_base_through_round=707,
        fresh_tests=dict(run=3,failures=0,errors=0),cumulative_numbered_tests=3326,
        cumulative_numbered_scientific_files=1436,unchanged_prior_evidence_files=2900,
        cumulative_unique_protected_evidence_files=2915,new_file_hashes=new,preserved_draft_hashes=preserved,
        text_checks=text,local_links_checked=links,broken_links=0,saved_results_reproduced=True,
        previous_results_unchanged=True,full_historical_science_rerun=False,
        previous_round_receipt_and_all_evidence_hashes_checked=True,relevant_recent_dependencies_reproduced=True,
        primary_code_and_note_review_completed=True,independent_final_code_and_draft_review_completed=False,
        inherited_results_not_claimed_as_new_theorems=True,visual_checks_performed=False,
        active_goal_unchanged=True,scope=result['scope'],all_reported_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-checks',action='store_true');args=p.parse_args()
    result=verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==core.read(TARGET)
    print(json.dumps({k:result[k] for k in ('round','cumulative_numbered_tests','all_reported_checks_passed')}))
