"""Verify699 exact offdiagonal full-Haar and sphere integral and strict negative Gram."""
import argparse
import ast
import json
import re
from pathlib import Path
import verify_interaction_rounds as core
import joint_offdiagonal_haar_certificate as model
HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_699_checks.json'


def verify():
    base=core.read(HERE/'research_round_698_checks.json')
    assert base['all_reported_checks_passed']
    assert (base['round'],base['cumulative_numbered_tests'],base['cumulative_numbered_scientific_files'])==(698,3303,1406)
    historical=dict(core.read(HERE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
    for n in range(584,699):
        receipt=core.read(HERE/f'research_round_{n}_checks.json')
        for key in ('new_file_hashes','preserved_draft_hashes'):
            for name,digest in receipt[key].items():
                assert name not in historical or historical[name]==digest
                historical[name]=digest
    supplement=core.read(HERE/'cognitive_foundation_bridge_605_navigation.json')
    for name,digest in supplement['supplementary_artifact_hashes'].items():
        assert name not in historical or historical[name]==digest
        historical[name]=digest
    assert len(historical)==2761
    for name,digest in historical.items():assert core.digest(HERE/name)==digest,name
    import joint_static_haar_majorant as prior_model
    assert prior_model.run()==core.read(prior_model.TARGET)
    result=core.read(model.TARGET)
    assert result==model.run()
    assert (result['tests_run'],result['failures'],result['errors'])==(2,0,0)
    for name,digest in result['dependency_hashes'].items():assert core.digest(HERE/name)==digest,name
    new={name:core.digest(HERE/name) for name in ('research_note_699.md','joint_offdiagonal_haar_certificate.py','joint_offdiagonal_haar_certificate_results.json')}
    for name,digest in core.read(HERE/'round699_drafts/entry_checks.json')['artifact_hashes'].items():
        assert core.digest(HERE/'round699_drafts'/name)==digest,name
    names=('unified_physics_condition_ledger_699.md', 'round699_drafts/research_note_699_draft.md', 'round699_drafts/final_review.txt', 'round700_drafts/STATUS.md', 'round699_drafts/free_domination_entry.py', 'round699_drafts/free_domination_entry_results.json', 'round699_drafts/free_domination_entry.md', 'round699_drafts/check_and_publish_entry.py', 'round699_drafts/entry_checks.json', 'round699_drafts/prepare_entry_publication.py', 'round699_drafts/clifford_moment_probe.py', 'round699_drafts/clifford_moment_probe_results.json', 'round699_drafts/clifford_moment_certificate.py', 'round699_drafts/offdiagonal_exact_sphere.py', 'round699_drafts/offdiagonal_exact_sphere_results.json', 'round699_drafts/offdiagonal_exact_group.py', 'round699_drafts/offdiagonal_exact_group_results.json')
    preserved={name:core.digest(HERE/name) for name in names}
    assert not set(new)&set(historical) and not set(preserved)&set(historical)
    assert len(historical|new|preserved)==2781
    assert (HERE/'research_note_699.md').read_bytes()==(HERE/'round699_drafts/research_note_699_draft.md').read_bytes()
    review=(HERE/'round699_drafts/final_review.txt').read_text('utf8')
    for name in (*new,'unified_physics_condition_ledger_699.md'):assert core.digest(HERE/name) in review,name
    text=core.text_checks(HERE/'research_note_699.md');assert text['display_formulas']==20
    assert not re.search(r'(?<!\\)\b(?:qquad|quad|longrightarrow|operatorname|stackrel)\b', (HERE/'research_note_699.md').read_text('utf8'))
    assert all(c>=32 or c in (10,13,9) for c in (HERE/'research_note_699.md').read_bytes())
    assert not re.search(r'\\[ \t]+(?:stackrel|frac|sqrt|begin|end|tag|operatorname)\b',(HERE/'research_note_699.md').read_text('utf8'))
    links=0
    for name in ('research_note_699.md','unified_physics_condition_ledger_699.md','round700_drafts/STATUS.md','round687_drafts/attribution_correction_653.md'):
        for link in core.link_parser()((HERE/name).read_text('utf8')):
            path=(HERE/name).parent.joinpath(link).resolve()
            assert path.exists() or path==TARGET.resolve(),(name,link)
            links+=1
    for name in ('joint_offdiagonal_haar_certificate.py',Path(__file__).name):ast.parse((HERE/name).read_text('utf8'))
    if TARGET.exists():
        old=core.read(TARGET)
        assert old['new_file_hashes']==new and old['preserved_draft_hashes']==preserved
    return dict(date='2026-10-02',round=699,scientific_base_through_round=698,
        fresh_tests=dict(run=2,failures=0,errors=0),cumulative_numbered_tests=3305,
        cumulative_numbered_scientific_files=1409,unchanged_prior_evidence_files=2761,
        cumulative_unique_protected_evidence_files=2781,new_file_hashes=new,preserved_draft_hashes=preserved,
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
