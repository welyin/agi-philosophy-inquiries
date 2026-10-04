"""Verify original spatial auxiliary measure, shared propagation and response."""
import argparse
import ast
import json
import re
from pathlib import Path
import verify_interaction_rounds as core
import joint_spatial_auxiliary_geometry as model
HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_656_checks.json'


def verify():
    base=core.read(HERE/'research_round_655_checks.json')
    assert base['all_reported_checks_passed']
    assert (base['round'],base['cumulative_numbered_tests'],base['cumulative_numbered_scientific_files'])==(655,3206,1277)
    historical=dict(core.read(HERE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
    for n in range(584,656):
        receipt=core.read(HERE/f'research_round_{n}_checks.json')
        for key in ('new_file_hashes','preserved_draft_hashes'):
            for name,digest in receipt[key].items():
                assert name not in historical or historical[name]==digest
                historical[name]=digest
    supplement=core.read(HERE/'cognitive_foundation_bridge_605_navigation.json')
    for name,digest in supplement['supplementary_artifact_hashes'].items():
        assert name not in historical or historical[name]==digest
        historical[name]=digest
    assert len(historical)==2218
    for name,digest in historical.items():assert core.digest(HERE/name)==digest,name
    import importlib.util
    entry_spec=importlib.util.spec_from_file_location('entry656',HERE/'round656_drafts/spatial_auxiliary_probe.py')
    entry=importlib.util.module_from_spec(entry_spec);entry_spec.loader.exec_module(entry)
    assert entry.run()==core.read(HERE/'round656_drafts/spatial_auxiliary_probe_results.json')
    result=core.read(model.TARGET)
    assert result==model.run()
    assert (result['tests_run'],result['failures'],result['errors'])==(4,0,0)
    for name,digest in result['dependency_hashes'].items():assert core.digest(HERE/name)==digest,name
    new={name:core.digest(HERE/name) for name in ('research_note_656.md','joint_spatial_auxiliary_geometry.py','joint_spatial_auxiliary_geometry_results.json')}
    names=('unified_physics_condition_ledger_656.md','round656_drafts/research_note_656_draft.md',
           'round656_drafts/final_review.txt','round657_drafts/STATUS.md',
           'round656_drafts/actual_spatial_measure_entry.md','round656_drafts/spatial_auxiliary_probe.py',
           'round656_drafts/spatial_auxiliary_probe_results.json',
           'round656_drafts/initial_planar_results.json')
    preserved={name:core.digest(HERE/name) for name in names}
    assert not set(new)&set(historical) and not set(preserved)&set(historical)
    assert len(historical|new|preserved)==2229
    assert (HERE/'research_note_656.md').read_bytes()==(HERE/'round656_drafts/research_note_656_draft.md').read_bytes()
    review=(HERE/'round656_drafts/final_review.txt').read_text('utf8')
    for name in (*new,'unified_physics_condition_ledger_656.md'):assert core.digest(HERE/name) in review,name
    text=core.text_checks(HERE/'research_note_656.md');assert text['display_formulas']==20
    assert not re.search(r'(?<!\\)\b(?:qquad|quad|longrightarrow|operatorname|stackrel)\b', (HERE/'research_note_656.md').read_text('utf8'))
    assert all(c>=32 or c in (10,9) for c in (HERE/'research_note_656.md').read_bytes())
    assert not re.search(r'\\[ \t]+(?:stackrel|frac|sqrt|begin|end|tag|operatorname)\b',(HERE/'research_note_656.md').read_text('utf8'))
    links=0
    for name in ('research_note_656.md','unified_physics_condition_ledger_656.md','round657_drafts/STATUS.md',
                 'round656_drafts/actual_spatial_measure_entry.md'):
        for link in core.link_parser()((HERE/name).read_text('utf8')):
            path=(HERE/name).parent.joinpath(link).resolve()
            assert path.exists() or path==TARGET.resolve(),(name,link)
            links+=1
    for name in ('joint_spatial_auxiliary_geometry.py',Path(__file__).name):ast.parse((HERE/name).read_text('utf8'))
    if TARGET.exists():
        old=core.read(TARGET)
        assert old['new_file_hashes']==new and old['preserved_draft_hashes']==preserved
    return dict(date='2026-10-02',round=656,scientific_base_through_round=655,
        fresh_tests=dict(run=4,failures=0,errors=0),cumulative_numbered_tests=3210,
        cumulative_numbered_scientific_files=1280,unchanged_prior_evidence_files=2218,
        cumulative_unique_protected_evidence_files=2229,new_file_hashes=new,preserved_draft_hashes=preserved,
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
