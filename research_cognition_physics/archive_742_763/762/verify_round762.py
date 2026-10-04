"""Reproduce 762 and verify all preceding frozen evidence and entry diagnostics."""
import argparse
import ast
import json
from pathlib import Path
import sys
import verify_interaction_rounds as core
import joint_gauss_first_order_process as model
HERE=Path(__file__).resolve().parent
TARGET=HERE/'research_round_762_checks.json'

def verify():
    base=core.read(HERE/'research_round_761_checks.json')
    assert base['all_reported_checks_passed']
    history=dict(core.read(HERE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
    for number in range(584,762):
        receipt=core.read(HERE/f'research_round_{number}_checks.json')
        for key in ('new_file_hashes','preserved_draft_hashes'):
            for name,digest in receipt[key].items():
                assert name not in history or history[name]==digest,name
                history[name]=digest
    history.update(core.read(HERE/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
    assert len(history)==3637
    for name,digest in history.items():assert core.digest(HERE/name)==digest,name
    sys.path.insert(0,str(HERE/'round762_drafts'))
    import verify_resource_scale_entry as entry
    assert entry.verify()==core.read(entry.TARGET)
    result=model.run()
    assert result==core.read(model.TARGET)
    assert (result['tests_run'],result['failures'],result['errors'])==(3,0,0)
    for name,digest in result['dependency_hashes'].items():assert core.digest(HERE/name)==digest
    main=('research_note_762.md','joint_gauss_first_order_process.py','joint_gauss_first_order_process_results.json')
    names=('unified_physics_condition_ledger_762.md','round762_drafts/research_note_762_draft.md',
           'round762_drafts/final_review.txt','round762_drafts/literature_scope_audit.json',
           'round762_drafts/scope_and_dedup_review.md','round763_drafts/STATUS.md',
           'round762_drafts/research_note_762_working.md','round762_drafts/resource_scale_entry.py',
           'round762_drafts/resource_scale_entry_results.json','round762_drafts/verify_resource_scale_entry.py',
           'round762_drafts/publish_resource_scale_entry.py','round762_drafts/resource_scale_entry_checks.json',
           'round762_drafts/resource_scale_entry_navigation_checks.json')
    new={n:core.digest(HERE/n) for n in main}
    preserved={n:core.digest(HERE/n) for n in names}
    assert not (set(new)|set(preserved))&set(history)
    assert len(history|new|preserved)==3653
    assert (HERE/main[0]).read_bytes()==(HERE/names[1]).read_bytes()
    review=(HERE/'round762_drafts/final_review.txt').read_text('utf8')
    for name in (*main,names[0]):assert core.digest(HERE/name) in review
    checks=core.text_checks(HERE/main[0])
    assert checks['display_formulas']==15
    links=0
    for name in (main[0],names[0],'round763_drafts/STATUS.md'):
        for link in core.link_parser()((HERE/name).read_text('utf8')):
            p=(HERE/name).parent.joinpath(link).resolve()
            assert p.exists() or p==TARGET.resolve(),(name,link)
            links+=1
    for name in (main[1],Path(__file__).name,'prepare_round762.py',
                 'publish_round762.py','postcheck_round762.py'):
        ast.parse((HERE/name).read_text('utf8'))
    return dict(date='2026-10-04',round=762,scientific_base_through_round=761,
        fresh_tests=dict(run=3,failures=0,errors=0),cumulative_numbered_tests=3473,
        cumulative_numbered_scientific_files=1598,unchanged_prior_evidence_files=3637,
        cumulative_unique_protected_evidence_files=3653,new_file_hashes=new,preserved_draft_hashes=preserved,
        text_checks=checks,local_links_checked=links,broken_links=0,saved_results_reproduced=True,
        previous_results_unchanged=True,full_historical_science_rerun=False,
        previous_round_receipt_and_all_evidence_hashes_checked=True,
        working_entry_reproduced_without_double_counting=True,
        gauss_packet_first_jet_and_matrix_order_checked=True,
        original_scalar_instrument_and_conditional_source_checked=True,
        fixed_graph_weak_expansion_only=True,first_order_trace_norm_state_theorem_proven=False,
        full_continuum_process_equivalence_proven=False,
        primary_code_and_note_review_completed=True,independent_final_code_and_draft_review_completed=False,
        inherited_results_not_claimed_as_new_theorems=True,visual_checks_performed=False,
        active_goal_unchanged=True,scope=result['scope'],all_reported_checks_passed=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args();result=verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==core.read(TARGET)
    print(json.dumps({k:result[k] for k in ('round','cumulative_numbered_tests','all_reported_checks_passed')}))

