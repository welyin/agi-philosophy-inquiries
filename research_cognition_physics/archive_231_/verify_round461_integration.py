"""Read-only integration of round 461; historical science remains immutable."""
import argparse
import ast
import json
from pathlib import Path
import re
import verify_round460_integration as previous
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
TARGET=HERE/'round461_integration_checks.json'

def verify(pending=False):
    old_target=previous.TARGET
    previous.TARGET=TARGET
    try:
        result=previous.verify(pending)
    finally:
        previous.TARGET=old_target
    checked=core.read(HERE/'research_round_461_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round']==459
    assert checked['fresh_tests']==dict(run=6,failures=0,errors=0)
    assert len(checked['new_file_hashes'])==3
    assert len(checked['preserved_draft_hashes'])==1
    for group in ['new_file_hashes','preserved_draft_hashes']:
        for name,sha in checked[group].items():
            assert core.digest(HERE/name)==sha,name
    saved=core.read(HERE/'independent_subject_admission_audit_results.json')
    assert saved['scope']==checked['scope']
    assert not saved['scope']['full_GR_goal_completed']
    assert not saved['scope']['phase_closure_triggered']
    assert (saved['tests_run'],saved['failures'],saved['errors'])==(6,0,0)
    note=HERE/'research_note_461.md'
    assert core.text_checks(note)==checked['text_checks']
    assert checked['text_checks']['display_formulas']==14
    for name in ['independent_subject_admission_audit.py','verify_independent_subject_admission_round.py',Path(__file__).name]:
        ast.parse((HERE/name).read_text(encoding='utf-8'))
    links=0
    for doc in [note,HERE/'recursive_functional_closure_proposal.md']:
        for link in core.link_parser()(doc.read_text(encoding='utf-8')):
            dest=(doc.parent/link).resolve()
            assert dest.exists() or (pending and dest==TARGET),link
            links+=1
    assert result.pop('science_hashes_verified_231_460')==691
    assert result.pop('unchanged_prior_science_hashes_231_459')==688
    assert result['stage_saved_tests']==2193
    assert result['total_protected_evidence_hashes']==727
    direction=(HERE.parent/'research_direction.md').read_text(encoding='utf-8')
    latest=re.search(r'最新科学轮次与检查数为(\d+)／(\d+)',direction)
    assert latest and int(latest[1])>=461 and int(latest[2])>=2199
    assert '递归功能闭合' in direction
    assert '## 89. 第461轮' in (HERE/'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24',rounds=[461],
        scientific_base_through_round_by_round={461:459},
        fresh_tests_by_round={461:6},fresh_tests=6,
        stage_saved_tests=2199,science_hashes_verified_231_461=694,
        unchanged_prior_science_hashes_231_460=691,
        total_protected_evidence_hashes=731,unchanged_prior_evidence_hashes=727,
        newly_preserved_draft_files=1,new_scientific_display_formulas_checked=14,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        independent_rounds_460_461_share_science_baseline_459=True,
        latest_round_scope=saved['scope'],
        recursive_functional_closure_recorded_as_user_candidate=True,
        cognitive_axioms_silently_strengthened=False,
        spark_project_files_modified=False,physical_positions_or_dimension_generated=False,
        full_cognition_to_gr_refuted=False,phase_closure_triggered=False,
        old_scientific_experiments_rerun=False,
        independent_final_agent_review_completed=False,
        parent_analytic_and_reproduction_review_completed=True,
        all_reported_checks_passed=True)
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args()
    result=verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x',encoding='utf-8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['rounds','stage_saved_tests','total_protected_evidence_hashes','all_reported_checks_passed']},ensure_ascii=False))
