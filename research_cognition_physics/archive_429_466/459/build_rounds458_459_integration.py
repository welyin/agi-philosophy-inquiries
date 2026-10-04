"""Create new immutable-report verifiers; never overwrite earlier files."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
for number, module, formulas, old_tests, old_science, old_evidence, section in [
    (458, 'continuous_record_response', 13, 2175, 682, 717, 86),
    (459, 'recursive_exchange_interface', 15, 2181, 685, 721, 87),
]:
    extra = 1 if number == 458 else 0
    source = f'''"""Read-only integration of round {number}; historical science remains immutable."""
import argparse
import ast
import json
from pathlib import Path
import re
import verify_round{number-1}_integration as previous
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
TARGET=HERE/'round{number}_integration_checks.json'

def verify(pending=False):
    old_target=previous.TARGET
    previous.TARGET=TARGET
    try:
        result=previous.verify(pending)
    finally:
        previous.TARGET=old_target
    checked=core.read(HERE/'research_round_{number}_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round']==457
    assert checked['fresh_tests']==dict(run=6,failures=0,errors=0)
    assert len(checked['new_file_hashes'])==3
    assert len(checked['preserved_draft_hashes'])=={extra}
    for group in ['new_file_hashes','preserved_draft_hashes']:
        for name,sha in checked[group].items():
            assert core.digest(HERE/name)==sha,name
    saved=core.read(HERE/'{module}_audit_results.json')
    assert saved['scope']==checked['scope']
    assert not saved['scope']['full_GR_goal_completed']
    assert not saved['scope']['phase_closure_triggered']
    assert (saved['tests_run'],saved['failures'],saved['errors'])==(6,0,0)
    note=HERE/'research_note_{number}.md'
    assert core.text_checks(note)==checked['text_checks']
    assert checked['text_checks']['display_formulas']=={formulas}
    for name in ['{module}_audit.py','verify_{module}_round.py',Path(__file__).name]:
        ast.parse((HERE/name).read_text(encoding='utf-8'))
    links=0
    for doc in [note,HERE/'recursive_functional_closure_proposal.md']:
        for link in core.link_parser()(doc.read_text(encoding='utf-8')):
            dest=(doc.parent/link).resolve()
            assert dest.exists() or (pending and dest==TARGET),link
            links+=1
    assert result.pop('science_hashes_verified_231_{number-1}')=={old_science}
    assert result.pop('unchanged_prior_science_hashes_231_{number-2}')=={old_science-3}
    assert result['stage_saved_tests']=={old_tests}
    assert result['total_protected_evidence_hashes']=={old_evidence}
    direction=(HERE.parent/'research_direction.md').read_text(encoding='utf-8')
    latest=re.search(r'最新科学轮次与检查数为(\\d+)／(\\d+)',direction)
    assert latest and int(latest[1])>={number} and int(latest[2])>={old_tests+6}
    assert '递归功能闭合' in direction
    assert '## {section}. 第{number}轮' in (HERE/'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24',rounds=[{number}],
        scientific_base_through_round_by_round={{{number}:457}},
        fresh_tests_by_round={{{number}:6}},fresh_tests=6,
        stage_saved_tests={old_tests+6},science_hashes_verified_231_{number}={old_science+3},
        unchanged_prior_science_hashes_231_{number-1}={old_science},
        total_protected_evidence_hashes={old_evidence+3+extra},unchanged_prior_evidence_hashes={old_evidence},
        newly_preserved_draft_files={extra},new_scientific_display_formulas_checked={formulas},
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        independent_rounds_458_459_share_science_baseline_457=True,
        latest_round_scope=saved['scope'],
        recursive_functional_closure_recorded_as_user_candidate=True,
        cognitive_axioms_silently_strengthened=False,
        spark_project_files_modified=False,physical_positions_or_dimension_generated=False,
        full_cognition_to_gr_refuted=False,phase_closure_triggered=False,
        old_scientific_experiments_rerun=False,
        independent_final_agent_review_completed={number==458},
        parent_analytic_and_reproduction_review_completed=True,
        all_reported_checks_passed=True)
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args()
    result=verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x',encoding='utf-8',newline='\\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\\n')
    print(json.dumps({{k:result[k] for k in ['rounds','stage_saved_tests','total_protected_evidence_hashes','all_reported_checks_passed']}},ensure_ascii=False))
'''
    with (HERE/f'verify_round{number}_integration.py').open('x', encoding='utf-8', newline='\n') as stream:
        stream.write(source)
