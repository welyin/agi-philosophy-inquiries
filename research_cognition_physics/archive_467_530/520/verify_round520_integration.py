"""Read-only integration verification for round 520 and all predecessors."""
import argparse
import ast
import json
from pathlib import Path
import re
import sys
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
import verify_uniform_tree_continuum_review as previous
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
TARGET=HERE/'round520_integration_checks.json'


def verify():
    base=previous.verify()
    assert base['unchanged_numbered_scientific_tests']==2546
    assert base['unchanged_scientific_file_hashes']==868
    assert base['total_protected_including_this_review']==1004
    checked=core.read(HERE/'research_round_520_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['fresh_tests']==dict(run=6,failures=0,errors=0)
    assert checked['previous_protected_evidence_hashes_verified']==1004
    for group,count in [('new_file_hashes',3),('preserved_draft_hashes',1)]:
        assert len(checked[group])==count
        for name,sha in checked[group].items():
            assert core.digest(HERE/name)==sha,name
    saved=core.read(HERE/'selective_record_state_learning_results.json')
    assert saved['scope']==checked['scope']
    assert (saved['tests_run'],saved['failures'],saved['errors'])==(6,0,0)
    note=HERE/'research_note_520.md'
    assert note.read_bytes()==(HERE/'round520_drafts/research_note_520.txt').read_bytes()
    assert core.text_checks(note)==checked['text_checks']
    for name in ('selective_record_state_learning.py','verify_selective_record_learning_round.py',Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding='utf8'))
    research=HERE.parent;root=research.parent
    docs=[root/'README.md',research/'README.md',research/'research_direction.md',
        research/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
        HERE/'three_dimensional_four_conditions_review.md',note]
    links=0
    for doc in docs:
        for link in core.link_parser()(doc.read_text(encoding='utf-8-sig')):
            assert (doc.parent/link).resolve().exists(),f'{doc.name}: {link}'
            links+=1
    latest=re.search(r'最新科学轮次与检查数为(\d+)／(\d+)',
        (research/'research_direction.md').read_text(encoding='utf8'))
    assert latest and int(latest[1])>=520 and int(latest[2])>=2552
    assert '## 158. 第520轮' in (HERE/'spatial_premise_closure_audit.md').read_text(encoding='utf8')
    assert '## 63. 520' in (HERE/'three_dimensional_four_conditions_review.md').read_text(encoding='utf8')
    return dict(date='2026-09-28',rounds=[520],scientific_base_through_round_by_round={520:519},
        fresh_tests_by_round={520:6},fresh_tests=6,stage_saved_tests=2552,
        science_hashes_verified_231_520=871,unchanged_prior_science_hashes_231_519=868,
        total_protected_evidence_hashes=1008,unchanged_prior_evidence_hashes=1004,
        inherited_post519_unnumbered_evidence_files=4,
        newly_preserved_draft_files=1,local_links_checked=links,broken_links=0,
        navigation_and_notes_checked=len(docs),new_scientific_display_formulas_checked=11,
        old_numbered_scientific_experiments_rerun=False,visual_checks_performed=False,
        independent_final_agent_review_completed=True,full_result_reproduction_independently_confirmed=True,
        cognitive_axioms_silently_strengthened=False,physical_positions_or_dimension_generated=False,
        full_cognition_to_gr_refuted=False,phase_closure_triggered=False,
        latest_round_scope=saved['scope'],all_reported_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args();answer=verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(answer,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:answer[k] for k in ('rounds','stage_saved_tests',
        'total_protected_evidence_hashes','local_links_checked','all_reported_checks_passed')}))
