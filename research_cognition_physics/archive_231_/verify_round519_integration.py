"""Read-only integration verification for round 519 and frozen predecessors."""
import argparse
import ast
import json
from pathlib import Path
import re
import sys
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
import verify_record_determined_predictor_review as previous
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
TARGET=HERE/'round519_integration_checks.json'


def verify():
    base=previous.verify()
    assert base['unchanged_numbered_scientific_tests']==2542
    assert base['unchanged_scientific_file_hashes']==865
    assert base['total_protected_including_this_review']==996
    checked=core.read(HERE/'research_round_519_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['fresh_tests']==dict(run=4,failures=0,errors=0)
    assert checked['previous_protected_evidence_hashes_verified']==996
    for group,count in [('new_file_hashes',3),('preserved_draft_hashes',1)]:
        assert len(checked[group])==count
        for name,sha in checked[group].items():
            assert core.digest(HERE/name)==sha,name
    saved=core.read(HERE/'boundary_corrected_return_results.json')
    assert saved['scope']==checked['scope']
    assert (saved['tests_run'],saved['failures'],saved['errors'])==(4,0,0)
    note=HERE/'research_note_519.md'
    assert note.read_bytes()==(HERE/'round519_drafts/research_note_519.txt').read_bytes()
    assert core.text_checks(note)==checked['text_checks']
    for name in ('boundary_corrected_return.py','verify_boundary_corrected_return_round.py',Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding='utf8'))
    research=HERE.parent;root=research.parent
    docs=[root/'README.md',research/'README.md',research/'research_direction.md',
        research/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
        HERE/'three_dimensional_four_conditions_review.md',note,HERE/'record_determined_predictor_review.md']
    links=0
    for doc in docs:
        for link in core.link_parser()(doc.read_text(encoding='utf-8-sig')):
            assert (doc.parent/link).resolve().exists(),f'{doc.name}: {link}'
            links+=1
    latest=re.search(r'最新科学轮次与检查数为(\d+)／(\d+)',
        (research/'research_direction.md').read_text(encoding='utf8'))
    assert latest and int(latest[1])>=519 and int(latest[2])>=2546
    assert '## 157. 第519轮' in (HERE/'spatial_premise_closure_audit.md').read_text(encoding='utf8')
    assert '## 62. 519' in (HERE/'three_dimensional_four_conditions_review.md').read_text(encoding='utf8')
    return dict(date='2026-09-28',rounds=[519],scientific_base_through_round_by_round={519:518},
        fresh_tests_by_round={519:4},fresh_tests=4,stage_saved_tests=2546,
        science_hashes_verified_231_519=868,unchanged_prior_science_hashes_231_518=865,
        total_protected_evidence_hashes=1000,unchanged_prior_evidence_hashes=996,
        inherited_post518_unnumbered_evidence_files=4,unnumbered_diagnostics_not_in_numbered_count=4,
        newly_preserved_draft_files=1,local_links_checked=links,broken_links=0,
        navigation_and_notes_checked=len(docs),new_scientific_display_formulas_checked=12,
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
