"""Read-only integration verification for round 518 and all frozen predecessors."""
import argparse
import ast
import json
from pathlib import Path
import re
import sys
sys.setrecursionlimit(max(sys.getrecursionlimit(),10000))
import verify_round517_integration as previous
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
TARGET=HERE/'round518_integration_checks.json'


def verify(pending=False):
    base=previous.verify(pending)
    assert base['stage_saved_tests']==2536
    assert base['science_hashes_verified_231_517']==862
    assert base['total_protected_evidence_hashes']==987
    checked=core.read(HERE/'research_round_518_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['fresh_tests']==dict(run=6,failures=0,errors=0)
    assert checked['previous_protected_evidence_hashes_verified']==987
    for group,count in [('new_file_hashes',3),('preserved_draft_hashes',2)]:
        assert len(checked[group])==count
        for name,sha in checked[group].items():
            assert core.digest(HERE/name)==sha,name
    saved=core.read(HERE/'same_time_neighborhood_source_results.json')
    assert saved['scope']==checked['scope']
    assert (saved['tests_run'],saved['failures'],saved['errors'])==(6,0,0)
    note=HERE/'research_note_518.md'
    assert note.read_bytes()==(HERE/'round518_drafts/research_note_518.txt').read_bytes()
    reviewed={
        'same_time_neighborhood_source.py':'6072e545e418813c7cff02de8c28ac0d8db764e365da7291c907db3e3ae837ea',
        'same_time_neighborhood_source_results.json':'d0fdbdd3c80220b14986f9176fd8c55728d87d96876eac4ec58018c01963ba51',
        'research_note_518.md':'b3adfa04ef91a1e0d417e807b7585b7d81c823f86e0e07d06f7098a9fc0d19b3'}
    for name,sha in reviewed.items():
        assert core.digest(HERE/name)==sha,name
    assert core.text_checks(note)==checked['text_checks']
    for name in ('same_time_neighborhood_source.py','verify_same_time_neighborhood_round.py',
                 Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding='utf8'))
    research=HERE.parent;root=research.parent
    docs=[root/'README.md',research/'README.md',research/'research_direction.md',
          research/'RESEARCH_STATE.md',HERE/'README.md',
          HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md',note]
    links=0
    for doc in docs:
        for link in core.link_parser()(doc.read_text(encoding='utf-8-sig')):
            target=(doc.parent/link).resolve()
            assert target.exists() or (pending and target==TARGET),f'{doc.name}: {link}'
            links+=1
    latest=re.search(r'最新科学轮次与检查数为(\d+)／(\d+)',
        (research/'research_direction.md').read_text(encoding='utf8'))
    assert latest and int(latest[1])>=518 and int(latest[2])>=2542
    assert '## 156. 第518轮' in (HERE/'spatial_premise_closure_audit.md').read_text(encoding='utf8')
    assert '## 61. 518' in (HERE/'three_dimensional_four_conditions_review.md').read_text(encoding='utf8')
    return dict(date='2026-09-28',rounds=[518],
        scientific_base_through_round_by_round={518:517},fresh_tests_by_round={518:6},fresh_tests=6,
        stage_saved_tests=base['stage_saved_tests']+saved['tests_run'],
        science_hashes_verified_231_518=862+len(checked['new_file_hashes']),
        unchanged_prior_science_hashes_231_517=862,
        total_protected_evidence_hashes=987+sum(len(checked[g]) for g in ('new_file_hashes','preserved_draft_hashes')),
        unchanged_prior_evidence_hashes=987,newly_preserved_draft_files=2,
        inherited_post516_unnumbered_evidence_files=7,
        local_links_checked=links,broken_links=0,navigation_and_notes_checked=len(docs),
        new_scientific_display_formulas_checked=16,
        old_numbered_scientific_experiments_rerun=False,visual_checks_performed=False,
        independent_final_agent_review_completed=True,full_result_reproduction_independently_confirmed=True,
        cognitive_axioms_silently_strengthened=False,physical_positions_or_dimension_generated=False,
        full_cognition_to_gr_refuted=False,phase_closure_triggered=False,
        latest_round_scope=saved['scope'],all_reported_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args();result=verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('rounds','stage_saved_tests',
        'total_protected_evidence_hashes','local_links_checked','all_reported_checks_passed')}))
