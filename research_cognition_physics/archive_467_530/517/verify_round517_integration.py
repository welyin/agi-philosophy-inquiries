"""Read-only round 517 integration including both post-516 scope audits."""
import argparse
import ast
import json
from pathlib import Path
import re
import sys
sys.setrecursionlimit(max(sys.getrecursionlimit(),10000))
import verify_qca_finite_region_scope as previous
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
TARGET=HERE/'round517_integration_checks.json'


def verify(pending=False):
    base=previous.verify()
    assert base['stage_saved_tests']==2531
    assert base['science_hashes_verified_231_516']==859
    assert base['total_protected_evidence_hashes']==982
    checked=core.read(HERE/'research_round_517_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['fresh_tests']==dict(run=5,failures=0,errors=0)
    assert checked['previous_protected_evidence_hashes_verified']==982
    for group,count in [('new_file_hashes',3),('preserved_draft_hashes',2)]:
        assert len(checked[group])==count
        for name,sha in checked[group].items():
            assert core.digest(HERE/name)==sha,name
    saved=core.read(HERE/'all_scale_monitored_graph_source_results.json')
    assert saved['scope']==checked['scope']
    assert (saved['tests_run'],saved['failures'],saved['errors'])==(5,0,0)
    note=HERE/'research_note_517.md'
    assert note.read_bytes()==(HERE/'round517_drafts/research_note_517.txt').read_bytes()
    assert core.digest(note)=='d60a182c3910d57933d9574db07c2ffe4cd7b2056dde784d5e7d49557780b28c'
    assert core.text_checks(note)==checked['text_checks']
    for name in ('all_scale_monitored_graph_source.py','verify_all_scale_graph_source_round.py',
                 Path(__file__).name,'verify_qca_finite_region_scope.py'):
        ast.parse((HERE/name).read_text(encoding='utf8'))
    research=HERE.parent;root=research.parent
    docs=[root/'README.md',research/'README.md',research/'research_direction.md',
          research/'RESEARCH_STATE.md',HERE/'README.md',
          HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md',
          note,HERE/'qca_finite_region_scope_addendum.md']
    links=0
    for doc in docs:
        body=doc.read_text(encoding='utf-8-sig')
        for link in core.link_parser()(body):
            target=(doc.parent/link).resolve()
            assert target.exists() or (pending and target==TARGET),f'{doc.name}: {link}'
            links+=1
    latest=re.search(r'最新科学轮次与检查数为(\d+)／(\d+)',
        (research/'research_direction.md').read_text(encoding='utf8'))
    assert latest and int(latest[1])>=517 and int(latest[2])>=2536
    assert '## 155. 第517轮' in (HERE/'spatial_premise_closure_audit.md').read_text(encoding='utf8')
    assert '## 60. 517' in (HERE/'three_dimensional_four_conditions_review.md').read_text(encoding='utf8')
    count=sum(len(checked[g]) for g in ('new_file_hashes','preserved_draft_hashes'))
    return dict(date='2026-09-28',rounds=[517],
        scientific_base_through_round_by_round={517:516},fresh_tests_by_round={517:5},fresh_tests=5,
        stage_saved_tests=base['stage_saved_tests']+saved['tests_run'],
        science_hashes_verified_231_517=859+len(checked['new_file_hashes']),
        unchanged_prior_science_hashes_231_516=859,
        total_protected_evidence_hashes=base['total_protected_evidence_hashes']+count,
        unchanged_prior_evidence_hashes=982,newly_preserved_draft_files=2,
        inherited_post516_unnumbered_evidence_files=7,
        local_links_checked=links,broken_links=0,navigation_and_notes_checked=len(docs),
        new_scientific_display_formulas_checked=12,
        old_numbered_scientific_experiments_rerun=False,visual_checks_performed=False,
        independent_final_agent_review_completed=True,
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
