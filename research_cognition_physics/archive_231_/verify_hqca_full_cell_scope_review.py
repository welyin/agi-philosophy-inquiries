"""Verify unnumbered whole-cell HQCA audit without rerunning legacy science."""
import argparse
import ast
import json
from pathlib import Path
import sys

sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
import verify_round516_integration as previous
import verify_interaction_rounds as core
import hqca_full_cell_scope_probe as probe

HERE=Path(__file__).resolve().parent
TARGET=HERE/'hqca_full_cell_scope_checks.json'
FILES=['hqca_full_cell_scope_probe.py','hqca_full_cell_scope_probe_results.json',
       'hqca_full_cell_scope_review.md','hqca_full_cell_scope_review_draft.txt',
       'hqca_full_cell_scope_probe_before_reference_check.txt']


def verify():
    base=previous.verify()
    assert base['stage_saved_tests']==2531
    assert base['science_hashes_verified_231_516']==859
    assert base['total_protected_evidence_hashes']==975
    if TARGET.exists():
        for name,sha in core.read(TARGET)['new_unnumbered_evidence_hashes'].items():
            assert core.digest(HERE/name)==sha,name
    assert (HERE/FILES[2]).read_bytes()==(HERE/FILES[3]).read_bytes()
    assert core.digest(HERE/FILES[0])=='c1bb343eecd2d895da818f05f70b457dd15a656ddb41c13667acf181cf84a16a'
    assert core.digest(HERE/FILES[2])=='6e83b048999a557b50894a46ec7b8ff521afa8c69eb9d8f3dabff7215a840d88'
    for name in (FILES[0],Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding='utf8'))
    saved=core.read(HERE/FILES[1])
    assert saved==json.loads(json.dumps(probe.run()))
    assert (saved['diagnostic_tests'],saved['failures'],saved['errors'])==(3,0,0)
    assert saved['scientific_baseline_round']==516
    assert saved['numbered_round_created'] is False
    assert saved['numbered_scientific_test_increment']==0
    for name,sha in saved['dependency_sha256'].items():
        assert core.digest(HERE/name)==sha,name
    checked=core.text_checks(HERE/FILES[2])
    assert checked['display_formulas']==8
    links=0
    for link in core.link_parser()((HERE/FILES[2]).read_text(encoding='utf8')):
        assert (HERE/link).resolve().exists(),link
        links+=1
    assert saved['scope']['all_cognitive_models_refuted'] is False
    assert saved['scope']['full_GR_goal_completed'] is False
    return dict(date='2026-09-28',scientific_base_through_round=516,
        unchanged_numbered_scientific_tests=2531,unchanged_scientific_file_hashes=859,
        unchanged_previously_protected_evidence_hashes=975,
        numbered_round_created=False,numbered_test_increment=0,
        independent_diagnostic_checks=dict(run=3,failures=0,errors=0),
        new_unnumbered_evidence_files=len(FILES),
        total_protected_including_this_review=975+len(FILES),
        new_unnumbered_evidence_hashes={name:core.digest(HERE/name) for name in FILES},
        results_reproduced=True,scientific_results_rewritten=False,
        independent_analytic_code_and_final_draft_review=True,
        visual_checks_performed=False,legacy_science_experiments_rerun=False,
        text_checks=checked,local_links_checked=links,broken_links=0,
        scope=saved['scope'],all_reported_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args();result=verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('scientific_base_through_round',
        'unchanged_numbered_scientific_tests','independent_diagnostic_checks',
        'total_protected_including_this_review','all_reported_checks_passed')}))
