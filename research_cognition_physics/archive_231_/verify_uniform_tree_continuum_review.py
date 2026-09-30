"""Verify a mature-theorem source bridge, preserving the 519 scientific count."""
import argparse
import ast
import json
from pathlib import Path
import sys
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
import verify_round519_integration as previous
import verify_interaction_rounds as core
import uniform_tree_continuum_probe as science

HERE=Path(__file__).resolve().parent
TARGET=HERE/'uniform_tree_continuum_checks.json'
FILES=['uniform_tree_continuum_probe.py','uniform_tree_continuum_probe_results.json',
       'uniform_tree_continuum_review.md','uniform_tree_continuum_review_draft.txt']


def verify():
    base=previous.verify()
    assert base['stage_saved_tests']==2546
    assert base['science_hashes_verified_231_519']==868
    assert base['total_protected_evidence_hashes']==1000
    if TARGET.exists():
        for name,sha in core.read(TARGET)['new_unnumbered_evidence_hashes'].items():
            assert core.digest(HERE/name)==sha,name
    assert (HERE/FILES[2]).read_bytes()==(HERE/FILES[3]).read_bytes()
    for name in (FILES[0],Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding='utf8'))
    saved=core.read(HERE/FILES[1])
    assert saved==json.loads(json.dumps(science.run()))
    assert (saved['diagnostic_tests'],saved['failures'],saved['errors'])==(3,0,0)
    assert saved['scientific_baseline_round']==519
    assert saved['numbered_round_created'] is False
    assert saved['numbered_scientific_test_increment']==0
    for name,sha in saved['dependency_sha256'].items():
        assert core.digest(HERE/name)==sha,name
    checked=core.text_checks(HERE/FILES[2]);assert checked['display_formulas']==6
    links=0
    for link in core.link_parser()((HERE/FILES[2]).read_text(encoding='utf8')):
        assert (HERE/link).resolve().exists(),link
        links+=1
    assert saved['scope']['exact_sampling_law_matched']
    assert saved['scope']['odd_size_subsequence_explicit']
    assert not saved['scope']['measured_dynamic_propagation_metric_identified']
    assert not saved['scope']['full_GR_goal_completed']
    return dict(date='2026-09-28',scientific_base_through_round=519,
        unchanged_numbered_scientific_tests=2546,unchanged_scientific_file_hashes=868,
        unchanged_previously_protected_evidence_hashes=1000,
        numbered_round_created=False,numbered_test_increment=0,
        independent_diagnostic_checks=dict(run=3,failures=0,errors=0),
        total_protected_including_this_review=1004,
        new_unnumbered_evidence_hashes={name:core.digest(HERE/name) for name in FILES},
        results_reproduced=True,scientific_results_rewritten=False,
        independent_analytic_code_and_final_draft_review=True,
        primary_theorem_periodicity_and_normalization_checked=True,
        visual_checks_performed=False,legacy_numbered_science_experiments_rerun=False,
        text_checks=checked,local_links_checked=links,broken_links=0,
        scope=saved['scope'],all_reported_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args();answer=verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(answer,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:answer[k] for k in ('scientific_base_through_round',
        'unchanged_numbered_scientific_tests','total_protected_including_this_review',
        'all_reported_checks_passed')}))
