"""Verify the unnumbered 518 source correction without changing old evidence."""
import argparse
import ast
import json
from pathlib import Path
import sys
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
import verify_round518_integration as previous
import verify_interaction_rounds as core
import record_determined_predictor_probe as science

HERE=Path(__file__).resolve().parent
TARGET=HERE/'record_determined_predictor_checks.json'
FILES=['record_determined_predictor_probe.py','record_determined_predictor_probe_results.json',
       'record_determined_predictor_review.md','record_determined_predictor_review_draft.txt']


def verify():
    base=previous.verify()
    assert base['stage_saved_tests']==2542
    assert base['science_hashes_verified_231_518']==865
    assert base['total_protected_evidence_hashes']==992
    if TARGET.exists():
        for name,sha in core.read(TARGET)['new_unnumbered_evidence_hashes'].items():
            assert core.digest(HERE/name)==sha,name
    assert (HERE/FILES[2]).read_bytes()==(HERE/FILES[3]).read_bytes()
    for name in (FILES[0],Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding='utf8'))
    saved=core.read(HERE/FILES[1])
    assert saved==json.loads(json.dumps(science.run()))
    assert (saved['diagnostic_tests'],saved['failures'],saved['errors'])==(4,0,0)
    assert saved['scientific_baseline_round']==518
    assert saved['numbered_round_created'] is False
    assert saved['numbered_scientific_test_increment']==0
    for name,sha in saved['dependency_sha256'].items():
        assert core.digest(HERE/name)==sha,name
    checked=core.text_checks(HERE/FILES[2])
    assert checked['display_formulas']==9
    links=0
    for link in core.link_parser()((HERE/FILES[2]).read_text(encoding='utf8')):
        assert (HERE/link).resolve().exists(),link
        links+=1
    assert saved['scope']['quantum_snapshot_not_required_for_this_forecast']
    assert not saved['scope']['certified_nontrivial_spatial_resolution']
    assert not saved['scope']['physical_three_dimensional_space_generated']
    return dict(date='2026-09-28',scientific_base_through_round=518,
        unchanged_numbered_scientific_tests=2542,unchanged_scientific_file_hashes=865,
        unchanged_previously_protected_evidence_hashes=992,
        numbered_round_created=False,numbered_test_increment=0,
        independent_diagnostic_checks=dict(run=4,failures=0,errors=0),
        total_protected_including_this_review=996,
        new_unnumbered_evidence_hashes={name:core.digest(HERE/name) for name in FILES},
        results_reproduced=True,scientific_results_rewritten=False,
        independent_analytic_code_and_final_draft_review=True,
        visual_checks_performed=False,legacy_science_experiments_rerun=False,
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
