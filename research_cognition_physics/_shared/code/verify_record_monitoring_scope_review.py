"""Read-only unnumbered follow-up audit; numbered science remains 514 / 2519."""
import argparse
import ast
import importlib.util
import json
from pathlib import Path
import sys
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
import verify_round514_integration as previous
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
TARGET=HERE/'record_monitoring_scope_checks.json'
FILES=['record_monitoring_scope_probe.py','record_monitoring_scope_probe_results.json',
       'record_monitoring_scope_review.md','record_monitoring_scope_review_draft.txt',
       'record_monitoring_scope_probe_before_record_budget_review.txt']


def verify():
    base=previous.verify()
    assert base['stage_saved_tests']==2519
    assert base['science_hashes_verified_231_514']==853
    assert base['total_protected_evidence_hashes']==960
    if TARGET.exists():
        for name,sha in core.read(TARGET)['new_unnumbered_evidence_hashes'].items():
            assert core.digest(HERE/name)==sha,name
    assert (HERE/FILES[2]).read_bytes()==(HERE/FILES[3]).read_bytes()
    assert core.digest(HERE/FILES[0])=='ae77dd031b594477e6965c0403535cb90f22c825f59eacda658decaf86fab774'
    assert core.digest(HERE/FILES[2])=='c78f503ca194700e482e979d89868dc0f13f00e0831aafe724008fdd3e8b4cce'
    for name in [FILES[0],Path(__file__).name]:ast.parse((HERE/name).read_text(encoding='utf8'))
    spec=importlib.util.spec_from_file_location('record_monitoring_scope_checked',HERE/FILES[0])
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module
    spec.loader.exec_module(module)
    saved=core.read(HERE/FILES[1])
    assert saved==json.loads(json.dumps(module.run()))
    assert (saved['diagnostic_tests'],saved['failures'],saved['errors'])==(4,0,0)
    assert saved['scientific_baseline_round']==514
    assert saved['numbered_round_created'] is False
    assert saved['numbered_scientific_test_increment']==0
    for name,sha in saved['dependency_sha256'].items():assert core.digest(HERE/name)==sha,name
    text_checked=core.text_checks(HERE/FILES[2]);assert text_checked['display_formulas']==9
    links=0
    for link in core.link_parser()((HERE/FILES[2]).read_text(encoding='utf8')):
        assert (HERE/link).resolve().exists(),link
        links+=1
    assert saved['scope']['fine_source_to_symmetric_vector_code_boundary_proved'] is True
    assert saved['scope']['shared_spatial_endpoint_generated'] is False
    assert saved['scope']['dimension_three_generated'] is False
    assert saved['scope']['full_GR_goal_completed'] is False
    return dict(date='2026-09-28',scientific_base_through_round=514,
        unchanged_numbered_scientific_tests=2519,unchanged_scientific_file_hashes=853,
        unchanged_previously_protected_evidence_hashes=960,
        numbered_round_created=False,numbered_test_increment=0,
        independent_diagnostic_checks=dict(run=4,failures=0,errors=0),
        new_unnumbered_evidence_files=len(FILES),
        total_protected_including_this_review=960+len(FILES),
        new_unnumbered_evidence_hashes={name:core.digest(HERE/name) for name in FILES},
        results_reproduced=True,scientific_results_rewritten=False,
        independent_analytic_code_and_final_draft_review=True,
        visual_checks_performed=False,legacy_science_experiments_rerun=False,
        text_checks=text_checked,local_links_checked=links,broken_links=0,
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
        'total_protected_including_this_review','all_reported_checks_passed')},ensure_ascii=False))
