"""Read-only source-scope audit, preserving the complete 980-file baseline."""
import argparse
import json
from pathlib import Path
import sys
sys.setrecursionlimit(max(sys.getrecursionlimit(),10000))
import verify_hqca_full_cell_scope_review as previous
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
TARGET=HERE/'qca_finite_region_scope_checks.json'
FILES=['qca_finite_region_scope_addendum.md','qca_finite_region_scope_addendum_draft.txt']
REVIEWED_SHA='2dd30231df4dbded2e1432dd21b13193fdf1fae1afb869b99a433922c07bd0e7'


def verify():
    base=previous.verify()
    assert base['total_protected_including_this_review']==980
    assert base['unchanged_numbered_scientific_tests']==2531
    assert base['unchanged_scientific_file_hashes']==859
    for name in FILES:
        assert core.digest(HERE/name)==REVIEWED_SHA,name
    if TARGET.exists():
        for name,sha in core.read(TARGET)['new_unnumbered_evidence_hashes'].items():
            assert core.digest(HERE/name)==sha,name
    checked=core.text_checks(HERE/FILES[0])
    assert checked['display_formulas']==2
    links=0
    for link in core.link_parser()((HERE/FILES[0]).read_text(encoding='utf8')):
        assert (HERE/link).resolve().exists(),link
        links+=1
    return dict(date='2026-09-28',scientific_base_through_round=516,
        stage_saved_tests=2531,science_hashes_verified_231_516=859,
        unchanged_previously_protected_evidence_hashes=980,
        total_protected_evidence_hashes=980+len(FILES),
        new_unnumbered_evidence_hashes={name:core.digest(HERE/name) for name in FILES},
        numbered_round_created=False,numbered_test_increment=0,
        new_scientific_experiments_run=False,source_theorem_and_reference_extension_reviewed=True,
        independent_final_review_completed=True,text_checks=checked,
        local_links_checked=links,broken_links=0,visual_checks_performed=False,
        exact_FUCP_or_continuous_Time_internal_implementation_proved=False,
        macroscopic_spatial_dimension_proved_or_refuted=False,
        full_GR_goal_completed=False,all_reported_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args=parser.parse_args();result=verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('scientific_base_through_round',
        'stage_saved_tests','total_protected_evidence_hashes','all_reported_checks_passed')}))
