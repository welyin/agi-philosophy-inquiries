"""Read-only verification of an unnumbered, explicitly conditional interface."""
import argparse
import ast
import json
from pathlib import Path
import sys
sys.setrecursionlimit(max(10000, sys.getrecursionlimit()))
import verify_joint_event_conjecture_review as previous
import verify_interaction_rounds as core
import correlated_remote_witness as science

HERE=Path(__file__).resolve().parent
TARGET=HERE/'remote_witness_interface_checks.json'
FILES=['correlated_remote_witness.py','correlated_remote_witness_results.json',
       'remote_witness_interface_review.md','remote_witness_interface_review_draft.txt']

def verify():
    base=previous.verify()
    assert base['unchanged_numbered_scientific_tests']==2552
    assert base['unchanged_scientific_file_hashes']==871
    assert base['total_protected_including_this_review']==1013
    if TARGET.exists():
        for name,sha in core.read(TARGET)['new_unnumbered_evidence_hashes'].items():
            assert core.digest(HERE/name)==sha,name
    assert (HERE/FILES[2]).read_bytes()==(HERE/FILES[3]).read_bytes()
    saved=core.read(HERE/FILES[1])
    assert saved==json.loads(json.dumps(science.run()))
    assert (saved['diagnostic_tests'],saved['failures'],saved['errors'])==(4,0,0)
    assert not saved['numbered_round_created']
    for name,sha in saved['dependency_sha256'].items():
        assert core.digest(HERE/name)==sha,name
    # Direct generator dependency is already protected by the inherited chain.
    direct='all_scale_monitored_graph_source.py'
    prior=core.read(HERE/'research_round_517_checks.json')
    assert core.digest(HERE/direct)==prior['new_file_hashes'][direct]
    for name in (FILES[0],Path(__file__).name):
        ast.parse((HERE/name).read_text('utf8'))
    checked=core.text_checks(HERE/FILES[2])
    assert checked['display_formulas']==7
    links=0
    for link in core.link_parser()((HERE/FILES[2]).read_text('utf8')):
        assert (HERE/link).resolve().exists(),link
        links+=1
    assert saved['scope']['explicit_added_local_coupling']
    assert not saved['scope']['spatial_dimension_or_GR_derived']
    return dict(date='2026-09-30',scientific_base_through_round=520,
        unchanged_numbered_scientific_tests=2552,unchanged_scientific_file_hashes=871,
        unchanged_previously_protected_evidence_hashes=1013,
        total_protected_including_this_review=1017,
        numbered_round_created=False,numbered_test_increment=0,
        independent_diagnostic_checks=dict(run=4,failures=0,errors=0),
        new_unnumbered_evidence_hashes={f:core.digest(HERE/f) for f in FILES},
        inherited_direct_dependency_verified={direct:core.digest(HERE/direct)},
        results_reproduced=True,independent_final_code_and_draft_review_completed=True,
        text_checks=checked,local_links_checked=links,broken_links=0,
        old_science_rewritten=False,visual_checks_performed=False,
        scope=saved['scope'],all_reported_checks_passed=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--write-checks',action='store_true')
    args=p.parse_args();answer=verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(answer,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:answer[k] for k in ('scientific_base_through_round',
        'unchanged_numbered_scientific_tests','total_protected_including_this_review',
        'all_reported_checks_passed')}))
