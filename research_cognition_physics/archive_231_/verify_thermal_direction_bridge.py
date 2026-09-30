"""Reproduce the conditional common field/record/qubit-direction model."""
import argparse
import ast
import json
from pathlib import Path
import sys
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
import verify_reference_acquisition_precision_interfaces as previous
import verify_interaction_rounds as core
import thermal_direction_bridge_model as science

HERE = Path(__file__).resolve().parent
TARGET = HERE/'thermal_direction_bridge_checks.json'
FILES = ['thermal_direction_bridge_model.py','thermal_direction_bridge_results.json',
         'thermal_direction_bridge_review.md','thermal_direction_bridge_review_draft.txt',
         'thermal_direction_bridge_review_before_final_review.txt']


def verify():
    base = previous.verify()
    assert base['total_protected_including_this_review']==1047
    assert base['unchanged_numbered_scientific_tests']==2552
    assert base['unchanged_scientific_file_hashes']==871
    if TARGET.exists():
        for name,digest in core.read(TARGET)['new_unnumbered_evidence_hashes'].items():
            assert core.digest(HERE/name)==digest,name
    saved = core.read(science.TARGET)
    assert saved==json.loads(json.dumps(science.run()))
    assert (saved['diagnostic_groups'],saved['failures'],saved['errors'])==(6,0,0)
    note = HERE/'thermal_direction_bridge_review.md'
    assert note.read_bytes()==(HERE/'thermal_direction_bridge_review_draft.txt').read_bytes()
    checked = core.text_checks(note)
    assert checked['display_formulas']==13
    links = 0
    for link in core.link_parser()(note.read_text('utf8')):
        assert (HERE/link).resolve().exists(),link
        links += 1
    for path in (HERE/'thermal_direction_bridge_model.py',Path(__file__)):
        ast.parse(path.read_text('utf8'))
    assert saved['scope']['conditional_three_dimensional_joint_model']
    assert not saved['scope']['fu_cp_unique_derivation']
    assert not saved['scope']['same_lattice_continuum_limit_proved']
    assert saved['scope']['absolute_precision_separate_requirement']
    return dict(date='2026-09-30',scientific_base_through_round=520,
        unchanged_numbered_scientific_tests=2552,unchanged_scientific_file_hashes=871,
        unchanged_previously_protected_evidence_hashes=1047,
        total_protected_including_this_review=1047+len(FILES),
        numbered_round_created=False,numbered_test_increment=0,
        diagnostic_groups=6,failures=0,errors=0,results_reproduced=True,
        text_checks=checked,local_links_checked=links,broken_links=0,
        new_unnumbered_evidence_hashes={f:core.digest(HERE/f) for f in FILES},
        independent_cross_review_completed=True,scope=saved['scope'],
        old_science_rewritten=False,visual_checks_performed=False,
        stage_complete=False,all_reported_checks_passed=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args = parser.parse_args()
    result = verify()
    if args.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('total_protected_including_this_review',
        'diagnostic_groups','all_reported_checks_passed')}))
