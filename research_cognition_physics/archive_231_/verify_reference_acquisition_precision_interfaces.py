"""Reproduce positive-rule acquisition and thermal coordinate precision reports."""
import argparse
import ast
import json
from pathlib import Path
import sys
sys.setrecursionlimit(max(10000, sys.getrecursionlimit()))
import verify_spatial_reference_interfaces as previous
import verify_interaction_rounds as core
import positive_reference_acquisition as acquisition
import thermal_reference_dimension_model as thermal

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'reference_acquisition_precision_checks.json'
GROUPS = {
    'positive_reference_acquisition': ('positive_reference_acquisition.py', acquisition, 'diagnostic_groups', 9),
    'thermal_reference_dimension': ('thermal_reference_dimension_model.py', thermal, 'diagnostic_tests', 12),
}
FILES = []
for stem, (code, _, _, _) in GROUPS.items():
    FILES.extend([code, stem+'_results.json', stem+'_review.md', stem+'_review_draft.txt'])
FILES.append('thermal_reference_dimension_review_before_final_review.txt')


def verify():
    base = previous.verify()
    assert base['total_protected_including_this_review'] == 1038
    assert base['unchanged_numbered_scientific_tests'] == 2552
    assert base['unchanged_scientific_file_hashes'] == 871
    if TARGET.exists():
        for name, digest in core.read(TARGET)['new_unnumbered_evidence_hashes'].items():
            assert core.digest(HERE/name) == digest, name
    reports, links = {}, 0
    for stem, (code, science, count_key, formulas) in GROUPS.items():
        saved = core.read(HERE/(stem+'_results.json'))
        assert saved == json.loads(json.dumps(science.run()))
        assert (saved[count_key], saved['failures'], saved['errors']) == (6, 0, 0)
        ast.parse((HERE/code).read_text('utf8'))
        note = HERE/(stem+'_review.md')
        assert note.read_bytes() == (HERE/(stem+'_review_draft.txt')).read_bytes()
        checked = core.text_checks(note)
        assert checked['display_formulas'] == formulas
        for link in core.link_parser()(note.read_text('utf8')):
            assert (HERE/link).resolve().exists(), link
            links += 1
        reports[stem] = dict(diagnostic_groups=6, failures=0, errors=0,
            results_reproduced=True, text_checks=checked, scope=saved['scope'])
    ast.parse(Path(__file__).read_text('utf8'))
    return dict(date='2026-09-30', scientific_base_through_round=520,
        unchanged_numbered_scientific_tests=2552, unchanged_scientific_file_hashes=871,
        unchanged_previously_protected_evidence_hashes=1038,
        total_protected_including_this_review=1038+len(FILES),
        numbered_round_created=False, numbered_test_increment=0, interfaces=reports,
        new_unnumbered_evidence_hashes={f:core.digest(HERE/f) for f in FILES},
        local_links_checked=links, broken_links=0, independent_cross_review_completed=True,
        old_science_rewritten=False, visual_checks_performed=False,
        independent_probe_acquisition_to_thermal_source_proved=False,
        spatial_dimension_selected=False, stage_complete=False,
        all_reported_checks_passed=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify()
    if args.write_checks:
        with TARGET.open('x', encoding='utf8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('scientific_base_through_round',
        'total_protected_including_this_review', 'all_reported_checks_passed')}))
