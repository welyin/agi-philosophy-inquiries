"""Reproduce an unnumbered conjecture audit; preserve round 520 and its evidence."""
import argparse
import ast
import json
from pathlib import Path
import re
import sys

sys.setrecursionlimit(max(10000, sys.getrecursionlimit()))
import verify_round520_integration as previous
import verify_interaction_rounds as core
import joint_event_conjecture_probe as science

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_event_conjecture_checks.json'
FILES = ['joint_event_conjecture_probe.py', 'joint_event_conjecture_probe_results.json',
         'joint_event_conjecture_review.md', 'joint_event_conjecture_review_draft.txt',
         'joint_event_conjecture_review_draft_v2.txt']


def verify():
    base = previous.verify()
    assert base['stage_saved_tests'] == 2552
    assert base['science_hashes_verified_231_520'] == 871
    assert base['total_protected_evidence_hashes'] == 1008
    if TARGET.exists():
        for name, sha in core.read(TARGET)['new_unnumbered_evidence_hashes'].items():
            assert core.digest(HERE/name) == sha, name
    assert (HERE/FILES[2]).read_bytes() == (HERE/FILES[4]).read_bytes()
    for name in (FILES[0], Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding='utf8'))
    saved = core.read(HERE/FILES[1])
    assert saved == json.loads(json.dumps(science.run()))
    assert (saved['diagnostic_tests'], saved['failures'], saved['errors']) == (3, 0, 0)
    assert saved['numbered_round_created'] is False
    assert saved['numbered_scientific_test_increment'] == 0
    for name, sha in saved['dependency_sha256'].items():
        assert core.digest(HERE.parents[1]/name) == sha, name
    # This is an unnumbered review, not a numbered research note: its three
    # display equations deliberately have no equation tags.
    content = (HERE/FILES[2]).read_text(encoding='utf8')
    assert not any(ord(c) < 32 and c not in '\r\n\t' for c in content)
    delimiters = re.findall(r'^\$\$\s*$', content, re.MULTILINE)
    assert len(delimiters) == 6
    assert '\\tag{' not in content
    text = dict(display_formulas=3, formulas_intentionally_unnumbered=True,
                delimiters_paired=True, sha256=core.digest(HERE/FILES[2]))
    links = 0
    for link in core.link_parser()((HERE/FILES[2]).read_text(encoding='utf8')):
        assert (HERE/link).resolve().exists(), link
        links += 1
    assert not saved['scope']['full_FUCP_reconstruction_refuted']
    assert not saved['scope']['spacetime_dimension_or_GR_or_standard_model_derived']
    return dict(date='2026-09-30', scientific_base_through_round=520,
        unchanged_numbered_scientific_tests=2552, unchanged_scientific_file_hashes=871,
        unchanged_previously_protected_evidence_hashes=1008,
        numbered_round_created=False, numbered_test_increment=0,
        independent_diagnostic_checks=dict(run=3, failures=0, errors=0),
        total_protected_including_this_review=1013,
        new_unnumbered_evidence_hashes={name:core.digest(HERE/name) for name in FILES},
        results_reproduced=True, scientific_results_rewritten=False,
        user_conjecture_unchanged=True,
        text_checks=text, local_links_checked=links, broken_links=0,
        visual_checks_performed=False, legacy_numbered_experiments_rerun=False,
        scope=saved['scope'], all_reported_checks_passed=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    answer = verify()
    if args.write_checks:
        with TARGET.open('x', encoding='utf8', newline='\n') as f:
            f.write(json.dumps(answer, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({k:answer[k] for k in ('scientific_base_through_round',
        'unchanged_numbered_scientific_tests', 'total_protected_including_this_review',
        'all_reported_checks_passed')}))
