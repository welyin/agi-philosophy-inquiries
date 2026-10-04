"""Reproduce round 530 and preserve the complete frozen evidence chain."""
import argparse
import ast
import json
from pathlib import Path
import finite_patch_reference_source as model
import verify_interaction_rounds as core
import verify_round529 as previous

HERE = Path(__file__).resolve().parent
TARGET = HERE/'research_round_530_checks.json'


def verify():
    if TARGET.exists():
        base = previous.verify()
    else:
        import verify_round388_integration as anchor
        original = anchor.verify
        def pending_once(pending_report=False):
            old_target = anchor.TARGET
            anchor.TARGET = TARGET
            try:
                return original(True)
            finally:
                anchor.TARGET = old_target
        anchor.verify = pending_once
        try:
            base = previous.verify()
        finally:
            anchor.verify = original
    assert base['round'] == 529
    assert base['cumulative_numbered_tests'] == 2616
    assert base['cumulative_numbered_scientific_files'] == 899
    assert base['cumulative_unique_protected_evidence_files'] == 1106
    result = core.read(model.TARGET)
    assert result == json.loads(json.dumps(model.run()))
    assert (result['tests_run'], result['failures'], result['errors']) == (6, 0, 0)
    for name, sha in result['dependency_hashes'].items():
        assert core.digest(HERE/name) == sha, name
    note = HERE/'research_note_530.md'
    assert note.read_bytes() == (HERE/'round530_drafts/research_note_530.txt').read_bytes()
    for name in ('scope_proposal.txt', 'buffer_review.txt', 'final_review.txt',
                 'research_note_530_before_final_review.txt'):
        assert (HERE/'round530_drafts'/name).exists()
    text = core.text_checks(note)
    assert text['display_formulas'] == 12
    hashes = {name:core.digest(HERE/name) for name in (
        'finite_patch_reference_source.py', 'finite_patch_reference_source_results.json',
        'research_note_530.md')}
    drafts = {p.relative_to(HERE).as_posix():core.digest(p)
              for p in sorted((HERE/'round530_drafts').glob('*.txt'))}
    if TARGET.exists():
        old = core.read(TARGET)
        assert old['new_file_hashes'] == hashes
        assert old['preserved_draft_hashes'] == drafts
    links = 0
    for link in core.link_parser()(note.read_text('utf8')):
        path = (HERE/link).resolve()
        assert path.exists() or path == TARGET.resolve(), link
        links += 1
    for name in ('finite_patch_reference_source.py', Path(__file__).name):
        ast.parse((HERE/name).read_text('utf8'))
    return dict(date='2026-09-30', round=530, scientific_base_through_round=529,
        fresh_tests=dict(run=6, failures=0, errors=0),
        cumulative_numbered_tests=2622, cumulative_numbered_scientific_files=902,
        unchanged_prior_evidence_files=1106,
        cumulative_unique_protected_evidence_files=1106+len(hashes)+len(drafts),
        new_file_hashes=hashes, preserved_draft_hashes=drafts,
        text_checks=text, local_links_checked=links, broken_links=0,
        saved_results_reproduced=True, previous_results_unchanged=True,
        inherited_theorems_not_counted_as_new=True,
        independent_final_code_and_draft_review_completed=True,
        visual_checks_performed=False, active_goal_unchanged=True,
        scope=result['scope'], all_reported_checks_passed=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify()
    if args.write_checks:
        with TARGET.open('x', encoding='utf8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    elif TARGET.exists():
        assert core.read(TARGET) == result
    print(json.dumps({k:result[k] for k in ('round', 'cumulative_numbered_tests',
        'cumulative_unique_protected_evidence_files', 'all_reported_checks_passed')}))
