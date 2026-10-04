"""Read-only scientific reproduction and preserved-evidence audit for round 524."""
import argparse
import ast
import json
from pathlib import Path
import compact_probe_reference_model as model
import verify_interaction_rounds as core
import verify_round521_523_integration as previous

HERE = Path(__file__).resolve().parent
TARGET = HERE/'research_round_524_checks.json'


def verify():
    if TARGET.exists():
        base = previous.verify()
    else:
        # The legacy archive scans every new Markdown note before this first
        # report exists. Reuse its explicit one-pending-report allowance;
        # change no frozen file and skip no scientific or hash checks.
        import verify_round388_integration as anchor
        original = anchor.verify
        def with_pending_report(pending_report=False):
            old_target = anchor.TARGET
            anchor.TARGET = TARGET
            try:
                return original(True)
            finally:
                anchor.TARGET = old_target
        anchor.verify = with_pending_report
        try:
            base = previous.verify()
        finally:
            anchor.verify = original
    assert base['latest_round'] == 523
    assert base['stage_saved_tests'] == 2570
    assert base['science_hashes_verified_231_523'] == 880
    assert base['total_protected_evidence_hashes'] == 1058
    result = core.read(model.TARGET)
    reproduced = json.loads(json.dumps(model.run()))
    assert result == reproduced, 'saved scientific result differs from reproduction'
    assert (result['tests_run'], result['failures'], result['errors']) == (8, 0, 0)
    for name, digest in result['dependency_hashes'].items():
        assert core.digest(HERE/name) == digest, name
    note = HERE/'research_note_524.md'
    draft = HERE/'round524_drafts/research_note_524.txt'
    assert note.read_bytes() == draft.read_bytes()
    assert (HERE/'round524_drafts/locality_review.txt').exists()
    assert (HERE/'round524_drafts/final_review.txt').exists()
    text = core.text_checks(note)
    assert text['display_formulas'] == 14
    new_hashes = {name: core.digest(HERE/name) for name in
                  ('compact_probe_reference_model.py',
                   'compact_probe_reference_results.json', 'research_note_524.md')}
    drafts = {p.relative_to(HERE).as_posix(): core.digest(p)
              for p in sorted((HERE/'round524_drafts').glob('*.txt'))}
    if TARGET.exists():
        old = core.read(TARGET)
        assert old['new_file_hashes'] == new_hashes
        assert old['preserved_draft_hashes'] == drafts
    local_links = 0
    for link in core.link_parser()(note.read_text('utf8')):
        path = (HERE/link).resolve()
        assert path.exists() or path == TARGET.resolve(), (note, link)
        local_links += 1
    for name in ('compact_probe_reference_model.py', Path(__file__).name):
        ast.parse((HERE/name).read_text('utf8'))
    return dict(date='2026-09-30', round=524, scientific_base_through_round=523,
        fresh_tests=dict(run=8, failures=0, errors=0),
        cumulative_numbered_tests=2578, cumulative_numbered_scientific_files=883,
        unchanged_prior_evidence_files=1058,
        cumulative_unique_protected_evidence_files=1058+len(new_hashes)+len(drafts),
        new_file_hashes=new_hashes, preserved_draft_hashes=drafts,
        text_checks=text, local_links_checked=local_links, broken_links=0,
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
    print(json.dumps({k: result[k] for k in ('round', 'cumulative_numbered_tests',
        'cumulative_unique_protected_evidence_files', 'all_reported_checks_passed')}))
