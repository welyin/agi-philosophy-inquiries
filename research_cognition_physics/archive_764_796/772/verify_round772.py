"""Reproduce 772 and verify preceding evidence; checks are scoped to the note."""
import argparse
import ast
import json
from pathlib import Path
import verify_interaction_rounds as core
import joint_brst_quartet_matching as model

HERE = Path(__file__).resolve().parent
TARGET = HERE/'research_round_772_checks.json'


def verify():
    base = core.read(HERE/'research_round_771_checks.json')
    assert base['all_reported_checks_passed']
    history = dict(core.read(HERE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
    for n in range(584, 772):
        receipt = core.read(HERE/f'research_round_{n}_checks.json')
        for key in ('new_file_hashes', 'preserved_draft_hashes'):
            for name, digest in receipt[key].items():
                assert name not in history or history[name] == digest, name
                history[name] = digest
    history.update(core.read(HERE/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
    assert len(history) == 3773
    for name, digest in history.items():
        assert core.digest(HERE/name) == digest, name
    result = model.run()
    assert result == core.read(model.TARGET)
    assert (result['tests_run'], result['failures'], result['errors']) == (3, 0, 0)
    for name, digest in result['dependency_hashes'].items():
        assert core.digest(HERE/name) == digest
    main = ('research_note_772.md', 'joint_brst_quartet_matching.py',
            'joint_brst_quartet_matching_results.json')
    extra = ('unified_physics_condition_ledger_772.md',
             'round772_drafts/research_note_772_draft.md',
             'round772_drafts/final_review.txt',
             'round772_drafts/literature_scope_audit.json',
             'round772_drafts/scope_and_dedup_review.md',
             'round773_drafts/STATUS.md',
             'round772_drafts/research_note_772_working.md',
             'round772_drafts/first_response_entry.py',
             'round772_drafts/first_response_entry_results.json',
             'round772_drafts/entry_scope_review.json',
             'round772_drafts/entry_condition_ledger.md',
             'round772_drafts/first_response_entry_checks.json',
             'round772_drafts/entry_postpublication_checks.json')
    entry = core.read(HERE/'round772_drafts/first_response_entry_checks.json')
    for p, digest in entry['working_artifact_hashes'].items():
        assert core.digest(HERE/p) == digest, p
    new = {p: core.digest(HERE/p) for p in main}
    preserved = {p: core.digest(HERE/p) for p in extra}
    assert not (set(new)|set(preserved)) & set(history)
    assert len(history|new|preserved) == 3789
    assert (HERE/main[0]).read_bytes() == (HERE/extra[1]).read_bytes()
    review = (HERE/extra[2]).read_text('utf8')
    for p in (*main, extra[0]):
        assert core.digest(HERE/p) in review
    checks = core.text_checks(HERE/main[0])
    assert checks['display_formulas'] == 16
    links = 0
    for p in (main[0], extra[0], extra[5]):
        for link in core.link_parser()((HERE/p).read_text('utf8')):
            target = ((HERE/p).parent/link).resolve()
            assert target.exists() or target == TARGET.resolve(), (p, link)
            links += 1
    for p in (main[1], 'prepare_round772.py', 'verify_round772.py', 'publish_round772.py'):
        ast.parse((HERE/p).read_text('utf8'))
    return dict(date='2026-10-04', round=772, scientific_base_through_round=771,
                fresh_tests=dict(run=3, failures=0, errors=0),
                cumulative_numbered_tests=3502, cumulative_numbered_scientific_files=1628,
                unchanged_prior_evidence_files=3773,
                cumulative_unique_protected_evidence_files=3789,
                new_file_hashes=new, preserved_draft_hashes=preserved,
                text_checks=checks, local_links_checked=links, broken_links=0,
                saved_results_reproduced=True, previous_results_unchanged=True,
                full_historical_science_rerun=False,
                original_physical_free_state_preserved=True,
                unphysical_extension_changed=True,
                algebraic_free_BRST_quartet_positivity_proven=True,
                source_and_mean_dictionary_matched=True,
                original_single_source_Ward_inherited=True,
                full_Wick_operator_domain_proven=False,
                interacting_QME_proven=False,
                finite_strength_self_consistency_proven=False,
                actual_instrument_or_Q_continuum_equivalence_proven=False,
                primary_code_and_note_review_completed=True,
                independent_final_code_and_draft_review_completed=False,
                inherited_results_not_claimed_as_new_theorems=True,
                visual_checks_performed=False, active_goal_unchanged=True,
                scope=result['scope'], all_reported_checks_passed=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify()
    payload = json.dumps(result, ensure_ascii=False, indent=2)+'\n'
    if args.write_checks:
        with TARGET.open('x', encoding='utf8', newline='\n') as f:
            f.write(payload)
    else:
        assert result == core.read(TARGET)
    print(json.dumps({k: result[k] for k in ('round', 'cumulative_numbered_tests', 'all_reported_checks_passed')}))
