"""Integrate round 430 while preserving all earlier scientific evidence."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round429_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round430_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_430_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 429
    assert checked['fresh_tests'] == dict(run=6, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    assert checked['preserved_draft_hashes'] == {}
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'exchange_relation_audit_results.json')
    assert (saved['tests_run'], saved['failures'], saved['errors']) == (6, 0, 0)
    for key in ('exact_programmable_processor_required', 'common_identification_and_couplings_derived',
                'internal_readout_apparatus_completed', 'full_cognitive_implementation_completed',
                'three_dimensional_space_unconditionally_derived', 'full_GR_goal_completed', 'phase_closure_triggered'):
        assert not saved['scope'][key], key
    for key in ('known_exchange_encoding_explicitly_attributed',
                'three_is_minimum_for_noncommuting_qubit_relations',
                'identical_complete_pair_marginals_can_evolve_differently', 'fixed_continuous_hamiltonian_used'):
        assert saved['scope'][key], key
    note = HERE / 'research_note_430.md'
    assert core.text_checks(note) == checked['text_checks']
    assert checked['text_checks']['display_formulas'] == 10
    for name in ('exchange_relation_audit.py', 'verify_exchange_relation_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    links = 0
    for link in core.link_parser()(note.read_text(encoding='utf-8')):
        dest = (HERE / link).resolve()
        assert dest.exists() or (pending and dest == TARGET), link
        links += 1
    assert result.pop('science_hashes_verified_231_429') == 598
    assert result.pop('unchanged_prior_science_hashes_231_428') == 595
    assert result['stage_saved_tests'] == 2002
    assert result['total_protected_evidence_hashes'] == 633
    direction = (HERE.parent / 'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 430 and int(latest[2]) >= 2008
    assert '统一内生演化' in direction and '连续自然演化' in direction
    assert '## 58. 第430轮' in (HERE / 'spatial_premise_closure_audit.md').read_text(encoding='utf-8')
    result.update(date='2026-09-24', rounds=[430],
        execution_mode='one complete continuous exchange-relation round with analytic and matrix verification',
        scientific_base_through_round_by_round={430: 429}, fresh_tests_by_round={430: 6}, fresh_tests=6,
        stage_saved_tests=2008, science_hashes_verified_231_430=601,
        unchanged_prior_science_hashes_231_429=598, total_protected_evidence_hashes=636,
        unchanged_prior_evidence_hashes=633, newly_preserved_draft_files=0,
        new_scientific_display_formulas_checked=10,
        local_links_checked=result['local_links_checked']+links,
        navigation_and_notes_checked=result['navigation_and_notes_checked']+1,
        three_body_exchange_relation_block_verified=True,
        established_exchange_encoding_attributed=True,
        identical_static_pair_marginals_and_different_future_reads_verified=True,
        continuous_natural_evolution_distinguished_from_exact_processor=True,
        original_229_theorem_not_rewritten=True,
        physical_positions_or_dimension_generated=False, full_cognition_to_gr_refuted=False,
        phase_closure_triggered=False, old_scientific_experiments_rerun=False,
        independent_final_agent_review_completed=False, all_reported_checks_passed=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
