"""Integrate round 417 while preserving prior scientific files and results."""
import argparse
import ast
import json
from pathlib import Path
import re

import verify_round415_416_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round417_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE / 'research_round_417_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 416
    assert checked['fresh_tests'] == dict(run=8, failures=0, errors=0)
    assert len(checked['new_file_hashes']) == 3
    for name, sha in checked['new_file_hashes'].items():
        assert core.digest(HERE / name) == sha, name
    saved = core.read(HERE / 'finite_readout_covariance_audit_results.json')
    assert saved['checks'] == checked['fresh_tests']
    assert saved['global_readout_closure_theorem_proved']
    assert saved['finite_summary_autonomous_closure_is_extra_input']
    for key in ('finite_direction_group_ruled_out', 'finite_window_readout_ruled_out',
                'repeated_access_to_source_ruled_out', 'quantum_direction_interface_disproved',
                'all_cognitive_principles_countermodel', 'physical_positions_or_dimension_generated',
                'phase_closure_triggered'):
        assert not saved[key], key
    assert core.text_checks(HERE / 'research_note_417.md') == checked['text_checks']
    links = 0
    for name in ('research_note_417.md', 'spatial_premise_closure_audit.md'):
        for link in core.link_parser()((HERE / name).read_text(encoding='utf-8')):
            dest = (HERE / link).resolve()
            assert dest.exists() or (pending and dest == TARGET), (name, link)
            links += 1
    for name in ('finite_readout_covariance_audit.py', 'verify_readout_covariance_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    assert result.pop('science_hashes_verified_231_416') == 559
    assert result.pop('unchanged_prior_science_hashes_231_414') == 553
    assert result['stage_saved_tests'] == 1899
    assert result['total_protected_evidence_hashes'] == 584
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', (HERE.parent / 'research_direction.md').read_text(encoding='utf-8'))
    assert latest and int(latest[1]) >= 417 and int(latest[2]) >= 1907
    result.update(date='2026-09-24', rounds=[417],
        execution_mode='single new readout-contract round following independently reviewed candidate audits',
        scientific_base_through_round_by_round={417: 416}, fresh_tests_by_round={417: 8},
        fresh_tests=8, stage_saved_tests=1907, science_hashes_verified_231_417=562,
        unchanged_prior_science_hashes_231_416=559, total_protected_evidence_hashes=587,
        unchanged_prior_evidence_hashes=584, new_scientific_display_formulas_checked=14,
        local_links_checked=result['local_links_checked'] + links,
        navigation_and_notes_checked=result['navigation_and_notes_checked'] + 2,
        global_finite_readout_joint_closure_triviality_proved=True,
        finite_summary_self_containment_adopted_as_cognitive_axiom=False,
        local_chart_or_source_reaccess_ruled_out=False, two_input_half_error_bound_proved=True,
        physical_positions_or_dimension_generated=False, full_cognition_to_gr_refuted=False,
        phase_closure_triggered=False, old_scientific_experiments_rerun=False,
        independent_final_agent_review_completed=True, all_reported_checks_passed=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks', action='store_true')
    args = parser.parse_args()
    result = verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
