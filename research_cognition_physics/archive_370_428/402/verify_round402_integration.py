"""Integrate round 402 while keeping the round-401 snapshot unchanged.

Round 401's historical verifier asserts a literal then-current navigation title.
Use round 400's cumulative audit and explicitly audit all round-401/402 science,
notes and metadata below; do not mutate that old verifier or historical reports.
"""
import argparse
import ast
import json
from pathlib import Path
import re
import verify_round400_integration as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'round402_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    snapshot = core.read(HERE / 'round401_integration_checks.json')
    assert snapshot['all_reported_checks_passed']
    assert snapshot['stage_saved_tests'] == 1789
    assert snapshot['total_protected_evidence_hashes'] == 535
    assert snapshot['science_hashes_verified_231_401'] == 514
    stems = {401: 'persistent_prefix_history_audit', 402: 'recoverable_region_intersection_audit'}
    tests = {401: 8, 402: 10}
    formulas = {}
    for number, stem in stems.items():
        checked = core.read(HERE / f'research_round_{number}_checks.json')
        assert checked['all_reported_checks_passed']
        assert checked['scientific_base_through_round'] == number-1
        assert checked['batch_scientific_dependencies'] == []
        assert len(checked['new_file_hashes']) == 3
        for name, sha in checked['new_file_hashes'].items():
            assert core.digest(HERE / name) == sha, name
        saved = core.read(HERE / (stem+'_results.json'))
        assert saved['checks'] == {'run': tests[number], 'failures': 0, 'errors': 0}
        assert not saved['spatial_dimension_generated']
        formulas[number] = core.text_checks(HERE / f'research_note_{number}.md')['display_formulas']
    assert formulas == {401: 11, 402: 14}
    old = core.read(HERE / 'persistent_prefix_history_audit_results.json')
    assert old['same_fixed_hamiltonian_for_all_finite_prefixes']
    assert old['infinite_internal_workspace_assumed']
    assert old['frozen_output_contract_required']
    for key in ('postselection_used', 'finite_time_exact_completion_claimed',
                'arbitrary_repeated_readout_guaranteed', 'whole_state_convergence_claimed',
                'full_sustained_cognition_model_constructed'):
        assert not old[key], key
    saved = core.read(HERE / 'recoverable_region_intersection_audit_results.json')
    assert saved['full_unknown_reference_preserved_by_authorized_recovery']
    assert saved['alternative_overlapping_coalitions_not_independent_receivers']
    assert saved['fixed_full_domain_dependency_intersection_proved']
    for key in ('unique_least_recovery_region_exists', 'recovery_regions_identified_as_point_neighborhoods',
                'autonomous_geometric_local_implementation_derived',
                'intervening_control_permissions_automatically_enlarged', 'full_cognition_to_gr_refuted'):
        assert not saved[key], key
    links = 0
    for path in (HERE/'research_note_401.md', HERE/'research_note_402.md',
                 HERE/'spatial_premise_closure_audit.md'):
        for link in core.link_parser()(path.read_text(encoding='utf-8')):
            destination = (path.parent / link).resolve()
            assert destination.exists() or (pending and destination == TARGET), (path, link)
            links += 1
    for name in ('persistent_prefix_history_audit.py', 'verify_persistent_prefix_round.py',
                 'verify_round401_integration.py', 'recoverable_region_intersection_audit.py',
                 'verify_recoverable_region_round.py', Path(__file__).name):
        ast.parse((HERE / name).read_text(encoding='utf-8'))
    assert result.pop('science_hashes_verified_231_400') == 511
    assert result.pop('unchanged_prior_science_hashes_231_399') == 508
    assert result['stage_saved_tests'] == 1781
    assert result['total_protected_evidence_hashes'] == 532
    direction = (HERE.parent/'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)', direction)
    assert latest and int(latest[1]) >= 402 and int(latest[2]) >= 1799
    result.update(rounds=[402], execution_mode='recovery versus full-domain dependency intersections',
                  scientific_base_through_round_by_round={402: 401},
                  fresh_tests_by_round={402: 10}, fresh_tests=10,
                  stage_saved_tests=1799, science_hashes_verified_231_402=517,
                  unchanged_prior_science_hashes_231_401=514,
                  total_protected_evidence_hashes=538, unchanged_prior_evidence_hashes=535,
                  new_scientific_display_formulas_checked=formulas[402],
                  local_links_checked=result['local_links_checked']+links,
                  navigation_and_notes_checked=result['navigation_and_notes_checked']+3,
                  round401_science_explicitly_reaudited=True,
                  historical_round401_navigation_assertion_replaced_by_current_lower_bound=True,
                  frozen_round401_integration_snapshot_sha256=core.digest(HERE/'round401_integration_checks.json'),
                  fixed_full_domain_dependency_intersection_proved=True,
                  autonomous_geometric_local_implementation_derived=False,
                  intervening_control_permissions_automatically_enlarged=False,
                  full_cognition_to_gr_refuted=False, spatial_dimension_generated=False,
                  independent_final_agent_review_completed=False,
                  all_reported_checks_passed=True)
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
