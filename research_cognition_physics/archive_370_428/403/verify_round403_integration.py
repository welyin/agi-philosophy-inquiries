"""Integrate round 403 without rerunning historical science or rewriting evidence."""
import argparse
import ast
import json
from pathlib import Path
import re
import verify_round402_integration as previous
import verify_passive_control_round as science
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
TARGET = HERE/'round403_integration_checks.json'


def verify(pending=False):
    old_target = previous.TARGET
    previous.TARGET = TARGET
    try:
        result = previous.verify(pending)
    finally:
        previous.TARGET = old_target
    checked = core.read(HERE/'research_round_403_checks.json')
    assert checked['all_reported_checks_passed']
    assert checked['scientific_base_through_round'] == 402
    assert checked['batch_scientific_dependencies'] == []
    assert len(checked['new_file_hashes']) == 3
    for name,sha in checked['new_file_hashes'].items():
        assert core.digest(HERE/name) == sha,name
    extra = science.correction_checks()
    assert extra == checked['round399_text_maintenance_hashes']
    saved = core.read(HERE/'passive_control_locality_audit_results.json')
    assert saved['checks'] == {'run':8,'failures':0,'errors':0}
    for flag in ('exact_zero_equivalence_proved','drift_norm_uniformly_bounded',
                 'control_amplitude_uniform_in_small_coupling','constant_control_exactly_solved'):
        assert saved[flag],flag
    for flag in ('uniform_passive_smallness_sufficient_for_unrestricted_control',
                 'propagation_time_uniform_in_small_coupling',
                 'integrated_control_action_uniform_in_small_coupling','outside_control_cost_ignored',
                 'full_autonomous_controller_derived','spatial_dimension_generated','full_cognition_to_gr_refuted'):
        assert not saved[flag],flag
    formulas = core.text_checks(HERE/'research_note_403.md')['display_formulas']
    assert formulas == 13
    links = 0
    for path in (HERE/'research_note_403.md',HERE/'research_note_399_text_v2.md',
                 HERE/'spatial_premise_closure_audit.md'):
        for link in core.link_parser()(path.read_text(encoding='utf-8')):
            destination = (path.parent/link).resolve()
            assert destination.exists() or (pending and destination == TARGET),(path,link)
            links += 1
    for name in ('passive_control_locality_audit.py','verify_passive_control_round.py',Path(__file__).name):
        ast.parse((HERE/name).read_text(encoding='utf-8'))
    assert result.pop('science_hashes_verified_231_402') == 517
    assert result.pop('unchanged_prior_science_hashes_231_401') == 514
    assert result['stage_saved_tests'] == 1799
    assert result['total_protected_evidence_hashes'] == 538
    assert result['protected_text_maintenance_files'] == 7
    direction = (HERE.parent/'research_direction.md').read_text(encoding='utf-8')
    latest = re.search(r'最新科学轮次与检查数为(\d+)／(\d+)',direction)
    assert latest and int(latest[1]) >= 403 and int(latest[2]) >= 1807
    result.update(rounds=[403],execution_mode='exact locality and resource-dependent intervention bounds',
                  scientific_base_through_round_by_round={403:402},
                  fresh_tests_by_round={403:8},fresh_tests=8,stage_saved_tests=1807,
                  science_hashes_verified_231_403=520,unchanged_prior_science_hashes_231_402=517,
                  total_protected_evidence_hashes=543,unchanged_prior_evidence_hashes=538,
                  protected_text_maintenance_files=9,round399_text_maintenance_hashes=extra,
                  round399_original_science_preserved=True,
                  new_scientific_display_formulas_checked=formulas,
                  local_links_checked=result['local_links_checked']+links,
                  navigation_and_notes_checked=result['navigation_and_notes_checked']+3,
                  exact_zero_intervention_equivalence_proved=True,
                  passive_error_alone_insufficient_even_with_bounded_control_amplitude=True,
                  uniform_time_or_integrated_control_cost_not_claimed=True,
                  full_autonomous_controller_derived=False,spatial_dimension_generated=False,
                  full_cognition_to_gr_refuted=False,independent_final_agent_review_completed=False,
                  all_reported_checks_passed=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checks',action='store_true')
    args = parser.parse_args()
    result = verify(args.write_checks)
    if args.write_checks:
        with TARGET.open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
