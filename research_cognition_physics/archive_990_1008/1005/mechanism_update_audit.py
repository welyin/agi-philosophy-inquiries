"""Reproduce the evidence and coverage audit; this is not a radiation solver."""
from pathlib import Path
import argparse
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]
TARGET = HERE / 'mechanism_update_results.json'


def read(p):
    return json.loads(p.read_text('utf-8-sig'))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run():
    base = read(STAGE / '1003/mechanism_synthesis_results.json')
    previous = read(STAGE / '1004/research_round_1004_checks.json')
    assert previous['all_delivery_checks_passed']
    assert previous['cumulative_numbered_test_groups_from_1003'] == 3786
    prose = (HERE / 'overall_operation_hypothesis_v2_1.md').read_text('utf-8')
    section = prose.split('## 3. 十六类现象', 1)[1].split('## 4.', 1)[0]
    chunks = re.split(r'^### P(\d\d) (.+)\n', section, flags=re.M)
    assert len(chunks) == 49
    keys = ('参与者与过程', '中间机制', '资源、记录与反作用', '输入与范围')
    coverage = []
    for i, row in enumerate(base['coverage']):
        num, title, body = chunks[1+3*i:4+3*i]
        assert num == f'{i+1:02d}' and title == row['phenomena']
        entries = {}
        for key in keys:
            found = re.findall(r'^- \*\*' + re.escape(key) + r'：\*\* (.+)$', body, re.M)
            assert len(found) == 1 and found[0].strip(), (num, key)
            entries[key] = found[0]
        # These are declared evidence classifications, not automatic physics judgments.
        coverage.append(dict(id=row['id'], section='P'+num, phenomena=title,
            conditions=row['conditions'], obligations=entries,
            complete_physical_recovery=False))
    assert {c for row in coverage for c in row['conditions']} == {f'C{i:02d}' for i in range(1, 28)}
    horizon = read(STAGE / '1004/finite_horizon_response_results.json')
    assert horizon['all_scientific_checks_passed']
    assert not horizon['same_background_is_two_self_consistent_backreaction_solutions']
    assert not horizon['numerical_calibration_is_four_dimensional']
    for term in ('386或425', '整体目标未完成', '独立物理输入尚未被证明减少'):
        assert term in prose, term
    evidence = [Path(__file__), HERE/'overall_operation_hypothesis_v2_1.md',
        HERE/'radiation_formation_adoption_v1.md', HERE/'drafts/adoption_decision.md',
        HERE/'drafts/review_notes.md', STAGE/'research_note_1005.md',
        STAGE/'1006/drafts/STATUS.md', STAGE/'1003/mechanism_synthesis_results.json',
        STAGE/'1003/overall_operation_hypothesis_v2.md',
        STAGE/'1004/research_round_1004_checks.json',
        STAGE/'1004/drafts/resource_adoption_checks.json',
        STAGE/'1004/resource_preparation_adoption_v1.md',
        STAGE/'1004/strong_field_adoption_v1.md',
        STAGE/'1004/finite_horizon_response_results.json',
        STAGE/'998/formation_resource_adoption_v1.md',
        STAGE/'1000/material_phase_transport_adoption_v1.md',
        STAGE/'1001/nuclear_material_resource_adoption_v1.md',
        STAGE.parent/'archive_342_369/research_note_355.md',
        STAGE.parent/'archive_531_553/research_note_539.md']
    return dict(round=1005, date='2026-10-08', kind='mechanism_adoption_and_coverage_audit',
        all_audit_checks_passed=True, formal_reports=1005, fresh_physical_test_groups=0,
        cumulative_test_groups=3786, mechanism_version='v2.1',
        conditions=base['conditions'], coverage=coverage,
        adoption_changes=['paired_matter_radiation_sources',
            'emission_net_energy_loss_and_temperature_distinguished',
            'scattering_momentum_and_heat_distinguished',
            'material_state_opacity_and_radiation_feedback',
            'task_resource_and_strong_field_updates_integrated'],
        adopted_analytic_examples=['nonzero_emission_with_zero_net_exchange_in_equilibrium',
            'coherent_elastic_rest_scattering_with_force_but_zero_heat'],
        newly_closed_cognitive_axiom_gaps=[], radiation_transport_solution_computed=False,
        same_cosmic_material_formation_history_certified=False,
        equilibrium_cooling_counterexample_is_new_physical_test=False,
        overall_mechanism_map_is_complete_physical_recovery=False,
        mathematical_axioms_957_changed=False, app_goal_changed=False,
        full_goal_completed=False,
        next_priority='conventional_superconductivity_as_a_macroscopic_quantum_phase',
        deferred_next_mechanism='neutrino_production_propagation_detection',
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in evidence})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write:
        with TARGET.open('x', encoding='utf-8') as dest:
            json.dump(result, dest, ensure_ascii=False, indent=2)
            dest.write('\n')
    else:
        assert result == read(TARGET)
    print(json.dumps({k:v for k,v in result.items() if k not in
        ('coverage','conditions','source_hashes')}, ensure_ascii=False, indent=2))
