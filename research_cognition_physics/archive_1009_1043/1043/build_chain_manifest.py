"""Resolve historical dependency deltas without editing them.

This is a provenance index, not a theorem prover. --write may update this
round's draft index until its final receipt exists; after that it is read-only.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
BASE = STAGE.parent
OUT = HERE / 'chain_manifest.json'
RECEIPT = HERE / 'research_round_1043_checks.json'

GROUPS = {
    'same_leading_action_chain': [
        'gauge_self_vertices', 'Higgs_gauge_vertices', 'portal_potential',
        'scalar_source_to_parent', 'gravity_ideal_to_parent', 'charged_source',
        'gravity_parent_bridge', 'lift_equivalence', 'parent_coupling_freedom'],
    'alternative_charge_mass_global_branches': [
        'hypercharge_narrow', 'mass_menu_broad', 'Higgs_family',
        'spectral_current', 'center_matter_kernel', 'electric_completion',
        'twisted_spin_mass'],
    'fixed_material_response': ['material_response'],
    'conditional_tools_not_automatically_same_parent': [
        'statistics_theorem', 'thermal_reference', 'chiral_species',
        'forward_dispersion', 'finite_internal_geometry', 'speed_attraction',
        'finite_scattering', 'classical_channel',
        'finite_geometry_effective_mass_transport', 'local_record_energy_tradeoff',
        'moving_direction_joint_statistics', 'fusion_pairing_quadratic_mass',
        'nonlinear_electromagnetic_common_cone',
        'whole_family_relay_time_selection',
        'complementary_recovery_central_entropy'],
    'unestablished_full_parent_transport': ['regional_CAR_witness', 'quotient_process'],
}


def read(path):
    return json.loads(path.read_text('utf8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path):
    return path.resolve().relative_to(BASE).as_posix()


def resolve_rows():
    paths = [STAGE / '1035/bridge_ledger_v0_2.json',
             STAGE / '1040/dependency_delta_1036_1040.json',
             STAGE / '1041/dependency_delta_1041.json',
             STAGE / '1042/dependency_delta_1042.json']
    docs = [read(p) for p in paths]
    rows = docs[0]['rows'].copy()
    for i, doc in enumerate(docs[1:], 1):
        key = 'base_ledger' if i == 1 else 'base_delta'
        digest_key = 'base_sha256' if i == 1 else 'base_delta_sha256'
        assert BASE / doc[key] == paths[i - 1]
        assert sha(paths[i - 1]) == doc[digest_key]
        assert len(rows) == doc.get('base_row_count', doc.get('base_resolved_rows'))
        rows.extend(doc['rows_append'])
        assert len(rows) == doc['resolved_row_count']
    assert len(rows) == len({r['id'] for r in rows}) == 34
    return rows, paths


def build():
    rows, ledger_paths = resolve_rows()
    memberships = {key: group for group, keys in GROUPS.items() for key in keys}
    assert sum(map(len, GROUPS.values())) == len(memberships) == len(rows)
    assert set(memberships) == {r['id'] for r in rows}
    input_ledger = STAGE / '1009/input_dependency_ledger_v0_1.md'
    original = BASE / 'archive_990_1008/1008/overall_operation_hypothesis_v2_2.md'
    inputs = re.findall(r'^\|(I\d\d) ([^|]+)\|', input_ledger.read_text('utf8'), re.M)
    phenomena = re.findall(r'^### (P\d\d) (.+)$', original.read_text('utf8'), re.M)
    note = STAGE / 'research_note_1043.md'
    body = note.read_text('utf8')
    assert [i for i, _ in inputs] == [f'I{i:02d}' for i in range(1, 15)]
    assert [i for i, _ in phenomena] == [f'P{i:02d}' for i in range(1, 17)]
    assert set(re.findall(r'^\|(I\d\d) ', body, re.M)) == {i for i, _ in inputs}
    assert set(re.findall(r'^\|(P\d\d) ', body, re.M)) == {i for i, _ in phenomena}
    source_paths = set(ledger_paths) | {input_ledger, original, note,
        STAGE/'1009/goal_start.json', HERE/'completion_audit.md'}
    source_paths.update(STAGE / f'research_note_{n}.md' for n in range(1009, 1043))
    for row in rows:
        source_paths.update(BASE / p for p in row['source_paths'])
    source_paths.update(STAGE/'1042'/name for name in (
        'joint_premise_history_audit.md', 'common_parent_coverage_audit.md',
        'generative_chain_scope_audit.md', 'representation_quotient_clarification.md'))
    historical_rounds = {
        'archive_217_222': [222], 'archive_223_230': [230],
        'archive_370_428': [382, 383, 384, 386, 425],
        'archive_467_530': [522, 523], 'archive_923_934': [929],
        'archive_935_955': [947, 951, 955],
        'archive_956_989': [956, 957, 958, 977, 981],
    }
    for folder, rounds in historical_rounds.items():
        source_paths.update(BASE / folder / f'research_note_{n}.md' for n in rounds)
    source_paths.add(BASE/'archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md')
    return {
        'schema': 'conditional_generation_provenance_v1', 'round': 1043,
        'claim_level': 'conditional_chains_and_stratified_remaining_freedom',
        'new_scientific_calibration_groups': 0, 'cumulative': 3819,
        'new_adopted_cognitive_axioms': 0,
        'full_cognitive_derivation_or_common_E_nat_claimed': False,
        'all_remaining_inputs_proved_independent': False,
        'groups_are_not_one_joint_model': True,
        'dependency_count': len(rows),
        'dependencies': [{
            'id': r['id'], 'group': memberships[r['id']],
            'original_status': r['status'], 'object': r['object'],
            'requires': r['requires'], 'not_proved': r['not_proved'],
            'source_rounds': r['source_rounds'], 'source_paths': r['source_paths']}
            for r in rows],
        'coverage_not_semantic_proof': True,
        'original_input_categories': dict(inputs),
        'original_phenomena': dict(phenomena),
        'mutually_incompatible_or_unmapped_uses': [
            'P981_fixed_nonzero_C5_and_exact_SpinZ4_without_changing_mass_menu',
            'pure_U1_dim4_full_mass_theorem_as_SM_dim5_species_classification',
            'same_g3_classical_family_as_full_A1_A7_quantum_countermodels',
            '929_both_updates_as_already_certified_947_field_extensions',
            '956_record_difference_as_1042_complementary_recovery',
            'same_marginals_as_1038_full_moving_direction_joint_data',
            '1042_first_second_label_swap_as_quotient_invariant_freedom'],
        'source_sha256': {rel(p): sha(p) for p in sorted(source_paths)},
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = build()
    if args.write:
        assert not RECEIPT.exists(), 'Final receipt already freezes this index'
        OUT.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', 'utf8')
    else:
        assert read(OUT) == result, 'Index or cited source changed'
    print(json.dumps({'passed': True, 'dependency_rows': len(result['dependencies']),
                      'sources': len(result['source_sha256']),
                      'semantic_proof_by_manifest': False}, ensure_ascii=False))
