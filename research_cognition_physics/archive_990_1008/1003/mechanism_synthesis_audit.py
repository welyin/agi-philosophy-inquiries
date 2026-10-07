"""Audit the declared mechanism map and inherited evidence, not physical recovery."""
from pathlib import Path
import argparse
import hashlib
import json
import re

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
TARGET=HERE/'mechanism_synthesis_results.json'


def read(path):return json.loads(path.read_text('utf-8-sig'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    base=read(STAGE/'999/mechanism_adoption_results.json')
    inventory=read(STAGE/'990/mechanism_inventory.json')
    previous=read(STAGE/'1002/research_round_1002_checks.json')
    assert previous['all_delivery_checks_passed'] and not previous['full_goal_completed']
    assert len(base['coverage'])==16
    prose=(HERE/'overall_operation_hypothesis_v2.md').read_text('utf-8')
    section=prose.split('## 3. 十六类现象',1)[1].split('## 4.',1)[0]
    chunks=re.split(r'^### P(\d\d) (.+)\n',section,flags=re.M)
    assert len(chunks)==49
    obligations=('参与者与过程','中间机制','资源、记录与反作用','输入与范围')
    coverage=[]
    # Classifications are declared editorial/evidence judgments, not theorems.
    status_updates={
        'records':'same_thermal_output_join_with_new_resources',
        'ordinary_matter':'formation_preservation_and_phase_mechanisms_adopted',
        'nuclear_weak':'reaction_and_environment_mechanisms_adopted',
        'phases_transport':'common_thermomechanical_example_with_phase_input'}
    gaps={
        'quantum':['interpretation','joint_implementation'],
        'records':['supply_mechanism','joint_implementation'],
        'spacetime':['physical_input_not_derived','joint_implementation'],
        'radiation':['matching'], 'motion_mass':['physical_input_not_derived','matching'],
        'ordinary_matter':['matching','formation_history'],
        'nuclear_weak':['parameters','formation_history'],
        'phases_transport':['parameters','noise_and_lifetime'],
        'thermal_arrow':['supply_mechanism','initial_boundary'],
        'gravity':['physical_input_not_derived','joint_implementation'],
        'horizons':['finite_task_mechanism','beyond_effective_domain'],
        'cosmic_formation':['initial_boundary','formation_history'],
        'asymmetry_neutrino':['actual_source','parameters','formation_history'],
        'dark_matter':['identity','parameters','formation_history'],
        'acceleration':['parameters','microscopic_origin'],
        'cross_scale':['joint_implementation','empirical_recovery']}
    for i,row in enumerate(base['coverage']):
        number,title,body=chunks[1+3*i:4+3*i]
        assert number==f'{i+1:02d}' and title==row['phenomena']
        entries={}
        for key in obligations:
            matches=re.findall(r'^- \*\*'+re.escape(key)+r'：\*\* (.+)$',body,re.M)
            assert len(matches)==1 and matches[0].strip(),(number,key)
            entries[key]=matches[0]
        coverage.append(dict(id=row['id'],section=f'P{number}',phenomena=title,
            conditions=row['conditions'],obligations=entries,
            previous_status=row['status'],status=status_updates.get(row['id'],row['status']),
            gaps=gaps[row['id']],complete_physical_recovery=False))
    assert {c for row in coverage for c in row['conditions']}=={f'C{i:02d}' for i in range(1,28)}
    assert [c['id'] for c in inventory['conditions']]==[f'C{i:02d}' for i in range(1,28)]
    for term in ('386或425','独立物理输入尚未被证明减少','不同实验允许不同准备',
                 '整体目标未完成','正熵产','持续QND','比较热态','场态','相位参考'):
        assert term in prose,term
    thermal=read(STAGE/'980/finite_thermal_records_results.json')
    reuse=read(STAGE/'1002/thermal_record_reuse_results.json')
    for key in ('thermal_populations','energies'):
        assert thermal['parameters'][key]==reuse['parameters'][key],key
    assert reuse['auxiliary_initial_pure_energy_state_is_input']
    assert reuse['engineered_energy_matched_interaction_is_input']
    for key in ('original_charge_interaction_implements_writer',
        'population_effect_is_qnd_under_full_active_H','cyclic_auxiliary_restoration_certified',
        'joint_clock_communication_lifecycle_certified','full_goal_completed'):
        assert not reuse[key],key
    joins=[
        dict(id='thermal_output_to_population_record',rounds=[980,1002],
             state_connection_verified=True,requires_new_resources=True,
             arbitrary_phase_lifecycle_verified=False),
        dict(id='phase_and_thermal_tasks',rounds=[976,980,1002],state_connection_verified=False),
        dict(id='nuclei_to_material_to_function',rounds=[1001,1000,1002],state_connection_verified=False),
        dict(id='dark_preparations',rounds=[991,995,996],state_connection_verified=False),
        dict(id='cosmic_windows',rounds=[992,994,996,997],state_connection_verified=False),
        dict(id='common_action_to_complete_process',rounds=[993],state_connection_verified=False)]
    changes=[
        'separate_existence_formation_preservation_operation',
        'adopt_phase_transport_without_claiming_phase_formation',
        'mixed_memory_can_record_with_additional_resources',
        'energy_auxiliary_also_receives_entropy_and_records',
        'separate_common_theory_from_common_history']
    evidence=[Path(__file__),HERE/'overall_operation_hypothesis_v2.md',
        HERE/'drafts/adoption_decision.md',HERE/'drafts/review_notes.md',
        STAGE/'research_note_1003.md',STAGE/'1004/drafts/STATUS.md',
        STAGE/'999/mechanism_adoption_results.json',STAGE/'999/overall_operation_hypothesis_v1.md',
        STAGE/'990/mechanism_inventory.json',STAGE/'957/drafts/unified_operation_hypotheses_v0_2.md',
        STAGE/'993/common_candidate_v1.md',STAGE/'997/cosmological_adoption_v1.md',
        STAGE/'998/formation_resource_adoption_v1.md',
        STAGE/'1000/material_phase_transport_adoption_v1.md',
        STAGE/'1001/nuclear_material_resource_adoption_v1.md',
        STAGE/'1002/organization_lifecycle_adoption_v1.md',
        STAGE/'980/finite_thermal_records_results.json',
        STAGE/'1002/thermal_record_reuse_results.json']
    evidence += [STAGE/f'research_note_{n}.md' for n in range(991,1003)]
    evidence += [STAGE/f'{n}/research_round_{n}_checks.json' for n in (976,980,999,1000,1001,1002)]
    for n in (976,980,999,1000,1001,1002):
        assert read(STAGE/f'{n}/research_round_{n}_checks.json')['all_delivery_checks_passed']
    return dict(round=1003,date='2026-10-07',kind='mechanism_synthesis_and_adoption_audit',
        all_audit_checks_passed=True,formal_reports=1003,fresh_physical_test_groups=0,
        cumulative_test_groups=previous['cumulative_numbered_test_groups_from_1001'],
        mechanism_version='v2',groups=inventory['groups'],
        conditions=[dict(id=c['id'],title=c['title'],newly_proved_closed=False)
                    for c in inventory['conditions']],
        coverage=coverage,adoption_changes=changes,connections=joins,
        next_priority='finite nonequilibrium supply to task preparation and relational references',
        mathematical_axioms_957_changed=False,application_goal_changed=False,
        new_physics_theorem=False,joint_lifecycle_certified=False,
        all_physical_phenomena_explained=False,full_goal_completed=False,
        coverage_is_physical_proof=False,old_numerical_models_rerun=False,
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in evidence})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as dest:
            json.dump(out,dest,ensure_ascii=False,indent=2);dest.write('\n')
    else:assert out==read(TARGET)
    print(json.dumps({k:v for k,v in out.items() if k not in
        ('groups','conditions','coverage','adoption_changes','connections','source_hashes')},
        ensure_ascii=False,indent=2))
