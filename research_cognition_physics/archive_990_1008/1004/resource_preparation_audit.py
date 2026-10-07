"""Classify saved resources and task obligations; no new physics experiment."""
from pathlib import Path
import argparse
import hashlib
import json

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
TARGET=HERE/'resource_preparation_audit_results.json'


def read(path):return json.loads(path.read_text('utf-8-sig'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    prev=read(STAGE/'1003/research_round_1003_checks.json')
    resource=read(STAGE/'992/cooling_availability_results.json')
    reuse=read(STAGE/'1002/thermal_record_reuse_results.json')
    assert resource['all_scientific_checks_passed'] and prev['all_delivery_checks_passed']
    rows=resource['rows'];p0=rows[0]['actual_populations']
    diagnostic=[]
    for row in rows:
        assert row['actual_populations']==p0
        assert abs(row['actual_material_entropy']-resource['initial_material_entropy'])<1e-14
        if row['a']!=1:assert row['availability']>0
        diagnostic.append(dict(a=row['a'],availability=row['availability'],
            actual_populations_unchanged=True,reference_is_prepared_bath=False,
            work_or_coherent_reference_generated=False))
    assert reuse['auxiliary_initial_pure_energy_state_is_input']
    assert reuse['engineered_energy_matched_interaction_is_input']
    contracts=[
        dict(id='population_and_finite_records',resource='population nonequilibrium plus actual receiver',
            inherited=[992,980,1002],natural_supply_certified=False,
            universal_pure_memory_required=False),
        dict(id='relational_encoding_and_joint_reference',resource='prepared accessible joint relation',
            inherited=[447,457,709,925],natural_supply_certified=False,
            external_absolute_phase_required=False),
        dict(id='bare_non_covariant_operation',resource='accounted asymmetry or changed relational task',
            inherited=[421],natural_supply_certified=False,
            supplied_by_free_energy_alone=False)]
    files=[Path(__file__),HERE/'drafts/STATUS.md',HERE/'drafts/NEXT.md',
        HERE/'drafts/adoption_decision.md',HERE/'resource_preparation_adoption_v1.md',
        STAGE/'research_note_1004_working.md',STAGE/'1003/research_round_1003_checks.json',
        STAGE/'1003/overall_operation_hypothesis_v2.md',STAGE/'992/cooling_availability_results.json',
        STAGE/'1002/thermal_record_reuse_results.json',STAGE/'research_note_925.md',
        STAGE/'research_note_992.md',STAGE/'_shared/notes/hqca_full_cell_scope_review.md',
        STAGE.parent/'archive_370_428/research_note_421.md',
        STAGE.parent/'archive_429_466/research_note_447.md',
        STAGE.parent/'archive_429_466/research_note_457.md',
        STAGE.parent/'archive_702_741/research_note_709.md',
        STAGE.parent/'archive_467_530/research_note_522.md',
        STAGE.parent/'archive_467_530/research_note_523.md']
    return dict(kind='working_resource_adoption_audit',date='2026-10-07',
        all_audit_checks_passed=True,formal_reports=1003,new_scientific_test_groups=0,
        cumulative_numbered_test_groups=prev['cumulative_numbered_test_groups_from_1002'],
        resource_classification=contracts,inherited_992_diagnostics=diagnostic,
        covariance_requires_additive_energy_and_invariant_joint_input=True,
        full_energy_conservation_alone_implies_bare_covariance=False,
        spatial_reference_522_523_generates_time_asymmetry=False,
        nonstationary_cosmos_is_global_thermal_operation=False,
        explicit_supply_mechanisms_are_adopted_not_constructed=True,
        natural_preparation_of_1002_auxiliary_certified=False,
        joint_lifecycle_certified=False,full_goal_completed=False,app_goal_changed=False,
        next_priority='finite horizon task: geometry, field state, detector and shared sources',
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as dest:
            json.dump(out,dest,ensure_ascii=False,indent=2);dest.write('\n')
    else:assert out==read(TARGET)
    print(json.dumps({k:v for k,v in out.items() if k not in
        ('source_hashes','resource_classification','inherited_992_diagnostics')},ensure_ascii=False,indent=2))
