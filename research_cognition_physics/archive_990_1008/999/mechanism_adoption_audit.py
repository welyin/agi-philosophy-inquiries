"""Reproducible mechanism/evidence audit, not a numerical physics experiment."""
from pathlib import Path
import argparse
import hashlib
import json

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
TARGET=HERE/'mechanism_adoption_results.json'


def read(path):return json.loads(path.read_text('utf-8-sig'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    base=read(STAGE/'990/mechanism_inventory.json')
    previous=read(STAGE/'998/research_round_998_checks.json')
    assert previous['all_delivery_checks_passed'] and not previous['full_goal_completed']
    groups=base['groups']
    assert sorted(n for values in groups.values() for n in values)==list(range(1,28))
    assert [r['id'] for r in base['conditions']]==[f'C{i:02d}' for i in range(1,28)]
    conditions=[dict(id=r['id'],title=r['title'],group=r['group'],
        newly_proved_closed=False) for r in base['conditions']]
    # This is a declared scope map. Coverage is not a proof of physical recovery.
    rows=[
        ('quantum','干涉、纠缠、测量反作用',[1,2,3,4],'inherited_conditional_interfaces'),
        ('records','稳定记录与宏观事实',[3,4,19,23],'inherited_selective_interaction'),
        ('spacetime','空间、方向与时间比较',[5,6,7,8,9,10],'conditional_geometry_and_physical_input'),
        ('radiation','光、辐射与光学响应',[14,17,18,21],'physical_input_and_partial_interfaces'),
        ('motion_mass','惯性、运动与质量',[10,15,17,18,22],'physical_input_and_partial_interfaces'),
        ('ordinary_matter','原子、化学与普通物质',[15,17,21],'external_stability_theorem_adopted'),
        ('nuclear_weak','核、弱作用与粒子反应',[14,15,16,17,21],'standard_model_input_not_new_derivation'),
        ('phases_transport','宏观相、弹性、声与输运',[17,19,23,26],'priority_mechanism_gap'),
        ('thermal_arrow','热化、耗散与时间箭头',[19,23,26],'conditional_finite_resource_mechanisms'),
        ('gravity','重力、钟效应、透镜与波',[10,11,12,18,22],'einstein_effective_input_and_partial_interfaces'),
        ('horizons','黑洞外部、视界与半经典效应',[11,12,13,25],'mature_effective_route_not_project_recovery'),
        ('cosmic_formation','膨胀、红移与结构形成',[9,17,19,24],'conditional_cosmic_and_formation_routes'),
        ('asymmetry_neutrino','物质不对称、中微子质量与混合',[15,16,17,21,24],'explicit_extension_not_actual_yield'),
        ('dark_matter','暗物质现象',[17,18,24],'optional_candidate_not_identification'),
        ('acceleration','加速膨胀',[12,18,24],'effective_lambda_input'),
        ('cross_scale','跨尺度共同预测',[20,27],'joint_certificate_open')]
    coverage=[dict(id=key,phenomena=title,conditions=[f'C{n:02d}' for n in ids],
        status=status,complete_physical_recovery=False) for key,title,ids,status in rows]
    assert {c for row in coverage for c in row['conditions']}=={r['id'] for r in conditions}
    summary=(HERE/'overall_operation_hypothesis_v1.md').read_text('utf-8')
    for row in coverage:assert '|'+row['phenomena']+'|' in summary,row['id']
    changes=[
        dict(id='selective_open_organization',conditions=['C03','C04','C17','C19','C23'],
             kind='adoption_of_existing_mechanisms',new_cognitive_axiom=False,
             dependencies=[972,975,976,980,988,929],
             claim='protected task information can coexist with selective physical exchanges',
             not_claimed='a joint material lifecycle has already been certified'),
        dict(id='ordinary_matter_stability',conditions=['C15','C17','C21'],
             kind='external_theorem_with_explicit_working_regime',new_cognitive_axiom=False,
             dependencies=[981,993],
             claim='fermionic nonrelativistic Coulomb matter has the stated extensive lower bound',
             not_claimed='all real material phases, formation rates or memory lifetimes follow'),
        dict(id='strong_gravity_scope',conditions=['C11','C12','C13','C25'],
             kind='separation_of_effective_routes_and_uv_questions',new_cognitive_axiom=False,
             dependencies=[990,993,997],
             claim='classical and semiclassical horizon tasks are not all ultraviolet completion',
             not_claimed='the project has solved complete evaporation or information loss')]
    joins=[
        dict(left=976,right=980,issue='different eta, reference, contact and history',certified_common_history=False),
        dict(left=975,right=980,issue='comparison Gibbs reference is not a prepared physical bath',certified_common_history=False),
        dict(left=980,right=988,issue='mixed thermal output is not the required pure dark blank',certified_common_history=False),
        dict(left=929,right=976,issue='engineered protocol Hamiltonian lacks native-material matching',certified_common_history=False),
        dict(left=998,right=980,issue='radiation receiver is not the prepared same-material copies',certified_common_history=False)]
    evidence=[STAGE/'990/mechanism_inventory.json',STAGE/'998/research_round_998_checks.json',
        STAGE/'999/drafts/STATUS.md',STAGE/'957/drafts/unified_operation_hypotheses_v0_2.md',
        STAGE/'981/drafts/common_parent_contract_v1.md',STAGE/'993/common_candidate_v1.md',
        STAGE/'998/formation_resource_adoption_v1.md',
        STAGE.parent/'archive_100_124/research_note_105.md',
        STAGE.parent/'archive_100_124/research_note_106.md',
        STAGE.parent/'archive_370_428/research_note_391.md']
    evidence += [STAGE/f'research_note_{n}.md' for n in (929,972,975,976,980,988,990,992,998)]
    evidence += [STAGE/f'{n}/research_round_{n}_checks.json' for n in (929,972,975,976,980,988)]
    evidence += [STAGE/'research_note_999.md',HERE/'overall_operation_hypothesis_v1.md',
        HERE/'drafts/adoption_decision.md',HERE/'drafts/review_notes.md',Path(__file__)]
    for n in (929,972,975,976,980,988):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        assert old['all_delivery_checks_passed'],n
        assert not old['full_goal_completed'],n
    sources=[
        dict(url='https://arxiv.org/pdf/math-ph/0209034',extent='author review sections 2 and 3',
             role='nonrelativistic Coulomb stability assumptions and extensive lower bound'),
        dict(url='https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.35.687',extent='publisher page',
             role='original fermion kinetic inequality and stability paper'),
        dict(url='https://arxiv.org/html/quant-ph/0408125',extent='sections VI and VII, read-only reviewer',
             role='pointer stability is insufficient for redundant accessible records'),
        dict(url='https://arxiv.org/html/0904.0418',extent='equations 3 and 10 with conditions, read-only reviewer',
             role='specific environment information relation, not a general heat-flow law'),
        dict(url='https://arxiv.org/html/1306.4352',extent='inherited theorem 3 and reviewer scope check',
             role='finite receiver and correlation accounting'),
        dict(url='https://arxiv.org/html/gr-qc/9912119',extent='sections 2 and 3',
             role='classical and semiclassical horizon scope, not a present experimental claim')]
    return dict(round=999,date='2026-10-07',kind='mechanism_adoption_and_scope_audit',
        all_audit_checks_passed=True,formal_reports=999,fresh_physical_test_groups=0,
        cumulative_test_groups=previous['cumulative_numbered_test_groups_from_997'],
        mechanism_version='v1',groups=groups,conditions=conditions,coverage=coverage,
        adoption_changes=changes,non_joinable_evidence=joins,primary_sources=sources,
        next_priority='material phases, collective modes and macroscopic transport',
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
        ('groups','conditions','coverage','adoption_changes','non_joinable_evidence','primary_sources','source_hashes')},
        ensure_ascii=False,indent=2))
