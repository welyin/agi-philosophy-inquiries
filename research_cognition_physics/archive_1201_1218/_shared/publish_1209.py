"""Accept the internal-budget theorem without granting physical control permissions."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util
import json
import re
import sys

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'research_cognition_physics'
STAGE=BASE/'archive_1201_'
HERE=STAGE/'_shared'
ROUND=STAGE/'1209'
OLD=BASE/'archive_1161_'
PREV=HERE/'progress_checks_1208.json'
OUT=HERE/'progress_checks_1209.json'
spec=importlib.util.spec_from_file_location('thermo_frozen_tools',HERE/'publish_1206.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
sha,read,write,run,links=m.sha,m.read,m.write,m.run,m.links


def verify():
    d=read(OUT);ns,nb=m.frozen(d)
    m.m.validate_map(d['management_sha256'])
    assert sha(PREV)==d['prior_progress_sha256']
    assert run(OLD/'_shared/update_navigation.py','--check')['changed_files']==0
    assert run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation')['status']=='PASS'
    ledger=read(HERE/'dependency_ledger.json')
    assert ledger['latest_accepted_round']==1209 and ledger['stage_check_groups']==18
    assert ledger['latest_management_receipt']==OUT.name
    assert ledger['goal_objective_sha256']==d['goal_objective_sha256']
    print(json.dumps(dict(status='PASS',latest_round=1209,stage_rounds=9,stage_check_groups=18,
        old_GR_science_files_verified=ns,old_GR_acceptance_bindings_verified=nb,
        current_stage_science_files_verified=len(d['science_files_sha256']),
        stage_complete=False,goal_changed=False,goal_status='active')))


def main():
    if '--verify' in sys.argv:return verify()
    assert not OUT.exists() and not (ROUND/'acceptance.json').exists(),'Already published'
    prev=read(PREV);ns,nb=m.frozen(prev);m.m.validate_map(prev['management_sha256'])
    reviewed=['proof.md','check.py','results.json','independent_check.py','independent_results.json',
              'NEXT.md','../research_note_1209.md']
    review=(ROUND/'independent_review.md').read_text(encoding='utf-8')
    assert '限定签核通过' in review
    for name in reviewed:assert sha((ROUND/name).resolve()) in review,name
    aout=run(ROUND/'check.py');iout=run(ROUND/'independent_check.py')
    assert aout['all_checks_passed'] and aout['results_match'] and iout['status']=='PASS'
    a=read(ROUND/'results.json');i=read(ROUND/'independent_results.json')
    assert a['all_checks_passed'] and i['all_checks_passed']
    for data in [a,i]:
        for name,h in data['source_sha256'].items():assert sha(STAGE/name)==h,name
    assert i['independent_code_sha256']==sha(ROUND/'independent_check.py')
    ac=a['native_record_inclusive_certificate'];ic=i['native_finite_bounds']['complete_record_first_hit']
    assert ac['blocks']==ic['blocks'] and ac['total_materials']==ic['material_count']
    assert int(ac['total_materials'])==200*int(ac['blocks'])
    assert ac['x_occupation']==ic['word_x_counts'] and ac['y_occupation']==ic['word_y_counts']
    assert F(ic['target_distance_strict_lower'])>F('.6313')
    assert F(ic['absolute_energy_change_strict_upper'])<F('.00339')
    assert F(ic['absolute_entropy_change_strict_upper'])<F('.00423')
    assert F(ic['auxiliary_overlap_loss_upper'])<F('.0212')**2
    assert not ac['original_QQ_realization_certified'] and not ac['physical_preparation_certified']
    assert iout['finite_cases']==15 and len(a['exact_small_models'])==8
    files=reviewed+['independent_review.md','../_shared/goal.md']
    hashes={name:sha((ROUND/name).resolve()) for name in files}
    write(ROUND/'acceptance.json',dict(schema='operational_thermodynamics_round_acceptance_v1',
        round=1209,date='2026-10-14',status='ACCEPTED_WITH_SCOPE',check_groups=2,
        scope='固定有限物理H、满秩平稳参考、独立副本、全部有限保能联合操作与独立辅助下，能源/双侧熵/状态三预算共同趋零的跨组合稳健性恰选择任实beta的Gibbs族；反向有限能梯与保留二值记录共同满足预算。',
        certified=dict(three_budget_classification=True,
            gibbs_trace_distance_bound='min(1,sqrt((abs(beta)*delta+kappa)/2))',
            non_gibbs_supremum_for_all_strictly_positive_budgets=1,
            complete_auxiliary_record_retained=True,initial_SA_independence_required=True,
            auxiliary_preparation_resources_finite_per_instance=True,
            uniform_initial_resource_bound=False,native_material_distance_lower='.6313',
            native_energy_change_upper='.00339',native_auxiliary_entropy_upper='.00423',
            native_auxiliary_trace_distance_upper='.0212',
            exact_zero_budget_classification_claimed=False,
            native_QQ_control_permission_proved=False,actual_clock_instrument_implementation_proved=False),
        remaining_inputs=['reference role adoption and qualified comparison family',
            'independent reference and auxiliary preparation; potentially unbounded initial resources across instances',
            'real charges, actual clock and region/flow identity'],
        limits=['the auxiliary energy budget is reference-specific, not the old all-input operator budget',
            'Hermitian-log autonomous representation uses a freely propagated comparison and is not a native local implementation',
            'finite effective references remain valid in their original task budgets'],
        new_adopted_axioms=0,new_empirical_calibration=0,GR_hard_deliverables_completed=0,
        stage_complete=False,root_executed_both_codes=True,files_sha256=hashes))
    status=('最新验收1209；新阶段9轮、18组。共同参考的稳健分类已提升到严格保能联合过程：'
        '辅助能源、双侧熵及状态返回三预算共同趋零，恰选择Gibbs族；有限反向构造保留末记录和全部辅助。'
        '原材料数学见证漂移大于0.6313，初始资源与整体操作权限公开。候选未采用，原QQ实现及几何共同身份仍待接。')
    ledger=read(HERE/'dependency_ledger.json');req=ledger['stage_requirement_status']
    req.update(thermodynamic_structural_basis='1209: Gibbs iff uniform stability in the finite energy-conserving joint-operation menu with auxiliary energy, absolute entropy, and state-return budgets; complete terminal record included',
        actual_comparison_permissions='Finite lift and all-resource budgets are mathematically certified; native QQ, local compilation and actual clock/readout permission remain unproved; no exact-zero-budget extrapolation',
        stage_synthesis='Reference-role candidate has complete finite-process selection force; cognitive adoption and actual same-object qualification remain open',
        stage_deliverables_all_completed=False)
    ledger.update(latest_accepted_round=1209,accepted_stage_rounds=list(range(1201,1210)),stage_check_groups=18,
        current_direction_document='../1209/NEXT.md',current_work=status,latest_management_receipt=OUT.name)
    write(HERE/'dependency_ledger.json',ledger)
    p=STAGE/'README.md';parts=p.read_text(encoding='utf-8').split('\n\n');parts[1]=status
    text='\n\n'.join(parts).replace('_shared/progress_checks_1208.json','_shared/'+OUT.name)
    text=text.replace('[当前下一项](1208/NEXT.md)','[1208原安排](1208/NEXT.md)')
    text=text.replace('## 接续与边界','- [1209报告](research_note_1209.md)：内部三预算、完整记录与热参考选择；[验收](1209/acceptance.json)、[当前下一项](1209/NEXT.md)。\n\n## 接续与边界',1)
    text=text.replace('[1208原安排](1208/NEXT.md)接续1208','[当前下一项](1209/NEXT.md)接续1209')
    p.write_text(text,encoding='utf-8')
    p=STAGE/'文件索引.md';text=p.read_text(encoding='utf-8')
    text=text.replace('最新验收1208；8轮、16组','最新验收1209；9轮、18组')
    text=text.replace('## 管理与接续','|1209|[内部预算与完整记录](research_note_1209.md)|[证明](1209/proof.md)、[代码](1209/check.py)、[结果](1209/results.json)、[独立核验](1209/independent_check.py)、[独审](1209/independent_review.md)、[验收](1209/acceptance.json)|\n\n## 管理与接续',1)
    text=re.sub(r'(\|[^\n]*\|)\n\n(?=\|[^\n]*\|)',r'\1\n',text)
    text=text.replace('[当前下一项](1208/NEXT.md)','[当前下一项](1209/NEXT.md)').replace('_shared/progress_checks_1208.json','_shared/'+OUT.name)
    text+='\n- [1208后控制与共同参考角色准入](_shared/after1208_reference_contract_admission.md)：未编号审计，1209另有正式证明。\n'
    p.write_text(text,encoding='utf-8')
    oldledger=read(OLD/'_shared/dependency_ledger.json')
    oldledger['latest_management_receipt']='../../archive_1201_/_shared/'+OUT.name
    write(OLD/'_shared/dependency_ledger.json',oldledger)
    p=OLD/'_shared/navigation_notice.md';prefix,rest=p.read_text(encoding='utf-8').split('以下保留前阶段的历史状态',1)
    prefix=re.sub(r'\*\*当前状态：\*\*.*?(?=\n\n)','**当前状态：** '+status+
        ' [1209报告](archive_1201_/research_note_1209.md) · [下一项](archive_1201_/1209/NEXT.md) · '
        '[当前核验](archive_1201_/_shared/'+OUT.name+')。旧GR阶段40轮、130组保留；目标继续进行。',prefix,flags=re.S)
    p.write_text(prefix+'以下保留前阶段的历史状态'+rest,encoding='utf-8')
    navapply=run(OLD/'_shared/update_navigation.py','--apply')
    navcheck=run(OLD/'_shared/update_navigation.py','--check');assert navcheck['changed_files']==0
    navverify=run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation');assert navverify['status']=='PASS'
    write(OUT,dict(status='PUBLISHING'))
    docs=[STAGE/'README.md',STAGE/'文件索引.md',STAGE/'research_note_1209.md',*ROUND.glob('*.md'),
          HERE/'after1208_reference_contract_admission.md',
          *[BASE/n for n in ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md']]]
    nlinks=links(docs)
    for p in docs:
        text=p.read_text(encoding='utf-8')
        assert not any(ord(c)<32 and c not in '\n\r\t' for c in text)
        assert text.count('$$')%2==0,p
    science=dict(prev['science_files_sha256'])
    science.update({(ROUND/name).resolve().relative_to(ROOT).as_posix():h for name,h in hashes.items()})
    audit=dict(prev['audit_files_sha256'])
    screen=HERE/'after1208_reference_contract_admission.md';audit[screen.relative_to(ROOT).as_posix()]=sha(screen)
    management={name:sha(ROOT/name) for name in prev['management_sha256']}
    management[Path(__file__).relative_to(ROOT).as_posix()]=sha(Path(__file__))
    accs=dict(prev['acceptance_sha256']);accs['1209']=sha(ROUND/'acceptance.json')
    write(OUT,dict(schema='thermodynamics_progress_1209_v1',date='2026-10-14',status='PASS',
        accepted_rounds=list(range(1201,1210)),stage_check_groups=18,
        application_goal_status='active',stage_complete=False,goal_changed=False,
        goal_objective_sha256=ledger['goal_objective_sha256'],acceptance_sha256=accs,
        science_files_sha256=science,management_sha256=management,
        screen_files_sha256=prev['screen_files_sha256'],audit_files_sha256=audit,
        prior_progress_path=PREV.relative_to(ROOT).as_posix(),prior_progress_sha256=sha(PREV),
        previous_GR_science_files_verified=ns,previous_GR_acceptance_bindings_verified=nb,
        navigation_apply=navapply,navigation_check=navcheck,navigation_verify=navverify,
        markdown_files_checked=len(docs),local_links_checked=nlinks,missing_links=[],
        stage_requirement_status=req,new_adopted_axioms=0,new_empirical_calibration=0,
        GR_hard_deliverables_completed=0,science_program_outputs=dict(author=aout,independent=iout),
        previous_goal_turn_classification='no new accepted round: user clarification; revalidated 1208 and existing audits',
        this_goal_turn_result='progress: internal-budget reference theorem, record-inclusive finite lift, deletion counterexample, independent checks and scoped acceptance'))
    verify()


if __name__=='__main__':main()
