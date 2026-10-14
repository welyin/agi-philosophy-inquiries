"""Publish scoped permission classification; retain all accepted scientific files."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT/'research_cognition_physics'
STAGE = BASE/'archive_1201_'
HERE = STAGE/'_shared'
ROUND = STAGE/'1208'
OLD = BASE/'archive_1161_'
PREV = HERE/'progress_checks_1207.json'
OUT = HERE/'progress_checks_1208.json'
spec=importlib.util.spec_from_file_location('thermo_frozen_tools',HERE/'publish_1206.py')
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
sha,read,write,run,links=m.sha,m.read,m.write,m.run,m.links


def verify():
    d=read(OUT); ns,nb=m.frozen(d)
    m.m.validate_map(d['management_sha256'])
    assert sha(PREV)==d['prior_progress_sha256']
    assert run(OLD/'_shared/update_navigation.py','--check')['changed_files']==0
    assert run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation')['status']=='PASS'
    ledger=read(HERE/'dependency_ledger.json')
    assert ledger['latest_accepted_round']==1208 and ledger['stage_check_groups']==16
    assert ledger['latest_management_receipt']==OUT.name
    assert ledger['goal_objective_sha256']==d['goal_objective_sha256']
    print(json.dumps({'status':'PASS','latest_round':1208,'stage_rounds':8,'stage_check_groups':16,
        'old_GR_science_files_verified':ns,'old_GR_acceptance_bindings_verified':nb,
        'current_stage_science_files_verified':len(d['science_files_sha256']),
        'goal_status':'active','goal_changed':False,'stage_complete':False}))


def main():
    if '--verify' in sys.argv:return verify()
    assert not OUT.exists() and not (ROUND/'acceptance.json').exists(),'Already published'
    prev=read(PREV); ns,nb=m.frozen(prev); m.m.validate_map(prev['management_sha256'])
    review=(ROUND/'independent_review.md').read_text(encoding='utf-8')
    assert '限定签核通过' in review
    for name in ['proof.md','check.py','results.json','independent_check.py','independent_results.json',
                 '../research_note_1208.md','NEXT.md']:
        assert sha((ROUND/name).resolve()) in review, 'Unbound review: '+name
    author=run(ROUND/'check.py'); independent=run(ROUND/'independent_check.py')
    assert author['all_checks_passed'] and author['results_match']
    assert independent['all_checks_passed'] and independent==read(ROUND/'independent_results.json')
    result=read(ROUND/'results.json')
    for data in [result,independent]:
        for name,digest in data['source_sha256'].items():assert sha(STAGE/name)==digest,name
    native=result['native_parity_witness']
    assert int(native['blocks'])==31894572775086737403576853338681
    assert int(native['total_materials'])==50*int(native['blocks'])
    assert F(native['trace_distance_strict_lower_exact'])>F('.47')
    assert F(native['energy_interval'][1])<F('.0008474425022')
    assert F(result['finite_charge_contrast']['52_block_distance_exact'])>F('.38')
    assert result['permission_lattice']['index_in_sum_zero_lattice']==2
    assert not native['physical_preparation_control_readout_certified']
    other=independent['native_doubled_cold_witness']
    assert native['blocks']==other['blocks'] and native['total_materials']==other['total_material_copies']
    assert native['occupation_difference']==other['new_occupation_difference']
    assert F(native['energy_interval'][0])<F(other['energy_error_lower'])
    assert F(other['energy_error_upper'])<F(native['energy_interval'][1])
    assert F(other['trace_distance_lower'])>F('.47')
    assert F(other['contrast_lower'])<=F(native['conditional_asymmetry_exact'])<=F(other['contrast_upper'])
    assert result['finite_charge_contrast']['52_block_distance_exact']==independent['qutrit_parity_witness']['trace_distance_main_exact']
    assert independent['native_QQ_algebra']['index_in_Lambda']==2
    assert not independent['native_QQ_algebra']['full_parity_block_unitary_reachability_claimed']
    assert result['finite_charge_return_counterexample']['WJ_equals_JT_exact_integer_columns']
    assert independent['finite_parity_catalyst_control']['qutrit_controlled_flip_W']['joint_isometry_WJ_equals_JT']
    assert independent['finite_parity_catalyst_control']['qutrit_controlled_flip_W']['P_T_P_is_neither_T_nor_minus_T']
    assert not independent['finite_parity_catalyst_control']['native_preparation_or_realization_claimed']
    files=['proof.md','check.py','results.json','independent_check.py','independent_results.json',
           'independent_review.md','permission_review.md','NEXT.md','../research_note_1208.md','../_shared/goal.md']
    hashes={name:sha((ROUND/name).resolve()) for name in files}
    scope=('固定谱、满秩对角参考、独立副本及人数模L的完整酉权限类下，跨组合稳健性恰等价于'
        'log p+beta E属于L实张成的正交补；有限指数限制不增加Gibbs自由，真正加性荷可以增加。'
        '原准备门保全材料人数，QQ破人数守奇偶；矩阵元生成不授予实际转换控制。'
        '有限相干辅助精确返回的反例表明整体奇偶守恒不能自动决定目标菜单。')
    write(ROUND/'acceptance.json',{'schema':'operational_thermodynamics_round_acceptance_v1',
        'round':1208,'date':'2026-10-14','status':'ACCEPTED_WITH_SCOPE','check_groups':2,'scope':scope,
        'certified':{'restricted_permission_classification':'exists real beta: log p+beta E annihilates L',
            'uniform_bound':'min(1,sqrt(abs(beta)*delta/2))',
            'nonqualifying_reference_positive_delta_supremum':1,
            'finite_index_subgroup_adds_no_reference_freedom':True,
            'actual_menu_requires_charge_preservation_and_witness_access':True,
            'native_preparation_count_preservation_includes_stocks':True,
            'native_QQ_matrix_element_lattice_index':2,
            'native_QQ_full_comparison_reachability_proven':False,
            'returned_finite_charge_reference_can_enable_target_parity_violation':True,
            'global_conservation_and_auxiliary_return_alone_imply_target_U_L':False,
            'qutrit_parity_example_52_block_distance_lower':'.38',
            'native_cold_parity_blocks':native['blocks'],
            'native_cold_parity_materials':native['total_materials'],
            'native_cold_distance_lower':'.47','native_cold_total_energy_error_upper':'.0008474425022',
            'physical_resources_for_large_witness_certified':False},
        'remaining_inputs':['fixed physical H, declared diagonal reference and iid preparations',
            'actual comparison primitives with both permission bounds, resources and timing',
            'reference role/stability adoption, conserved charges and physical temperature',
            'common clock, region, full flux and independent geometric response'],
        'limits':['integer generation is not controllability','finite-index symmetry is not a full permission certificate',
            'finite-basis checks do not prove the all-scale theorem; proof does',
            'finite effective references remain valid in their stated task budgets'],
        'new_adopted_axioms':0,'new_empirical_calibration':0,'GR_hard_deliverables_completed':0,
        'stage_complete':False,'root_executed_both_codes':True,'files_sha256':hashes})
    status=('最新验收1208；新阶段8轮、16组。受限比较的稳健参考分类已明确：有限模数如奇偶不增加'
        'Gibbs自由，真正加性荷保留广义参考。原准备门保人数，QQ矩阵元破人数守奇偶；实际跨规模'
        '转换权限尚未取得，不能靠扩布尔程序补足。相干辅助归还可绕过目标奇偶限制；目标菜单须单独认证。'
        '候选仍未采用，几何共同身份待接。')
    ledger=read(HERE/'dependency_ledger.json'); req=ledger['stage_requirement_status']
    req.update(thermodynamic_structural_basis='1208: restricted comparison classification; finite-index constraints preserve Gibbs selection, additive charges retain generalized references. Actual native comparison access and cognitive adoption remain open',
        actual_comparison_permissions='Existing preparation/copy/swap preserves all material counts; native QQ breaks counts but preserves parity. Matrix-element lattice index two is not controllability; an exactly returned coherent finite-charge reference can bypass target parity, so the target menu requires its own qualification',
        stage_synthesis='Conditional reference-selection structure refined; actual qualified comparisons and common geometric identities remain open',
        stage_deliverables_all_completed=False)
    ledger.update(latest_accepted_round=1208,accepted_stage_rounds=list(range(1201,1209)),stage_check_groups=16,
        current_direction_document='../1208/NEXT.md',current_work=status,latest_management_receipt=OUT.name,
        stage_requirement_status=req)
    write(HERE/'dependency_ledger.json',ledger)
    p=STAGE/'README.md'; parts=p.read_text(encoding='utf-8').split('\n\n'); parts[1]=status
    text='\n\n'.join(parts).replace('_shared/progress_checks_1207.json','_shared/'+OUT.name)
    text=text.replace('[当前下一项](1207/NEXT.md)','[1207原安排](1207/NEXT.md)')
    text=text.replace('## 接续与边界','- [1208报告](research_note_1208.md)：许可比较与守恒荷的稳健分类；[验收](1208/acceptance.json)、[当前下一项](1208/NEXT.md)。\n\n## 接续与边界',1)
    text=text.replace('[1207原安排](1207/NEXT.md)接续1207','[当前下一项](1208/NEXT.md)接续1208')
    p.write_text(text,encoding='utf-8')
    p=STAGE/'文件索引.md'; text=p.read_text(encoding='utf-8')
    text=text.replace('最新验收1207；7轮、14组','最新验收1208；8轮、16组')
    text=text.replace('## 管理与接续','|1208|[许可比较与守恒荷](research_note_1208.md)|[证明](1208/proof.md)、[代码](1208/check.py)、[结果](1208/results.json)、[独立核验](1208/independent_check.py)、[独审](1208/independent_review.md)、[权限](1208/permission_review.md)、[验收](1208/acceptance.json)|\n\n## 管理与接续',1)
    text=re.sub(r'(\|[^\n]*\|)\n\n(?=\|[^\n]*\|)',r'\1\n',text)
    text=text.replace('[当前下一项](1207/NEXT.md)','[当前下一项](1208/NEXT.md)').replace('_shared/progress_checks_1207.json','_shared/'+OUT.name)
    p.write_text(text,encoding='utf-8')
    oldledger=read(OLD/'_shared/dependency_ledger.json')
    oldledger['latest_management_receipt']='../../archive_1201_/_shared/'+OUT.name
    write(OLD/'_shared/dependency_ledger.json',oldledger)
    p=OLD/'_shared/navigation_notice.md'; prefix,rest=p.read_text(encoding='utf-8').split('以下保留前阶段的历史状态',1)
    prefix=re.sub(r'\*\*当前状态：\*\*.*?(?=\n\n)','**当前状态：** '+status+
        ' [1208报告](archive_1201_/research_note_1208.md) · [下一项](archive_1201_/1208/NEXT.md) · '
        '[当前核验](archive_1201_/_shared/'+OUT.name+')。旧GR阶段40轮、130组保留；目标继续进行。',prefix,flags=re.S)
    p.write_text(prefix+'以下保留前阶段的历史状态'+rest,encoding='utf-8')
    navapply=run(OLD/'_shared/update_navigation.py','--apply')
    navcheck=run(OLD/'_shared/update_navigation.py','--check'); assert navcheck['changed_files']==0
    navverify=run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation'); assert navverify['status']=='PASS'
    write(OUT,{'status':'PUBLISHING'})
    docs=[STAGE/'README.md',STAGE/'文件索引.md',STAGE/'research_note_1208.md',*ROUND.glob('*.md'),
          *[BASE/name for name in ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md']]]
    nlinks=links(docs)
    for p in docs:assert not any(ord(c)<32 and c not in '\n\r\t' for c in p.read_text(encoding='utf-8'))
    science=dict(prev['science_files_sha256'])
    science.update({(ROUND/name).resolve().relative_to(ROOT).as_posix():h for name,h in hashes.items()})
    management={name:sha(ROOT/name) for name in prev['management_sha256']}
    management[Path(__file__).relative_to(ROOT).as_posix()]=sha(Path(__file__))
    accs=dict(prev['acceptance_sha256']); accs['1208']=sha(ROUND/'acceptance.json')
    write(OUT,{'schema':'thermodynamics_progress_1208_v1','date':'2026-10-14','status':'PASS',
        'accepted_rounds':list(range(1201,1209)),'stage_check_groups':16,
        'application_goal_status':'active','stage_complete':False,'goal_changed':False,
        'goal_objective_sha256':ledger['goal_objective_sha256'],'acceptance_sha256':accs,
        'science_files_sha256':science,'management_sha256':management,
        'screen_files_sha256':prev['screen_files_sha256'],'audit_files_sha256':prev['audit_files_sha256'],
        'prior_progress_path':PREV.relative_to(ROOT).as_posix(),'prior_progress_sha256':sha(PREV),
        'previous_GR_science_files_verified':ns,'previous_GR_acceptance_bindings_verified':nb,
        'navigation_apply':navapply,'navigation_check':navcheck,'navigation_verify':navverify,
        'markdown_files_checked':len(docs),'local_links_checked':nlinks,'missing_links':[],
        'stage_requirement_status':req,'new_adopted_axioms':0,'new_empirical_calibration':0,
        'GR_hard_deliverables_completed':0,'science_program_outputs':{'author':author,
            'independent_all_checks_passed':True,'independent_result_sha256':sha(ROUND/'independent_results.json')},
        'previous_goal_turn_classification':'progress: 1207 accepted with conditional classification and next task',
        'this_goal_turn_result':'progress: restricted reference theorem; finite-index vs additive-charge distinction; native primitive audit localizes the genuine conversion-access gap'})
    verify()


if __name__=='__main__':main()
