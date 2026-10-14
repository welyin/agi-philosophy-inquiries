"""Accept scoped event/resource and shared-carrier results without altering the goal."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util
import json
import re
import sys

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'research_cognition_physics';STAGE=BASE/'archive_1201_'
HERE=STAGE/'_shared';OLD=BASE/'archive_1161_'
PREV=HERE/'progress_checks_1216.json';OUT=HERE/'progress_checks_1218.json'
spec=importlib.util.spec_from_file_location('frozen_tools',HERE/'publish_1206.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
sha,read,write,run,links=m.sha,m.read,m.write,m.run,m.links
ROUNDS=[1217,1218]
AUDIT=HERE/'after1216_four_laws_scope_screen.md'
STATUS=('最新验收1218；本阶段18轮、36组。1217把同一末事件校准接到完整输入的谱尾与可加反向相对熵预算，'
        '区分逐份有限协议不能精确达零与固定预算整族不能任意逼近；前者无需统一资源增长律。'
        '1218给同一载体跨段安全资格及联合缺陷1/4的删除，逐段原准备校准和单体均能不够。'
        '参考／工作角色及真实菜单资格仍待明确，未采用新公理，四律尚未全部完成。')

def verify():
    d=read(OUT);old,bindings=m.frozen(d);m.m.validate_map(d['management_sha256'])
    assert sha(PREV)==d['prior_progress_sha256']
    ledger=read(HERE/'dependency_ledger.json')
    assert ledger['latest_accepted_round']==1218 and ledger['stage_check_groups']==36
    assert ledger['latest_management_receipt']==OUT.name
    assert not ledger['stage_requirement_status']['stage_deliverables_all_completed']
    assert ledger['goal_objective_sha256']==d['goal_objective_sha256']
    assert run(OLD/'_shared/update_navigation.py','--check')['changed_files']==0
    assert run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation')['status']=='PASS'
    print(json.dumps(dict(status='PASS',latest_round=1218,stage_rounds=18,stage_check_groups=36,
        old_GR_science_files_verified=old,old_GR_acceptance_bindings_verified=bindings,
        current_stage_science_files_verified=len(d['science_files_sha256']),
        goal_changed=False,application_goal_status='active',stage_complete=False)))

def scientific_checks():
    a=read(STAGE/'1217/results.json');b=read(STAGE/'1217/independent_results.json')
    assert a['all_checks_passed'] and b['all_checks_passed']
    assert a['check_groups']==8 and b['check_groups']==6
    ag=a['groups'];bg={x['name']:x for x in b['groups']}
    assert all(x['pass'] for x in ag.values()) and all(x['passed'] for x in bg.values())
    native=ag['native_objects_and_quantifiers']
    assert native['D']==16 and native['bad_rank']==12
    assert not native['same_trial_equality_claim'] and not native['uses_1205_correlated_preparation_as_channel']
    assert not native['new_hardware_or_axioms']
    cert=ag['actual_A_and_K_spectral_bounds']
    assert F(cert['exact_unital_tail'])==F(bg['native_exact_joint_spectrum']['ideal_rank_12_tail'])
    independent={str(x['rectangular_contact_duration']):x for x in bg['same_native_actual_A_and_K_event_certificates']['events']}
    binary={x['contact']:x for x in ag['binary_and_additive_uniform_bounds']['native_contacts']}
    for x in cert['contacts']:
        y=independent[x['contact']];z=binary[x['contact']]
        assert F(x['a'])==F(y['q_lower']) and F(x['alpha'])==F(y['effect_trace_lower'])
        assert F(x['sharp_tail_exact'])==F(y['exact_fractional_spectral_tail'])
        assert F(x['sharp_tail_exact'])>F(x['strict_display_lower'])>F(x['old_1214_lower'])
        assert F(z['binary_root_bracket']['lower'])<=F(y['binary_root_interval'][0])
        assert F(y['binary_root_interval'][1])<=F(z['binary_root_bracket']['upper'])
        assert F(z['explicit_exponential_rational_lower'])==F(y['strict_explicit_lower'])
        assert set(x['event_certificates'])=={'A_excited','K_nonzero'}
    av=ag['reverse_relative_entropy_budget']['V_CH_product_sum']
    bv=bg['independent_rational_log_and_native_V0']['sum_of_marginal_V0_intervals']
    assert F(bv[0])<=F(av['lower'])<=F(av['upper'])<=F(bv[1])
    q=bg['fixed_dimension_resource_family_and_quantifier_deletion']
    assert q['Q1_pointwise_no_uniform_resource_growth_required']
    assert q['Q2_uniform_bound_not_implied_by_individual_faithfulness']
    assert q['Q3_time_budget_is_an_extra_sufficient_interface']
    assert ag['binary_and_additive_uniform_bounds']['uniform_time_budget_not_derived']
    assert ag['correlations_pure_inventory_and_terminal_record_scope']['no_universe_wide_faithfulness_assumption']
    assert not ag['three_quantifiers_and_stopping_line']['four_laws_complete']
    a=read(STAGE/'1218/results.json');b=read(STAGE/'1218/independent_results.json')
    assert a['all_checks_passed'] and b['all_checks_passed']
    assert a['check_groups']==8 and b['check_groups']==6
    ag={x['name']:x for x in a['groups']};bg=b['groups']
    assert all(x['passed'] for x in bg.values())
    deletion=ag['same finite carrier path deletion']
    assert deletion['normalized_unital_defect']=='1/4'
    for x in bg['same_preparation_individual_checks_do_not_compose']['examples']:
        assert x['joint_diagonal_exact']==deletion['exact_joint_probabilities']==['1/8','3/8','3/8','1/8']
        assert x['normalized_unital_defect_exact']=='1/4'
        assert x['mean_A_exact']==x['mean_B_exact']=='1/2'
        assert F(x['mean_C_exact'])==F(x['N'],2)
        assert x['both_single_original_preparations_unital'] and x['both_output_marginals_maxmixed']
    assert ag['safe cross-material closure']['same_carrier_not_reset']
    assert ag['positive common-phase representation']['same_phase_for_all_contacts']
    assert ag['actual middle nonfeedback QND']['all_records_retained']
    er=ag['reachable operator calibration and path error']
    assert er['uses_operator_supremum_not_average_path_probability']
    assert F(er['deleted_path_eta_squared'])==F(1,2)
    assert F(er['calibration_a_squared'])<F(er['strict_half_diamond_upper'])**2
    assert bg['full_reachable_subspace_error_budget']['endpoint_error']<=bg['full_reachable_subspace_error_budget']['summed_bound']
    assert not b['scope']['work_role_adopted'] and not b['scope']['actual_native_QQ_controller']
    assert not b['scope']['all_time_safe_reuse']
    assert b['scope']['no_universal_thermodynamic_or_third_law_failure_claim']

def main():
    if '--verify' in sys.argv:return verify()
    assert not OUT.exists()
    prev=read(PREV);old,bindings=m.frozen(prev);m.m.validate_map(prev['management_sha256'])
    outputs={};files={}
    for number in ROUNDS:
        folder=STAGE/str(number);assert not (folder/'acceptance.json').exists()
        names=['proof.md','check.py','results.json','independent_check.py',
               'independent_results.json','NEXT.md',f'../research_note_{number}.md']
        review=(folder/'independent_review.md').read_text(encoding='utf-8')
        assert '限定签核通过' in review
        for name in names:assert sha((folder/name).resolve()) in review,(number,name)
        for name in ['proof.md','NEXT.md',f'../research_note_{number}.md']:
            value=(folder/name).read_text(encoding='utf-8')
            assert not any(x in value for x in (r'\[',r'\]',r'\(',r'\)'))
        a=run(folder/'check.py');b=run(folder/'independent_check.py')
        for result in (a,b):
            assert result['status']=='PASS' and result['all_checks_passed'] and result['results_match']
        for result in (read(folder/'results.json'),read(folder/'independent_results.json')):
            for name,digest in result.get('source_sha256',{}).items():
                path=next((p for p in (BASE/name,STAGE/name,ROOT/name) if p.is_file()),None)
                assert path and sha(path)==digest,(number,name)
        names+=['independent_review.md','../_shared/goal_four_laws.md']
        files[number]={name:sha((folder/name).resolve()) for name in names}
        outputs[str(number)]={'author':a,'independent':b}
    scientific_checks()
    certified={
      1217:dict(same_terminal_event_calibration=True,complete_input_spectrum_required=True,
                reverse_relative_entropy_budget=True,pointwise_vs_uniform_quantifiers_separated=True,
                uniform_resource_growth_not_required_for_pointwise_finite_no_exact_cooling=True,
                physical_infinite_cost_from_infinite_vacancy_claimed=False,
                universal_finite_time_bound_claimed=False),
      1218:dict(same_carrier_without_reset=True,whole_reachable_input_subspace_required=True,
                individual_original_preparation_tests_insufficient=True,
                exact_joint_normalized_unital_defect='1/4',both_single_marginals_unchanged=True,
                all_three_mean_energies_unchanged=True,all_time_safe_reuse_claimed=False,
                work_role_adopted=False,actual_native_QQ_controller_claimed=False)}
    scopes={
      1217:'同一真实末失败事件的中性校准与完整输入谱给精确分数谱尾；反向相对熵给可加资源预算。逐份有限合格协议排精确零无需统一n(t)，固定时间整族界仍需独立资格；不禁止纯态或宣称其物理成本无穷。',
      1218:'同一有限工作载体全路径安全时整段保幺；原准备逐段合格仍可在连续复用时产生联合1/4缺陷，各单体与均能不变。保完整关联、记录及误差，未赋原QQ新权限或采用公理。'}
    for number in ROUNDS:
        write(STAGE/str(number)/'acceptance.json',dict(
            schema='operational_thermodynamics_round_acceptance_v1',round=number,date='2026-10-15',
            status='ACCEPTED_WITH_SCOPE',check_groups=2,scope=scopes[number],
            root_executed_both_codes=True,files_sha256=files[number],certified=certified[number],
            new_adopted_axioms=0,new_empirical_calibration=0,GR_hard_deliverables_completed=0,
            stage_complete=False))
    ledger=read(HERE/'dependency_ledger.json');req=ledger['stage_requirement_status']
    req.update(common_process='1217 checks terminal-event qualification in the existing actual contact; 1218 separates safe whole-process work-role composition from individually calibrated reuse',
        third_law='Pointwise finite faithful qualified protocols exclude exact zero without uniform n(t). 1217 gives terminal-event spectral and additive reverse-relative-entropy resource certificates; fixed-time family bounds and infinite-participant extensions retain independent conditions',
        stage_synthesis='Zeroth/first/second-law prior scoped proofs and conditional Nernst remain; roles, cyclic qualifications and real menus require explicit common-scope audit rather than repairing every candidate',
        stage_deliverables_all_completed=False)
    ledger.update(latest_accepted_round=1218,accepted_stage_rounds=list(range(1201,1219)),
        stage_check_groups=36,current_direction_document='../1217/NEXT.md',
        current_work=STATUS,latest_management_receipt=OUT.name)
    write(HERE/'dependency_ledger.json',ledger)
    p=STAGE/'README.md';value=p.read_text(encoding='utf-8');parts=value.split('\n\n');parts[1]=STATUS
    value='\n\n'.join(parts).replace('_shared/'+PREV.name,'_shared/'+OUT.name)
    rows='\n'.join(f'- [{n}报告](research_note_{n}.md)：{scopes[n]} [验收]({n}/acceptance.json)、[下一项]({n}/NEXT.md)。' for n in ROUNDS)
    value=value.replace('## 接续与边界',rows+'\n\n## 接续与边界',1)
    p.write_text(value,encoding='utf-8')
    p=STAGE/'文件索引.md';value=p.read_text(encoding='utf-8')
    value=value.replace('最新验收1216；16轮、32组','最新验收1218；18轮、36组')
    value=value.replace('_shared/'+PREV.name,'_shared/'+OUT.name)
    titles={1217:'末事件校准与有限冷却资源',1218:'同一载体跨段安全与联合校准'}
    rows='\n'.join(f'|{n}|[{titles[n]}](research_note_{n}.md)|[证明]({n}/proof.md)、[代码]({n}/check.py)、[结果]({n}/results.json)、[独审]({n}/independent_review.md)、[验收]({n}/acceptance.json)|' for n in ROUNDS)
    value=value.replace('\n## 管理与接续','\n'+rows+'\n\n## 管理与接续',1)
    value=re.sub(r'(\|[^\n]*\|)\n\n(?=\|[^\n]*\|)',r'\1\n',value)
    value+='\n[1216后共同范围审计](_shared/after1216_four_laws_scope_screen.md)：参考、工作、正温、循环分别核验；逐份不可达不以统一资源增长为前提。该审计不另计科研轮次。\n'
    p.write_text(value,encoding='utf-8')
    oldledger=read(OLD/'_shared/dependency_ledger.json')
    oldledger['latest_management_receipt']='../../archive_1201_/_shared/'+OUT.name
    write(OLD/'_shared/dependency_ledger.json',oldledger)
    p=OLD/'_shared/navigation_notice.md';prefix,history=p.read_text(encoding='utf-8').split('以下保留前阶段的历史状态',1)
    prefix=re.sub(r'\*\*当前状态：\*\*.*?(?=\n\n)','**当前状态：** '+STATUS+
        ' [1217报告](archive_1201_/research_note_1217.md) · [1218报告](archive_1201_/research_note_1218.md) · [当前核验](archive_1201_/_shared/'+OUT.name+')。旧目标均未追溯标记完成。',prefix,flags=re.S)
    p.write_text(prefix+'以下保留前阶段的历史状态'+history,encoding='utf-8')
    navapply=run(OLD/'_shared/update_navigation.py','--apply')
    navcheck=run(OLD/'_shared/update_navigation.py','--check');assert navcheck['changed_files']==0
    navverify=run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation');assert navverify['status']=='PASS'
    write(OUT,dict(status='PUBLISHING'))
    docs=[STAGE/'README.md',STAGE/'文件索引.md',AUDIT]
    for n in ROUNDS:docs += [STAGE/f'research_note_{n}.md',*(STAGE/str(n)).glob('*.md')]
    docs += [BASE/name for name in ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md']]
    nlinks=links(docs)
    for p in docs:assert not any(ord(c)<32 and c not in '\n\r\t' for c in p.read_text(encoding='utf-8'))
    management={name:sha(ROOT/name) for name in prev['management_sha256']}
    management[Path(__file__).relative_to(ROOT).as_posix()]=sha(Path(__file__))
    science=dict(prev['science_files_sha256'])
    for n in ROUNDS:
        science.update({(STAGE/str(n)/name).resolve().relative_to(ROOT).as_posix():digest for name,digest in files[n].items()})
    audits={**prev['audit_files_sha256'],AUDIT.relative_to(ROOT).as_posix():sha(AUDIT)}
    write(OUT,dict(schema='thermodynamics_progress_1218_v1',date='2026-10-15',status='PASS',
        accepted_rounds=list(range(1201,1219)),stage_check_groups=36,goal_changed=False,
        application_goal_status='active',stage_complete=False,goal_objective_sha256=prev['goal_objective_sha256'],
        acceptance_sha256={**prev['acceptance_sha256'],**{str(n):sha(STAGE/str(n)/'acceptance.json') for n in ROUNDS}},
        science_files_sha256=science,management_sha256=management,audit_files_sha256=audits,
        screen_files_sha256=prev['screen_files_sha256'],prior_progress_path=PREV.relative_to(ROOT).as_posix(),prior_progress_sha256=sha(PREV),
        old_GR_science_files_verified=old,old_GR_acceptance_bindings_verified=bindings,
        navigation_apply=navapply,navigation_check=navcheck,navigation_verify=navverify,
        local_links_checked=nlinks,missing_links=[],stage_requirement_status=req,
        new_adopted_axioms=0,new_empirical_calibration=0,GR_hard_deliverables_completed=0,
        science_program_outputs=outputs,prior_goal_turn_classification='progress: 1215 and 1216 accepted and verified',
        this_goal_turn_result='progress: event-calibrated resource bounds and whole-process shared-carrier qualification; pointwise/uniform third-law requirements separated'))
    verify()

if __name__=='__main__':main()
