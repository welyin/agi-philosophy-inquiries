"""Accept two scoped rounds; preserve the active four-law objective and history."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util
import json
import re
import sys

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'research_cognition_physics'; STAGE=BASE/'archive_1201_'
HERE=STAGE/'_shared'; OLD=BASE/'archive_1161_'
PREV=HERE/'progress_checks_1214.json'; OUT=HERE/'progress_checks_1216.json'
spec=importlib.util.spec_from_file_location('frozen_tools',HERE/'publish_1206.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
sha,read,write,run,links=m.sha,m.read,m.write,m.run,m.links
ROUNDS=[1215,1216]
AUDIT=HERE/'after1214_positive_temperature_admission.md'
STATUS=('最新验收1216；本阶段16轮、32组。1215以有限安全余额窗的完整振幅一致性，'
        '给出1214共同保幺角色的不含熵定义的充分实现，并核相位地址与纯记录反馈的删除。'
        '1216以同一原自主接触的两时刻严格证书，排除将目标A直接视为无记忆CP传播。'
        '工作角色未采用为公理，正温参考、实际权限及第三律资源增长仍保留条件；四律尚未全部完成。')

def verify():
    d=read(OUT);old,bindings=m.frozen(d);m.m.validate_map(d['management_sha256'])
    assert sha(PREV)==d['prior_progress_sha256']
    ledger=read(HERE/'dependency_ledger.json')
    assert ledger['latest_accepted_round']==1216 and ledger['stage_check_groups']==32
    assert ledger['latest_management_receipt']==OUT.name
    assert not ledger['stage_requirement_status']['stage_deliverables_all_completed']
    assert ledger['goal_objective_sha256']==d['goal_objective_sha256']
    assert run(OLD/'_shared/update_navigation.py','--check')['changed_files']==0
    assert run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation')['status']=='PASS'
    print(json.dumps(dict(status='PASS',latest_round=1216,stage_rounds=16,stage_check_groups=32,
        old_GR_science_files_verified=old,old_GR_acceptance_bindings_verified=bindings,
        current_stage_science_files_verified=len(d['science_files_sha256']),
        goal_changed=False,application_goal_status='active',stage_complete=False)))

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
    ir=read(STAGE/'1215/independent_results.json')
    assert ir['all_checks_passed'] and ir['check_groups']==7
    assert ir['groups']['probability_only_contract_fails']['normalized_unital_defect_exact']=='1/4'
    assert ir['groups']['balance_address_reset']['normalized_unital_defect_exact']=='1/2'
    assert ir['groups']['QND_record_then_feedback_scope']['complete_SK_channel_unital']
    assert ir['independent_code_sha256']==sha(STAGE/'1215/independent_check.py')
    ar15=read(STAGE/'1215/results.json')
    assert ar15['all_checks_passed'] and ar15['check_groups']==7
    assert ar15['code_sha256']==sha(STAGE/'1215/check.py')
    phase=next(g for g in ar15['groups'] if g['name']=='population calibration is insufficient')
    assert phase['normalized_unital_defect']=='3/8' and phase['late_energy_copy_defect']=='0'
    ar16=read(STAGE/'1216/results.json');g16=ar16['groups']
    assert ar16['all_checks_passed'] and ar16['check_groups']==8
    assert all(g['pass'] for g in g16.values())
    assert g16['fixed_original_process_and_history']['physical_times']==[301,601]
    assert not g16['fixed_original_process_and_history']['intermediate_environment_reset']
    assert F(g16['rational_normalized_Choi_witness']['exact_real_Rayleigh'])<F('-0.06')
    transport=g16['same_clock_transport_and_actual_non_CP']
    assert F(transport['normalized_Choi_actual_to_Y_bound'])<F('.018')
    assert F(transport['rational_actual_Rayleigh_upper'])<F(transport['certified_actual_Rayleigh_less_than'])==F('-0.042')
    assert not g16['scope_and_no_new_machine']['full_AB_non_CP_divisibility_claim']
    assert not g16['scope_and_no_new_machine']['third_law_or_thermodynamics_failure_claim']
    ir16=read(STAGE/'1216/independent_results.json')
    assert ir16['all_checks_passed'] and ir16['check_groups']==6
    ig16={g['name']:g for g in ir16['groups']}
    assert all(g['passed'] for g in ig16.values())
    assert ig16['independent_full_integer_reconstruction']['exact_ShThat_match']
    assert F(ig16['independent_integer_normalized_Choi_witness']['exact_Rayleigh'])==F(g16['rational_normalized_Choi_witness']['exact_real_Rayleigh'])
    it=ig16['same_clock_actual_invertibility_and_non_CP_transport']
    assert it['physical_times']==[301,601] and it['program_v']==1
    assert F(it['actual_Rayleigh_upper'])==F(transport['rational_actual_Rayleigh_upper'])
    assert F(it['normalized_Choi_perturbation'])==F(transport['normalized_Choi_actual_to_Y_bound'])
    assert not ir16['scope']['full_AB_non_CP_divisibility_claimed']
    drafts=STAGE/'1216/_drafts/pre_review_replay_fix'
    preserved={drafts/'check.py':'6cc66a9e5d5dc9293fb5fabf37c953fc30584777cf0715c0e2806ec9890780f9',
               drafts/'results.json':'916c62c373a8dd65c98c71c8a98e83cfe94c6329948cc172d45bc116e9d2c552'}
    for path,digest in preserved.items():assert sha(path)==digest
    certified={
      1215:dict(finite_positive_balance_window=True,complete_amplitude_including_phase_required=True,
                full_boundary_energy_conservation=True,terminal_QND_scope_explicit=True,
                feedback_purity_not_free=True,old_1209_witness_membership=True,
                sufficient_role_not_adopted_axiom=True,actual_native_controller_claimed=False),
      1216:dict(same_native_parent_and_preparation=True,physical_times=[301,601],
                finite_integer_negative_Choi_certificate=True,
                actual_normal_clock_error_included=True,target_A_CP_divisibility=False,
                full_AB_non_CP_divisibility_claimed=False,thermodynamics_failure_claimed=False,
                universal_finite_time_third_law_claimed=False)}
    scopes={
      1215:'有限安全余额窗完整振幅一致性与明确末QND记录给保幺材料角色；相位/余额地址和纯记录反馈删除说明概率表不足与完整材料域必要。角色为新充分合同，未采用，不授原QQ实现。',
      1216:'同一原v=1自主父H、同冷B准备，t=301/601的A约化通道唯一中间映射归一Choi有严格负方向，计完整钟误差；只否定目标A的CP可分接口，不否定完整热力学或全AB过程。'}
    for number in ROUNDS:
        write(STAGE/str(number)/'acceptance.json',dict(
            schema='operational_thermodynamics_round_acceptance_v1',round=number,date='2026-10-15',
            status='ACCEPTED_WITH_SCOPE',check_groups=2,scope=scopes[number],
            root_executed_both_codes=True,files_sha256=files[number],certified=certified[number],
            new_adopted_axioms=0,new_empirical_calibration=0,GR_hard_deliverables_completed=0,
            stage_complete=False))
    ledger=read(HERE/'dependency_ledger.json');req=ledger['stage_requirement_status']
    req.update(common_process='1215 supplies a finite balance-window sufficient implementation of the 1214 unital role; phase and memory scope explicit, actual native QQ permissions not inferred',
        first_law='D-P gives time/complete-energy identity; 1212 full energy ledger and 1215 finite energy-exchange work-role candidate distinguish energy accounting from purity supply',
        second_law='1214/1215 common-material Clausius role retained; 1212 cyclic necessary bounds remain, no automatic reusable engine or free record reset',
        third_law='Finite faithful-input bounds retained. 1216 excludes a direct positive-rate Markov shortcut for the actual target-A path; full resource/time growth and source conditions remain',
        stage_synthesis='Work-role sufficient mechanism and actual non-Markov qualification audited; positive-temperature normal-extension shortcut declined without new round',
        stage_deliverables_all_completed=False)
    ledger.update(latest_accepted_round=1216,accepted_stage_rounds=list(range(1201,1217)),
        stage_check_groups=32,current_direction_document='../1215/NEXT.md',
        current_work=STATUS,latest_management_receipt=OUT.name)
    write(HERE/'dependency_ledger.json',ledger)
    p=STAGE/'README.md';text=p.read_text(encoding='utf-8');parts=text.split('\n\n');parts[1]=STATUS
    text='\n\n'.join(parts).replace('_shared/'+PREV.name,'_shared/'+OUT.name)
    rows='\n'.join(f'- [{n}报告](research_note_{n}.md)：{scopes[n]} [验收]({n}/acceptance.json)、[下一项]({n}/NEXT.md)。' for n in ROUNDS)
    text=text.replace('## 接续与边界',rows+'\n\n## 接续与边界',1)
    p.write_text(text,encoding='utf-8')
    p=STAGE/'文件索引.md';text=p.read_text(encoding='utf-8')
    text=text.replace('最新验收1214；14轮、28组','最新验收1216；16轮、32组')
    text=text.replace('_shared/'+PREV.name,'_shared/'+OUT.name)
    titles={1215:'有限余额完整转换与工作角色',1216:'原实际接触的非CP可分性'}
    rows='\n'.join(f'|{n}|[{titles[n]}](research_note_{n}.md)|[证明]({n}/proof.md)、[代码]({n}/check.py)、[结果]({n}/results.json)、[独审]({n}/independent_review.md)、[验收]({n}/acceptance.json)|' for n in ROUNDS)
    text=text.replace('\n## 管理与接续','\n'+rows+'\n\n## 管理与接续',1)
    text=re.sub(r'(\|[^\n]*\|)\n\n(?=\|[^\n]*\|)',r'\1\n',text)
    text+='\n[1214后正温准入审计](_shared/after1214_positive_temperature_admission.md)：工作钟不是平衡参考；有限窗正常性不能代替统一紧性。本审计不另计科研轮次。\n'
    p.write_text(text,encoding='utf-8')
    oldledger=read(OLD/'_shared/dependency_ledger.json')
    oldledger['latest_management_receipt']='../../archive_1201_/_shared/'+OUT.name
    write(OLD/'_shared/dependency_ledger.json',oldledger)
    p=OLD/'_shared/navigation_notice.md';prefix,history=p.read_text(encoding='utf-8').split('以下保留前阶段的历史状态',1)
    prefix=re.sub(r'\*\*当前状态：\*\*.*?(?=\n\n)','**当前状态：** '+STATUS+
        ' [1215报告](archive_1201_/research_note_1215.md) · [1216报告](archive_1201_/research_note_1216.md) · [当前核验](archive_1201_/_shared/'+OUT.name+')。旧目标均未追溯标记完成。',prefix,flags=re.S)
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
    audits.update({path.relative_to(ROOT).as_posix():digest for path,digest in preserved.items()})
    write(OUT,dict(schema='thermodynamics_progress_1216_v1',date='2026-10-15',status='PASS',
        accepted_rounds=list(range(1201,1217)),stage_check_groups=32,goal_changed=False,
        application_goal_status='active',stage_complete=False,goal_objective_sha256=prev['goal_objective_sha256'],
        acceptance_sha256={**prev['acceptance_sha256'],**{str(n):sha(STAGE/str(n)/'acceptance.json') for n in ROUNDS}},
        science_files_sha256=science,management_sha256=management,audit_files_sha256=audits,
        screen_files_sha256=prev['screen_files_sha256'],
        prior_progress_path=PREV.relative_to(ROOT).as_posix(),prior_progress_sha256=sha(PREV),
        old_GR_science_files_verified=old,old_GR_acceptance_bindings_verified=bindings,
        navigation_apply=navapply,navigation_check=navcheck,navigation_verify=navverify,
        local_links_checked=nlinks,missing_links=[],stage_requirement_status=req,
        new_adopted_axioms=0,new_empirical_calibration=0,GR_hard_deliverables_completed=0,
        science_program_outputs=outputs,prior_goal_turn_classification='progress: 1214 accepted and verified',
        this_goal_turn_result='progress: finite work-role sufficient mechanism and actual non-CP-divisibility evidence; no extra axioms adopted'))
    verify()

if __name__=='__main__':main()
