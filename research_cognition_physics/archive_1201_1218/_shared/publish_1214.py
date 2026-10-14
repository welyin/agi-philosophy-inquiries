"""Accept the common operational-class bridge; leave the four-law goal active."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util
import json
import re
import sys

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'research_cognition_physics'; STAGE=BASE/'archive_1201_'
HERE=STAGE/'_shared'; OLD=BASE/'archive_1161_'; ROUND=STAGE/'1214'
PREV=HERE/'progress_checks_1213.json'; OUT=HERE/'progress_checks_1214.json'
spec=importlib.util.spec_from_file_location('frozen_tools',HERE/'publish_1206.py')
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
sha,read,write,run,links=m.sha,m.read,m.write,m.run,m.links
REVIEWED=['proof.md','check.py','results.json','independent_check.py',
          'independent_results.json','NEXT.md','../research_note_1214.md']


def verify():
    d=read(OUT); old,bindings=m.frozen(d); m.m.validate_map(d['management_sha256'])
    assert sha(PREV)==d['prior_progress_sha256']
    ledger=read(HERE/'dependency_ledger.json')
    assert ledger['latest_accepted_round']==1214 and ledger['stage_check_groups']==28
    assert ledger['latest_management_receipt']==OUT.name
    assert not ledger['stage_requirement_status']['stage_deliverables_all_completed']
    assert ledger['goal_objective_sha256']==d['goal_objective_sha256']
    assert run(OLD/'_shared/update_navigation.py','--check')['changed_files']==0
    assert run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation')['status']=='PASS'
    print(json.dumps(dict(status='PASS',latest_round=1214,stage_rounds=14,stage_check_groups=28,
        old_GR_science_files_verified=old,old_GR_acceptance_bindings_verified=bindings,
        current_stage_science_files_verified=len(d['science_files_sha256']),
        goal_changed=False,application_goal_status='active',stage_complete=False)))


def main():
    if '--verify' in sys.argv:return verify()
    assert not OUT.exists() and not (ROUND/'acceptance.json').exists()
    prev=read(PREV); old,bindings=m.frozen(prev); m.m.validate_map(prev['management_sha256'])
    review=(ROUND/'independent_review.md').read_text(encoding='utf-8')
    assert '限定签核通过' in review
    for name in REVIEWED:assert sha((ROUND/name).resolve()) in review,name
    for name in ['proof.md','NEXT.md','../research_note_1214.md']:
        text=(ROUND/name).read_text(encoding='utf-8')
        assert not any(x in text for x in (r'\[',r'\]',r'\(',r'\)'))
    author=run(ROUND/'check.py'); independent=run(ROUND/'independent_check.py')
    for obj in (author,independent):
        assert obj['status']=='PASS' and obj['all_checks_passed'] and obj['results_match']
    for obj in (read(ROUND/'results.json'),read(ROUND/'independent_results.json')):
        for name,digest in obj.get('source_sha256',{}).items():
            p=next((p for p in (STAGE/name,BASE/name,ROOT/name) if p.is_file()),None)
            assert p and sha(p)==digest,name
    ir=read(ROUND/'independent_results.json')
    assert ir['all_checks_passed'] and ir['independent_code_sha256']==sha(ROUND/'independent_check.py')
    ar=read(ROUND/'results.json'); ag={g['name']:g for g in ar['groups']}
    assert ar['all_checks_passed'] and ar['check_groups']==7 and ir['check_groups']==6
    assert all(g['passed'] for g in ar['groups'])
    native=ag['exact_native_population_and_clock_calibration']
    other=ir['groups']['native_independent_contact_bridge']
    assert F(native['joint_lambda_exact'])==F(other['lambda_exact'])
    for a,b in zip(native['bounds'],other['rows']):
        assert a['QQ_time']==b['QQ_duration']
        assert F(a['half_diamond_exact_upper'])==F(b['epsilon_exact'])
        assert F(a['strong_lower_using_exact_clock_formula'])==F(b['actual_failure_lower_exact'])
        assert F(a['strong_failure_exact_lower_using_outward_error'])>F(a['strict_decimal_lower'])>F('0.26933')
    assert ir['groups']['first_hit_lift_with_complete_pointer']['TP_and_unital_integer_Gram_exact']
    assert ir['groups']['first_hit_lift_with_complete_pointer']['unknown_SR_Gram_channel_residual']<1e-12
    assert not ar['scope']['new_adopted_axioms']
    assert not ar['scope']['general_finite_time_third_law_claimed']
    assert not ar['scope']['actual_1205_prepared_correlated_cut_has_fixed_AB_channel_claimed']
    names=REVIEWED+['independent_review.md','../_shared/goal_four_laws.md']
    hashes={name:sha((ROUND/name).resolve()) for name in names}
    scope=('在未采用的完整材料保幺角色子类中，1209有限能梯与末copy仍保留非Gibbs反向选择见证；'
           '同类连接共同温标、完整材料/热源Clausius与忠实输入冷却界，并运输到原独立正常钟接触。'
           '不将目标子系统单独要求保幺，不覆盖相关工厂准备，不由有限证书推出一般时间增长律。')
    write(ROUND/'acceptance.json',dict(schema='operational_thermodynamics_round_acceptance_v1',
        round=1214,date='2026-10-15',status='ACCEPTED_WITH_SCOPE',check_groups=2,
        scope=scope,root_executed_both_codes=True,files_sha256=hashes,
        new_adopted_axioms=0,new_empirical_calibration=0,GR_hard_deliverables_completed=0,
        stage_complete=False,certified=dict(common_unital_material_role_candidate=True,
            native_1209_reverse_witness_retained=True,full_material_not_reduced_target_unital=True,
            finite_faithful_input_bound=True,independent_1204_contact_bridge=True,
            related_1205_factory_channel_automatically_covered=False,
            positive_temperature_selected_by_unitality=False,
            separate_auxiliary_budgets_automatically_compose=False,
            full_finite_time_third_law=False,role_adopted_as_cognitive_axiom=False)))
    status=('最新验收1214；本阶段14轮、28组。1212—1213已判清有限参考循环损失、实际机器未归还及弱能源条件的零温边界；'
            '1214在同一完整材料保幺角色子类中保留Gibbs选择，并连接温标、Clausius与忠实输入冷却界及原独立正常钟接触。'
            '该角色尚未采用为公理，相关工厂权限、资源随时间增长及四律共同结项仍未完成。')
    ledger=read(HERE/'dependency_ledger.json'); req=ledger['stage_requirement_status']
    req.update(common_process='1214: conditional common unital whole-material class retains the 1209 witness and connects Gibbs/equilibrium, Clausius and faithful cooling; actual independent 1204 contact has a controlled bridge, correlated factory kept separate',
        zeroth_law='1210/1214: conditional common Gibbs parameter and equilibrium role survives the unital sub-menu; no new reference axiom adopted',
        second_law='1212/1214: full cyclic-work necessary bounds and whole-material/bath entropy role; actual heat-engine/reversible-cycle permissions remain open',
        third_law='1211/1213/1214: source/role distinctions, weak-energy endpoint counterexample, and robust finite faithful-input bound; no general resource/time growth or universal zero-temperature statement',
        stage_synthesis='Common-scope bridge advanced; role adoption and four-law completion audit remain pending',
        stage_deliverables_all_completed=False)
    ledger.update(latest_accepted_round=1214,accepted_stage_rounds=list(range(1201,1215)),
        stage_check_groups=28,current_direction_document='../1214/NEXT.md',
        current_work=status,latest_management_receipt=OUT.name)
    write(HERE/'dependency_ledger.json',ledger)
    p=STAGE/'README.md'; value=p.read_text(encoding='utf-8'); parts=value.split('\n\n');parts[1]=status
    value='\n\n'.join(parts).replace('_shared/'+PREV.name,'_shared/'+OUT.name)
    value=value.replace('## 接续与边界','- [1214报告](research_note_1214.md)：共同材料角色与四律接口；[验收](1214/acceptance.json)、[下一项](1214/NEXT.md)。\n\n## 接续与边界',1)
    p.write_text(value,encoding='utf-8')
    p=STAGE/'文件索引.md';value=p.read_text(encoding='utf-8')
    value=value.replace('最新验收1213；13轮、26组','最新验收1214；14轮、28组')
    value=value.replace('_shared/'+PREV.name,'_shared/'+OUT.name)
    row='|1214|[共同材料角色与四律接口](research_note_1214.md)|[证明](1214/proof.md)、[代码](1214/check.py)、[结果](1214/results.json)、[独审](1214/independent_review.md)、[验收](1214/acceptance.json)|\n'
    value=value.replace('\n## 管理与接续','\n'+row+'\n## 管理与接续',1)
    value=re.sub(r'(\|[^\n]*\|)\n\n(?=\|[^\n]*\|)',r'\1\n',value)
    p.write_text(value,encoding='utf-8')
    oldledger=read(OLD/'_shared/dependency_ledger.json');oldledger['latest_management_receipt']='../../archive_1201_/_shared/'+OUT.name
    write(OLD/'_shared/dependency_ledger.json',oldledger)
    p=OLD/'_shared/navigation_notice.md';prefix,history=p.read_text(encoding='utf-8').split('以下保留前阶段的历史状态',1)
    prefix=re.sub(r'\*\*当前状态：\*\*.*?(?=\n\n)','**当前状态：** '+status+
        ' [1214报告](archive_1201_/research_note_1214.md) · [当前核验](archive_1201_/_shared/'+OUT.name+')。旧目标均未追溯标记完成。',prefix,flags=re.S)
    p.write_text(prefix+'以下保留前阶段的历史状态'+history,encoding='utf-8')
    navapply=run(OLD/'_shared/update_navigation.py','--apply')
    navcheck=run(OLD/'_shared/update_navigation.py','--check');assert navcheck['changed_files']==0
    navverify=run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation');assert navverify['status']=='PASS'
    write(OUT,dict(status='PUBLISHING'))
    docs=[STAGE/'README.md',STAGE/'文件索引.md',STAGE/'research_note_1214.md',*ROUND.glob('*.md')]
    docs += [BASE/name for name in ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md']]
    nlinks=links(docs)
    for p in docs:
        assert not any(ord(c)<32 and c not in '\n\r\t' for c in p.read_text(encoding='utf-8'))
    management={name:sha(ROOT/name) for name in prev['management_sha256']}
    management[Path(__file__).relative_to(ROOT).as_posix()]=sha(Path(__file__))
    science={**prev['science_files_sha256'],**{(ROUND/name).resolve().relative_to(ROOT).as_posix():digest for name,digest in hashes.items()}}
    write(OUT,dict(schema='thermodynamics_progress_1214_v1',date='2026-10-15',status='PASS',
        accepted_rounds=list(range(1201,1215)),stage_check_groups=28,goal_changed=False,
        application_goal_status='active',stage_complete=False,goal_objective_sha256=prev['goal_objective_sha256'],
        acceptance_sha256={**prev['acceptance_sha256'],'1214':sha(ROUND/'acceptance.json')},
        science_files_sha256=science,management_sha256=management,
        audit_files_sha256=prev['audit_files_sha256'],screen_files_sha256=prev['screen_files_sha256'],
        prior_progress_path=PREV.relative_to(ROOT).as_posix(),prior_progress_sha256=sha(PREV),
        old_GR_science_files_verified=old,old_GR_acceptance_bindings_verified=bindings,
        navigation_apply=navapply,navigation_check=navcheck,navigation_verify=navverify,
        local_links_checked=nlinks,missing_links=[],stage_requirement_status=req,
        new_adopted_axioms=0,new_empirical_calibration=0,GR_hard_deliverables_completed=0,
        science_program_outputs={'1214':dict(author=author,independent=independent)},
        prior_goal_turn_classification='progress: cyclic-work and endpoint resource boundaries accepted',
        this_goal_turn_result='progress: common operational-class theorem and actual finite-error bridge, no new axiom adopted'))
    verify()


if __name__=='__main__':main()
