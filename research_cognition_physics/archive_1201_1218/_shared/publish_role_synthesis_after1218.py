"""Publish the reviewed role/theorem proposal without adding a scientific round."""
from pathlib import Path
import importlib.util
import json
import re
import sys

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'research_cognition_physics';STAGE=BASE/'archive_1201_'
HERE=STAGE/'_shared';OLD=BASE/'archive_1161_'
PREV=HERE/'progress_checks_1218.json'
OUT=HERE/'role_synthesis_checks_after1218.json'
spec=importlib.util.spec_from_file_location('frozen_tools',HERE/'publish_1206.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
sha,read,write,run,links=m.sha,m.read,m.write,m.run,m.links
DOCS=['four_laws_joint_theorem_after1218.md',
      'cognitive_thermal_roles_proposal_after1218.md',
      'after1218_role_scope_review.md',
      'after1218_common_realization_review.md',
      'after1218_joint_theorem_independent_review.md']
REPLAY_FILES=['after1218_common_prediction_check.py','after1218_common_prediction_results.json']
STATUS=('最新科学验收1218；本阶段18轮、36组。四律联合条件定理和热学角色提案已完成采用前审阅，'
        '分开数学角色、共同非空范围与原实际菜单资格。完整第二律不要求先把材料通道保幺；'
        '第三律逐协议不可达与统一限时界分别陈述。提案未采用为公理，待明确热学角色的理论定位后形成正式论文；'
        '不新增科学轮次，不宣布四律或GR完成。')

def verify():
    d=read(OUT);old,bindings=m.frozen(d);m.m.validate_map(d['management_sha256'])
    assert sha(PREV)==d['prior_progress_sha256']
    ledger=read(HERE/'dependency_ledger.json')
    assert ledger['latest_management_receipt']==OUT.name
    assert ledger['latest_accepted_round']==1218 and ledger['stage_check_groups']==36
    assert ledger['new_adopted_axioms']==[]
    assert not ledger['stage_requirement_status']['stage_deliverables_all_completed']
    assert run(OLD/'_shared/update_navigation.py','--check')['changed_files']==0
    assert run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation')['status']=='PASS'
    print(json.dumps(dict(status='PASS',latest_scientific_round=1218,stage_rounds=18,
        stage_check_groups=36,old_GR_science_files_verified=old,
        old_GR_acceptance_bindings_verified=bindings,
        current_stage_science_files_verified=len(d['science_files_sha256']),
        new_scientific_rounds=0,new_adopted_axioms=0,goal_changed=False,stage_complete=False)))

def main():
    if '--verify' in sys.argv:return verify()
    assert not OUT.exists(),'Already published'
    prev=read(PREV);old,bindings=m.frozen(prev);m.m.validate_map(prev['management_sha256'])
    for name in DOCS:assert (HERE/name).is_file(),name
    review=(HERE/DOCS[-1]).read_text(encoding='utf-8')
    assert '限定审阅通过' in review
    for name in DOCS[:2]:assert sha(HERE/name) in review,name
    for name in REPLAY_FILES:assert sha(HERE/name) in review,name
    replay=run(HERE/REPLAY_FILES[0])
    assert replay['status']=='PASS' and replay['all_checks_passed'] and replay['results_match']
    prediction=read(HERE/REPLAY_FILES[1])
    for name,digest in prediction['source_sha256'].items():assert sha(BASE/name)==digest,name
    assert prediction['groups']['native_factory_event_transport']['actual_1205_post_preparation_channel_unital_claimed'] is False
    assert prediction['groups']['native_factory_event_transport']['1203_parent_or_error_substituted_for_1205'] is False
    for name in DOCS:
        v=(HERE/name).read_text(encoding='utf-8')
        assert not any(x in v for x in (r'\[',r'\]',r'\(',r'\)')),name
        assert not re.search(r'^\s*\$\s*$',v,flags=re.M),name
        assert not any(ord(c)<32 and c not in '\n\r\t' for c in v),name
    ledger=read(HERE/'dependency_ledger.json')
    req=ledger['stage_requirement_status']
    req.update(stage_synthesis='Reviewed joint conditional theorem and cognitive thermal-role proposal. Definition, mathematical nonempty class, and actual operational qualification separated. Roles not adopted; final paper and requirement-by-requirement closure remain.',
        second_law='Complete unitary entropy ledger with references, correlations, machines and records does not need material unitality; unital work is a sufficient role interface only. Clausius-to-Kelvin conditional relation retained; reverse equivalence not claimed.',
        common_process='Existing finite original process qualifies with explicit non-Gibbs defects and complete preparation/clock/record budgets; its pure-stock whole preparation is not a faithful cold-resource input or a cyclic engine.',
        stage_deliverables_all_completed=False)
    ledger.update(current_direction_document='cognitive_thermal_roles_proposal_after1218.md',
        current_scope_audit='four_laws_joint_theorem_after1218.md',current_work=STATUS,
        latest_management_receipt=OUT.name)
    write(HERE/'dependency_ledger.json',ledger)
    oldledger=read(OLD/'_shared/dependency_ledger.json')
    oldledger['latest_management_receipt']='../../archive_1201_/_shared/'+OUT.name
    write(OLD/'_shared/dependency_ledger.json',oldledger)
    for name in ['README.md','文件索引.md']:
        p=STAGE/name;v=p.read_text(encoding='utf-8')
        if name=='README.md':
            parts=v.split('\n\n');parts[1]=STATUS;v='\n\n'.join(parts)
        v=v.replace('_shared/'+PREV.name,'_shared/'+OUT.name)
        v+='\n## 1218后共同理论与采用前审查\n\n'
        v+='- [四律联合条件定理](_shared/four_laws_joint_theorem_after1218.md)\n'
        v+='- [认知热学角色提案](_shared/cognitive_thermal_roles_proposal_after1218.md)\n'
        v+='- [角色范围审计](_shared/after1218_role_scope_review.md)、[共同实现审计](_shared/after1218_common_realization_review.md)、[独立综合审阅](_shared/after1218_joint_theorem_independent_review.md)\n'
        v+='\n本次为已有证明的综合与采用前审查，不新增研究编号或公理；目标继续有效，阶段未完成。\n'
        p.write_text(v,encoding='utf-8')
    p=OLD/'_shared/navigation_notice.md'
    prefix,history=p.read_text(encoding='utf-8').split('以下保留前阶段的历史状态',1)
    prefix=re.sub(r'\*\*当前状态：\*\*.*?(?=\n\n)','**当前状态：** '+STATUS+
        ' [共同定理](archive_1201_/_shared/four_laws_joint_theorem_after1218.md) · '
        '[角色提案](archive_1201_/_shared/cognitive_thermal_roles_proposal_after1218.md) · '
        '[当前核验](archive_1201_/_shared/'+OUT.name+')。旧目标均未追溯标记完成。',prefix,flags=re.S)
    p.write_text(prefix+'以下保留前阶段的历史状态'+history,encoding='utf-8')
    navapply=run(OLD/'_shared/update_navigation.py','--apply')
    navcheck=run(OLD/'_shared/update_navigation.py','--check');assert navcheck['changed_files']==0
    navverify=run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation');assert navverify['status']=='PASS'
    write(OUT,dict(status='PUBLISHING'))
    docs=[HERE/name for name in DOCS]+[STAGE/'README.md',STAGE/'文件索引.md']
    docs += [BASE/name for name in ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md']]
    nlinks=links(docs)
    management={name:sha(ROOT/name) for name in prev['management_sha256']}
    management[Path(__file__).relative_to(ROOT).as_posix()]=sha(Path(__file__))
    audits={**prev['audit_files_sha256'],**{(HERE/name).relative_to(ROOT).as_posix():sha(HERE/name) for name in DOCS+REPLAY_FILES}}
    write(OUT,dict(schema='four_laws_role_synthesis_after1218_v1',date='2026-10-15',status='PASS',
        accepted_rounds=prev['accepted_rounds'],stage_check_groups=36,goal_changed=False,
        application_goal_status='active',stage_complete=False,new_scientific_rounds=0,
        goal_objective_sha256=prev['goal_objective_sha256'],
        acceptance_sha256=prev['acceptance_sha256'],science_files_sha256=prev['science_files_sha256'],
        management_sha256=management,audit_files_sha256=audits,screen_files_sha256=prev['screen_files_sha256'],
        prior_progress_path=PREV.relative_to(ROOT).as_posix(),prior_progress_sha256=sha(PREV),
        old_GR_science_files_verified=old,old_GR_acceptance_bindings_verified=bindings,
        navigation_apply=navapply,navigation_check=navcheck,navigation_verify=navverify,
        local_links_checked=nlinks,missing_links=[],stage_requirement_status=req,
        new_adopted_axioms=0,new_empirical_calibration=0,GR_hard_deliverables_completed=0,
        supporting_existing_certificate_replay=replay,
        prior_goal_turn_classification='progress: accepted and verified 1217 and 1218',
        this_goal_turn_result='progress: reviewed common four-law theorem and concrete cognitive thermal-role proposal, removing inappropriate universal unitality/hardware requirements while preserving actual qualification gaps'))
    verify()

if __name__=='__main__':main()
