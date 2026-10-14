"""Publish scoped round1202, after validating previous frozen evidence."""
from pathlib import Path
import hashlib,importlib.util,json,re,subprocess,sys
ROOT=Path(__file__).resolve().parents[3];BASE=ROOT/'research_cognition_physics'
STAGE=BASE/'archive_1201_';HERE=STAGE/'_shared';ROUND=STAGE/'1202';OLD=BASE/'archive_1161_'
s=importlib.util.spec_from_file_location('publish1201_tools',HERE/'publish_1201.py')
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
sha,read,write,run,links,old_check,validate_map=m.sha,m.read,m.write,m.run,m.links,m.old_check,m.validate_map
OUT=HERE/'progress_checks_1202.json'
def verify():
    d=read(OUT);old,count=old_check(read(HERE/'startup_checks.json'))
    validate_map(d['science_files_sha256']);validate_map(d['management_sha256']);validate_map(d['screen_files_sha256'])
    for number,h in d['acceptance_sha256'].items():
        p=STAGE/number/'acceptance.json';assert sha(p)==h
        for name,digest in read(p)['files_sha256'].items():assert sha((p.parent/name).resolve())==digest
    assert run(OLD/'_shared/update_navigation.py','--check')['changed_files']==0
    assert run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation')['status']=='PASS'
    print(json.dumps({'status':'PASS','latest_round':1202,'new_stage_rounds':2,'stage_check_groups':4,
       'prior_GR_science_files_verified':len(old['science_files_sha256']),'prior_GR_bindings':count,
       'current_stage_science_files':len(d['science_files_sha256']),'goal_status':'active','stage_complete':False},ensure_ascii=False))
def main():
    if '--verify' in sys.argv:return verify()
    assert not OUT.exists() and not (ROUND/'acceptance.json').exists()
    prev=read(HERE/'progress_checks_1201.json');previous_hash=sha(HERE/'progress_checks_1201.json')
    validate_map(prev['science_files_sha256']);validate_map(prev['management_sha256']);validate_map(prev['screen_files_sha256'])
    assert sha(STAGE/'1201/acceptance.json')==prev['acceptance_sha256']
    old,count=old_check(read(HERE/'startup_checks.json'))
    assert '最终限定签核' in (ROUND/'resource_review.md').read_text(encoding='utf-8')
    a=run(ROUND/'check.py');b=run(ROUND/'independent_check.py')
    assert a['status']=='PASS_CONDITIONAL_PREPARATION_AND_USE' and b['all_checks_passed']
    r=read(ROUND/'results.json')
    assert r['preparation']['dimension']==b['exact_counting']['shell_dimension']
    for i in range(4):
        assert abs(r['one_marginal'][i]-b['source_marginals']['one'][i])<1e-14
        for j in range(4):
            assert abs(r['two_marginal'][i][j]-b['source_marginals']['two'][i][j])<1e-14
            assert abs(r['native_contacts']['joint_outputs'][i][j]-b['original_1201_connection']['joint_600_300'][i][j])<3e-12
    required=['check.py','results.json','proof.md','independent_check.py','independent_results.json',
       'independent_review.md','resource_review.md','NEXT.md','../research_note_1202.md','../_shared/goal.md']
    hashes={name:sha((ROUND/name).resolve()) for name in required}
    acceptance={'schema':'operational_thermodynamics_round_acceptance_v1','round':1202,'date':'2026-10-14',
       'status':'ACCEPTED_WITH_SCOPE','check_groups':2,
       'scope':'固定100材料总能窗、1730位独立种子、保记录unrank及基态目标/已知能级库存的受控交换合同，成功参考严格平稳且关联保留；接原QQ的单/双源概率、完整平均能源及有限Clausius缺陷。未签原QQ实现SWAP、自主控制、精确Gibbs或几何身份。',
       'certified':{'N':100,'total_bare_energy_window':['39/2','41/2'],
          'exact_shell_dimension':r['preparation']['dimension'],'seed_bits':1730,
          'all_rejected_probability_upper':'2.11e-7','conditional_success_uniformity':True,
          'target_initial_level':'e0, distinct from probe A center c',
          'material_counts':{'ground_targets':100,'known_energy_inventory':400},
          'inventory_supplies_target_energy_exact':True,'preparation_switch_net_work':0,
          'temporary_control_energy_and_clock_not_free':True,
          'effective_beta_comparison':'6/5, not an input to the preparation rule',
          'single_reference_trace_error_upper':'.01001','two_reference_initial_relative_entropy_upper':'.000912',
          'two_reference_clausius_lower':'-.000912, initially independent probes and kept passive source records',
          'native_joint_probability_entry_numerical_error_upper':'3e-10',
          'native_correlation_prediction_TV_lower':'.00069502',
          'correlation_experimental_resolution_certified':False},
       'inputs_retained':r['remaining_inputs'],'new_adopted_axioms':0,'new_empirical_calibration':0,
       'GR_hard_deliverables_completed':0,'stage_complete':False,
       'limits':r['not_claimed'],
       'review':'Separate integer/rank and native-contact reconstruction, independent relative entropy/resource audit; root reran both. Candidate reviews retained as original review history; this acceptance gives adopted scope.',
       'root_executed_both_codes':True,'files_sha256':hashes}
    write(ROUND/'acceptance.json',acceptance)
    status='最新验收1202；新阶段2轮、4组。固定能源条件已给保失败/库存/记录的有限参考准备，接原QQ并认证两源熵—能源缺陷小于.000912。均匀种子、已知能级库存、受控交换与程序仍为公开条件；下一步核实际内部钟和原接触的共同能源账。'
    req={'conditional_identification_and_prediction':'completed_for_declared_model',
        'independently_generated_stable_reference':'conditional_finite_preparation_constructed_in_1202; autonomous control inputs remain',
        'internal_control_clock_record_preparation_resource_closure':'incomplete; inventory and seed/record account explicit, autonomous timing/program incomplete',
        'actual_geometry_identity':'incomplete; reviewed old clock/region objects and identified mismatch',
        'stage_deliverables_all_completed':False}
    ledger=read(HERE/'dependency_ledger.json')
    ledger.update(latest_accepted_round=1202,accepted_stage_rounds=[1201,1202],stage_check_groups=4,
        current_direction_document='../1202/NEXT.md',current_work=status,latest_management_receipt='progress_checks_1202.json',
        stage_requirement_status=req)
    write(HERE/'dependency_ledger.json',ledger)
    text=(STAGE/'README.md').read_text(encoding='utf-8')
    paragraphs=text.split('\n\n');paragraphs[1]=status
    text='\n\n'.join(paragraphs).replace('1201只有条件连接，尚未完成全部交付','1201—1202完成条件连接，内部控制和区域身份仍未闭合')
    text=text.replace('## 接续与边界','- [1202报告](research_note_1202.md)：独立能源条件准备、库存反作用、保关联接触与有限热力学界；[验收](1202/acceptance.json)、[下一项](1202/NEXT.md)。\n\n## 接续与边界',1)
    text=text.replace('_shared/progress_checks_1201.json','_shared/progress_checks_1202.json')
    text+='\n[实际钟与区域复用预审](_shared/clock_geometry_reuse_screen.md)：旧材料实现不自动运输给本轮h与QQ。\n'
    (STAGE/'README.md').write_text(text,encoding='utf-8')
    index=(STAGE/'文件索引.md').read_text(encoding='utf-8')
    index=index.replace('最新验收1201；1轮、2组','最新验收1202；2轮、4组')
    index=index.replace('## 管理与接续','|1202|[独立能量条件准备与保关联的有限热力学](research_note_1202.md)|[证明](1202/proof.md)、[代码](1202/check.py)、[结果](1202/results.json)、[独立核验](1202/independent_check.py)、[资源审阅](1202/resource_review.md)、[验收](1202/acceptance.json)|\n\n## 管理与接续',1)
    index=index.replace('[下一项](1201/NEXT.md)','[下一项](1202/NEXT.md)').replace('_shared/progress_checks_1201.json','_shared/progress_checks_1202.json')
    index+='\n- [实际钟与区域复用预审](_shared/clock_geometry_reuse_screen.md)\n'
    (STAGE/'文件索引.md').write_text(index,encoding='utf-8')
    oldledger=read(OLD/'_shared/dependency_ledger.json')
    oldledger['latest_management_receipt']='../../archive_1201_/_shared/progress_checks_1202.json'
    write(OLD/'_shared/dependency_ledger.json',oldledger)
    notice=OLD/'_shared/navigation_notice.md';text=notice.read_text(encoding='utf-8')
    prefix,rest=text.split('以下保留前阶段的历史状态',1)
    prefix=re.sub(r'\*\*当前状态：\*\*.*?(?=\n\n)','**当前状态：** '+status+
       ' [1202报告](archive_1201_/research_note_1202.md) · [验收](archive_1201_/1202/acceptance.json) · [当前核验](archive_1201_/_shared/progress_checks_1202.json)。旧GR阶段40轮、130组与全部冻结稿保留，阶段整体尚未结项。',prefix,flags=re.S)
    notice.write_text(prefix+'以下保留前阶段的历史状态'+rest,encoding='utf-8')
    navapply=run(OLD/'_shared/update_navigation.py','--apply')
    navcheck=run(OLD/'_shared/update_navigation.py','--check');assert navcheck['changed_files']==0
    navverify=run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation');assert navverify['status']=='PASS'
    write(OUT,{'status':'PUBLISHING'})
    paths=[STAGE/'README.md',STAGE/'文件索引.md',STAGE/'research_note_1202.md',
       *ROUND.glob('*.md'),HERE/'clock_geometry_reuse_screen.md',
       *[BASE/n for n in ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md']]]
    nlinks=links(paths)
    science=dict(prev['science_files_sha256'])
    science.update({(ROUND/n).resolve().relative_to(ROOT).as_posix():h for n,h in hashes.items()})
    management={p:sha(ROOT/p) for p in prev['management_sha256']}
    management[Path(__file__).relative_to(ROOT).as_posix()]=sha(Path(__file__))
    screens=dict(prev['screen_files_sha256'])
    p=HERE/'clock_geometry_reuse_screen.md';screens[p.relative_to(ROOT).as_posix()]=sha(p)
    receipt={'schema':'thermodynamics_progress_1202_v1','date':'2026-10-14','status':'PASS',
       'accepted_rounds':[1201,1202],'stage_check_groups':4,'application_goal_status':'active','stage_complete':False,
       'goal_objective_sha256':ledger['goal_objective_sha256'],
       'acceptance_sha256':{'1201':sha(STAGE/'1201/acceptance.json'),'1202':sha(ROUND/'acceptance.json')},
       'science_files_sha256':science,'management_sha256':management,'screen_files_sha256':screens,
       'prior_progress_path':(HERE/'progress_checks_1201.json').relative_to(ROOT).as_posix(),'prior_progress_sha256':previous_hash,
       'previous_GR_science_files_verified':len(old['science_files_sha256']),'previous_GR_acceptances_verified':40,
       'previous_GR_bindings':count,'previous_stage_round_preserved':1201,
       'navigation_apply':navapply,'navigation_check':navcheck,'navigation_verify':navverify,
       'markdown_files_checked':len(paths),'local_links_checked':nlinks,'missing_links':[],
       'stage_requirement_status':req,'new_adopted_axioms':0,'new_empirical_calibration':0,
       'GR_hard_deliverables_completed':0,
       'this_goal_turn_result':'progress: finite independent energy-window reference preparation, kept inventory/records, native correlated contacts and finite Clausius defect; old clock/geometry reuse audited; autonomous-contact admission is next'}
    write(OUT,receipt);verify()
if __name__=='__main__':main()
