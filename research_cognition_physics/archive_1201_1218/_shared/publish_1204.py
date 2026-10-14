"""Publish the reviewed two-reference flow round, keeping old scientific files."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util,json,re
ROOT=Path(__file__).resolve().parents[3];BASE=ROOT/'research_cognition_physics'
STAGE=BASE/'archive_1201_';HERE=STAGE/'_shared';ROUND=STAGE/'1204';OLD=BASE/'archive_1161_'
s=importlib.util.spec_from_file_location('publish1201_tools',HERE/'publish_1201.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
sha,read,write,run,links,old_check,validate_map=m.sha,m.read,m.write,m.run,m.links,m.old_check,m.validate_map
OUT=HERE/'progress_checks_1204.json'

def verify():
    d=read(OUT);old,n=old_check(read(HERE/'startup_checks.json'))
    for key in ['science_files_sha256','management_sha256','screen_files_sha256']:validate_map(d[key])
    for number,h in d['acceptance_sha256'].items():
        p=STAGE/number/'acceptance.json';assert sha(p)==h
        for name,digest in read(p)['files_sha256'].items():assert sha((p.parent/name).resolve())==digest
    assert run(OLD/'_shared/update_navigation.py','--check')['changed_files']==0
    assert run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation')['status']=='PASS'
    print(json.dumps({'status':'PASS','latest_round':1204,'new_stage_rounds':4,'stage_check_groups':8,
        'old_GR_science_files_verified':len(old['science_files_sha256']),'old_GR_bindings_verified':n,
        'current_stage_science_files':len(d['science_files_sha256']),'goal_status':'active','stage_complete':False}))

def main():
    import sys
    if '--verify' in sys.argv:return verify()
    assert not OUT.exists() and not (ROUND/'acceptance.json').exists()
    pp=HERE/'progress_checks_1203.json';prev=read(pp)
    for key in ['science_files_sha256','management_sha256','screen_files_sha256']:validate_map(prev[key])
    for number,h in prev['acceptance_sha256'].items():assert sha(STAGE/number/'acceptance.json')==h
    old,n=old_check(read(HERE/'startup_checks.json'))
    assert '限定签核' in (ROUND/'independent_review.md').read_text(encoding='utf-8')
    a=run(ROUND/'check.py');b=run(ROUND/'independent_check.py')
    assert a['status']=='AUTHOR_DRAFT_CHECKS_PASS_NOT_ACCEPTANCE' and b['all_checks_passed']
    r=read(ROUND/'results.json')
    assert str(b['hot_shell']['dimension'])==r['new_reference']['dimension']
    assert b['hot_shell']['types']==r['new_reference']['selected_types']
    assert [F(q) for q in b['hot_shell']['one_rational']]==[F(q) for q in r['new_reference']['one_marginal_rational']]
    for t in ['600','300']:
        for i in range(16):
            for j in range(16):assert abs(float(r['full_native_contacts'][t]['full_transition_P'][i][j])-b['contacts'][t]['full_transition_P'][i][j])<2e-11
    assert F(r['actual_flow_predictions']['t300_lambda1/2']['Q_dir_actual_interval'][0])>F('.01373963')
    assert F(r['multi_reference_thermodynamics']['cold_hot_initial_relative_entropy_defect_upper'])<F('.000502766')
    files=['proof.md','check.py','results.json','independent_check.py','independent_results.json',
           'independent_review.md','NEXT.md','../research_note_1204.md','../_shared/goal.md']
    hashes={name:sha((ROUND/name).resolve()) for name in files}
    scope='固定两独立100材料能窗、保失败/标签的λ三点准备和原自主QQ/能源指针，给实际定向流、均温对照、有限响应/耗竭及完整多参考熵账。对换A读的统计代理与同次实际Q推断分开；未授普遍温度序、自动准备程序、区域几何或GR。'
    acc={'schema':'operational_thermodynamics_round_acceptance_v1','round':1204,'date':'2026-10-14',
        'status':'ACCEPTED_WITH_SCOPE','check_groups':2,'scope':scope,
        'certified':{'cold_window':['39/2','41/2'],'hot_window':['79/2','81/2'],
            'hot_shell_dimension':r['new_reference']['dimension'],'hot_selected_types':2391,
            'hot_seed_bits':1950,'hot_failure_upper':r['new_reference']['all_rejected_upper'],
            'conditional_mixture_weight_accounts_for_unequal_failures':True,
            'hot_beta_comparison':'1/3, after independent preparation',
            'hot_reference_trace_error_upper':r['temperature_comparisons']['hot']['local_trace_distance_upper'],
            'initial_two_reference_relative_entropy_defect_upper':r['multi_reference_thermodynamics']['cold_hot_initial_relative_entropy_defect_upper'],
            'four_nonzero_label_branches_actual_direction_positive':True,
            'weakest_half_preparation_actual_Q_lower':'.01373963',
            'actual_swap_symmetry_of_readout_parent_claimed':False,
            'future_measurement_inference_uses_two_model_error_allowances':True,
            'all_possible_observations_guaranteed_positive_confidence_lower':False,
            'mixed_AB_clock_rank_bound':1024,'Y_clock_record_rank_bound':256,
            'clock_work_direction_for_600_certified':False},
        'remaining_inputs':r['remaining_inputs'],'limits':r['not_claimed'],
        'new_adopted_axioms':0,'new_empirical_calibration':0,'GR_hard_deliverables_completed':0,
        'stage_complete':False,'root_executed_both_codes':True,
        'review':'Independent integer and native matrix reconstruction, fixed-budget inference and retained-resource/entropy audit; root read proofs and ran both programs. No physical experiment.',
        'files_sha256':hashes}
    write(ROUND/'acceptance.json',acc)
    status='最新验收1204；新阶段4轮、8组。两独立能窗参考在同一原QQ／正能钟／实体记录过程中给有限定向流、标签响应与耗竭；失败偏置和反向读数推断误差已计入。下一项补固定参考准备程序的有限空间、门数、控制能源及总误差，再作阶段综合。'
    req={'conditional_identification_and_prediction':'same_H_clock_contact_record_finite_prediction_completed_in_1203',
        'independently_generated_stable_reference':'fixed independent energy-window references; explicit stocks and failures; preparation program finite resource bound pending',
        'internal_control_clock_record_preparation_resource_closure':'contact/readout clock complete within declared model; fixed preparation work-space/gate/control/combined-error bounds remain',
        'physical_temperature_flow_and_finite_nonequilibrium':'1204 scoped two-reference actual flow, equal-reference comparison, labelled response and one-use depletion certified',
        'actual_geometry_identity':'object/access/full contact-energy flux table completed; inherited-space and physical area/Unruh identity still open',
        'stage_synthesis':'pending fixed-preparation resource bridge and consolidated conditional theorem',
        'stage_deliverables_all_completed':False}
    ledger=read(HERE/'dependency_ledger.json')
    ledger.update(latest_accepted_round=1204,accepted_stage_rounds=[1201,1202,1203,1204],stage_check_groups=8,
        current_direction_document='../1204/NEXT.md',current_work=status,latest_management_receipt='progress_checks_1204.json',stage_requirement_status=req)
    write(HERE/'dependency_ledger.json',ledger)
    p=STAGE/'README.md';txt=p.read_text(encoding='utf-8');para=txt.split('\n\n');para[1]=status;txt='\n\n'.join(para)
    txt=txt.replace('## 接续与边界','- [1204报告](research_note_1204.md)：两参考实际流、保失败标签响应与多参考热力学；[验收](1204/acceptance.json)、[下一项](1204/NEXT.md)。\n\n## 接续与边界',1)
    txt=txt.replace('_shared/progress_checks_1203.json','_shared/progress_checks_1204.json')
    txt=txt.replace('1201—1203已接原参考、内部接触与实体记录；温度实际流用途和区域身份继续核验','1201—1204已接参考、内部接触／记录及有限实际流；固定准备资源上界和阶段综合仍待完成')
    txt+='\n[同一过程的对象与几何接口](_shared/current_process_geometry_interfaces.md)：材料、参考、钟、记录、控制边及完整能流；不把内部X当作三维位置。\n'
    p.write_text(txt,encoding='utf-8')
    p=STAGE/'文件索引.md';txt=p.read_text(encoding='utf-8').replace('最新验收1203；3轮、6组','最新验收1204；4轮、8组')
    txt=txt.replace('## 管理与接续','|1204|[独立参考之间的能源流与有限非平衡响应](research_note_1204.md)|[证明](1204/proof.md)、[代码](1204/check.py)、[结果](1204/results.json)、[独立核验](1204/independent_check.py)、[审阅](1204/independent_review.md)、[验收](1204/acceptance.json)|\n\n## 管理与接续',1)
    # Old appended rows contained blank lines; keep one valid Markdown table.
    while re.search(r'(\|[^\n]*\|)\n\n(?=\|[^\n]*\|)',txt):txt=re.sub(r'(\|[^\n]*\|)\n\n(?=\|[^\n]*\|)',r'\1\n',txt)
    txt=txt.replace('[下一项](1203/NEXT.md)','[下一项](1204/NEXT.md)').replace('_shared/progress_checks_1203.json','_shared/progress_checks_1204.json')
    txt+='\n- [阶段要求审计](_shared/stage_completion_audit_after1203.md)：准备资源的最小后补与后续几何桥分开。\n'
    txt+='- [对象、访问与几何接口](_shared/current_process_geometry_interfaces.md)\n'
    p.write_text(txt,encoding='utf-8')
    d=read(OLD/'_shared/dependency_ledger.json');d['latest_management_receipt']='../../archive_1201_/_shared/progress_checks_1204.json';write(OLD/'_shared/dependency_ledger.json',d)
    p=OLD/'_shared/navigation_notice.md';txt=p.read_text(encoding='utf-8');prefix,rest=txt.split('以下保留前阶段的历史状态',1)
    prefix=re.sub(r'\*\*当前状态：\*\*.*?(?=\n\n)','**当前状态：** '+status+' [1204报告](archive_1201_/research_note_1204.md) · [验收](archive_1201_/1204/acceptance.json) · [当前核验](archive_1201_/_shared/progress_checks_1204.json)。旧GR阶段40轮、130组和冻结稿均保留；目标继续进行。',prefix,flags=re.S)
    p.write_text(prefix+'以下保留前阶段的历史状态'+rest,encoding='utf-8')
    navapply=run(OLD/'_shared/update_navigation.py','--apply');navcheck=run(OLD/'_shared/update_navigation.py','--check');assert navcheck['changed_files']==0
    navverify=run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation');assert navverify['status']=='PASS'
    write(OUT,{'status':'PUBLISHING'})
    paths=[STAGE/'README.md',STAGE/'文件索引.md',STAGE/'research_note_1204.md',*ROUND.glob('*.md'),
           HERE/'stage_completion_audit_after1203.md',HERE/'current_process_geometry_interfaces.md',*[BASE/n for n in ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md']]]
    nlinks=links(paths)
    science=dict(prev['science_files_sha256']);science.update({(ROUND/name).resolve().relative_to(ROOT).as_posix():h for name,h in hashes.items()})
    management={name:sha(ROOT/name) for name in prev['management_sha256']};management[Path(__file__).relative_to(ROOT).as_posix()]=sha(Path(__file__))
    screens=dict(prev['screen_files_sha256'])
    for name in ['stage_completion_audit_after1203.md','current_process_geometry_interfaces.md']:
        p=HERE/name;screens[p.relative_to(ROOT).as_posix()]=sha(p)
    accs=dict(prev['acceptance_sha256']);accs['1204']=sha(ROUND/'acceptance.json')
    write(OUT,{'schema':'thermodynamics_progress_1204_v1','date':'2026-10-14','status':'PASS',
        'accepted_rounds':[1201,1202,1203,1204],'stage_check_groups':8,'application_goal_status':'active','stage_complete':False,
        'goal_objective_sha256':ledger['goal_objective_sha256'],'acceptance_sha256':accs,
        'science_files_sha256':science,'management_sha256':management,'screen_files_sha256':screens,
        'prior_progress_path':pp.relative_to(ROOT).as_posix(),'prior_progress_sha256':sha(pp),
        'previous_GR_science_files_verified':len(old['science_files_sha256']),'previous_GR_bindings':n,
        'navigation_apply':navapply,'navigation_check':navcheck,'navigation_verify':navverify,
        'markdown_files_checked':len(paths),'local_links_checked':nlinks,'missing_links':[],
        'stage_requirement_status':req,'new_adopted_axioms':0,'new_empirical_calibration':0,'GR_hard_deliverables_completed':0,
        'this_goal_turn_result':'progress: 1203 autonomous clock/record/thermodynamics and 1204 two-reference physical flow/finite response; fixed preparation resource closure admitted next'})
    verify()

if __name__=='__main__':main()
