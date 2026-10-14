"""Publish only independently reviewed round1203; preserve historical bindings."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util,json,re,sys
ROOT=Path(__file__).resolve().parents[3]; BASE=ROOT/'research_cognition_physics'
STAGE=BASE/'archive_1201_'; HERE=STAGE/'_shared'; ROUND=STAGE/'1203'; OLD=BASE/'archive_1161_'
s=importlib.util.spec_from_file_location('publish1201_tools',HERE/'publish_1201.py')
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
sha,read,write,run,links,old_check,validate_map=m.sha,m.read,m.write,m.run,m.links,m.old_check,m.validate_map
OUT=HERE/'progress_checks_1203.json'

def verify():
    d=read(OUT);old,n=old_check(read(HERE/'startup_checks.json'))
    for key in ['science_files_sha256','management_sha256','screen_files_sha256']:validate_map(d[key])
    for number,h in d['acceptance_sha256'].items():
        p=STAGE/number/'acceptance.json';assert sha(p)==h
        for name,hsh in read(p)['files_sha256'].items():assert sha((p.parent/name).resolve())==hsh
    assert run(OLD/'_shared/update_navigation.py','--check')['changed_files']==0
    assert run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation')['status']=='PASS'
    print(json.dumps({'status':'PASS','latest_round':1203,'new_stage_rounds':3,'stage_check_groups':6,
        'old_GR_science_files_verified':len(old['science_files_sha256']),'old_GR_bindings_verified':n,
        'current_stage_science_files':len(d['science_files_sha256']),'goal_status':'active','stage_complete':False}))

def main():
    if '--verify' in sys.argv:return verify()
    assert not OUT.exists() and not (ROUND/'acceptance.json').exists()
    prevpath=HERE/'progress_checks_1202.json';prev=read(prevpath)
    for key in ['science_files_sha256','management_sha256','screen_files_sha256']:validate_map(prev[key])
    for number,h in prev['acceptance_sha256'].items():assert sha(STAGE/number/'acceptance.json')==h
    old,count=old_check(read(HERE/'startup_checks.json'))
    assert '最终限定签核' in (ROUND/'resource_review.md').read_text(encoding='utf-8')
    a=run(ROUND/'check.py');b=run(ROUND/'independent_check.py')
    assert a['status']=='PASS_SCOPED_CLOCK_RECORD_CONNECTION'
    # Independent schema is inspected by root; do not infer success from file existence.
    assert b.get('all_checks_passed') is True
    r=read(ROUND/'results.json')
    assert F(r['combined_isometry_error_exact_upper'])==F(b['joint']['isometry_bound_exact'])
    for i in range(4):
        for j in range(4):assert abs(r['pointer_joint_probabilities'][i][j]-b['joint']['pointer_table'][i][j])<2e-11
    assert F(r['thermodynamics']['AK_only_clausius_lower'])<=F(str(b['entropy']['actual_AK_plus_beta_delta_EB_lower_outward_15dp']))
    assert F(r['combined_isometry_error_upper'])<F('.00017')
    assert F(r['thermodynamics']['AK_only_clausius_lower'])>F('-.000913')
    files=['proof.md','check.py','results.json','independent_check.py','independent_results.json',
           'independent_review.md','resource_review.md','NEXT.md','../research_note_1203.md','../_shared/goal.md']
    hashes={name:sha((ROUND/name).resolve()) for name in files}
    acceptance={'schema':'operational_thermodynamics_round_acceptance_v1','round':1203,
        'date':'2026-10-14','status':'ACCEPTED_WITH_SCOPE','check_groups':2,
        'scope':'同一正能正常钟父H接原QQ、独立能窗参考和实体A能源记录；共同终时606的全未知输入/参考误差、完整钟供能及保钟/可见记录Clausius缺陷。参考准备程序、B能源指针、区域几何和实际误差标定未完成。',
        'certified':{'clock_H':'P^2/(2M)','M':'1e24','sigma':'1e-10','common_terminal_time':606,
            'original_contact_windows':[600,300],'same_H_program_labels':[1,2],
            'physical_A_energy_pointer_dimension':4,'postselect_clock':False,
            'combined_prepared_isometry_error_upper':r['combined_isometry_error_upper'],
            'joint_pointer_table_TV_error_upper':r['pointer_joint_table_TV_error_upper'],
            'retained_AKC_clausius_lower':'-.000912',
            'visible_AK_clausius_lower':r['thermodynamics']['AK_only_clausius_lower'],
            'each_clock_entropy_rank_bound':256,'clock_impurity_probability_is_squared_amplitude_bound':True,
            'energy_interaction_tails_include_contact_and_record':True,
            'native_switch_work_supplied_by_actual_clock_with_error':True,
            'both_clocks_mean_energy_losses_strictly_positive':True},
        'inputs_retained':r['remaining_inputs'],'limits':r['not_claimed'],
        'new_adopted_axioms':0,'new_empirical_calibration':0,'GR_hard_deliverables_completed':0,
        'stage_complete':False,'root_executed_both_codes':True,
        'review':'Independent rational clock budget and native material/pointer rebuild; separate resource and entropy review. Internal mathematical review, not empirical validation.',
        'files_sha256':hashes}
    write(ROUND/'acceptance.json',acceptance)
    status='最新验收1203；新阶段3轮、6组。独立参考、原QQ、正能内部钟与实体能源记录已在同一候选过程接续，联合误差小于.00017；钟实际承担原开关平均功，可见记录的Clausius下界大于−.000912389042。参考自主准备和几何身份仍有公开条件；下一项检验两参考之间的温度—流及有限非平衡响应。'
    req={'conditional_identification_and_prediction':'same_H_clock_contact_record_finite_prediction_completed_in_1203',
         'independently_generated_stable_reference':'1202 independent finite energy-window preparation; autonomous preparation remains conditional',
         'internal_control_clock_record_preparation_resource_closure':'contact and A record autonomously connected with clock backreaction; preparation, final interrogation and physical calibration still inputs',
         'physical_temperature_flow_and_finite_nonequilibrium':'next: independent references, flow versus clock work, finite depletion',
         'actual_geometry_identity':'incomplete; internal clock X not certified as old 3D position or region boundary',
         'stage_deliverables_all_completed':False}
    ledger=read(HERE/'dependency_ledger.json')
    ledger.update(latest_accepted_round=1203,accepted_stage_rounds=[1201,1202,1203],stage_check_groups=6,
        current_direction_document='../1203/NEXT.md',current_work=status,latest_management_receipt='progress_checks_1203.json',
        stage_requirement_status=req)
    write(HERE/'dependency_ledger.json',ledger)
    p=STAGE/'README.md';txt=p.read_text(encoding='utf-8');paragraphs=txt.split('\n\n');paragraphs[1]=status;txt='\n\n'.join(paragraphs)
    txt=txt.replace('## 接续与边界','- [1203报告](research_note_1203.md)：正能内部钟、实体能源记录、完整钟供能及有限热力学；[验收](1203/acceptance.json)、[下一项](1203/NEXT.md)。\n\n## 接续与边界',1)
    txt=txt.replace('_shared/progress_checks_1202.json','_shared/progress_checks_1203.json')
    txt=txt.replace('1201—1202完成条件连接，内部控制和区域身份仍未闭合','1201—1203已接原参考、内部接触与实体记录；温度实际流用途和区域身份继续核验')
    p.write_text(txt,encoding='utf-8')
    p=STAGE/'文件索引.md';txt=p.read_text(encoding='utf-8').replace('最新验收1202；2轮、4组','最新验收1203；3轮、6组')
    txt=txt.replace('## 管理与接续','|1203|[内部正能钟、实体记录与有限热力学](research_note_1203.md)|[证明](1203/proof.md)、[代码](1203/check.py)、[结果](1203/results.json)、[独立核验](1203/independent_check.py)、[资源审阅](1203/resource_review.md)、[验收](1203/acceptance.json)|\n\n## 管理与接续',1)
    txt=txt.replace('[下一项](1202/NEXT.md)','[下一项](1203/NEXT.md)').replace('_shared/progress_checks_1202.json','_shared/progress_checks_1203.json')
    txt+='\n- [1203后方向预审](_shared/after1203_scope_screen.md)：温度实际流与有限非平衡，未编号。\n'
    p.write_text(txt,encoding='utf-8')
    d=read(OLD/'_shared/dependency_ledger.json');d['latest_management_receipt']='../../archive_1201_/_shared/progress_checks_1203.json';write(OLD/'_shared/dependency_ledger.json',d)
    p=OLD/'_shared/navigation_notice.md';txt=p.read_text(encoding='utf-8');prefix,rest=txt.split('以下保留前阶段的历史状态',1)
    prefix=re.sub(r'\*\*当前状态：\*\*.*?(?=\n\n)','**当前状态：** '+status+
        ' [1203报告](archive_1201_/research_note_1203.md) · [验收](archive_1201_/1203/acceptance.json) · [当前核验](archive_1201_/_shared/progress_checks_1203.json)。旧GR阶段40轮、130组和全部冻结稿保留；阶段未结项。',prefix,flags=re.S)
    p.write_text(prefix+'以下保留前阶段的历史状态'+rest,encoding='utf-8')
    navapply=run(OLD/'_shared/update_navigation.py','--apply');navcheck=run(OLD/'_shared/update_navigation.py','--check');assert navcheck['changed_files']==0
    navverify=run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation');assert navverify['status']=='PASS'
    write(OUT,{'status':'PUBLISHING'})
    paths=[STAGE/'README.md',STAGE/'文件索引.md',STAGE/'research_note_1203.md',*ROUND.glob('*.md'),
        HERE/'autonomous_contact_screen.md',HERE/'after1203_scope_screen.md',
        *[BASE/n for n in ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md']]]
    nlinks=links(paths)
    science=dict(prev['science_files_sha256']);science.update({(ROUND/n).resolve().relative_to(ROOT).as_posix():h for n,h in hashes.items()})
    management={p:sha(ROOT/p) for p in prev['management_sha256']};management[Path(__file__).relative_to(ROOT).as_posix()]=sha(Path(__file__))
    screens=dict(prev['screen_files_sha256'])
    for n in ['autonomous_contact_screen.md','after1203_scope_screen.md']:
        p=HERE/n;screens[p.relative_to(ROOT).as_posix()]=sha(p)
    accs=dict(prev['acceptance_sha256']);accs['1203']=sha(ROUND/'acceptance.json')
    write(OUT,{'schema':'thermodynamics_progress_1203_v1','date':'2026-10-14','status':'PASS',
        'accepted_rounds':[1201,1202,1203],'stage_check_groups':6,'application_goal_status':'active','stage_complete':False,
        'goal_objective_sha256':ledger['goal_objective_sha256'],'acceptance_sha256':accs,
        'science_files_sha256':science,'management_sha256':management,'screen_files_sha256':screens,
        'prior_progress_path':prevpath.relative_to(ROOT).as_posix(),'prior_progress_sha256':sha(prevpath),
        'previous_GR_science_files_verified':len(old['science_files_sha256']),'previous_GR_bindings':count,
        'navigation_apply':navapply,'navigation_check':navcheck,'navigation_verify':navverify,
        'markdown_files_checked':len(paths),'local_links_checked':nlinks,'missing_links':[],
        'stage_requirement_status':req,'new_adopted_axioms':0,'new_empirical_calibration':0,
        'GR_hard_deliverables_completed':0,
        'this_goal_turn_result':'progress: same autonomous clock/contact/physical record and finite thermodynamics connected, clock backreaction and entropy certified; actual temperature-flow test admitted next'})
    verify()

if __name__=='__main__':main()
