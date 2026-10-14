"""Publish round1201 once; validate frozen science and reversible navigation."""
from pathlib import Path
import hashlib,json,re,subprocess,sys
from urllib.parse import unquote
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'research_cognition_physics'
STAGE=BASE/'archive_1201_'; HERE=STAGE/'_shared'; ROUND=STAGE/'1201'
OLD=BASE/'archive_1161_'; RECEIPT=HERE/'progress_checks_1201.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def write(p,obj):p.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def validate_map(m):
    for name,digest in m.items():assert sha(ROOT/name)==digest,'Changed: '+name
def old_check(startup):
    p=ROOT/startup['previous_stage_receipt'];assert sha(p)==startup['previous_stage_receipt_sha256']
    old=read(p);validate_map(old['science_files_sha256']);count=0
    for number,digest in old['acceptance_sha256'].items():
        p=OLD/number/'acceptance.json';assert sha(p)==digest
        for name,h in read(p)['files_sha256'].items():
            assert sha((p.parent/name).resolve())==h,(number,name);count+=1
    assert count==485
    return old,count
def run(*args):
    r=subprocess.run([sys.executable,'-B','-X','utf8',*map(str,args)],capture_output=True,text=True,encoding='utf-8')
    assert r.returncode==0,r.stderr+'\n'+r.stdout
    return json.loads(r.stdout)
def links(paths):
    count=0;fence=chr(96)*3
    for p in paths:
        text=p.read_text(encoding='utf-8-sig')
        text=re.sub(r'(?ms)^'+fence+r'.*?^'+fence+r'[^\n]*','',text)
        assert text.count('$$')%2==0,str(p)
        for target in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)',text):
            target=target.strip('<> ').split('#')[0]
            if not target or re.match(r'^[a-zA-Z]+:',target):continue
            assert (p.parent/unquote(target)).resolve().exists(),(p,target)
            count+=1
    return count
def verify():
    receipt=read(RECEIPT);old,count=old_check(read(HERE/'startup_checks.json'))
    validate_map(receipt['science_files_sha256']);validate_map(receipt['management_sha256'])
    validate_map(receipt['screen_files_sha256'])
    acc=ROUND/'acceptance.json';assert sha(acc)==receipt['acceptance_sha256']
    for name,h in read(acc)['files_sha256'].items():assert sha((ROUND/name).resolve())==h
    nav=run(OLD/'_shared/update_navigation.py','--check');assert nav['changed_files']==0
    actual=run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation');assert actual['status']=='PASS'
    print(json.dumps({'status':'PASS','round':1201,'previous_science_files':len(old['science_files_sha256']),
        'previous_acceptance_bindings':count,'new_science_files':len(receipt['science_files_sha256']),
        'goal_status':'active','stage_complete':False},ensure_ascii=False))
def main():
    if '--verify' in sys.argv:return verify()
    assert not RECEIPT.exists() and not (ROUND/'acceptance.json').exists(),'Already published'
    startup=read(HERE/'startup_checks.json')
    validate_map(startup['new_stage_files_sha256']);validate_map(startup['management_sha256'])
    old,count=old_check(startup)
    assert '主证明补充后的限定签核' in (ROUND/'certificate_review.md').read_text(encoding='utf-8'),'Final review pending'
    author=run(ROUND/'check.py');independent=run(ROUND/'independent_check_1201.py')
    assert author['status']=='PASS_FINITE_CONDITIONAL_CERTIFICATE' and independent['all_candidate_assertions_passed']
    result=read(ROUND/'results.json')
    required=['check.py','results.json','proof.md','independent_check_1201.py',
        'independent_results_1201.json','independent_candidate_review.md','protocol_resource_audit.md',
        'certificate_review.md','NEXT.md','../research_note_1201.md','../_shared/goal.md']
    hashes={s:sha((ROUND/s).resolve()) for s in required}
    acceptance={'schema':'operational_thermodynamics_round_acceptance_v1','round':1201,'date':'2026-10-14',
        'status':'ACCEPTED_WITH_SCOPE','check_groups':2,
        'scope':'原完整16维QQ有限接触在未知平稳参考族上的人口识别、全实温度族有限排除及另一固定窗能源预测；保开关均值能源、指针等距与任意BR。内部参考生成、实际钟/控制/读出与三维区域身份未完成。',
        'certified':{'g':'1/10','identification_time':600,'holdout_time':300,'matrix_entry_error_upper':'1e-10',
            'actual_inverse_inf_norm_upper':'3.14','population_error_upper':'.015701',
            'nonthermal_compatible_set_gibbs_inf_distance_lower':'.093598',
            'model_holdout_energy_errors_upper':['.009262','.009294','.000100'],
            'holdout_instrument_errors_must_be_added':True,'fresh_iid_trials_for_declared_statistics':250000,
            'statistical_failure_upper':8*__import__('math').exp(-8),
            'exact_deletion':'center/outer reference coherences invisible to this parity-respecting menu'},
        'inputs_retained':result['remaining_inputs'],
        'stage_requirement_status':{'conditional_identification_and_prediction':'completed_for_declared_model',
            'independently_generated_stable_reference':'incomplete',
            'internal_control_clock_record_preparation_resource_closure':'incomplete',
            'actual_geometry_identity':'incomplete','stage_deliverables_all_completed':False},
        'new_adopted_axioms':0,'new_empirical_calibration':0,'GR_hard_deliverables_completed':0,
        'review':'Internal independent native NumPy rebuild plus separate rational/rounding review; not blind experimental validation. Candidate-review pending language is historical; final scope is this acceptance.',
        'files_sha256':hashes,'root_executed_both_codes':True,
        'previous_scientific_files_verified':len(old['science_files_sha256']),
        'previous_acceptance_bindings_verified':count}
    write(ROUND/'acceptance.json',acceptance)
    state='最新验收为1201，新阶段完成1轮、2组核验。原完整QQ接触已给未知平稳参考的人口识别、非热排除及300窗能源留出，逆范数界3.14；独立参考准备、内部控制/钟/读出和三维区域身份仍未完成。当前接续独立能源条件准备，不继续优化测温矩阵。'
    (STAGE/'README.md').write_text('# 1201起：认知操作热力学\n\n'+state+'\n\n'
        '目标保持进行中：[完整目标](_shared/goal.md) · [启动回执](_shared/goal_start_receipt.json)。\n\n'
        '推进顺序：认知操作 → 操作热力学 → 热力学—几何对应 → 有效GR。阶段要求同一过程的有限预测及完整资源/误差说明；1201只有条件连接，尚未完成全部交付。\n\n'
        '## 已验收\n\n'
        '- [1201报告](research_note_1201.md)：实际16维接触、有限精度、任意BR、完整交互平均能源与记录等距；模型预测误差与实际仪器误差分开。\n'
        '- [验收](1201/acceptance.json) · [证明](1201/proof.md) · [下一项](1201/NEXT.md)。\n\n'
        '## 接续与边界\n\n'
        '[三体准备预筛](_shared/reference_preparation_screen.md)未得到稳定单温度参考，且原共振交换保各能级人数。此为未编号筛选，不计新轮，不推广为一般热化不可能。\n\n'
        '保留[旧GR阶段](../archive_1161_/README.md)1161—1200的40轮、130组及全部冻结成果，旧目标未被标成完成。面积熵、Unruh与动态几何仍须独立桥接；暗部门保留为可能共同结果。\n\n'
        '[轮次索引](文件索引.md) · [历史复用](_shared/history_reuse.md) · [当前核验](_shared/progress_checks_1201.json)。\n',encoding='utf-8')
    (STAGE/'文件索引.md').write_text('# 认知操作热力学阶段索引\n\n2026-10-14。最新验收1201；1轮、2组。目标启用、文献回查和预筛不计研究轮次。\n\n'
        '|轮次|报告与结论|代码和证据|\n|---|---|---|\n'
        '|1201|[未知参考的有限接触测温与能源留出](research_note_1201.md)|[证明](1201/proof.md)、[代码](1201/check.py)、[结果](1201/results.json)、[独立代码](1201/independent_check_1201.py)、[数学审阅](1201/certificate_review.md)、[验收](1201/acceptance.json)|\n\n'
        '## 管理与接续\n\n'
        '- [阶段说明](README.md)、[完整目标](_shared/goal.md)、[交接](_shared/handoff.md)、[启动回执](_shared/goal_start_receipt.json)\n'
        '- [历史复用](_shared/history_reuse.md)、[未编号三体准备预筛](_shared/reference_preparation_screen.md)\n'
        '- [下一项](1201/NEXT.md)、[当前核验](_shared/progress_checks_1201.json)、[前阶段1200](../archive_1161_/research_note_1200.md)\n',encoding='utf-8')
    ledger=read(HERE/'dependency_ledger.json')
    ledger.update(latest_accepted_round=1201,accepted_stage_rounds=[1201],stage_check_groups=2,
        current_direction_document='../1201/NEXT.md',current_work=state,
        latest_management_receipt='progress_checks_1201.json',
        stage_requirement_status=acceptance['stage_requirement_status'])
    write(HERE/'dependency_ledger.json',ledger)
    oldledger=read(OLD/'_shared/dependency_ledger.json')
    oldledger['latest_management_receipt']='../../archive_1201_/_shared/progress_checks_1201.json'
    write(OLD/'_shared/dependency_ledger.json',oldledger)
    notice=OLD/'_shared/navigation_notice.md'
    text=notice.read_text(encoding='utf-8-sig');before,rest=text.split('以下保留前阶段的历史状态',1)
    before=re.sub(r'\*\*当前状态：\*\*.*?(?=\n\n)', '**当前状态：** '+state+
        ' [1201报告](archive_1201_/research_note_1201.md) · [验收](archive_1201_/1201/acceptance.json) · [核验](archive_1201_/_shared/progress_checks_1201.json)。旧GR阶段40轮、130组保留，完整GR未结项。',before,flags=re.S)
    notice.write_text(before+'以下保留前阶段的历史状态'+rest,encoding='utf-8')
    navapply=run(OLD/'_shared/update_navigation.py','--apply')
    navcheck=run(OLD/'_shared/update_navigation.py','--check');assert navcheck['changed_files']==0
    navverify=run(ROOT/'scripts/run_research_gr_active.py','--verify-navigation');assert navverify['status']=='PASS'
    write(RECEIPT,{'status':'PUBLISHING'})
    markdown=list(STAGE.rglob('*.md'))+[BASE/n for n in ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md']]
    checked=links(markdown)
    science={(ROUND/s).resolve().relative_to(ROOT).as_posix():h for s,h in hashes.items()}
    mgmt={p:sha(ROOT/p) for p in startup['management_sha256']}
    for p in [STAGE/'README.md',STAGE/'文件索引.md',HERE/'dependency_ledger.json',Path(__file__)]:
        mgmt[p.relative_to(ROOT).as_posix()]=sha(p)
    screens={p.relative_to(ROOT).as_posix():sha(p) for p in HERE.glob('*screen.md')}
    receipt={'schema':'thermodynamics_progress_1201_v1','date':'2026-10-14','status':'PASS',
        'accepted_rounds':[1201],'stage_check_groups':2,'application_goal_status':'active','stage_complete':False,
        'goal_objective_sha256':ledger['goal_objective_sha256'],'acceptance_sha256':sha(ROUND/'acceptance.json'),
        'science_files_sha256':science,'management_sha256':mgmt,'screen_files_sha256':screens,
        'previous_stage_receipt':startup['previous_stage_receipt'],'previous_stage_receipt_sha256':startup['previous_stage_receipt_sha256'],
        'previous_science_files_verified':len(old['science_files_sha256']),'previous_acceptances_verified':40,
        'previous_acceptance_bindings':count,'navigation_apply':navapply,'navigation_check':navcheck,
        'navigation_verify':navverify,'markdown_files_checked':len(markdown),'local_links_checked':checked,'missing_links':[],
        'new_adopted_axioms':0,'new_empirical_calibration':0,'GR_hard_deliverables_completed':0,
        'this_goal_turn_result':'progress: accepted finite unknown-reference thermometry and independent-time energy prediction; stopped failed small-chain preparation; next reference source screened separately'}
    write(RECEIPT,receipt);verify()
if __name__=='__main__':main()
