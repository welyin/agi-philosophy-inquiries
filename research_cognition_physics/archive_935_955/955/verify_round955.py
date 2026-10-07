"""955 delivery audit, separate from recalculation of the scientific witness."""
from pathlib import Path
import argparse, ast, hashlib, json, re, sys
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_955_checks.json'

def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,955):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):
        add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'953/drafts/priority_reaudit_checks.json');add(au['frozen_evidence_hashes']);add(au['audit_file_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify();r=read(HERE/'joint_effective_window_results.json')
    assert r['round']==955 and r['all_scientific_checks_passed']
    for rel,digest in r['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    b=r['bounds'];scope=r['scope'];center=r['fixed_center']
    assert b['global_cq_instrument_distance']<.000622
    assert b['weak_higgs_total_error']<.000083
    assert b['higgs_probability_contrast_lower']>.000099
    assert b['record_probability_contrast_lower']>.0127
    assert b['gravity_probability_contrast_lower']>.0217
    assert b['using_global_bound_for_Higgs_lower']<0
    assert center['Choi_recovery_error']<1e-12 and center['instrument_Gram_minimum']>-1e-12
    assert r['recovery']['matter_relative_kernel_bound']<.003
    assert r['recovery']['record_W_Q1']>.002
    assert scope['actual_two_protocol_unknown_input_and_passive_reference_transported']
    assert scope['full_instrument_bound_not_confused_with_effect_bound']
    assert scope['new_quadratic_field_vacuum_declared']
    assert not scope['complete_covariant_stress_and_quantum_gravity_certified']
    assert not scope['all_access_instruments_autonomously_manufactured']
    assert not scope['full_goal_completed']
    note=STAGE/'research_note_955.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==16 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,9)]
    for phrase in ('整体目标未完成','Heisenberg','停止','不是实验约束','真实运动'):
        if phrase=='不是实验约束':
            assert '不是全部现实实验预算' in (HERE/'drafts/common_model_scope_v0_5.md').read_text('utf-8')
        else:assert phrase in prose,phrase
    newdocs=[note,HERE/'drafts/joint_window_decision.md',HERE/'drafts/common_model_scope_v0_5.md',STAGE/'956/drafts/STATUS.md']
    next_prose=newdocs[-1].read_text('utf-8')
    assert '先比较共同模型接口' in next_prose and '采用决定' in next_prose
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split('## 六条共同协议：全局缺口对应与检验优先级（截至955）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    links=0
    for doc in newdocs+nav:
        content=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',content):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if m and p.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=955)==list(range(1,956))
    assert '001—955轮共955份' in nav[0].read_text('utf-8-sig')
    assert '231—955的725份' in nav[5].read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '955：物质恢复点与有限协议共同采用' in s and '956/drafts/STATUS.md' in s
    oldfiles=[STAGE/'954/research_round_954_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'joint_effective_window.py',HERE/'joint_effective_window_results.json',HERE/'drafts/update_navigation955.py']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(oldfiles[0])
    out=dict(round=955,date='2026-10-07',all_delivery_checks_passed=True,
       fresh_test_groups=1,formal_reports=955,
       cumulative_numbered_test_groups_from_954=previous['cumulative_numbered_test_groups_from_953']+1,
       historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
       same_parameter_finite_instrument_and_matter_window_verified=True,
       observable_specific_error_not_full_instrument_bound=True,
       general_bound_negative_control_verified=True,
       new_vacuum_and_preexisting_resource_scope_declared=True,
       next_work_requires_candidate_adoption_comparison=True,
       no_full_SM_or_quantum_gravity_certificate_claimed=True,
       full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
       frozen_inputs={str(q.relative_to(ROOT)):sha(q) for q in oldfiles},
       new_scientific_and_entry_files={str(q.relative_to(ROOT)):sha(q) for q in newfiles})
    if not writing:
        before=read(TARGET)
        for k in ('frozen_inputs','new_scientific_and_entry_files'):assert before[k]==out[k],k
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    out=run(a.write)
    if a.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
