"""957 delivery audit, separate from recalculation of the scientific witness."""
from pathlib import Path
import argparse, ast, hashlib, json, re, sys
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_957_checks.json'

def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,957):
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
    layout=Layout().verify();r=read(HERE/'task_domain_budget_results.json')
    assert r['round']==957 and r['all_scientific_checks_passed']
    for rel,digest in r['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    a=r['actual_955_budget'];b=r['actual_956_budget'];scope=r['scope']
    assert a['sum_integrated_action_upper']<5214
    assert a['total_H_rms_upper']<13.248
    assert a['term_rms_upper']['source_field_interaction']<.773
    assert abs(b['all_state_integrated_action_upper']-47.1238898038469)<1e-10
    example=r['separating_example']
    assert example['full_interaction_operator_norm']=='infinite'
    rows=example['fixed_second_moment_rows']
    assert all(x['unrestricted_trace_distance']>.999999 for x in rows)
    assert all(x['task_trace_distance']<=x['task_trace_distance_upper'] for x in rows)
    assert scope['A1_D_is_declared_revision_not_old_A1_proved']
    assert scope['all_955_H_terms_controlled_during_saved_free_run']
    assert not scope['arbitrary_955_terminal_instrument_moment_closure_certified']
    assert not scope['old_propagation_or_dimension_theorems_automatically_inherited']
    assert not scope['complete_native_parent_model_certified']
    assert not scope['full_goal_completed']
    note=STAGE/'research_note_957.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==18 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,10)]
    for phrase in ('整体目标未完成','A1_D','Duhamel','停止','保留'):
        assert phrase in prose,phrase
    newdocs=[note,HERE/'drafts/common_contract_decision.md',
        HERE/'drafts/unified_operation_hypotheses_v0_2.md',STAGE/'958/drafts/STATUS.md']
    contract=newdocs[2].read_text('utf-8')
    assert '唯一实质修订是A1' in contract
    assert all(f'A{i}' in contract for i in range(2,8))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split('## 六条共同协议：全局缺口对应与检验优先级（截至957）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=957)==list(range(1,958))
    assert '001—957轮共957份' in nav[0].read_text('utf-8-sig')
    assert '231—957的727份' in nav[5].read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '957：任务域预算与共同假说v0.2' in s and '957/drafts/STATUS.md' in s
    oldfiles=[STAGE/'956/research_round_956_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'task_domain_budget.py',HERE/'task_domain_budget_results.json',HERE/'drafts/update_navigation957.py']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(oldfiles[0])
    out=dict(round=957,date='2026-10-07',all_delivery_checks_passed=True,
       fresh_test_groups=1,formal_reports=957,
       cumulative_numbered_test_groups_from_956=previous['cumulative_numbered_test_groups_from_955']+1,
       historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
       explicit_A1_working_revision_documented=True,
       actual_955_free_run_all_terms_budgeted=True,
       actual_956_full_material_strong_budget_reused=True,
       unrestricted_and_task_limited_norms_separated=True,
       instruments_and_sources_not_automatically_certified=True,
       no_old_dimension_or_propagation_premise_removed=True,
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
