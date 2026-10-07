"""943 delivery audit, separate from recalculation of the scientific witness."""
from pathlib import Path
import argparse, ast, hashlib, json, re, sys
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_943_checks.json'

def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,943):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):
        add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify();r=read(HERE/'material_content_access_results.json')
    assert r['round']==943 and r['all_scientific_checks_passed']
    for rel,digest in r['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    h=r['actual_929_history'];q=r['invariant_bulk_probe'];b=r['new_flavour_portal']
    assert h['full_declared_flavour_dimension']==4063232 and h['computed_history_shape']==[31,8192,2]
    assert h['actual_gate_weighted_intertwiner_error']<1e-12
    assert abs(h['event_probability_gap']-.6)<1e-13
    assert max(q['input_output_trace_distances'])<1e-12 and q['passive_reference_factorization_error']<1e-12
    assert q['extreme_mass_conditional_probe_trace_distance']>.1
    assert b['mass_commutator_error']<1e-13 and b['original_flavour_group_commutator_norm']>.5
    assert b['finite_pointer_probability_gap']>.7 and b['new_Yukawa_matrix_is_additional_physical_input']
    assert not b['actual_p_to_R_field_readout_proved'] and not q['full_Einstein_evolution_simulated']
    assert r['analytic_scope']['candidate_access_failure_not_whole_program_failure'] and not r['full_goal_completed']
    note=STAGE/'research_note_943.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==20 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,11)]
    for phrase in ('整体目标未完成','新增交互是物理输入','停止此编码','0.30','与质量对易','不是完整Einstein演化模拟'):assert phrase in prose,phrase
    newdocs=[note,HERE/'drafts/access_decision.md',STAGE/'944/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split('## 六条共同协议：全局缺口对应与检验优先级（截至943）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=943)==list(range(1,944))
    assert '001—943轮共943份' in nav[0].read_text('utf-8-sig')
    assert '231—943的713份' in nav[5].read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '943：内部记录能否进入共同物理交互' in s and '944/drafts/STATUS.md' in s
    oldfiles=[STAGE/'942/research_round_942_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'material_content_access.py',HERE/'material_content_access_results.json',HERE/'drafts/update_navigation943.py']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(oldfiles[0])
    out=dict(round=943,date='2026-10-07',all_delivery_checks_passed=True,
       fresh_test_groups=1,formal_reports=943,
       cumulative_numbered_test_groups_from_942=previous['cumulative_numbered_test_groups_from_941']+1,
       historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
       old_spectral_and_access_lemmas_reused=True,
       actual_942_action_and_929_logical_content_access_tested=True,
       finite_task_error_lower_bound_documented=True,
       new_interaction_declared_input_not_original_action_claim=True,
       candidate_access_limit_not_program_failure=True,
       candidate_refinement_stopped_after_decisive_test=True,
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
