"""947 delivery audit, separate from recalculation of the scientific witness."""
from pathlib import Path
import argparse, ast, hashlib, json, re, sys
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_947_checks.json'

def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,947):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):
        add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify();r=read(HERE/'protocol_field_transport_results.json')
    assert r['round']==947 and r['all_scientific_checks_passed']
    for rel,digest in r['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    q=r['same_physical_protocol'];j=r['joint_instrument'];b=r['finite_joint_control']
    assert q['event_coupling_h_commutator_error']<1e-12
    assert q['complete_event_intertwining_error']<1e-12
    assert q['internal_fault_sectors']==16 and q['repeats']==2
    assert j['choi_branch_recovery_error']<1e-12
    assert j['minimum_event_effect_eigenvalue']>.04
    assert b['full_cq_instrument_distance_bound']<.000635
    assert b['postevent_normalized_private_state_distance_upper']<.033
    assert b['record_probability_contrast_lower']>.012
    assert b['higgs_probability_contrast_lower']>.00138
    assert b['gravity_probability_contrast_lower']>.021
    assert r['scope']['original_929_finite_six_protocol_mathematical_contract_jointly_transported']
    assert not r['scope']['autonomous_construction_of_all_access_instruments_certified']
    assert not r['scope']['full_goal_completed']
    note=STAGE/'research_note_947.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==14 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,8)]
    for phrase in ('整体目标未完成','0.000633467','停止','原六协议有限数学合同','尚未签收'):assert phrase in prose,phrase
    newdocs=[note,HERE/'drafts/protocol_transport_decision.md',STAGE/'948/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split('## 六条共同协议：全局缺口对应与检验优先级（截至947）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=947)==list(range(1,948))
    assert '001—947轮共947份' in nav[0].read_text('utf-8-sig')
    assert '231—947的717份' in nav[5].read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '947：完整有限协议与同一场—引力过程的共同运输' in s and '948/drafts/STATUS.md' in s
    oldfiles=[STAGE/'946/research_round_946_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'protocol_field_transport.py',HERE/'protocol_field_transport_results.json',HERE/'drafts/update_navigation947.py']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(oldfiles[0])
    out=dict(round=947,date='2026-10-07',all_delivery_checks_passed=True,
       fresh_test_groups=1,formal_reports=947,
       cumulative_numbered_test_groups_from_946=previous['cumulative_numbered_test_groups_from_945']+1,
       historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
       original_protocol_transported_through_same_physical_H=True,
       actual_event_sensitive_coupling_intertwines_with_mass=True,
       unknown_private_state_and_reference_retained_conditionally=True,
       declared_faults_and_two_repeats_preserved=True,
       finite_mass_redshift_and_time_window_jointly_bounded=True,
       old_field_tables_reused_with_explicit_transport_bound=True,
       all_access_instruments_and_full_parent_matching_not_claimed=True,
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
