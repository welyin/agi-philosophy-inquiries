"""939 delivery and evidence checks, separate from scientific recalculation."""
from pathlib import Path
import argparse, ast, hashlib, json, re, sys
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_939_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,939):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):
        add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'receiver_relational_access_results.json')
    assert result['round']==939 and result['all_scientific_checks_passed']
    for rel,digest in result['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    assert result['full_CAR']['original_record_charge_commutator']>.9
    assert result['full_CAR']['relational_record_charge_commutator']<1e-13
    assert result['neutral_pointer_dilation']['unitary_error']<1e-13
    assert result['neutral_pointer_dilation']['total_flavour_conservation_error']<1e-13
    assert result['exact_visibility_ratio']=='1/8'
    assert result['exact_hbar_cubed_leading_coefficient']=='-9/15680'
    assert abs(result['relational_parity_content_difference'])>1e-5
    assert result['maximum_joint_probability_difference']>2e-6
    assert result['shared_reference_state_is_globally_invariant']
    assert result['old_853_formal_writing_theorem_preserved']
    for k in ('full_original_action_generates_readout_proved','continuum_finite_coupling_or_full_backreaction_computed','full_goal_completed'):assert not result[k]
    note=STAGE/'research_note_939.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==20 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,11)]
    for phrase in ('整体目标未完成','不是证明原S_rec已经产生U_read','没有证明所有内容均不可读','不继续本参考的精度、寿命或器件设计'):assert phrase in prose,phrase
    newdocs=[note,HERE/'drafts/common_model_identity_working.md',HERE/'drafts/common_model_recovery_map.md',STAGE/'940/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    ledger=nav[-1].read_text('utf-8-sig')
    mapping=ledger.split('## 六条共同协议：全局缺口对应与检验优先级（截至939）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    links=0
    for doc in newdocs+nav:
        s=doc.read_text('utf-8-sig')
        content=re.sub(r'\$\$.*?\$\$','',s,flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',content):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if m and p.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=939)==list(range(1,940))
    assert '001—939轮共939份' in nav[0].read_text('utf-8-sig')
    assert '231—939的709份' in nav[5].read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig')
        assert '939：共同恢复表与接收记录的有限内部核验' in s and '940/drafts/STATUS.md' in s
    oldfiles=[STAGE/'938/research_round_938_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'receiver_relational_access.py',HERE/'receiver_relational_access_results.json',HERE/'drafts/update_navigation939.py']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(oldfiles[0])
    out=dict(round=939,date='2026-10-07',all_delivery_checks_passed=True,
       fresh_test_groups=1,formal_reports=939,
       cumulative_numbered_test_groups_from_938=previous['cumulative_numbered_test_groups_from_937']+1,
       historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
       whole_model_recovery_map_completed=True,old_action_and_preparation_used=True,
       full_CAR_not_confused_with_one_particle_qubit=True,
       neutral_readout_restriction_and_relational_positive_witness_separated=True,
       original_source_distribution_reused=True,
       added_reference_preparation_and_readout_not_claimed_free=True,
       reference_device_refinement_stopped=True,
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
