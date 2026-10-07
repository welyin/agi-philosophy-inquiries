"""956 delivery audit, separate from recalculation of the scientific witness."""
from pathlib import Path
import argparse, ast, hashlib, json, re, sys
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_956_checks.json'

def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,956):
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
    layout=Layout().verify();r=read(HERE/'native_material_interface_results.json')
    assert r['round']==956 and r['all_scientific_checks_passed']
    for rel,digest in r['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    m=r['material'];v=r['verification'];scope=r['scope']
    assert .035 < m['charge_probability_singlet'] < .036
    assert m['charge_probability_triplet']==0
    assert v['exact_H_intertwining_error']<1e-12
    assert v['U_source_error']<1e-12 and v['hopping_source_error']<1e-12
    assert v['old_434_logical_readout_transport_error']<1e-12 and v['SU2_commutator']<1e-12
    assert m['static_U_curvature']<-.064 and m['bare_compressed_U_curvature']==0
    assert .068 < m['unread_charge_instrument_leakage_singlet'] < .070
    assert m['omitted_identity_source_error']>.025
    assert r['old_434_match']['all_14_weights_positive']
    assert not r['old_434_match']['full_array_high_order_and_long_time_transport_certified']
    assert scope['native_electronic_exchange_charge_and_static_sources_common']
    assert not scope['original_955_bounds_transferred']
    assert not scope['spin_only_postmeasurement_closure_claimed']
    assert not scope['all_six_protocols_autonomously_implemented']
    assert not scope['material_to_full_SM_and_GR_error_certified']
    assert not scope['full_goal_completed']
    note=STAGE/'research_note_956.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==20 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,11)]
    for phrase in ('整体目标未完成','Hubbard','停止','SU(2)','静态'):
        assert phrase in prose,phrase
    newdocs=[note,HERE/'drafts/adoption_comparison.md',STAGE/'957/drafts/STATUS.md']
    next_prose=newdocs[-1].read_text('utf-8')
    assert '统一的对象与合同' in next_prose and 'A1' in next_prose
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split('## 六条共同协议：全局缺口对应与检验优先级（截至956）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=956)==list(range(1,957))
    assert '001—956轮共956份' in nav[0].read_text('utf-8-sig')
    assert '231—956的726份' in nav[5].read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '956：原生电子材料的交换、记录与来源' in s and '957/drafts/STATUS.md' in s
    oldfiles=[STAGE/'955/research_round_955_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'native_material_interface.py',HERE/'native_material_interface_results.json',HERE/'drafts/update_navigation956.py']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(oldfiles[0])
    out=dict(round=956,date='2026-10-07',all_delivery_checks_passed=True,
       fresh_test_groups=1,formal_reports=956,
       cumulative_numbered_test_groups_from_955=previous['cumulative_numbered_test_groups_from_954']+1,
       historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
       common_material_exchange_charge_and_source_verified=True,
       old_434_relation_code_transported_to_charge_effect=True,
       actual_charge_poststate_boundary_verified=True,
       scalar_offset_and_virtual_response_negative_controls=True,
       no_full_array_protocol_or_parent_physics_claimed=True,
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
